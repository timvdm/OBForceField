/*********************************************************************
MMFF94TorsionTerm - MMFF94 force field torsion stratching term

Copyright (C) 2006-2008,2009 by Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include "torsion.h"
#include <OBFFType>
#include <OBParameterDB>
#include <OBFunction>
#include <OBFunctionTerm>

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
  namespace OBFFs {
 
    const std::string TorsionHarmonic::m_name = "Torsion Harmonic";

    TorsionHarmonic::TorsionHarmonic(OBFunction *function, std::string tableName)
      : OBFunctionTerm(function), m_tableName(tableName), m_value(999999.99), m_calcs(NULL), m_i(NULL), m_numTorsions(0) {}


    TorsionHarmonic::~TorsionHarmonic()
    {
      delete [] m_i;
      delete [] m_calcs;
    }

    // Conventions taken from Lammps potentials

    // E = K (1 + d * cos(n * phi)) 
    // K (energy/torsion^2)
    // phi0 (torsion)

    void TorsionHarmonic::Compute(OBFunction::Computation computation)
    {
      m_value = 0.0;
      double phi, delta, delta2, e, cosine;
	
      if (computation == OBFunction::Gradients) {
	size_t ia, ib, ic, id;
	Eigen::Vector3d Fa, Fb, Fc, Fd;
	double dE, sine;
	for (unsigned int i = 0; i < m_numTorsions; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  ic = m_i[i].iC;
	  id = m_i[i].iD;
	  phi = VectorTorsionDerivative(m_function->GetPositions()[ia], m_function->GetPositions()[ib], m_function->GetPositions()[ic], m_function->GetPositions()[id], Fa, Fb, Fc, Fd); 
	  if (!isfinite(phi))
	    phi = 0.0;
	  sine = sin(DEG_TO_RAD* m_calcs[i].n * phi);	  
	  dE = m_calcs[i].K * m_calcs[i].d * m_calcs[i].n * sine;
	  Fa *= dE;
	  Fb *= dE;
	  Fc *= dE;
	  Fd *= dE;
	  m_function->GetGradients()[ia] += Fa;
	  m_function->GetGradients()[ib] += Fb;
	  m_function->GetGradients()[ic] += Fc;
	  m_function->GetGradients()[id] += Fd;

	  cosine = cos(DEG_TO_RAD * m_calcs[i].n * phi);
	  e = m_calcs[i].K * (1.0 + m_calcs[i].d * cosine);
	  m_value += e;
	}
      } else {
	for (unsigned int i = 0; i < m_numTorsions; ++i) {
	  phi = VectorTorsion(m_function->GetPositions()[m_i[i].iA], m_function->GetPositions()[m_i[i].iB],
			    m_function->GetPositions()[m_i[i].iC], m_function->GetPositions()[m_i[i].iD]);
	  if (!isfinite(phi))
	    phi = 0.0;

	  cosine = cos(DEG_TO_RAD * m_calcs[i].n * phi);
	  e = m_calcs[i].K * (1.0 + m_calcs[i].d * cosine);
	  m_value += e;
	}
      }
    }
  
    bool TorsionHarmonic::Setup()
    {
      // combine the typing stored in obfftype with the parameters from the parameter database
      OBParameterDBTable * pTable = ((m_function->GetParameterDB())->GetTable(m_tableName));
      OBFFType * pOBFFType(m_function->GetOBFFType());
      vector<OBFFType::TorsionIdentifier> torsions(pOBFFType->GetTorsions());
      vector<OBParameterDBTable::Query> query;
      std::vector< std::vector<OBVariant> > rows;
      std::vector< std::vector<OBVariant> >::const_iterator itr2;
      Index i;
      std::vector< Index > v_i;
      std::vector< Parameter > v_calcs;
      Parameter parameter;
      multimap<string,Parameter> parameters;
      multimap<string,Parameter>::const_iterator itr;
      pair< multimap<string,Parameter>::const_iterator , multimap<string,Parameter>::const_iterator > ret;

      if ( (pTable==NULL) || (pOBFFType==NULL) )
	return false;

      //The torsion potential can be a sum of terms, i.e., more than one entry per dihedral
      //Every term in such a sum will be a separate entry into m_i and m_calcs
      //Therefore we use multimap, FindRows, and do not know m_numTorsions a priori

      v_i.reserve(torsions.size());
      v_calcs.reserve(torsions.size());
      for(unsigned int j=0;j != torsions.size();++j){
	ret=parameters.equal_range(torsions[j].name);
	i.iA = torsions[j].iA;
	i.iB = torsions[j].iB;
	i.iC = torsions[j].iC;
	i.iD = torsions[j].iD;
	if (ret.first==ret.second){
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(torsions[j].name)));
	  rows = pTable->FindRows(query);
	  for(itr2 = rows.begin(); itr2 != rows.end(); ++itr2){
	    parameter.K = itr2->at(5).AsDouble();
	    parameter.d = itr2->at(6).AsDouble();
	    parameter.n = itr2->at(7).AsDouble();
	    parameters.insert(pair<string,Parameter>(torsions[j].name,parameter));
	    v_i.push_back(i);
	    v_calcs.push_back(parameter);
	  }
	}
	else {
	  for(itr = ret.first; itr != ret.second; ++itr){
	    parameter=itr->second;
	    v_i.push_back(i);
	    v_calcs.push_back(parameter);
	  }
	}
      }
      m_numTorsions = v_i.size();
      delete [] m_i;
      delete [] m_calcs;
      m_i = new Index [m_numTorsions];
      m_calcs = new Parameter [m_numTorsions];
      for (unsigned int i=0; i< m_numTorsions; ++i){
	m_i[i]=v_i[i];
	m_calcs[i]=v_calcs[i];
      }
      return true;
    }  
  }
} // end namespace OpenBabel

