/*********************************************************************
MMFF94BondTerm - MMFF94 force field bond stratching term

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

#include "bondcubicharmonic.h"
#include <OBFFType>
#include <OBParameterDB>
#include <OBFunction>
#include <OBFunctionTerm>

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

#include <map>

using namespace std;

namespace OpenBabel {
  namespace OBFFs {
 
    BondCubicHarmonicTerm::BondCubicHarmonicTerm(OBFunction *function, double prefactor, double cs, 
        double cs2, const std::string &tableName, int forceConstantColumn, int bondLengthColumn)
      : OBFunctionTerm(function), m_tableName(tableName), m_forceConstantColumn(forceConstantColumn),
      m_bondLengthColumn(bondLengthColumn), m_prefactor(prefactor), m_cs(cs), m_cs2(cs2), 
      m_value(999999.99), m_calcs(0), m_i(0), m_numBonds(0)
    {
    }

    BondCubicHarmonicTerm::~BondCubicHarmonicTerm() 
    {
      delete [] m_i;
      delete [] m_calcs;
    }

    void BondCubicHarmonicTerm::Compute(OBFunction::Computation computation)
    {
      m_value = 0.0;
      unsigned int ia, ib;
      double rab, delta, delta2, e;
      Eigen::Vector3d Fa, Fb;

      if (computation == OBFunction::Gradients) {
	double dE;
	for (unsigned int i = 0; i < m_numBonds; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  rab = VectorBondDerivative(m_function->GetPositions()[ia], m_function->GetPositions()[ib], Fa, Fb);
	  delta = rab - m_calcs[i].r0;
	  delta2 = delta * delta;
	  dE = m_prefactor * m_calcs[i].K * delta * (1.0 + 3.0 * m_cs * delta + 4.0 * m_cs2 * delta2);
	  Fa *= dE;
	  Fb *= dE;
	  m_function->GetGradients()[ia] += Fa;
	  m_function->GetGradients()[ib] += Fb;
	  e = m_prefactor * m_calcs[i].K * delta2 * (1.0 + m_cs * delta + m_cs2 * delta2);
	  m_value += e;
	}
      } else {
	for (unsigned int i = 0; i < m_numBonds; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  const Eigen::Vector3d ab = m_function->GetPositions()[ia] - m_function->GetPositions()[ib];
	  rab = ab.norm();
	  delta = rab - m_calcs[i].r0;
	  delta2 = delta * delta;
	  e = m_prefactor * m_calcs[i].K * delta2 * (1.0 + m_cs * delta + m_cs2 * delta2);
	  m_value += e;      
	}
      }
    }
 
    bool BondCubicHarmonicTerm::Setup()
    {
      // combine the typing stored in obfftype with the parameters from the parameter database
      OBParameterDBTable *table = m_function->GetParameterDB()->GetTable(m_tableName);
      if (!table)
        return false;

      OBFFType *obfftype = m_function->GetOBFFType();
      if (!obfftype)
        return false;

      vector<OBFFType::BondIdentifier> bonds(obfftype->GetBonds());
      vector<OBParameterDBTable::Query> query;
      vector<OBVariant> row;
      Parameter parameter;
      map<string,Parameter> parameters;
      map<string,Parameter>::iterator itr;

      m_numBonds = bonds.size();
      delete [] m_i;
      delete [] m_calcs;
      m_i = new Index [m_numBonds];
      m_calcs = new Parameter [m_numBonds];
      for (unsigned int i = 0; i != m_numBonds; ++i) {
	itr = parameters.find(bonds[i].name);
	if (itr == parameters.end()){
          // create the parameter
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(bonds[i].name)));
	  row = table->FindRow(query);
          if (!row.size()) {
            OBLogFile *logFile = m_function->GetLogFile();
            std::stringstream ss;
            ss << "Could not find parameters for bond with name: " << bonds[i].name << endl;
            logFile->Write(ss.str());
            return false;
          }
	  parameter.K = row.at(m_forceConstantColumn).AsDouble();
	  parameter.r0 = row.at(m_bondLengthColumn).AsDouble();
	  parameters.insert(pair<string,Parameter>(bonds[i].name,parameter));
	} else {
          // this parameter already exists
	  parameter = itr->second;
	}
          
        
        // store the calculation
	m_i[i].iA = bonds[i].iA;
	m_i[i].iB  = bonds[i].iB;
	m_calcs[i] = parameter;
      }
      return true;
    }
 
  }
} // end namespace OpenBabel

