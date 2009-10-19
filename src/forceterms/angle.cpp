/*********************************************************************
MMFF94AngleTerm - MMFF94 force field angle stratching term

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

#include "angle.h"
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
 
    const std::string AngleHarmonic::m_name = "Angle Harmonic";

    AngleHarmonic::AngleHarmonic(OBFunction *function, std::string tableName)
      : OBFunctionTerm(function), m_tableName(tableName), m_i(0), m_calcs(0),
      m_value(999999.99) 
    {
    }


    AngleHarmonic::~AngleHarmonic()
    {
      if (m_i)
        delete [] m_i;
      if (m_calcs)
        delete [] m_calcs;
    }

    // Conventions taken from Lammps potentials

    // E = K (theta - theta0)^2
    // K (energy/angle^2)
    // theta0 (angle)

    void AngleHarmonic::Compute(OBFunction::Computation computation)
    {
      m_value = 0.0;
      double theta, delta, delta2, e;
	
      if (computation == OBFunction::Gradients) {
	size_t ia, ib, ic;
	Eigen::Vector3d Fa, Fb, Fc;
	double dE;
	for (unsigned int i = 0; i < m_numAngles; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  ic = m_i[i].iC;
	  theta = VectorAngleDerivative(m_function->GetPositions()[ia], m_function->GetPositions()[ib], m_function->GetPositions()[ic], Fa, Fb, Fc); 
	  delta = DEG_TO_RAD * (theta - m_calcs[i].theta0);
	  if (!isfinite(theta))
	    theta = 0.0;
	  dE = 2.0 * m_calcs[i].K * delta;
	  Fa *= dE;
	  Fb *= dE;
	  Fc *= dE;
	  m_function->GetGradients()[ia] += Fa;
	  m_function->GetGradients()[ib] += Fb;
	  m_function->GetGradients()[ic] += Fc;
	  delta2 = delta * delta;
	  e = m_calcs[i].K * delta2;
	  m_value += e;
	}
      } else {
	Eigen::Vector3d ab, bc;
	for (unsigned int i = 0; i < m_numAngles; ++i) {
	  ab = m_function->GetPositions()[m_i[i].iA] - m_function->GetPositions()[m_i[i].iB];
	  bc = m_function->GetPositions()[m_i[i].iC] - m_function->GetPositions()[m_i[i].iB];
	  theta = VectorAngle(ab, bc);
	  if (!isfinite(theta))
	    theta = 0.0;
	  delta = DEG_TO_RAD * (theta - m_calcs[i].theta0);
	  delta2 = delta * delta;
	  e = m_calcs[i].K * delta2;
	  m_value += e;
	}
      }
    }
  
    bool AngleHarmonic::Setup()
    {
      // combine the typing stored in obfftype with the parameters from the parameter database
      OBParameterDBTable * pTable = ((m_function->GetParameterDB())->GetTable(m_tableName));
      OBFFType * pOBFFType(m_function->GetOBFFType());
      vector<OBFFType::AngleIdentifier> angles(pOBFFType->GetAngles());
      vector<OBParameterDBTable::Query> query;
      std::vector<OBVariant> row;
      Parameter parameter;
      map<string,Parameter> parameters;
      map<string,Parameter>::iterator itr;

      if ( (pTable==NULL) || (pOBFFType==NULL) )
	return false;
      m_numAngles=angles.size();
      if (m_i)
        delete [] m_i;
      if (m_calcs)
        delete [] m_calcs;
      m_i = new Index [m_numAngles];
      m_calcs = new Parameter [m_numAngles];
      for(unsigned int i=0;i != m_numAngles;++i){
	itr=parameters.find(angles[i].name);
	if (itr==parameters.end()){
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(angles[i].name)));
	  row = pTable->FindRow(query);
	  parameter.K = row.at(4).AsDouble();
	  parameter.theta0 = row.at(5).AsDouble();
	  parameters.insert(pair<string,Parameter>(angles[i].name,parameter));
	}
	else {
	  parameter=itr->second;
	}
	m_i[i].iA = angles[i].iA;
	m_i[i].iB  = angles[i].iB;
	m_i[i].iC  = angles[i].iC;
	m_calcs[i] = parameter;
      }
      return true;
    }
  }
} // end namespace OpenBabel

