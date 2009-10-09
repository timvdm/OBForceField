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

#include "bond.h"
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
 
    const std::string BondHarmonic::m_name = "Bond Harmonic";

    BondHarmonic::BondHarmonic(OBFunction *function, std::string tableName)
      : OBFunctionTerm(function), m_tableName(tableName), m_value(999999.99), m_calcs(NULL), m_i(NULL), m_numBonds(0) {}

    BondHarmonic::~BondHarmonic() 
    {
      delete [] m_i;
      delete [] m_calcs;
    }


    // Conventions taken from Lammps potentials

    // E = K * (r-r0)^2
    // K (energy/distance^2)
    // r0 (distance)

    void BondHarmonic::Compute(OBFunction::Computation computation)
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
	  dE = 2.0 * m_calcs[i].K * delta;
	  Fa *= dE;
	  Fb *= dE;
	  m_function->GetGradients()[ia] += Fa;
	  m_function->GetGradients()[ib] += Fb;
	  e = m_calcs[i].K * delta2;
	  m_value += e;
	}
      }      
      else {
	for (unsigned int i = 0; i < m_numBonds; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  const Eigen::Vector3d ab = m_function->GetPositions()[ia] - m_function->GetPositions()[ib];
	  rab = ab.norm();
	  delta = rab - m_calcs[i].r0;
	  delta2 = delta * delta;
	  e = m_calcs[i].K * delta2;
	  m_value += e;      
	}
      }
    }
  
    bool BondHarmonic::Setup()
    {
      // combine the typing stored in obfftype with the parameters from the parameter database
      OBParameterDBTable * pTable = ((m_function->GetParameterDB())->GetTable(m_tableName));
      OBFFType * pOBFFType(m_function->GetOBFFType());
      vector<OBFFType::BondIdentifier> bonds(pOBFFType->GetBonds());
      vector<OBParameterDBTable::Query> query;
      vector<OBVariant> row;
      Parameter parameter;
      map<string,Parameter> parameters;
      map<string,Parameter>::iterator itr;

      if ( (pTable==NULL) || (pOBFFType==NULL) )
	return false;
      m_numBonds=bonds.size();
      delete [] m_i;
      delete [] m_calcs;
      m_i = new Index [m_numBonds];
      m_calcs = new Parameter [m_numBonds];
      for(unsigned int i=0;i != m_numBonds;++i){
	itr=parameters.find(bonds[i].name);
	if (itr==parameters.end()){
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(bonds[i].name)));
	  row = pTable->FindRow(query);
	  parameter.K = row.at(3).AsDouble();
	  parameter.r0 = row.at(4).AsDouble();
	  parameters.insert(pair<string,Parameter>(bonds[i].name,parameter));
	}
	else {
	  parameter=itr->second;
	}
	m_i[i].iA = bonds[i].iA;
	m_i[i].iB  = bonds[i].iB;
	m_calcs[i] = parameter;
      }
      return true;
    }

    const std::string BondClass2::m_name = "Bond Class 2";

    BondClass2::BondClass2(OBFunction *function, std::string tableName)
      : OBFunctionTerm(function), m_tableName(tableName), m_value(999999.99), m_calcs(NULL), m_i(NULL), m_numBonds(0) {}


    BondClass2::~BondClass2() 
    {
      delete [] m_i;
      delete [] m_calcs;
    }

    // Conventions taken from Lammps potentials

    // E = K * (r-r0)^2
    // K (energy/distance^2)
    // r0 (distance)

    void BondClass2::Compute(OBFunction::Computation computation)
    {
      m_value = 0.0;
      unsigned int ia, ib;
      double rab, delta, delta2, e, dE;
      Eigen::Vector3d Fa, Fb;

      if (computation == OBFunction::Gradients) {
	for (unsigned int i = 0; i < m_numBonds; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  rab = VectorBondDerivative(m_function->GetPositions()[ia], m_function->GetPositions()[ib], Fa, Fb);
	  delta = rab - m_calcs[i].r0;
	  delta2 = delta * delta;
	  dE = delta * (2.0 * m_calcs[i].K2 + 3.0 * m_calcs[i].K3 * delta + 4.0 * m_calcs[i].K4 * delta2);
	  Fa *= dE;
	  Fb *= dE;
	  m_function->GetGradients()[ia] += Fa;
	  m_function->GetGradients()[ib] += Fb;
	  e = delta2 * (m_calcs[i].K2 + m_calcs[i].K3 * delta + m_calcs[i].K4 * delta2);
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
	  e = delta2 * (m_calcs[i].K2 + m_calcs[i].K3 * delta + m_calcs[i].K4 * delta2);
	  m_value += e;      
	}
      }
    }
  
    bool BondClass2::Setup()
    {
      // combine the typing stored in obfftype with the parameters from the parameter database
      OBParameterDBTable * pTable = ((m_function->GetParameterDB())->GetTable(m_tableName));
      OBFFType * pOBFFType(m_function->GetOBFFType());
      vector<OBFFType::BondIdentifier> bonds(pOBFFType->GetBonds());
      vector<OBParameterDBTable::Query> query;
      std::vector<OBVariant> row;
      Parameter parameter;
      map<string,Parameter> parameters;
      map<string,Parameter>::iterator itr;

      if ( (pTable==NULL) || (pOBFFType==NULL) )
	return false;
      m_numBonds=bonds.size();
      delete [] m_i;
      delete [] m_calcs;
      m_i = new Index [m_numBonds];
      m_calcs = new Parameter [m_numBonds];
      for(unsigned int i=0;i != m_numBonds;++i){
	itr=parameters.find(bonds[i].name);
	if (itr==parameters.end()){
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(bonds[i].name)));
	  row = pTable->FindRow(query);
	  parameter.K2 = row.at(3).AsDouble();
	  parameter.K3 = row.at(4).AsDouble();
	  parameter.K4 = row.at(5).AsDouble();
	  parameter.r0 = row.at(6).AsDouble();
	  parameters.insert(pair<string,Parameter>(bonds[i].name,parameter));
	}
	else {
	  parameter=itr->second;
	}
	m_i[i].iA = bonds[i].iA;
	m_i[i].iB  = bonds[i].iB;
	m_calcs[i] = parameter;
      }
      return true;
    }
 
  }
} // end namespace OpenBabel

