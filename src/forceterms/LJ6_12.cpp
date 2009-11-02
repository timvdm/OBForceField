/*********************************************************************
Non-bonded Lennard-Jones Term 

Copyright (C) 2009 by Frank Peters
 
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

#include "LJ6_12.h"
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
 
    const std::string LJ6_12::m_name = "Lennard-Jones 6-12";

    template<> void LJ6_12::Mix<LJ6_12::geometric>(double & sigma, double & epsilon, const double & sigma_1,  const double & epsilon_1,  const double & sigma_2,  const double & epsilon_2)
    {
      sigma = sqrt(sigma_1 * sigma_2);
      epsilon = sqrt(epsilon_1 * epsilon_2);
    }

    template<> void LJ6_12::Mix<LJ6_12::arithmetic>(double & sigma, double & epsilon, const double & sigma_1,  const double & epsilon_1,  const double & sigma_2,  const double & epsilon_2)
    {
      sigma = 0.5*(sigma_1 + sigma_2);
      epsilon = sqrt(epsilon_1 * epsilon_2);
    }

    template<> void LJ6_12::Mix<LJ6_12::sixthpower>(double & sigma, double & epsilon, const double & sigma_1,  const double & epsilon_1,  const double & sigma_2,  const double & epsilon_2)
    {
      sigma = pow(0.5*(pow(sigma_1, 6.0) + pow(sigma_2, 6.0)), 1./6.);
      epsilon = (2 * sqrt(epsilon_1 * epsilon_2) * pow(sigma_1, 3.0) * pow(sigma_2, 3.0))/(pow(sigma_1, 6.0) + pow(sigma_2, 6.0));
    }

    LJ6_12::LJ6_12(OBFunction *function, const double factorOneFour, const LJ6_12::MixingRule rule, const std::string tableName)
      : OBFunctionTerm(function), m_tableName(tableName), m_value(999999.99), m_calcs(NULL), m_i(NULL), m_numPairs(0), m_factorOneFour(factorOneFour) 
    {
      switch (rule)
	{
	case geometric: m_Mix = & LJ6_12::Mix<geometric>; break; 
	case arithmetic: m_Mix = & LJ6_12::Mix<arithmetic>; break; 
	case sixthpower: m_Mix = & LJ6_12::Mix<sixthpower>; break; 
	}
    }

    LJ6_12::~LJ6_12() 
    {
      delete [] m_i;
      delete [] m_calcs;
    }

    // Conventions taken from Lammps potentials
    // E = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 )
    // epsilon (energy)
    // sigma (distance)

    void LJ6_12::Compute(OBFunction::Computation computation)
    {
      m_value = 0.0;
      double rab, term, term3, term6, term12, e;
      Eigen::Vector3d Fa, Fb;

      if (computation == OBFunction::Gradients) {
	size_t ia, ib;
	double dE;
	for (unsigned int i = 0; i < m_numPairs; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  rab = VectorBondDerivative(m_function->GetPositions()[ia], m_function->GetPositions()[ib], Fa, Fb);
	  term = m_calcs[i].sigma / rab;
	  term3 = term * term * term;
	  term6 = term3 * term3;
	  term12 = term6 * term6;
	  dE = 24.* m_calcs[i].epsilon * (-2.0*term12 + term6)/rab;
	  Fa *= dE;
	  Fb *= dE;
	  m_function->GetGradients()[ia] += Fa;
	  m_function->GetGradients()[ib] += Fb;
	  e = 4.0 * m_calcs[i].epsilon * (term12-term6);
	  m_value += e;
	}
      }      
      else {
	for (unsigned int i = 0; i < m_numPairs; ++i) {
	  const Eigen::Vector3d ab = m_function->GetPositions()[m_i[i].iA] - m_function->GetPositions()[m_i[i].iB];
	  rab = ab.norm();
	  term = m_calcs[i].sigma / rab;
	  term3 = term*term*term;
	  term6 = term3*term3;
	  term12 = term6*term6;
	  e = 4.0 * m_calcs[i].epsilon * (term12-term6);
	  m_value += e;      
	}
      }
    }
  
    bool LJ6_12::Setup()
    {
      // combine the typing stored in obfftype with the parameters from the parameter database
      OBParameterDBTable * pTable = ((m_function->GetParameterDB())->GetTable(m_tableName));
      OBFFType * pOBFFType(m_function->GetOBFFType());
      vector<OBFFType::AtomIdentifier> atoms(pOBFFType->GetAtoms());
      vector<OBParameterDBTable::Query> query;
      vector<OBVariant> row;
      Parameter parameter;
      string name;
      map<string,Parameter> parameters;
      map<string,Parameter>::iterator itr;
      Index i;
      vector <Index> v_i;
      vector <Parameter> v_calcs;
      OBAtom *a, *b;

      if ( (pTable==NULL) || (pOBFFType==NULL) )
	return false;

      double sigma_j, sigma_k, epsilon_j, epsilon_k;
      for(unsigned int j=0; j != atoms.size(); ++j){
	for(unsigned int k= j+1; k != atoms.size(); ++k){
	  if (atoms[j]<atoms[k])
	    name = atoms[j] + "-" + atoms[k];
	  else
	    name = atoms[k] + "-" + atoms[j];  
	  itr=parameters.find(name);
	  if (itr==parameters.end()){
	    query.clear();
	    query.push_back( OBParameterDBTable::Query(0, OBVariant(atoms[j])));
	    row = pTable->FindRow(query);
	    sigma_j = row.at(1).AsDouble();
	    epsilon_j = row.at(2).AsDouble();
	    query.clear();
	    query.push_back( OBParameterDBTable::Query(0, OBVariant(atoms[k])));
	    row = pTable->FindRow(query);
	    sigma_k = row.at(1).AsDouble();
	    epsilon_k = row.at(2).AsDouble();
	    (*m_Mix)(parameter.sigma, parameter.epsilon, sigma_j, epsilon_j, sigma_k, epsilon_k);
	    parameters.insert(pair<string,Parameter>(name,parameter));
	  }
	  else
	    parameter=itr->second;
	  i.iA = j;
	  i.iB = k;
	  if (pOBFFType->IsConnected(i.iA, i.iB))
	    continue;
	  if (pOBFFType->IsOneThree(i.iA, i.iB))
	    continue;
	  if (pOBFFType->IsOneFour(i.iA, i.iB))
	    parameter.epsilon *= m_factorOneFour;
	  v_i.push_back(i);
	  v_calcs.push_back(parameter);
	}
      }

      m_numPairs = v_i.size();
      delete [] m_i;
      delete [] m_calcs;
      m_i = new Index [m_numPairs];
      m_calcs = new Parameter [m_numPairs];
      for (size_t i=0; i< m_numPairs; ++i){
	m_i[i] = v_i[i];
	m_calcs[i] = v_calcs[i];
      }
      return true;
    }
  }
} // end namespace OpenBabel

