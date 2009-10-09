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

#include "Coulomb.h"
#include <OBFFType>
#include <OBParameterDB>
#include <OBChargeMethod>
#include <OBFunction>
#include <OBFunctionTerm>

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
  namespace OBFFs {
 
    const std::string Coulomb::m_name = "Coulomb";

    Coulomb::Coulomb(OBFunction *function, const double factorOneFour, const double relativePermittivity)
      : OBFunctionTerm(function), m_value(999999.99), m_calcs(NULL), m_i(NULL), m_numPairs(0), m_factorOneFour(factorOneFour), m_relativePermittivity(relativePermittivity) {}

    Coulomb::~Coulomb() 
    {
      delete [] m_i;
      delete [] m_calcs;
    }

    // Conventions taken from Lammps potentials
    // E = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 )
    // epsilon (energy)
    // sigma (distance)

    void Coulomb::Compute(OBFunction::Computation computation)
    {
      m_value = 0.0;
      unsigned int ia, ib;
      double rab, term, e;
      Eigen::Vector3d Fa, Fb;

      if (computation == OBFunction::Gradients) {
	double dE;
	for (unsigned int i = 0; i < m_numPairs; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  rab = VectorBondDerivative(m_function->GetPositions()[ia], m_function->GetPositions()[ib], Fa, Fb);
	  term = 1.0 / rab;
	  e = m_calcs[i].qq * term;
	  dE = - e * term;
	  Fa *= dE;
	  Fb *= dE;
	  m_function->GetGradients()[ia] += Fa;
	  m_function->GetGradients()[ib] += Fb;
	  m_value += e;
	}
      }      
      else {
	for (unsigned int i = 0; i < m_numPairs; ++i) {
	  ia = m_i[i].iA;
	  ib = m_i[i].iB;
	  const Eigen::Vector3d ab = m_function->GetPositions()[ia] - m_function->GetPositions()[ib];
	  rab = ab.norm();
	  e =  m_calcs[i].qq / rab;
	  m_value +=  e;
	}
      }
    }
  
    bool Coulomb::Setup()
    {
      OBChargeMethod * pOBChargeMethod(m_function->GetOBChargeMethod());
      OBFFType * pOBFFType(m_function->GetOBFFType());
      Parameter parameter;
      Index i;
      vector <Index> v_i;
      vector <Parameter> v_calcs;
      const double factor = 332.0716 / m_relativePermittivity; // energy scale: kcal/mol

      if ( (pOBFFType==NULL) || (pOBChargeMethod==NULL))
	return false;

      const vector<double> & partialCharge = (pOBChargeMethod->GetPartialCharges());

      for(unsigned int j=0; j != partialCharge.size();++j){
	for(unsigned int k= j+1 ;k != partialCharge.size();++k){
	  i.iA = j;
	  i.iB = k;
	  if (pOBFFType->IsConnected(i.iA, i.iB))
	    continue;
	  if (pOBFFType->IsOneThree(i.iA, i.iB))
	    continue;
	  parameter.qq = factor * partialCharge[j] * partialCharge[k];
	  if (pOBFFType->IsOneFour(i.iA, i.iB))
	    parameter.qq *= m_factorOneFour;
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

