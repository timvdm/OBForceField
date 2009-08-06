/*********************************************************************
MMFF94ElectroTerm - MMFF94 force field bond stratching term

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

#include "electro.h"
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { Class, TypeA, TypeB, TypeC, Ka, Theta0 };
  };
 
  MMFF94ElectroTerm::MMFF94ElectroTerm(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
   
  void MMFF94ElectroTerm::Compute(OBFunction::Computation computation)
  {
    cout << "MMFF94ElectroTerm::Compute" << endl;
    cout << "m_calcs.size = " << m_calcs.size() << endl;
    m_value = 0.0;
 
    unsigned int idxA, idxB;
    double rab, e;
    Eigen::Vector3d Fa, Fb;
    for (int i = 0; i < (int)m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;
  
      if (computation == OBFunction::Gradients) {
        rab = VectorDistanceDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], Fa, Fb);
        rab += 0.05; // ??
        const double rab2 = rab * rab;
        const double dE = - (m_calcs[i].qq / rab2);
        Fa *= dE;
        Fb *= dE;
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
      } else {
        const Eigen::Vector3d ab = m_function->GetPositions()[idxA] - m_function->GetPositions()[idxB];
        rab = ab.norm();
        rab += 0.05; // ??
      }

      e = m_calcs[i].qq / rab;
      m_value += e;
     
     /* 
      if (m_log->IsHigh()) {
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d   %8.3f  %8.3f  %8.3f  %8.3f\n", 
                atoi(m_calcs[i].a->GetType()), atoi(m_calcs[i].b->GetType()), 
                rab, m_calcs[i].a->GetPartialCharge(), 
                m_calcs[i].b->GetPartialCharge(), e);
        m_log->Log(_logbuf);
      }
      */
    }
    /*
    if (m_log->IsMedium()) {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL ELECTROSTATIC ENERGY = %8.5f %s\n", energy, GetUnit().c_str());
      m_log->Log(_logbuf);
    }
    */
 
    cout << "E_ele = " << m_value << endl;
 
  }
  
  bool MMFF94ElectroTerm::Setup(/*const*/ OBMol &mol)
  {
    OBAtom *a, *b, *c;

    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();
 
    // 
    // Electrostatic Calculations
    //
    if (logFile->IsLow())
      logFile->Write("SETTING UP ELECTROSTATIC CALCULATIONS...\n");
 
    Calculation elecalc;
    m_calcs.clear();
    
    FOR_PAIRS_OF_MOL(p, mol) {
      a = mol.GetAtom((*p)[0]);
      b = mol.GetAtom((*p)[1]);
      
      elecalc.qq = 332.0716 * m_common->GetPartialCharge(a->GetIdx()-1) * m_common->GetPartialCharge(b->GetIdx()-1);
      
      if (elecalc.qq) {
        elecalc.idx1 = a->GetIdx() - 1;
        elecalc.idx2 = b->GetIdx() - 1;
        
        // 1-4 scaling
        if (a->IsOneFour(b))
          elecalc.qq *= 0.75;
	  
        m_calcs.push_back(elecalc);
      }
    }

    return true;


  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
