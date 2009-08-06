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
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { Class, TypeA, TypeB, Kb, R0 };
  };
  struct EmpiricalRule {
    enum { Species1, Species2, R0Ref, KbRef };
  };
 
  MMFF94BondTerm::MMFF94BondTerm(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
      
  // 
  // MMFF part I - page 494
  //      
  //                   kb_ij                              7
  // EB_ij = 143.9325 ------- /\r_ij^2 (1 + cs /\_rij + ---- cs^2 r_ij^2)
  //                     2                               12
  //
  // kb_ij	force constant (md/A)
  //
  // /\r_ij 	r_ij - r0_ij (A)
  //
  // cs		cubic stretch constant = -2 A^(-1)
  //
  void MMFF94BondTerm::Compute(OBFunction::Computation computation)
  {
    m_value = 0.0;
    
    unsigned int idxA, idxB;
    double rab, delta, delta2, e;
    Eigen::Vector3d Fa, Fb;
    for (unsigned int i = 0; i < m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;

      if (computation == OBFunction::Gradients) {
        rab = VectorBondDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], Fa, Fb);
        delta = rab - m_calcs[i].r0;
        delta2 = delta * delta;
      
        const double dE = 143.9325 * m_calcs[i].kb * delta * (1.0 - 3.0 * delta + 14.0/3.0 * delta2);
      
        Fa *= dE;
        Fb *= dE;
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
      } else {
        const Eigen::Vector3d ab = m_function->GetPositions()[idxA] - m_function->GetPositions()[idxB];
        rab = ab.norm();
        delta = rab - m_calcs[i].r0;
        delta2 = delta * delta;
      }
    
      e = m_calcs[i].kb * delta2 * (1.0 - 2.0 * delta + 7.0/3.0 * delta2);
      m_value += e;      
    }
   
    m_value *= 143.9325 * 0.5;
  }
  
  bool MMFF94BondTerm::Setup(/*const*/ OBMol &mol)
  {
    OBAtom *a, *b;
    bool found;

    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();
    
    if (logFile->IsLow())
      logFile->Write("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");
    
    // 
    // Bond Calculations
    //
    // no "step-down" procedure
    // MMFF part V - page 625 (empirical rule)
    //
    if (logFile->IsLow())
      logFile->Write("SETTING UP BOND CALCULATIONS...\n");
 
    Calculation bondcalc;
    m_calcs.clear();
    
    FOR_BONDS_OF_MOL(bond, mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();
     
      int bondType = m_common->GetBondType(a, b);
      int type_a = m_common->GetCachedType(a);
      int type_b = m_common->GetCachedType(b);
     
      std::vector<OBParameterDB::Query> query;
      query.push_back( OBParameterDB::Query(Parameter::Class, OBVariant(bondType)) );
      query.push_back( OBParameterDB::Query(Parameter::TypeA, OBVariant(type_a), true) );
      query.push_back( OBParameterDB::Query(Parameter::TypeB, OBVariant(type_b), true) );
      std::vector<OBVariant> row = database->FindRow(MMFF94SimpleParameterDB::BondParameters, query);
      if (row.size()) {
        bondcalc.idx1 = a->GetIdx() - 1;
        bondcalc.idx2 = b->GetIdx() - 1;
        bondcalc.kb = row.at(Parameter::Kb).AsDouble();
        bondcalc.r0 = row.at(Parameter::R0).AsDouble();
        bondcalc.bondType = bondType;

        m_calcs.push_back(bondcalc);
      } else {
        query.clear();
        query.push_back( OBParameterDB::Query(Parameter::TypeA, OBVariant(static_cast<int>(a->GetAtomicNum())), true) );
        query.push_back( OBParameterDB::Query(Parameter::TypeB, OBVariant(static_cast<int>(b->GetAtomicNum())), true) );
        row = database->FindRow(MMFF94SimpleParameterDB::EmpiricalBondRules, query);
        if (row.empty()) { 
          if (logFile->IsLow()) {
            // This should never happen
            std::stringstream ss;
            ss << "    COULD NOT FIND PARAMETERS FOR BOND " << m_common->GetCachedType(a) << "-" << m_common->GetCachedType(b) << " (MMFF94 Atom Types)..." << endl;
            logFile->Write(ss.str());
          }
          return false;
        } else {
          if (logFile->IsLow()) {
            std::stringstream ss;
            ss << "    USING EMPIRICAL RULE FOR BOND STRETCHING " << m_common->GetCachedType(a) << "-" << m_common->GetCachedType(b) << " (MMFF94 Atom Types)..." << endl;
            logFile->Write(ss.str());
          }

          double rr, rr2, rr4, rr6;
          bondcalc.idx1 = a->GetIdx() - 1;
          bondcalc.idx2 = b->GetIdx() - 1;
          bondcalc.r0 = m_common->GetRuleBondLength(a, b); 
	  
          rr = row.at(EmpiricalRule::R0Ref).AsDouble() / bondcalc.r0;
          rr2 = rr * rr;
          rr4 = rr2 * rr2;
          rr6 = rr4 * rr2;
	  
          bondcalc.kb = row.at(EmpiricalRule::KbRef).AsDouble() * rr6;
          bondcalc.bondType = bondType;
          
          m_calcs.push_back(bondcalc);
        }
      }
    }

  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
