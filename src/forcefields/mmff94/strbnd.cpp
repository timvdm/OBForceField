/*********************************************************************
MMFF94StrBndTerm - MMFF94 force field bond stratching term

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

#include "strbnd.h"
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { Class, TypeA, TypeB, TypeC, Kba123, Kba321 };
  };
  struct EmpiricalRule {
    enum { Row1, Row2, Row3, F123, F321 };
  };
 
  MMFF94StrBndTerm::MMFF94StrBndTerm(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
   
  void MMFF94StrBndTerm::Compute(OBFunction::Computation computation)
  {
    m_value = 0.0;
  
    unsigned int idxA, idxB, idxC;
    double theta, rab, rbc, delta_theta, delta_rab, delta_rbc, e;
    Eigen::Vector3d Fa, Fb, Fc, Fabc_a, Fabc_c;
    for (unsigned int i = 0; i < m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;
      idxC = m_calcs[i].idx3;

      if (computation == OBFunction::Gradients) {
        theta = VectorAngleDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], 
            m_function->GetPositions()[idxC], Fabc_a, Fb, Fabc_c);
        rab = VectorDistanceDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], Fa, Fb);
        rbc = VectorDistanceDerivative(m_function->GetPositions()[idxB], m_function->GetPositions()[idxC], Fb, Fc);
        if (!isfinite(theta))
          theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

        delta_theta = theta - m_calcs[i].theta0;
        delta_rab = rab - m_calcs[i].rab0;
        delta_rbc = rbc - m_calcs[i].rbc0;
        const double factor = RAD_TO_DEG * (m_calcs[i].kba123 * delta_rab + 
            m_calcs[i].kba321 * delta_rbc);

        e = DEG_TO_RAD * factor * delta_theta;
        Fa = 2.51210 * (Fa * m_calcs[i].kba123 * delta_theta + Fabc_a * factor);
        Fc = 2.51210 * (Fc * m_calcs[i].kba321 * delta_theta + Fabc_c * factor);
        Fb = - (Fa + Fc);
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
        m_function->GetGradients()[idxC] += Fc;
      } else {
        const Eigen::Vector3d ab = m_function->GetPositions()[idxA] - m_function->GetPositions()[idxB];
        const Eigen::Vector3d bc = m_function->GetPositions()[idxC] - m_function->GetPositions()[idxB];
        theta = VectorAngle(ab, bc);
        rab = ab.norm();
        rbc = bc.norm();
       
        if (!isfinite(theta))
        theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

        delta_theta = theta - m_calcs[i].theta0;
        delta_rab = rab - m_calcs[i].rab0;
        delta_rbc = rbc - m_calcs[i].rbc0;
        /*
        cout << "rab0 = " << m_calcs[i].rab0 << endl;
        cout << "rab = " << rab << endl;
        cout << "rbc0 = " << m_calcs[i].rbc0 << endl;
        cout << "rbc = " << rbc << endl;
        */

        const double factor = RAD_TO_DEG * (m_calcs[i].kba123 * delta_rab + 
            m_calcs[i].kba321 * delta_rbc);

        e = DEG_TO_RAD * factor * delta_theta;
      }

      m_value += 2.51210 * e;
/*
      char _logbuf[1000];
      snprintf(_logbuf, 1000, "%2d   %2d   %2d     %2d   %8.3f   %8.3f   %8.3f   %8.3f   %8.3f\n", 
                m_calcs[i].idx1, m_calcs[i].idx2, m_calcs[i].idx3, m_calcs[i].strbndType, 
                theta, delta_theta, 
                m_calcs[i].kba123, m_calcs[i].kba321, 
                2.51210 * e);
      m_function->GetLogFile()->Write(_logbuf);
*/
    }
	
    cout << "E{strbnd} = " << m_value << endl;
  }
  
  bool MMFF94StrBndTerm::Setup(/*const*/ OBMol &mol)
  {
    OBAtom *a, *b, *c;

    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();
 
    //
    // Angle Calculations
    //
    // MMFF part I - page 513 ("step-down" prodedure)
    // MMFF part I - page 519 (reference 68 is actually a footnote)
    // MMFF part V - page 627 (empirical rule)
    //
    // First try and find an exact match, if this fails, step down using the equivalences from mmffdef.par
    // five-stage protocol: 1-1-1, 2-2-2, 3-2-3, 4-2-4, 5-2-5
    // If this fails, use empirical rules
    // Since 1-1-1 = 2-2-2, we will only try 1-1-1 before going to 3-2-3
    // 
    // Stretch-Bend Calculations
    //
    if (logFile->IsLow())
      logFile->Write("SETTING UP STRETCH-BEND CALCULATIONS...\n");
 
    Calculation strbndcalc;
 
    m_calcs.clear();
    
    FOR_ANGLES_OF_MOL(angle, mol) {
      b = mol.GetAtom((*angle)[0] + 1);
      a = mol.GetAtom((*angle)[1] + 1);
      c = mol.GetAtom((*angle)[2] + 1);

      int type_a = m_common->GetCachedType(a);
      int type_b = m_common->GetCachedType(b);
      int type_c = m_common->GetCachedType(c);
      
      int angleType = m_common->GetAngleType(a, b, c);
      int strbndType = m_common->GetStrBndType(a, b, c);
      
      if (m_common->HasLinSet(type_b))
        continue;

      // theta from angle parameters...
      std::vector<OBParameterDB::Query> query;
      query.push_back( OBParameterDB::Query(Parameter::Class, OBVariant(angleType)) );
      query.push_back( OBParameterDB::Query(Parameter::TypeA, OBVariant(type_a), true) );
      query.push_back( OBParameterDB::Query(Parameter::TypeB, OBVariant(type_b), true) );
      query.push_back( OBParameterDB::Query(Parameter::TypeC, OBVariant(type_c), true) );
      std::vector<OBVariant> row = database->FindRow(MMFF94SimpleParameterDB::AngleParameters, query);
      if (row.empty()) {
        // try 3-2-3
        query[1] = OBParameterDB::Query(Parameter::TypeA, OBVariant(m_common->EqLvl3(type_a)), true);
        query[3] = OBParameterDB::Query(Parameter::TypeC, OBVariant(m_common->EqLvl3(type_c)), true);
        row = database->FindRow(MMFF94SimpleParameterDB::AngleParameters, query);
      }
      if (row.empty()) {
        // try 4-2-4
        query[1] = OBParameterDB::Query(Parameter::TypeA, OBVariant(m_common->EqLvl4(type_a)), true);
        query[3] = OBParameterDB::Query(Parameter::TypeC, OBVariant(m_common->EqLvl4(type_c)), true);
        row = database->FindRow(MMFF94SimpleParameterDB::AngleParameters, query);
      }
      if (row.empty()) {
        // try 5-2-5
        query[1] = OBParameterDB::Query(Parameter::TypeA, OBVariant(m_common->EqLvl5(type_a)), true);
        query[3] = OBParameterDB::Query(Parameter::TypeC, OBVariant(m_common->EqLvl5(type_c)), true);
        row = database->FindRow(MMFF94SimpleParameterDB::AngleParameters, query);
      }
 
      if (!row.empty()) {
        strbndcalc.theta0 = row.at(5).AsDouble();
      } else {
        if (logFile->IsLow()) {
          std::stringstream ss;
          ss << "   USING DEFAULT ANGLE FOR " << m_common->GetCachedType(a) << "-" 
             << m_common->GetCachedType(b) << "-" << m_common->GetCachedType(c)
             << " (MMFF94 Atom Types)..." << endl;
          ss << "   USING EMPIRICAL RULE FOR ANGLE BENDING " << m_common->GetCachedType(a) << "-" 
             << m_common->GetCachedType(b) << "-" << m_common->GetCachedType(c)
             << " (MMFF94 Atom Types)..." << endl;
          logFile->Write(ss.str());
        }

        strbndcalc.theta0 = 120.0;

        if (m_common->GetCrd(type_b) == 4)
          strbndcalc.theta0 = 109.45;
        
        if ((m_common->GetCrd(type_b) == 2) && b->IsOxygen())
          strbndcalc.theta0 = 105.0;
	
        if (b->GetAtomicNum() > 10)
          strbndcalc.theta0 = 95.0;
	
        if (m_common->HasLinSet(type_b))
          strbndcalc.theta0 = 180.0;
	
        if ((m_common->GetCrd(type_b) == 3) && (m_common->GetVal(type_b) == 3) && !m_common->GetMltb(type_b)) {
          if (b->IsNitrogen()) {
            strbndcalc.theta0 = 107.0;
          } else {
            strbndcalc.theta0 = 92.0;
          }
        }

        if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && m_common->IsInSameRing(a, c))
          strbndcalc.theta0 = 60.0;
	
        if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && m_common->IsInSameRing(a, c))
          strbndcalc.theta0 = 90.0;
      }

      bool swapped;
      query[0] = OBParameterDB::Query(Parameter::Class, OBVariant(strbndType));
      query[1] = OBParameterDB::Query(Parameter::TypeA, OBVariant(type_a), true);
      query[3] = OBParameterDB::Query(Parameter::TypeC, OBVariant(type_c), true);
      row = database->FindRow(MMFF94SimpleParameterDB::StretchBendParameters, query, &swapped);

      if (row.size()) {
        if (!swapped) {
          strbndcalc.kba123 = row.at(Parameter::Kba123).AsDouble();
          strbndcalc.kba321 = row.at(Parameter::Kba321).AsDouble();
        } else {
          strbndcalc.kba123 = row.at(Parameter::Kba321).AsDouble();
          strbndcalc.kba321 = row.at(Parameter::Kba123).AsDouble();
        }
      } else {
        // This is not a real empirical rule...
        //if (m_log->IsLow()) {
        //  snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR STRETCH-BENDING FORCE CONSTANT %d-%d-%d (IDX)...\n", a->GetIdx(), b->GetIdx(), c->GetIdx());
        //  Log(_logbuf);
        //}
    
        int row_a = m_common->GetElementRow(a);
        int row_b = m_common->GetElementRow(b);
        int row_c = m_common->GetElementRow(c);
       
        query.clear();
        query.push_back( OBParameterDB::Query(EmpiricalRule::Row1, OBVariant(row_a), true) );
        query.push_back( OBParameterDB::Query(EmpiricalRule::Row2, OBVariant(row_b), true) );
        query.push_back( OBParameterDB::Query(EmpiricalRule::Row3, OBVariant(row_c), true) );
        row = database->FindRow(MMFF94SimpleParameterDB::EmpiricalStretchBendRules, query, &swapped);
 
        if (row.empty()) {
          // This should never happen
          if (logFile->IsLow()) {
            std::stringstream ss;
            ss << "    COULD NOT FIND PARAMETERS FOR STRETCH-BEND " << m_common->GetCachedType(a) << "-" 
               << m_common->GetCachedType(b) << "-" << m_common->GetCachedType(c) << " (MMFF94 Atom Types)..." << endl;
            logFile->Write(ss.str());
          }
          return false;
        }

        if (!swapped) {
          strbndcalc.kba123 = row.at(EmpiricalRule::F123).AsDouble();
          strbndcalc.kba321 = row.at(EmpiricalRule::F321).AsDouble();
        } else {
          strbndcalc.kba123 = row.at(EmpiricalRule::F321).AsDouble();
          strbndcalc.kba321 = row.at(EmpiricalRule::F123).AsDouble();
        }
      }
      
      strbndcalc.idx1 = a->GetIdx() - 1;
      strbndcalc.idx2 = b->GetIdx() - 1;
      strbndcalc.idx3 = c->GetIdx() - 1;
      strbndcalc.rab0 = m_common->GetBondLength(a, b);
      strbndcalc.rbc0 = m_common->GetBondLength(b ,c); 
      strbndcalc.strbndType = strbndType;
      
      m_calcs.push_back(strbndcalc);
    }

    return true;
  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
