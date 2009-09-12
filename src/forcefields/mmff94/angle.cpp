/*********************************************************************
MMFF94AngleTerm - MMFF94 force field bond stratching term

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
 
  MMFF94AngleTerm::MMFF94AngleTerm(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
   
  // 
  // MMFF part I - page 495
  //      
  //                       ka_ijk                       
  // EA_ijk = 0.438449325 -------- /\0_ijk^2 (1 + cs /\0_ijk)
  //                         2                          
  //
  // ka_ijk	force constant (md A/rad^2)
  //
  // /\0_ijk 	0_ijk - 00_ijk (degrees)
  //
  // cs		cubic bend constant = -0.007 deg^-1 = -0.4 rad^-1
  //    
  void MMFF94AngleTerm::Compute(OBFunction::Computation computation)
  {
    m_value = 0.0;
   
    unsigned int idxA, idxB, idxC;
    double theta, e, delta, dE;
    Eigen::Vector3d Fa, Fb, Fc;
    for (unsigned int i = 0; i < m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;
      idxC = m_calcs[i].idx3;
   
      if (computation == OBFunction::Gradients) {
        theta = VectorAngleDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], 
            m_function->GetPositions()[idxC], Fa, Fb, Fc);
      
        if (!isfinite(theta))
          theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;
      
        delta = theta - m_calcs[i].theta0;
       
        if (m_calcs[i].linear) {
          e  = 143.9325 * m_calcs[i].ka * (1.0 + cos((theta) * DEG_TO_RAD));
          dE = -sin((theta) * DEG_TO_RAD) * 143.9325 * m_calcs[i].ka;
        } else {
          const double delta2 = delta * delta;
          e  = 0.043844 * 0.5 * m_calcs[i].ka * delta2 * (1.0 - 0.007 * delta);
          dE = RAD_TO_DEG * 0.043844 * m_calcs[i].ka * delta * (1.0 - 1.5 * 0.007 * delta);
        }
      
        Fa *= dE;
        Fb *= dE;
        Fc *= dE;
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
        m_function->GetGradients()[idxC] += Fc;
      } else {
        const Eigen::Vector3d ab = m_function->GetPositions()[idxA] - m_function->GetPositions()[idxB];
        const Eigen::Vector3d bc = m_function->GetPositions()[idxC] - m_function->GetPositions()[idxB];
        theta = VectorAngle(ab, bc);
      
        if (!isfinite(theta))
          theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;
      
        delta = theta - m_calcs[i].theta0;
      
        if (m_calcs[i].linear) {
          e = 143.9325 * m_calcs[i].ka * (1.0 + cos(theta * DEG_TO_RAD));
        } else {
          const double delta2 = delta * delta;
          e = 0.043844 * 0.5 * m_calcs[i].ka * delta2 * (1.0 - 0.007 * delta);
        }
      }

      m_value += e;
    }
 
  }
  
  bool MMFF94AngleTerm::Setup(/*const*/ OBMol &mol)
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
      logFile->Write("SETTING UP ANGLE BENDING CALCULATIONS...\n");
 
    Calculation anglecalc;
 
    m_calcs.clear();
    
    FOR_ANGLES_OF_MOL(angle, mol) {
      b = mol.GetAtom((*angle)[0] + 1);
      a = mol.GetAtom((*angle)[1] + 1);
      c = mol.GetAtom((*angle)[2] + 1);

      int type_a = m_common->GetCachedType(a);
      int type_b = m_common->GetCachedType(b);
      int type_c = m_common->GetCachedType(c);
      
      int angleType = m_common->GetAngleType(a, b, c);
      
      if (m_common->HasLinSet(type_b)) {
        anglecalc.linear = true;
      } else {
        anglecalc.linear = false;
      }

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
        anglecalc.ka = row.at(Parameter::Ka).AsDouble();
        anglecalc.theta0 = row.at(Parameter::Theta0).AsDouble();
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

        anglecalc.ka = 0.0;
        anglecalc.theta0 = 120.0;

        if (m_common->GetCrd(type_b) == 4)
          anglecalc.theta0 = 109.45;
        
        if ((m_common->GetCrd(type_b) == 2) && b->IsOxygen())
          anglecalc.theta0 = 105.0;
	
        if (b->GetAtomicNum() > 10)
          anglecalc.theta0 = 95.0;
	
        if (m_common->HasLinSet(type_b))
          anglecalc.theta0 = 180.0;
	
        if ((m_common->GetCrd(type_b) == 3) && (m_common->GetVal(type_b) == 3) && !m_common->GetMltb(type_b)) {
          if (b->IsNitrogen()) {
            anglecalc.theta0 = 107.0;
          } else {
            anglecalc.theta0 = 92.0;
          }
        }

        if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && m_common->IsInSameRing(a, c))
          anglecalc.theta0 = 60.0;
	
        if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && m_common->IsInSameRing(a, c))
          anglecalc.theta0 = 90.0;
      }
      
      // empirical rule for 0-b-0 and standard angles
      if (anglecalc.ka == 0.0) {
        if (logFile->IsLow()) {
          std::stringstream ss;
          ss << "   USING EMPIRICAL RULE FOR ANGLE BENDING FORCE CONSTANT " << m_common->GetCachedType(a) << "-" 
             << m_common->GetCachedType(b) << "-" << m_common->GetCachedType(c)
             << " (MMFF94 Atom Types)..." << endl;
          logFile->Write(ss.str());
        }

        double beta, Za, Zc, Cb, r0ab, r0bc, theta, theta2, D, rr, rr2;
        Za = m_common->GetZParam(a);
        Cb = m_common->GetCParam(b); // Fixed typo -- PR#2741658
        Zc = m_common->GetZParam(c);
	
        r0ab = m_common->GetBondLength(a, b);
        r0bc = m_common->GetBondLength(b, c);
        rr = r0ab + r0bc;
        rr2 = rr * rr;
        D = (r0ab - r0bc) / rr2;

        theta = anglecalc.theta0;
        theta2 = theta * theta;

        beta = 1.75;
        if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && m_common->IsInSameRing(a, c))
          beta = 0.85 * beta;
        if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && m_common->IsInSameRing(a, c))
          beta = 0.05 * beta;
        
        anglecalc.ka = (beta * Za * Cb * Zc * exp(-2 * D)) / (rr * theta2 * DEG_TO_RAD * DEG_TO_RAD); // PR#2741669
      }
      
      anglecalc.idx1 = a->GetIdx() - 1;
      anglecalc.idx2 = b->GetIdx() - 1;
      anglecalc.idx3 = c->GetIdx() - 1;
      anglecalc.angleType = angleType;
      
      m_calcs.push_back(anglecalc);
    }

  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
