/*********************************************************************
MMFF94TorsionTerm - MMFF94 force field bond stratching term

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
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { TorsionType = 0, TypeA, TypeB, TypeC, TypeD, V1, V2, V3 };
  };
 
  MMFF94TorsionTerm::MMFF94TorsionTerm(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
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
  void MMFF94TorsionTerm::Compute(OBFunction::Computation computation)
  {
    m_value = 0.0;
  
    unsigned int idxA, idxB, idxC, idxD;
    double tor, e;
    Eigen::Vector3d Fa, Fb, Fc, Fd;
    for (int i = 0; i < (int)m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;
      idxC = m_calcs[i].idx3;
      idxD = m_calcs[i].idx4;
      
      if (computation == OBFunction::Gradients) {
        tor = VectorTorsionDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], 
            m_function->GetPositions()[idxC], m_function->GetPositions()[idxD], Fa, Fb, Fc, Fd);
        if (!isfinite(tor))
          tor = 1.0e-3;
      
        const double sine1 = sin(DEG_TO_RAD * tor);
        const double sine2 = sin(2.0 * DEG_TO_RAD * tor);
        const double sine3 = sin(3.0 * DEG_TO_RAD * tor);
        const double dE = 0.5 * (m_calcs[i].v1 * sine1 - 
                          2.0 *  m_calcs[i].v2 * sine2 + 
                          3.0 *  m_calcs[i].v3 * sine3); // MMFF
        Fa *= dE; 
        Fb *= dE; 
        Fc *= dE; 
        Fd *= dE; 
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
        m_function->GetGradients()[idxC] += Fc;
        m_function->GetGradients()[idxD] += Fd;
      } else {
        tor = VectorTorsion(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], 
            m_function->GetPositions()[idxC], m_function->GetPositions()[idxD]);
        if (!isfinite(tor))
          tor = 1.0e-3;
      }
    
      const double cosine1 = cos(DEG_TO_RAD * tor);
      const double cosine2 = cos(DEG_TO_RAD * 2 * tor);
      const double cosine3 = cos(DEG_TO_RAD * 3 * tor);
      
      const double phi1 = 1.0 + cosine1;
      const double phi2 = 1.0 - cosine2;
      const double phi3 = 1.0 + cosine3;
    
      e = (m_calcs[i].v1 * phi1 + m_calcs[i].v2 * phi2 + m_calcs[i].v3 * phi3);
      m_value += e;
    }
 
    m_value *= 0.5;
  }
  
  bool MMFF94TorsionTerm::Setup(/*const*/ OBMol &mol)
  {
    OBAtom *a, *b, *c, *d;

    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();

    //
    // Torsion Calculations
    //
    // MMFF part I - page 513 ("step-down" prodedure)
    // MMFF part I - page 519 (reference 68 is actually a footnote)
    // MMFF part IV - page 631 (empirical rule)
    //
    // First try and find an exact match, if this fails, step down using the equivalences from mmffdef.par
    // five-stage protocol: 1-1-1-1, 2-2-2-2, 3-2-2-5, 5-2-2-3, 5-2-2-5
    // If this fails, use empirical rules
    // Since 1-1-1-1 = 2-2-2-2, we will only try 1-1-1-1 before going to 3-2-2-5
    //
    if (logFile->IsLow())
      logFile->Write("SETTING UP TORSION CALCULATIONS...\n");
 
    Calculation torsioncalc;
    m_calcs.clear();
 
    FOR_TORSIONS_OF_MOL(t, mol) {
      a = mol.GetAtom((*t)[0] + 1);
      b = mol.GetAtom((*t)[1] + 1);
      c = mol.GetAtom((*t)[2] + 1);
      d = mol.GetAtom((*t)[3] + 1);
      
      int type_a = m_common->GetCachedType(a);
      int type_b = m_common->GetCachedType(b);
      int type_c = m_common->GetCachedType(c);
      int type_d = m_common->GetCachedType(d);
      
      int torsionType = m_common->GetTorsionType(a, b, c, d);
      // CXT = MC*(J*MA**3 + K*MA**2 + I*MA + L) + TTijkl  MC = 6, MA = 136
      int order = (type_c*2515456 + type_b*18496 + type_d*136 + type_a) 
        - (type_b*2515456 + type_c*18496 + type_a*136 + type_d);

      std::vector<OBParameterDB::Query> query;
      query.push_back( OBParameterDB::Query(0, OBVariant(torsionType)) );
      query.push_back( OBParameterDB::Query(1, OBVariant(type_a), true) );
      query.push_back( OBParameterDB::Query(2, OBVariant(type_b), true) );
      query.push_back( OBParameterDB::Query(3, OBVariant(type_c), true) );
      query.push_back( OBParameterDB::Query(4, OBVariant(type_d), true) );
 
      /*
      query.push_back( OBParameterDB::Query(Parameter::TorsionType, OBVariant(torsionType)) );
      query.push_back( OBParameterDB::Query(Parameter::TypeA, OBVariant(type_a), true) );
      query.push_back( OBParameterDB::Query(Parameter::TypeB, OBVariant(type_b), true) );
      query.push_back( OBParameterDB::Query(Parameter::TypeC, OBVariant(type_c), true) );
      query.push_back( OBParameterDB::Query(Parameter::TypeD, OBVariant(type_d), true) );
      */
      std::vector<OBVariant> row = database->FindRow(MMFF94SimpleParameterDB::TorsionParameters, query);
      if (row.empty()) {
        // try 3-2-2-5
        query[1] = OBParameterDB::Query(Parameter::TypeA, OBVariant(m_common->EqLvl3(type_a)), true);
        query[4] = OBParameterDB::Query(Parameter::TypeD, OBVariant(m_common->EqLvl5(type_d)), true);
        row = database->FindRow(MMFF94SimpleParameterDB::TorsionParameters, query);
      }
      if (row.empty()) {
        // try 5-2-2-3
        query[1] = OBParameterDB::Query(Parameter::TypeA, OBVariant(m_common->EqLvl5(type_a)), true);
        query[4] = OBParameterDB::Query(Parameter::TypeD, OBVariant(m_common->EqLvl3(type_d)), true);
        row = database->FindRow(MMFF94SimpleParameterDB::TorsionParameters, query);
      }
      if (row.empty()) {
        // try 5-2-2-5
        query[1] = OBParameterDB::Query(Parameter::TypeA, OBVariant(m_common->EqLvl5(type_a)), true);
        query[4] = OBParameterDB::Query(Parameter::TypeD, OBVariant(m_common->EqLvl5(type_d)), true);
        row = database->FindRow(MMFF94SimpleParameterDB::TorsionParameters, query);
      }
 
      if (!row.empty()) {
        torsioncalc.v1 = row.at(Parameter::V1).AsDouble();
        torsioncalc.v2 = row.at(Parameter::V2).AsDouble();
        torsioncalc.v3 = row.at(Parameter::V3).AsDouble();
      } else {
        bool found_rule = false;
	
        //if (m_log->IsLow()) {
        //  snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR TORSION FORCE CONSTANT %d-%d-%d-%d (IDX)...\n", 
        //    a->GetIdx(), b->GetIdx(), c->GetIdx(), d->GetIdx());
        //  Log(_logbuf);
        //}

        // rule (a) page 631
        if (m_common->HasLinSet(type_b) || m_common->HasLinSet(type_c))
          continue;
        
        // rule (b) page 631
        if (b->GetBond(c)->IsAromatic()) {
          double Ub, Uc, pi_bc, beta;
          Ub = m_common->GetUParam(b);
          Uc = m_common->GetUParam(c);

          if (!m_common->HasPilpSet(type_b) && !m_common->HasPilpSet(type_c))
            pi_bc = 0.5;
          else
            pi_bc = 0.3;

          if (((m_common->GetVal(type_b) == 3) && (m_common->GetVal(type_c) == 4)) || 
              ((m_common->GetVal(type_b) == 4) && (m_common->GetVal(type_c) == 3)))
            beta = 3.0;
          else
            beta = 6.0;
	  
          torsioncalc.v1 = 0.0;
          torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
          torsioncalc.v3 = 0.0;
          found_rule = true;
        } else {
          // rule (c) page 631	
       	  double Ub, Uc, pi_bc, beta;
          Ub = m_common->GetUParam(b);
          Uc = m_common->GetUParam(c);

          if (((m_common->GetMltb(type_b) == 2) && (m_common->GetMltb(type_c) == 2)) && a->GetBond(b)->IsDouble())
            pi_bc = 1.0;
          else
            pi_bc = 0.4;

          beta = 6.0;
          torsioncalc.v1 = 0.0;
          torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
          torsioncalc.v3 = 0.0;
          found_rule = true;
        }

        // rule (d) page 632
        if (!found_rule)
          if (((m_common->GetCrd(type_b) == 4) && (m_common->GetCrd(type_c) == 4))) {
            double Vb, Vc;
            Vb = m_common->GetVParam(b);
            Vc = m_common->GetVParam(c);

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = 0.0;
            torsioncalc.v3 = sqrt(Vb * Vc) / 9.0;
            found_rule = true;
          }
	
        // rule (e) page 632
        if (!found_rule)
          if (((m_common->GetCrd(type_b) == 4) && (m_common->GetCrd(type_c) != 4))) {
            if (m_common->GetCrd(type_c) == 3) // case (1)
              if ((m_common->GetVal(type_c) == 4) || (m_common->GetVal(type_c) == 34) || (m_common->GetMltb(type_c) != 0))
                continue;
	    
            if (m_common->GetCrd(type_c) == 2) // case (2)
              if ((m_common->GetVal(type_c) == 3) || (m_common->GetMltb(type_c) != 0))
                continue;
	    
            // case (3) saturated bonds -- see rule (h)
          }
	
        // rule (f) page 632
        if (!found_rule)
          if (((m_common->GetCrd(type_b) != 4) && (m_common->GetCrd(type_c) == 4))) {
            if (m_common->GetCrd(type_b) == 3) // case (1)
              if ((m_common->GetVal(type_b) == 4) || (m_common->GetVal(type_b) == 34) || (m_common->GetMltb(type_b) != 0))
                continue;
	    
            if (m_common->GetCrd(type_b) == 2) // case (2)
              if ((m_common->GetVal(type_b) == 3) || (m_common->GetMltb(type_b) != 0))
                continue;
	    
            // case (3) saturated bonds
          }
	
        // rule (g) page 632
        if (!found_rule)
          if (b->GetBond(c)->IsSingle() && (
                                            (m_common->GetMltb(type_b) && m_common->GetMltb(type_c)) ||
                                            (m_common->GetMltb(type_b) && m_common->HasPilpSet(type_c)) ||
                                            (m_common->GetMltb(type_c) && m_common->HasPilpSet(type_b))  )) {
            if (m_common->HasPilpSet(type_b) && m_common->HasPilpSet(type_c)) // case (1)
              continue;
	    
            double Ub, Uc, pi_bc, beta;
            Ub = m_common->GetUParam(b);
            Uc = m_common->GetUParam(c);
            beta = 6.0;
  
            if (m_common->HasPilpSet(type_b) && m_common->GetMltb(type_c)) { // case (2)
              if (m_common->GetMltb(type_c) == 1)
                pi_bc = 0.5;
              else if ((m_common->GetElementRow(b) == 1) && (m_common->GetElementRow(c) == 1))
                pi_bc = 0.3;
              else
                pi_bc = 0.15;
              found_rule = true;
            }
	    
            if (m_common->HasPilpSet(type_c) && m_common->GetMltb(type_b)) { // case (3)
              if (m_common->GetMltb(type_b) == 1)
                pi_bc = 0.5;
              else if ((m_common->GetElementRow(b) == 1) && (m_common->GetElementRow(c) == 1))
                pi_bc = 0.3;
              else
                pi_bc = 0.15;
              found_rule = true;
            }
            
            if (!found_rule)
              if (((m_common->GetMltb(type_b) == 1) || (m_common->GetMltb(type_c) == 1)) && (!b->IsCarbon() || !c->IsCarbon())) {
                pi_bc = 0.4;
                found_rule = true;
              }
	    
            if (!found_rule)
              pi_bc = 0.15;
	    
            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = beta * pi_bc * sqrt(Ub * Uc);
            torsioncalc.v3 = 0.0;
            found_rule = true;
          }
	
        // rule (h) page 632
        if (!found_rule)
          if ((b->IsOxygen() || b->IsSulfur()) && (c->IsOxygen() || c->IsSulfur())) {
            double Wb, Wc;

            if (b->IsOxygen()) {
              Wb = 2.0;
            }
            else {
              Wb = 8.0;
            }
	    
            if (c->IsOxygen()) {
              Wc = 2.0;
            }
            else {
              Wc = 8.0;
            }

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = -sqrt(Wb * Wc);
            torsioncalc.v3 = 0.0;
          } else {
            double Vb, Vc, Nbc;
            Vb = m_common->GetVParam(b);
            Vc = m_common->GetVParam(c);

            if (logFile->IsLow()) {
              std::stringstream ss;
              ss << "   USING EMPIRICAL RULE FOR TORSION FORCE CONSTANT " 
                 << m_common->GetCachedType(a) << "-" << m_common->GetCachedType(b) << "-"
                 << m_common->GetCachedType(c) << "-" << m_common->GetCachedType(d)
                 << " (NNFF94 Atom Types)..." << endl;
              logFile->Write(ss.str());
            }
            
            Nbc = m_common->GetCrd(type_b) * m_common->GetCrd(type_c);

            torsioncalc.v1 = 0.0;
            torsioncalc.v2 = 0.0;
            torsioncalc.v3 = sqrt(Vb * Vc) / Nbc;
          }
      }
      
      torsioncalc.idx1 = a->GetIdx() - 1;
      torsioncalc.idx2 = b->GetIdx() - 1;
      torsioncalc.idx3 = c->GetIdx() - 1;
      torsioncalc.idx4 = d->GetIdx() - 1;
      torsioncalc.torsionType = torsionType;

      m_calcs.push_back(torsioncalc);
    }

  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
