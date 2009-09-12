/*********************************************************************
MMFF94VDWTerm - MMFF94 force field bond stratching term

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

#include "vdw.h"
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { Type, Alpha_i, N_i, A_i, G_i, DA };
  };
 
  MMFF94VDWTerm::MMFF94VDWTerm(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
   
  void MMFF94VDWTerm::Compute(OBFunction::Computation computation)
  {
    m_value = 0.0;
   
    unsigned int idxA, idxB;
    double rab, e;
    Eigen::Vector3d Fa, Fb;
    for (int i = 0; i < (int)m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;
  
      if (computation == OBFunction::Gradients) {
        rab = VectorDistanceDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], Fa, Fb);
      } else {
        const Eigen::Vector3d ab = m_function->GetPositions()[idxA] - m_function->GetPositions()[idxB];
        rab = ab.norm();
      }
    
      const double rab7 = rab*rab*rab*rab*rab*rab*rab;
      double erep = (1.07 * m_calcs[i].R_AB) / (rab + 0.07 * m_calcs[i].R_AB); //***
      double erep7 = erep*erep*erep*erep*erep*erep*erep;
      double eattr = (((1.12 * m_calcs[i].R_AB7) / (rab7 + 0.12 * m_calcs[i].R_AB7)) - 2.0);
      
      e = m_calcs[i].epsilon * erep7 * eattr;
    
      if (computation == OBFunction::Gradients) { 
        const double q = rab / m_calcs[i].R_AB;
        const double q6 = q*q*q*q*q*q;
        const double q7 = q6 * q;
        erep = 1.07 / (q + 0.07); 
        erep7 = erep*erep*erep*erep*erep*erep*erep;
        const double term = q7 + 0.12;
        const double term2 = term * term;
        eattr = (-7.84 * q6) / term2 + ((-7.84 / term) + 14) / (q + 0.07);
        const double dE = (m_calcs[i].epsilon / m_calcs[i].R_AB) * erep7 * eattr;
        Fa *= dE;
        Fb *= dE;
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
      }

      m_value += e;
    }
   
    OBLogFile *logFile = m_function->GetLogFile();
    if (logFile->IsMedium()) {
      std::stringstream ss;
      ss << "     TOTAL VAN DER WAALS ENERGY = " << m_value << " " << m_function->GetUnit() << std::endl;
      logFile->Write(ss.str());
    }
    
  }
  
  bool MMFF94VDWTerm::Setup(/*const*/ OBMol &mol)
  {
    OBAtom *a, *b, *c;

    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();

    // 
    // VDW Calculations
    //  
    if (logFile->IsLow())
      logFile->Write("SETTING UP VAN DER WAALS CALCULATIONS...\n");
 
    Calculation vdwcalc;
    double alpha_a, Na, Aa, Ga;
    double alpha_b, Nb, Ab, Gb;
    double R_AA, R_BB, R_AB6, g_AB, g_AB2;
    double R_AB2, R_AB4, sqrt_a, sqrt_b;
    int aDA, bDA;
 
    m_calcs.clear();
 
    FOR_PAIRS_OF_MOL(p, mol) {
      a = mol.GetAtom((*p)[0]);
      b = mol.GetAtom((*p)[1]);
      
      std::vector<OBParameterDB::Query> query;
      query.push_back( OBParameterDB::Query(Parameter::Type, OBVariant(m_common->GetCachedType(a))) );
      std::vector<OBVariant> row_a = database->FindRow(MMFF94SimpleParameterDB::VanDerWaalsParameters, query);
      query[0] = OBParameterDB::Query(Parameter::Type, OBVariant(m_common->GetCachedType(b)));
      std::vector<OBVariant> row_b = database->FindRow(MMFF94SimpleParameterDB::VanDerWaalsParameters, query);
 
      if ((row_a.empty()) || (row_b.empty())) {
        if (logFile->IsLow()) {
          std::stringstream ss;
          ss << "   COULD NOT FIND VAN DER WAALS PARAMETERS FOR " << m_common->GetCachedType(a) << "-"
             << m_common->GetCachedType(b) << " (MMFF94 Atom Types)..." << endl;
          logFile->Write(ss.str());
        }

        return false;
      }
      
      vdwcalc.idx1 = a->GetIdx() - 1;
      vdwcalc.idx2 = b->GetIdx() - 1;
      /*
      vdwcalc.alpha_a = parameter_a->_dpar[0];
      vdwcalc.Na = parameter_a->_dpar[1];
      vdwcalc.Aa = parameter_a->_dpar[2];
      vdwcalc.Ga = parameter_a->_dpar[3];
      vdwcalc.aDA = parameter_a->_ipar[0];
      
      vdwcalc.alpha_b = parameter_b->_dpar[0];
      vdwcalc.Nb = parameter_b->_dpar[1];
      vdwcalc.Ab = parameter_b->_dpar[2];
      vdwcalc.Gb = parameter_b->_dpar[3];
      vdwcalc.bDA = parameter_b->_ipar[0];
      */
      alpha_a = row_a.at(Parameter::Alpha_i).AsDouble();
      Na = row_a.at(Parameter::N_i).AsDouble();
      Aa = row_a.at(Parameter::A_i).AsDouble();
      Ga = row_a.at(Parameter::G_i).AsDouble();
      aDA = row_a.at(Parameter::DA).AsInt();
 
      alpha_b = row_b.at(Parameter::Alpha_i).AsDouble();
      Nb = row_b.at(Parameter::N_i).AsDouble();
      Ab = row_b.at(Parameter::A_i).AsDouble();
      Gb = row_b.at(Parameter::G_i).AsDouble();
      bDA = row_b.at(Parameter::DA).AsInt();
      
      //these calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      R_AA = Aa * pow(alpha_a, 0.25);
      R_BB = Ab * pow(alpha_b, 0.25);
      sqrt_a = sqrt(alpha_a / Na);
      sqrt_b = sqrt(alpha_b / Nb);
      
      if (aDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5 * (R_AA + R_BB);
        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        
        if (bDA == 2) { // hydrogen bond acceptor
          vdwcalc.epsilon = 0.5 * (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
          // R_AB is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled. 
          vdwcalc.R_AB *= 0.8;
        } else
          vdwcalc.epsilon = (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
	
        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;
      } else if (bDA == 1) { // hydrogen bond donor
        vdwcalc.R_AB = 0.5 * (R_AA + R_BB);
       	R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        
        if (aDA == 2) { // hydrogen bond acceptor
          vdwcalc.epsilon = 0.5 * (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
          // R_AB is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled. 
          vdwcalc.R_AB *= 0.8;
        } else
          vdwcalc.epsilon = (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
	
        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;
      } else {
        g_AB = (R_AA - R_BB) / ( R_AA + R_BB);
        g_AB2 = g_AB * g_AB;
        vdwcalc.R_AB =  0.5 * (R_AA + R_BB) * (1.0 + 0.2 * (1.0 - exp(-12.0 * g_AB2)));
        R_AB2 = vdwcalc.R_AB * vdwcalc.R_AB;
        R_AB4 = R_AB2 * R_AB2;
        R_AB6 = R_AB4 * R_AB2;
        vdwcalc.R_AB7 = R_AB6 * vdwcalc.R_AB;
        vdwcalc.epsilon = (181.16 * Ga * Gb * alpha_a * alpha_b) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
      }
      
      m_calcs.push_back(vdwcalc);
    }

  } 

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
