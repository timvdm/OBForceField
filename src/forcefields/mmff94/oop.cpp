/*********************************************************************
MMFF94OutOfPlaneTerm - MMFF94 force field bond stratching term

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

#include "oop.h"
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { TypeA, TypeB, TypeC, TypeD, Koop };
  };
 
  MMFF94OutOfPlaneTerm::MMFF94OutOfPlaneTerm(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
   
  void MMFF94OutOfPlaneTerm::Compute(OBFunction::Computation computation)
  {
    m_value = 0.0;
 
    unsigned int idxA, idxB, idxC, idxD;
    double angle, e;
    Eigen::Vector3d Fa, Fb, Fc, Fd;
    for (int i = 0; i < (int)m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;
      idxC = m_calcs[i].idx3;
      idxD = m_calcs[i].idx4;
 
      if (computation == OBFunction::Gradients) {
        angle = VectorOOPDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], 
            m_function->GetPositions()[idxC], m_function->GetPositions()[idxD], Fa, Fb, Fc, Fd);
      
        const double dE =  (-1.0 * RAD_TO_DEG * 0.043844 * angle * m_calcs[i].koop) / cos(angle * DEG_TO_RAD);
      
        Fa *= dE;
        Fb *= dE;
        Fc *= dE;
        Fd *= dE;
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
        m_function->GetGradients()[idxC] += Fc;
        m_function->GetGradients()[idxD] += Fd;
      } else {
        angle = VectorOOP(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], 
            m_function->GetPositions()[idxC], m_function->GetPositions()[idxD]); 
      }
     
      if (!isfinite(angle))
        angle = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;
    
      const double angle2 = angle * angle;
      e = m_calcs[i].koop * angle2;
      m_value += e;
     
     /* 
        snprintf(_logbuf, BUFF_SIZE, "%2d   %2d   %2d   %2d      0   %8.3f   %8.3f     %8.3f\n", 
                atoi(m_calcs[i].a->GetType()), atoi(m_calcs[i].b->GetType()), 
                atoi(m_calcs[i].c->GetType()), atoi(m_calcs[i].d->GetType()), 
                angle, m_calcs[i].koop, 
                0.043844 * 0.5 * e);
        m_log->Log(_logbuf);
        */
    }
 
    /*
    if (m_log->IsMedium()) {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL OUT-OF-PLANE BENDING ENERGY = %8.5f %s\n", 
          0.043844 * 0.5 * energy, GetUnit().c_str());
      m_log->Log(_logbuf);
    }
    */

    m_value *= 0.043844 * 0.5;

  }
  
  bool MMFF94OutOfPlaneTerm::Setup(/*const*/ OBMol &mol)
  {
    OBAtom *a, *b, *c, *d;

    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();

    //
    // Out-Of-Plane Calculations
    //
    if (logFile->IsLow())
      logFile->Write("SETTING UP OOP CALCULATIONS...\n");
 
    Calculation oopcalc;
    m_calcs.clear();
 
    FOR_ATOMS_OF_MOL(atom, mol) {
      b = (OBAtom*) &*atom;

      bool found = false;
      int type_b = m_common->GetCachedType(b);
     
      const std::vector<std::vector<OBVariant> > &rows = database->GetAllRows(MMFF94SimpleParameterDB::OutOfPlaneParameters);
      for (unsigned int idx = 0; idx < rows.size(); idx++) {
        const std::vector<OBVariant> &row = rows.at(idx);
        if (type_b == row.at(Parameter::TypeB).AsInt()) {
          a = c = d = 0;

          FOR_NBORS_OF_ATOM(nbr, b) {
            if (!a)
              a = (OBAtom*) &*nbr;
            else if (!c)
              c = (OBAtom*) &*nbr;
            else
              d = (OBAtom*) &*nbr;
          }
	  
          if (!a || !c || !d)
            break;
          
          int type_a = m_common->GetCachedType(a);
          int type_c = m_common->GetCachedType(c);
          int type_d = m_common->GetCachedType(d);

          int row_a = row.at(Parameter::TypeA).AsInt();
          int row_c = row.at(Parameter::TypeC).AsInt();
          int row_d = row.at(Parameter::TypeD).AsInt();
 
          if (((type_a == row_a) && (type_c == row_c) && (type_d == row_d)) ||
              ((type_c == row_a) && (type_a == row_c) && (type_d == row_d)) ||
              ((type_c == row_a) && (type_d == row_c) && (type_a == row_d)) ||
              ((type_d == row_a) && (type_c == row_c) && (type_a == row_d)) ||
              ((type_a == row_a) && (type_d == row_c) && (type_c == row_d)) ||
              ((type_d == row_a) && (type_a == row_c) && (type_c == row_d)))
          {
            found = true;

            oopcalc.koop = row.at(Parameter::Koop).AsDouble();

            // A-B-CD || C-B-AD  PLANE = ABC
            oopcalc.idx1 = a->GetIdx() - 1;
            oopcalc.idx2 = b->GetIdx() - 1;
            oopcalc.idx3 = c->GetIdx() - 1;
            oopcalc.idx4 = d->GetIdx() - 1;

            m_calcs.push_back(oopcalc);

            // C-B-DA || D-B-CA  PLANE BCD
            oopcalc.idx1 = d->GetIdx() - 1;
            oopcalc.idx4 = a->GetIdx() - 1;

            m_calcs.push_back(oopcalc);

            // A-B-DC || D-B-AC  PLANE ABD
            oopcalc.idx1 = a->GetIdx() - 1;
            oopcalc.idx3 = d->GetIdx() - 1;
            oopcalc.idx4 = c->GetIdx() - 1;

            m_calcs.push_back(oopcalc);
          }

            // *-XX-*-*
          if ((row_a == 0) && (row_c == 0) && (row_d == 0) && !found) {
            oopcalc.koop = row.at(Parameter::Koop).AsDouble();

            // A-B-CD || C-B-AD  PLANE = ABC
            oopcalc.idx1 = a->GetIdx() - 1;
            oopcalc.idx2 = b->GetIdx() - 1;
            oopcalc.idx3 = c->GetIdx() - 1;
            oopcalc.idx4 = d->GetIdx() - 1;

            m_calcs.push_back(oopcalc);

            // C-B-DA || D-B-CA  PLANE BCD
            oopcalc.idx1 = d->GetIdx() - 1;
            oopcalc.idx4 = a->GetIdx() - 1;

            m_calcs.push_back(oopcalc);

            // A-B-DC || D-B-AC  PLANE ABD
            oopcalc.idx1 = a->GetIdx() - 1;
            oopcalc.idx3 = d->GetIdx() - 1;
            oopcalc.idx4 = c->GetIdx() - 1;

            m_calcs.push_back(oopcalc);
          }
        }
      }
    }

    return true;
  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
