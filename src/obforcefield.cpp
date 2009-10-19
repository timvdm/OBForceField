/**********************************************************************
forcefield.cpp - Handle OBForceField class.

Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
Some portions Copyright (C) 2007-2008 by Geoffrey Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

      development status:
      - src/forcefield.cpp
      - LineSearch(): finished
      - SteepestDescent(): finished
      - ConjugateGradients(): finished
      - SystematicRotorSearch(): finished
      - RandomRotorSearch(): finished
      - WeightedRotorSearch: finished
      - DistanceGeometry(): needs matrix operations (Eigen)

      Constraints:
      - Fix Atom: working 
      - Fix Atom X: working 
      - Fix Atom Y: working 
      - Fix Atom Z: working 
      - Distance: working
      - Angle: working
      - Torsion: working 
      - Chirality: TODO

      src/forcefields/forcefieldghemical.cpp
      - Atom typing: finished
      - Charges: finished
      - Energy terms: finished
      - Analytical gradients: finished
      - Validation: finished

      src/forcefields/forcefieldmmff94.cpp
      - Atom typing: done.
      - Charges: done.
      - Energy terms: finished (small problems with SSSR 
        algorithm not finding all bridged rings)
      - Analytical gradients: finished 
      - Validation: http://home.scarlet.be/timvdm/MMFF94_validation_output.gz

      src/forcefields/forcefielduff.cpp
      - Energy terms: finished
      - OOP: needs validation
      - Gradients: need OOP gradient
      - Validation in progress...


      forcefield.cpp layout:

        1) Constraints                XXXX01
        2) Methods for logging        XXXX02
        3) General                    XXXX03
        4) Structure generation       XXXX04
        5) Cut-off                    XXXX05
        6) Interaction groups         XXXX06
        8) Molecular Dynamics         XXXX08
        9) Misc.                      XXXX09

***********************************************************************/


#include <OBForceField>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/rotamer.h>
#include <openbabel/rotor.h>

#include <OBCutOffGrid>
#include <OBMinimize>
#include <OBVectorMath>


#include <openbabel/babelconfig.h>

using namespace std;

namespace OpenBabel {
namespace OBFFs {

  class OBForceFieldPrivate 
  {
    public:
    // general variables
    bool 	validSetup; //!< was the last call to Setup succesfull
  };

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Methods for logging  (XXXX02)
  //
  //////////////////////////////////////////////////////////////////////////////////
  
//   void OBForceField::PrintAtomTypes()
//   {
//     OBLogFile *logFile = GetLogFile();
//     if (logFile->IsLow()) {
//       logFile->Write("\nA T O M   T Y P E S\n\n");
//       logFile->Write("IDX\tTYPE\n");
      
//       std::stringstream ss;
//       std::vector<std::string> types = GetAtomTypes();
//       for (unsigned int i = 0; i < types.size(); ++i) {
//         ss << i << "\t" << types.at(i) << std::endl;
//         logFile->Write(ss.str());
//         ss.str("");
//       }
//     }
//   }
    
//   void OBForceField::PrintFormalCharges()
//   {
//     OBLogFile *logFile = GetLogFile();
//     if (logFile->IsLow()) {
//       logFile->Write("\nF O R M A L   C H A R G E S\n\n");
//       logFile->Write("IDX\tCHARGE\n");
 
//       std::stringstream ss;
//       std::vector<double> charges = GetFormalCharges();
//       for (unsigned int i = 0; i < charges.size(); ++i) {
//         ss << i << "\t" << charges.at(i) << std::endl;
//         logFile->Write(ss.str());
//         ss.str("");
//       }
//     }
//   }

//   void OBForceField::PrintPartialCharges()
//   {
//     OBLogFile *logFile = GetLogFile();
//     if (logFile->IsLow()) {
//       logFile->Write("\nP A R T I A L   C H A R G E S\n\n");
//       logFile->Write("IDX\tCHARGE\n");
 
//       std::stringstream ss;
//       std::vector<double> charges = GetPartialCharges();
//       for (unsigned int i = 0; i < charges.size(); ++i) {
//         ss << i << "\t" << charges.at(i) << std::endl;
//         logFile->Write(ss.str());
//         ss.str("");
//       }
//     }
//   }

  void OBForceField::ProcessOptions(std::vector<Option> &options)
  {
  
  }


} // end namespace OpenBabel
}

//! \file forcefield.cpp
//! \brief Handle OBForceField class
