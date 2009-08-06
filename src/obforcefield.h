/**********************************************************************
forcefield.h - Handle OBForceField class.
 
Copyright (C) 2006-2007 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#ifndef OB_FORCEFIELD_H
#define OB_FORCEFIELD_H

#include <vector>
#include <string>

#include <OBFunction>       // full header to derive classes

namespace OpenBabel {

  class OBMol;
  class OBAtom;

namespace OBFFs {
  
  //class OBForceFieldPrivate;

  class OBForceField : public OBFunction
  {
    protected:
      // other data stored using d pointer
//      OBForceFieldPrivate * const d; //!< the d-pointer

  public:
    /** 
     * @brief Destructor.
     */
//    virtual ~OBForceField();
    /** 
     * @brief Get the force field atom types. 
     */
    virtual std::vector<std::string> GetAtomTypes() const = 0;
    /** 
     * @brief Get the force field formal charges.
     */
    virtual std::vector<double> GetFormalCharges() const = 0;
    /** 
     * @brief Get the force field partial charges.
     */
    virtual std::vector<double> GetPartialCharges() const = 0;
 
    /////////////////////////////////////////////////////////////////////////
    // Logging                                                             //
    /////////////////////////////////////////////////////////////////////////
      
    //! \name Methods for logging
    //@{
    /** 
     * @brief Print the atom types to the log.
     */
    void PrintAtomTypes();
    /** 
     * @brief Print the formal charges to the log (atom.GetPartialCharge(), 
     *  MMFF94 FC's are not always int).
     */ 
    void PrintFormalCharges();
    /** 
     * @brief Print the partial charges to the log.
     */ 
    void PrintPartialCharges();
    /** 
     * @brief Print the velocities to the log.
     */ 
    void PrintVelocities();
    //@}

    protected:
      virtual void ProcessOptions(std::vector<Option> &options);

  }; // class OBForceField

}// namespace OpenBabel
}

#endif   // OB_FORCEFIELD_H

//! @file forcefield.h
//! @brief Handle forcefields
