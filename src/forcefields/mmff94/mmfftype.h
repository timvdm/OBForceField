/**********************************************************************
forcefieldmmff94.h - MMFF94
 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#include <OBFFType>
#include <vector>
#include <string>

#include <openbabel/atom.h>
#include <openbabel/bond.h>

namespace OpenBabel {

  class OBRing;
  class OBMol;

namespace OBFFs {
  
  class OBParameterDB;
  
  // Class OBForceFieldMMFF94
  // class introduction in forcefieldmmff94.cpp
  class MMFF94Type : public OBFFType
  {
    public:
      /** 
       * Constructor
       */
      MMFF94Type();
      //! Destructor
      virtual ~MMFF94Type();
      
      //! Get the unit in which the energy is expressed
      std::string GetUnit() 
      { 
        return std::string("kcal/mol"); 
      }
 
      //! detect which rings are aromatic
      bool PerceiveAromaticity(/*const*/ OBMol &mol);
      //! \return Get the MMFF94 atom type for atom
      int GetType(OBAtom *atom);
      //! \return Sets atomtypes to MMFF94 in _mol
      bool SetTypes(const OBMol &mol);
      //!  Sets formal charges
      bool SetFormalCharges(/*const*/ OBMol &mol);
      //!  Sets partial charges
      bool SetPartialCharges(/*const*/ OBMol &mol);
      //! \return The row of the element atom in the periodic table
      int GetElementRow(OBAtom *atom);
      //! \return The bond type (BTIJ)
      int GetBondType(OBAtom* a, OBAtom* b);
      //! \return The angle type (ATIJK)
      int GetAngleType(OBAtom* a, OBAtom* b, OBAtom *c);
      //! \return The stretch-bend type (SBTIJK)
      int GetStrBndType(OBAtom* a, OBAtom* b, OBAtom *c);
      //! \return The torsion type (TTIJKL)
      int GetTorsionType(OBAtom* a, OBAtom* b, OBAtom *c, OBAtom *d);
      //! \return true if atomtype has sbmb set in mmffprop.par
      bool HasSbmbSet(int atomtype);
      //! \return true if atomtype has pilp set in mmffprop.par
      bool HasPilpSet(int atomtype);
      //! \return true if atomtype has arom set in mmffprop.par
      bool HasAromSet(int atomtype);
      //! \return true if atomtype has lin set in mmffprop.par
      bool HasLinSet(int atomtype);
      //! \return the crd value for the atomtype in mmffprop.par
      int GetCrd(int atomtype);
      //! \return the val value for the atomtype in mmffprop.par
      int GetVal(int atomtype);
      //! \return the mltb value for the atomtype in mmffprop.par
      int GetMltb(int atomtype);
      //! \return the level 2 equivalent atom type for type (mmffdef.par)
      int EqLvl2(int type);
      //! \return the level 3 equivalent atom type for type (mmffdef.par)
      int EqLvl3(int type);
      //! \return the level 4 equivalent atom type for type (mmffdef.par)
      int EqLvl4(int type);
      //! \return the level 5 equivalent atom type for type (mmffdef.par)
      int EqLvl5(int type);
      //! \return the canonical bond index
      unsigned int GetCXB(int type, int a, int b);
      //! \return the canonical angle index
      unsigned int GetCXA(int type, int a, int b, int c);
      //! \return the canonical stretch-bend index
      unsigned int GetCXS(int type, int a, int b, int c);
      //! \return the canonical out-of-plane index
      unsigned int GetCXO(int a, int b, int c, int d);
      //! \return the canonical torsion index
      unsigned int GetCXT(int type, int a, int b, int c, int d);
      //! \return the canonical bond-charge-increment index
      unsigned int GetCXQ(int type, int a, int b);
      //! \return the U value for the atom from table X page 631
      double GetUParam(OBAtom* atom);
      //! \return the Z value for the atom from table VI page 628
      double GetZParam(OBAtom* atom);
      //! \return the C value for the atom from table VI page 628
      double GetCParam(OBAtom* atom);
      //! \return the V value for the atom from table X page 631
      double GetVParam(OBAtom* atom);
      //! return the covalent radius from Blom and Haaland, value from etab if not available
      double GetCovalentRadius(OBAtom* a);
      //! return the bond length calculated with a modified version of the Schomaker-Stevenson rule
      double GetRuleBondLength(OBAtom* a, OBAtom* b);
      //! return the bond length from mmffbond.par, if not found, one is calculated with a modified version of the Schomaker-Stevenson rule
      double GetBondLength(OBAtom* a, OBAtom* b);
      /** 
       * @brief Check if two atoms are in the same ring. [NOTE: this function uses SSSR, 
       * this means that not all rings are found for bridged rings. This causes 
       * some problems with the MMFF94 validation.]
       *  
       * @param a Atom a.
       * @param b Atom b.
       * 
       * @return True if atom a and b are in the same ring
       */
      bool IsInSameRing(OBAtom* a, OBAtom* b);

      bool IsAromatic(OBAtom *atom) const { return m_aromAtoms.at(atom->GetIdx()-1); }
      bool IsAromatic(OBBond *bond) const { return m_aromBonds.at(bond->GetIdx()); }
      bool IsAromatic(OBRing *ring) const;


      int GetCachedType(OBAtom *atom) const { return m_types.at(atom->GetIdx()-1); }
      int GetCachedType(unsigned int idx) const { return m_types.at(idx); }
      double GetPartialCharge(unsigned int idx) { return m_pCharges.at(idx); }
      OBParameterDB *GetParameterDB() { return m_database; }

      std::string MakeBondName(const OBMol &mol, unsigned int iA, unsigned int iB);
      std::string MakeAngleName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC);
      std::string MakeStrBndName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC);
      std::string MakeTorsionName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD);
      std::string MakeOOPName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD);


      OBParameterDB             *m_database;
      
      //OBMol m_mol;
      std::vector<int> m_types;
      std::vector<double> m_fCharges;
      std::vector<double> m_pCharges;
      std::vector<int> m_aromAtoms;
      std::vector<int> m_aromBonds;
 
  }; // class OBForceFieldMM2

}
}// namespace OpenBabel

//! \file forcefieldmmff94.h
//! \brief MMFF94 force field

