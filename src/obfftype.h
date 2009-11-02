/**********************************************************************

 
Copyright (C) 2009 by Frank Peters <e.a.j.f.peters@gmail.com>
 
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

#include <vector>
#include <string>
#include <set>

namespace OpenBabel {

  class OBMol;
  class OBAtom;
  class OBFFParameterDB;

  namespace OBFFs {
  
    class OBFFType
    {
    public:
      typedef std::string AtomIdentifier;

      /**
       * Bond identifier
       */
      struct BondIdentifier
      {
	std::string name; //!< The name for the bond (e.g. "C-O", "1:C-O", "2-6").
	unsigned int iA; //!< The atom index for the bond's begin atom (indexed from 0 to N-1).
	unsigned int iB; //!< The atom index for the bond's end atom (indexed from 0 to N-1).
      };

      /**
       * Angle identifier.
       */
      struct AngleIdentifier
      {
	std::string name; //!< The name for the angle (e.g. "C-C-O", "1:C-C-O", "2-2-6").
	unsigned int iA; //!< The atom index for the angle's first terminal atom (indexed from 0 to N-1).
	unsigned int iB; //!< The atom index for the angle's vertex atom (indexed from 0 to N-1).
	unsigned int iC; //!< The atom index for the angle's second terminal atom (indexed from 0 to N-1).
      };

      /**
       * Torsion identifier.
       */
      struct TorsionIdentifier
      {
	std::string name; //!< The name for the torsion (e.g. "C-C-C-O", "1:C-C-C-O", "2-2-2-6").
	unsigned int iA; //!< The atom index for the torsion's first atom (indexed from 0 to N-1).
	unsigned int iB; //!< The atom index for the torsion's second atom (indexed from 0 to N-1).
	unsigned int iC; //!< The atom index for the torsion's third atom (indexed from 0 to N-1).
	unsigned int iD; //!< The atom index for the torsion's fourth atom (indexed from 0 to N-1).
      };

      /**
       * Out-of-plane identifier.
       */
      struct OOPIdentifier
      {
	std::string name; //!< The name for the out-of plane angle (e.g. "C-C-C-O", "1:C-C-C-O", "2-2-2-6").
	unsigned int iA; //!< The atom index for the oop's first atom (indexed from 0 to N-1).
	unsigned int iB; //!< The atom index for the oop's second atom (indexed from 0 to N-1).
	unsigned int iC; //!< The atom index for the oop's third atom (indexed from 0 to N-1).
	unsigned int iD; //!< The atom index for the oop's fourth atom (indexed from 0 to N-1).
      };

      bool Setup(const OBMol &mol);
      
      /**
       * Get the atom for atom with index @p index.
       */
      virtual const std::string & GetAtomType(unsigned int index) const
      {
        if (index < m_atoms.size())
          return m_atoms.at(index);
        return m_nullType;
      }
      /**
       * Get the atoms. These are the force field atom types as std::string.
       * @sa AtomIdentifier
       */
      virtual const std::vector<AtomIdentifier> & GetAtoms() const
      {
        return m_atoms;
      }
      /**
       * Get the bonds. These are BondIdentifier structs containing the name 
       * (e.g. "2-5", "C-O", "C-X") for the bond and the atom I's.
       *
       * @sa I() Idx() BondIdentifier
       */
      virtual const std::vector<BondIdentifier> & GetBonds() const
      {
        return m_bonds;
      }
      /**
       * Get the angles. These are AngleIdentifier structs containing the name 
       * (e.g. "2-5-2", "C-O-C", "X-C-X") for the angle and the atom I's.
       *
       * @sa I() Idx() AngleIdentifier
       */
      virtual const std::vector<AngleIdentifier> & GetAngles() const
      {
        return m_angles;
      }
      /**
       * Get the torsions. These are TorsionIdentifier structs containing the 
       * name (e.g. "2-5-6-2", "C-O-C-C", "X-C-C-X") for the angle and the atom I's.
       *
       * @sa I() Idx() TorsionIdentifier
       */
      virtual const std::vector<TorsionIdentifier> & GetTorsions() const
      {
        return m_torsions;
      }
      /**
       * Get the out-of-plane angles. These are OOPIdentifier structs containing the 
       * name (e.g. "2-5-6-2", "C-O-C-C", "X-C-C-X") for the oop angle and the atom I's.
       *
       * @sa I() Idx() TorsionIdentifier
       */

      virtual const std::vector<OOPIdentifier> & GetOOPs() const
      {
        return m_oops;
      }
      /**
       * @return True if atoms with index iA & iB are connected.
       */
      virtual bool IsConnected(unsigned int iA, unsigned int iB) const;
      /**
       * @return True if atoms with index iA & iB are in a 1-3 relation.
       */
      virtual bool IsOneThree(unsigned int iA, unsigned int iB) const;
      /**
       * @return True if atoms with index iA & iB are in a 1-4 relation.
       */
      virtual bool IsOneFour(unsigned int iA, unsigned int iB) const;
    protected:
      /**
       * Find atom types and initialize atom identifiers. 
       */
      virtual bool SetTypes(const OBMol & mol) = 0;

      /**
       * Initialize m_connected, m_oneThree and m_oneFour
       */
      void InitOneX(const OBMol &mol);
      
     
      /**
       * Initialize the term identifiers (e.g. BondIdentifier, AngleIdentifier).
       */
      void InitIdentifiers(const OBMol &mol);
      /**
       * Subclasses must implement this function to return the bond name for 
       * the bond between atoms with index iA and iB. These bond names are used
       * to find the correct parameters for a given bond. In the simplest case
       * these bond names can consist of the two atom types separated by a '-'
       * character (e.g. "C-O", "1-3"). 
       *
       * The actual format of the bond name doesn't matter as long as it is 
       * compatible with the OBFFParameterDB implementation used. For example,
       * the MMFF94 force field also has a "bond class" and possible bond
       * names look like "0:3-3" or "1:3-3". These two bond names are correctly
       * recognised by the MMFFParameterDB as a bond between atom types 3 and
       * bond class 0 and 1 respectively.
       */
      virtual std::string MakeBondName(const OBMol &mol, unsigned int iA, unsigned int iB) = 0;
      /**
       * Subclasses must implement this function to return the angle name for 
       * the angle between atoms with index iA, iB and iC. These angle names are 
       * used to find the correct parameters for a given angle. In the simplest case
       * these angle names can consist of the three atom types separated by a '-'
       * character (e.g. "C-C-O", "1-1-3"). 
       *
       * The actual format of the angle name doesn't matter as long as it is 
       * compatible with the OBFFParameterDB implementation used. For example,
       * the MMFF94 force field also has a "angle class" and possible angle
       * names look like "0:3-3-3" or "4:3-3-3". These two angle names are correctly
       * recognised by the MMFFParameterDB as a angle between atom types 3 and
       * angle class 0 and 4 respectively.
       */
      virtual std::string MakeAngleName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC) = 0;
      /**
       * Subclasses must implement this function to return the stretch-bend name for 
       * the angle between atoms with index iA, iB and iC. These stretch-bend names are 
       * used to find the correct parameters for a given stretch-bend interactions. In 
       * the simplest case these stretch-bend names can consist of the three atom types separated 
       * by a '-' character (e.g. "C-C-O", "1-1-3"). When a force field doesn't have a 
       * stretch-bend term, this function can just return an empty std::string.
       *
       * The actual format of the stretch-bend name doesn't matter as long as it is 
       * compatible with the OBFFParameterDB implementation used. For example,
       * the MMFF94 force field also has a "strbnd class" and possible stretch-bend
       * names look like "0:3-3-3" or "4:3-3-3". These two stretch-bend names are correctly
       * recognised by the MMFFParameterDB as a stretch-bend interaction between the atom types 
       * 3 and stretch-bend class 0 and 4 respectively.
       */
      virtual std::string MakeStrBndName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC) = 0;
      /**
       * Subclasses must implement this function to return the torsion name for 
       * the torsion between atoms with index iA, iB and iC. These torsion names are 
       * used to find the correct parameters for a given torsion. In the simplest case
       * these torsion names can consist of the three atom types separated by a '-'
       * character (e.g. "C-C-O", "1-1-3"). 
       *
       * The actual format of the torsion name doesn't matter as long as it is 
       * compatible with the OBFFParameterDB implementation used. For example,
       * the MMFF94 force field also has a "torsion class" and possible torsion
       * names look like "0:3-3-3" or "4:3-3-3". These two torsion names are correctly
       * recognised by the MMFFParameterDB as a torsion between atom types 3 and
       * torsion class 0 and 4 respectively.
       */
      virtual std::string MakeTorsionName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD) = 0;
      /**
       * Subclasses must implement this function to return the out-of-plane (oop) name for 
       * the oop between atoms with index iA, iB and iC. These oop names are 
       * used to find the correct parameters for a given oop. In the simplest case
       * these oop names can consist of the three atom types separated by a '-'
       * character (e.g. "C-C-O", "1-1-3"). 
       *
       * The actual format of the oop name doesn't matter as long as it is 
       * compatible with the OBFFParameterDB implementation used. For example,
       * the MMFF94 force field also has a "oop class" and possible oop
       * names look like "0:3-3-3" or "4:3-3-3". These two oop names are correctly
       * recognised by the MMFFParameterDB as a oop between atom types 3 and
       * oop class 0 and 4 respectively.
       */
      virtual std::string MakeOOPName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD) = 0;


      std::vector<AtomIdentifier>    m_atoms;
      std::vector<BondIdentifier>    m_bonds;
      std::vector<AngleIdentifier>   m_angles;
      std::vector<TorsionIdentifier> m_torsions;
      std::vector<OOPIdentifier> m_oops;
      std::string m_nullType;

      unsigned int m_numAtoms;
      std::set<unsigned long int> m_Connected;
      std::set<unsigned long int> m_OneThree;
      std::set<unsigned long int> m_OneFour;
    }; 
  }
}// namespace OpenBabel

//! \brief OBFFType force field atom types

