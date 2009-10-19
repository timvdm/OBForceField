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
	std::string name;
	size_t iA, iB;
      };

      /**
       * Angle identifier.
       */
      struct AngleIdentifier
      {
	std::string name;
	size_t iA, iB, iC;
      };

      /**
       * Torsion identifier.
       */
      struct TorsionIdentifier
      {
	std::string name;
	size_t iA, iB, iC, iD;
      };

      /**
       * Out-of-plane identifier.
       */
      struct OOPIdentifier
      {
	std::string name;
	size_t iA, iB, iC, iD;
      };
      
      /**
       * Find atom types and initialize atom, bond, torsion and out-of-plane
       * identifiers. This function will also initialize the internal cache
       * for IsConnected, IsOneThree, IsOneFour.
       */
      virtual bool SetTypes(const OBMol & mol) = 0;

      virtual const std::string & GetAtomType(const size_t & i) const =0;
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
      virtual bool IsConnected(const size_t & iA, const size_t & iB) const =0;
      /**
       * @return True if atoms with index iA & iB are in a 1-3 relation.
       */
      virtual bool IsOneThree(const size_t & iA, const size_t & iB) const =0;
      /**
       * @return True if atoms with index iA & iB are in a 1-4 relation.
       */
      virtual bool IsOneFour(const size_t & iA, const size_t & iB) const =0;
    protected:
      std::vector<AtomIdentifier>    m_atoms;
      std::vector<BondIdentifier>    m_bonds;
      std::vector<AngleIdentifier>   m_angles;
      std::vector<TorsionIdentifier> m_torsions;
      std::vector<OOPIdentifier> m_oops;
    }; 
  }
}// namespace OpenBabel

//! \brief OBFFType force field atom types

