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

      struct BondIdentifier
      {
	std::string name;
	size_t iA, iB;
      };

      struct AngleIdentifier
      {

	std::string name;
	size_t iA, iB, iC;
      };

      struct TorsionIdentifier
      {
	std::string name;
	size_t iA, iB, iC, iD;
      };

      struct OOPIdentifier
      {
	std::string name;
	size_t iA, iB, iC, iD;
      };
      virtual bool SetTypes(const OBMol & mol) =0;
      virtual unsigned int Idx(const size_t & i) const =0;
      virtual size_t I(const unsigned int & idx) const =0;
      virtual const std::string & GetAtomType(const size_t & i) const =0;
      virtual const std::vector<AtomIdentifier> & GetAtoms() const =0;
      virtual const std::vector<BondIdentifier> & GetBonds() const =0;
      virtual const std::vector<AngleIdentifier> & GetAngles() const =0;
      virtual const std::vector<TorsionIdentifier> & GetTorsions() const =0;
      virtual const std::vector<OOPIdentifier> & GetOOPs() const =0;
      virtual bool IsConnected(const size_t & iA, const size_t & iB) const =0;
      virtual bool IsOneThree(const size_t & iA, const size_t & iB) const =0;
      virtual bool IsOneFour(const size_t & iA, const size_t & iB) const =0;
    }; 
  }
}// namespace OpenBabel

//! \brief OBFFType force field atom types

