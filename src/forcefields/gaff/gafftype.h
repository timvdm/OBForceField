/**********************************************************************
forcefieldmmff94.h - GAFF
 
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
#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/typer.h>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <OBFFType>

namespace OpenBabel {
  class OBMol;
  namespace OBFFs {

    class GAFFTypeRules;
    class GAFFParameterDB;

    class GAFFType: public OBFFType
    {
    public:
      GAFFType(GAFFTypeRules * ptyperules=NULL) : p_typerules(ptyperules) {}
      bool IsInitialized();
      GAFFTypeRules * GetGAFFTypeRules() const {return p_typerules;}
      void SetGAFFTypeRules(GAFFTypeRules * ptyperules) {p_typerules=ptyperules;}
      bool SetTypes(const OBMol &mol);
      bool ValidateTypes(GAFFParameterDB * pdatabase); //Check if types are in database. If not check for default patterns. If found change name. If still not found remove interaction from list.
//      const std::string & GetAtomType(unsigned int idx) const;
      //const std::vector<AtomIdentifier> & GetAtoms() const;
      //const std::vector<BondIdentifier> & GetBonds() const;
      //const std::vector<AngleIdentifier> & GetAngles() const;
      //const std::vector<TorsionIdentifier> & GetTorsions() const;
      //const std::vector<OOPIdentifier> & GetOOPs() const;
      //bool IsConnected(const size_t & idxA, const size_t & idxB) const;
      //bool IsOneThree(const size_t & idxA, const size_t & idxB) const;
      //bool IsOneFour(const size_t & idxA, const size_t & idxB) const;
    protected:
      static std::vector<std::string> MakeAlternativeAtomNames(std::string aName);
      static std::string MakeBondName(std::string aName, std::string bName);
      static std::vector<std::string> MakeAlternativeBondNames(std::string aName, std::string bName);
      static std::string MakeAngleName(std::string aName, std::string bName, std::string cName);
      static std::vector<std::string> MakeAlternativeAngleNames(std::string aName, std::string bName, std::string cName);
      static std::string MakeTorsionName(std::string aName, std::string bName, std::string cName, std::string dName);
      static std::vector<std::string> MakeAlternativeTorsionNames(std::string aName, std::string bName, std::string cName, std::string dName);
      static std::string MakeOOPName(std::string aName, std::string bName, std::string cName, std::string dName);
      static std::vector<std::string> MakeAlternativeOOPNames(std::string aName, std::string bName, std::string cName, std::string dName);

      std::string MakeBondName(const OBMol &mol, unsigned int iA, unsigned int iB);
      std::string MakeAngleName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC);
      std::string MakeStrBndName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC) { return ""; }
      std::string MakeTorsionName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD);
      std::string MakeOOPName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD);


    private:
      GAFFTypeRules * p_typerules;
      unsigned int m_numAtoms;
      /* now defined in OBFFType
      std::vector<AtomIdentifier> m_atoms;
      std::vector<BondIdentifier> m_bonds;
      std::vector<AngleIdentifier> m_angles;
      std::vector<TorsionIdentifier> m_torsions;
      std::vector<OOPIdentifier> m_oops;
      */
      std::vector<double> m_masses;
      std::vector<double> m_charges;
      //std::set<unsigned long int> m_Connected;
      //std::set<unsigned long int> m_OneThree;
      //std::set<unsigned long int> m_OneFour;
      friend class GAFFParameterDB;
    };

    class GAFFTypeRules
    {
    public:
      GAFFTypeRules(const std::string & filename);
      bool IsInitialized();
    private:
      std::string m_filename;
      bool ParseParamFile();
      std::vector<std::pair<OBSmartsPattern*,std::string> > m_vexttyp; // external atom type rules
      bool m_initialized;
      friend bool GAFFType::SetTypes(const OBMol &mol);
    };

  }
}// namespace OpenBabel

//! \file forcefieldmmff94.h
//! \brief GAFF force field
