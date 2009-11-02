/*********************************************************************
OBFF parameter database

Copyright (C) 2006-2009 by Tim Vandermeersch
 
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

#include "obfftype.h"

#include <openbabel/mol.h>

using namespace std;

namespace OpenBabel {
  namespace OBFFs {

    bool OBFFType::Setup(const OBMol &mol)
    {
      cout << "OBFFType::Setup()" << endl;
      if (!SetTypes(mol))
        return false;
      InitIdentifiers(mol);
      InitOneX(mol);
      return true;
    }
   
    void OBFFType::InitIdentifiers(const OBMol &mol)
    {
      cout << "OBFFType::InitIdentifiers()" << endl;
      m_bonds.clear();
      m_bonds.reserve(mol.NumBonds());
      OBFFType::BondIdentifier bondID;
      FOR_BONDS_OF_MOL(bond, const_cast<OBMol&>(mol)){
	bondID.iA = bond->GetBeginAtom()->GetIdx() - 1;
	bondID.iB = bond->GetEndAtom()->GetIdx() - 1;
	bondID.name = MakeBondName(mol, bondID.iA, bondID.iB);
        cout << "    " << bondID.name << endl;
	m_bonds.push_back(bondID);
      }

      m_angles.clear();
      OBFFType::AngleIdentifier angleID;
      FOR_ANGLES_OF_MOL(angle,const_cast<OBMol&>(mol)){
	angleID.iB = (*angle)[0];
	angleID.iA = (*angle)[1];
	angleID.iC = (*angle)[2];
	angleID.name = MakeAngleName(mol, angleID.iA, angleID.iB, angleID.iC);
	m_angles.push_back(angleID);
      }

      m_torsions.clear();
      OBFFType::TorsionIdentifier torsionID;
      FOR_TORSIONS_OF_MOL(t,const_cast<OBMol&>(mol)) {
	torsionID.iA = (*t)[0];
	torsionID.iB = (*t)[1];
	torsionID.iC = (*t)[2];
	torsionID.iD = (*t)[3];
	torsionID.name = MakeTorsionName(mol, torsionID.iA, torsionID.iB, torsionID.iC, torsionID.iD);
	m_torsions.push_back(torsionID);
      }

      m_oops.clear();
      OBFFType::OOPIdentifier oopID;
      FOR_ATOMS_OF_MOL(atom, const_cast<OBMol&>(mol)) {
	OBAtom *b = (OBAtom*) &*atom;
	oopID.iB = b->GetIdx() - 1;
	for( OBAtomAtomIter nbr1(b); nbr1; ++nbr1 ) {
	  oopID.iA = nbr1->GetIdx() -1;
	  OBAtomAtomIter nbr2=nbr1;
	  ++nbr2;
	  for( ; nbr2; ++nbr2 ) {
	    oopID.iC = nbr2->GetIdx() - 1;
	    OBAtomAtomIter nbr3=nbr2;
	    ++nbr3;
	    for( ; nbr3; ++nbr3 ) {
	      oopID.iD = nbr3->GetIdx() - 1;
	      oopID.name = MakeOOPName(mol, oopID.iA, oopID.iB, oopID.iC, oopID.iD);
	      m_oops.push_back(oopID);
	    }
	  }
	}
      }
//      m_oops.clear();
     
    }


    void OBFFType::InitOneX(const OBMol &mol)
    {
      m_numAtoms = mol.NumAtoms();
      unsigned int ia, ib, ic, id;
      OBAtom *a, *b, *c, *d;
      OBBond *bond1, *bond2, *bond3;
      OBBondIterator itr3, itr4, itr5;
      FOR_ATOMS_OF_MOL(atom, const_cast<OBMol&>(mol)) {
	a = (OBAtom*) &*atom;
	ia = a->GetIdx() - 1;
	for (bond1 = a->BeginBond(itr3);bond1;bond1 = a->NextBond(itr3)){
	  if (bond1->GetBeginAtom() == a)
	    b = bond1->GetEndAtom();
	  else
	    b = bond1->GetBeginAtom();
	  ib = b->GetIdx() - 1;
	  m_Connected.insert(ia+ib*m_numAtoms);
	  for (bond2 = b->BeginBond(itr4);bond2;bond2 = b->NextBond(itr4)){
	    if (bond2->GetBeginAtom() == b)
	      c = bond2->GetEndAtom();
	    else
	      c = bond2->GetBeginAtom();
	    if (c==a) continue;
	    ic = c->GetIdx() - 1;
	    m_OneThree.insert(ia+ic*m_numAtoms);
	    for (bond3 = c->BeginBond(itr5);bond3;bond3 = c->NextBond(itr5)){
	      if (bond3->GetBeginAtom() == c)
		d = bond3->GetEndAtom();
	      else
		d = bond3->GetBeginAtom();
	      if (d==b||d==a) continue;
	      id = d->GetIdx() - 1;
	      m_OneFour.insert(ia+id*m_numAtoms);
	    }
	  }
	}
      }

    }

    bool OBFFType::IsConnected(unsigned int idxA, unsigned int idxB) const
    {
      return (m_Connected.find((idxA)+m_numAtoms*(idxB)) != m_Connected.end());
    }

    bool OBFFType::IsOneThree(unsigned int idxA, unsigned int idxB) const
    {
      return (m_OneThree.find((idxA)+m_numAtoms*(idxB)) != m_OneThree.end());
    }

    bool OBFFType::IsOneFour(unsigned int idxA, unsigned int idxB) const
    {
      return (m_OneFour.find((idxA)+m_numAtoms*(idxB))!=m_OneFour.end());
    }


  } // end namespace OpenBabel
}
//! \brief OBFF parameter database
