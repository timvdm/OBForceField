/**********************************************************************
obchargemethod.cpp - Provides charges for electrostatic term of a forcefield
 
Some portions Copyright (C) 2009 by Frank Peters
 
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
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <OBChargeMethod>

using namespace std;

namespace OpenBabel {
  namespace OBFFs {

    void OBChargeMethod::CopyFromMol(OBMol & mol)
    {
      OBAtom *atom;
      vector<OBAtom*>::iterator itr;
      m_partialCharges.clear();
      m_partialCharges.reserve(mol.NumAtoms());
      m_formalCharges.clear();
      m_formalCharges.reserve(mol.NumAtoms());

      for (atom = mol.BeginAtom(itr);atom;atom = mol.NextAtom(itr))
      {
	m_partialCharges.push_back((double)atom->GetPartialCharge());
	m_formalCharges.push_back((double)atom->GetFormalCharge());
      }
    }

    bool OBChargeMethod::CopyToMol(OBMol & mol) const
    {
      OBAtom *atom;
      vector<OBAtom*>::iterator itr;
      vector<double>::const_iterator itr2, itr3;

      if ( m_partialCharges.size() != mol.NumAtoms() )
	return false;
      
      for (atom = mol.BeginAtom(itr), itr2 = m_partialCharges.begin(), itr3 = m_formalCharges.begin();
	   atom; atom = mol.NextAtom(itr), ++itr2, ++itr3) {
	atom->SetPartialCharge(*itr2);
	atom->SetFormalCharge(*itr3);
      }

      mol.SetPartialChargesPerceived();
	  
    }
      
    const std::vector<double> & OBChargeMethod::GetPartialCharges() const
    {
      return m_partialCharges;
    }

    const std::vector<double> & OBChargeMethod::GetFormalCharges() const
    {
      return m_formalCharges;
    }

    bool OBChargeMethod::ComputeCharges(OBMol & mol)
    {
      return false;
    }

  }
} // end namespace OpenBabel

//! \file molchrg.cpp
//! \brief Assign charges
