#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include "obgasteiger.h"

using namespace std;

namespace OpenBabel {
  namespace OBFFs {

    bool OBGasteiger::ComputeCharges(OBMol & mol)
    {
      if (mol.AutomaticFormalCharge() && mol.AutomaticPartialCharge() && mol.HasPartialChargesPerceived()) //Gasteiger charges are already stored in mol 
	this->CopyFromMol(mol);
      else {
	OBAtom *atom;
	vector<OBAtom*>::iterator itr;
	vector<double>::iterator itr2;
	vector<short>::iterator itr3;
	vector<double> partialCharges;
	vector<short> formalCharges;

	// backup the charges stored in the molecule
	bool isSetPerceived=mol.HasPartialChargesPerceived();
	bool isAutomaticFormal(mol.AutomaticFormalCharge());
	bool isAutomaticPartial(mol.AutomaticPartialCharge());
	for (atom = mol.BeginAtom(itr); atom; atom = mol.NextAtom(itr))
	  {
	    partialCharges.push_back((double)atom->GetPartialCharge());
	    formalCharges.push_back(atom->GetFormalCharge());
	  }

	// new Gasteiger charges are computed and copied
	mol.SetAutomaticFormalCharge(true);
	mol.SetAutomaticPartialCharge(true);
	mol.UnsetPartialChargesPerceived();
	m_partialCharges.clear();
	m_partialCharges.reserve(mol.NumAtoms());
	m_formalCharges.clear();
	m_formalCharges.reserve(mol.NumAtoms());
	for (atom = mol.BeginAtom(itr);atom;atom = mol.NextAtom(itr)) {
	  m_partialCharges.push_back(atom->GetPartialCharge());
	  m_formalCharges.push_back((double)atom->GetFormalCharge());
	}
	// restore charges of molecule
	for (atom = mol.BeginAtom(itr), itr2=partialCharges.begin(), itr3= formalCharges.begin() ;atom;atom = mol.NextAtom(itr), ++itr2, ++itr3) {
	  atom->SetPartialCharge(*itr2);
	  atom->SetFormalCharge(*itr3);
	}
	mol.SetAutomaticFormalCharge(isAutomaticFormal);
	mol.SetAutomaticPartialCharge(isAutomaticPartial);
	(isSetPerceived)? mol.SetPartialChargesPerceived() : mol.UnsetPartialChargesPerceived();
      }
      return true;
    }
  }
} // end namespace OpenBabel
