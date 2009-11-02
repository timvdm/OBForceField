#include "OBVariant"
#include "gafftype.h"
#include "gaffparameter.h"
#include <openbabel/babelconfig.h>
#include <openbabel/oberror.h>
#include <openbabel/locale.h>
#include <openbabel/mol.h>
#include <openbabel/typer.h>
#include <algorithm>

using namespace std;

namespace OpenBabel {
  namespace OBFFs {

    GAFFTypeRules::GAFFTypeRules(const string &filename)
      : m_filename(filename) {
      m_initialized = ParseParamFile();
    }

    bool GAFFTypeRules::IsInitialized()
    {
      return m_initialized;
    }

    bool GAFFTypeRules::ParseParamFile()
    {
      // Set the locale for number parsing to avoid locale issues: PR#1785463
      obLocale.SetLocale();
  
      vector<string> vs;
      char buffer[256];
      bool valid_line=true;
      OBSmartsPattern * sp;

      // open data/_parFile
      ifstream ifs;
      if (OpenDatafile(ifs, m_filename).length() == 0) {
	obErrorLog.ThrowError(__FUNCTION__, "Cannot open type definitions file", obError);
	return false;
      }

      while (ifs.getline(buffer, 255)) {
	if (EQn(buffer, "atom", 4)) {
	  tokenize(vs, buffer);
	  sp = new OBSmartsPattern;
	  if (sp->Init(vs[1])){
	    m_vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
	  }
	  else {
	    delete sp;
	    sp = NULL;
	    obErrorLog.ThrowError(__FUNCTION__, " Could not parse atom type table from gaff.prm", obError);
	    return false;
	  }
	}
      }
      // return the locale to the original one
      obLocale.RestoreLocale();
      
      return true;
    }
    
    bool GAFFType::IsInitialized()
    {
      if (this->p_typerules==NULL)
	return false;
      else
	return this->p_typerules->IsInitialized();
    }

    bool GAFFType::SetTypes(const OBMol &mol)
    {
      if (!IsInitialized()){
	obErrorLog.ThrowError(__FUNCTION__, "Force field types are not initialized", obError);
      }

      m_numAtoms = mol.NumAtoms();

      OBAtom *atom;
      vector<OBAtom*>::iterator itr;
      unsigned int idx;
      size_t i;

      vector<vector<int> > mlist;
      vector<pair<OBSmartsPattern*,string> >::const_iterator itr1;
      vector<vector<int> >::const_iterator itr2;

      m_atoms.clear();
      m_atoms.resize(m_numAtoms);
      for (itr1 = p_typerules -> m_vexttyp.begin();itr1 != p_typerules -> m_vexttyp.end();++itr1) {
	if (itr1->first->Match(const_cast<OBMol&>(mol))) {
	  mlist = itr1->first->GetMapList();
	  for (itr2 = mlist.begin();itr2 != mlist.end();++itr2) {
	    m_atoms[ (*itr2)[0] - 1 ]=(itr1->second).c_str();
	  }
	}
      }
      
      // Implementation of a special feature of GAFF concerning conjugated bonds
      // In a conjugated ring system cc-cc, and cd-cd are single conjugated bonds, cc-cd are double ones
      // All these bonds are initially of type cc-cc.
      // Some cc's are changed to cd.
      // Similar for ce and cf, cp and cq, nc and nd, ne and nf, pc and pd, pe and pf
      // Initially types cc, ce, cp, nc, ne , pc and pe are assigned.
      // Some of these are changed into cd, cf, cq, nd, nf, pd and pf
      OBBitVec visited;
      OBAtom *a, *b;
      unsigned int BO;
      size_t ia, ib;
      visited.Resize(m_numAtoms);
      visited.Clear();
      for (i=0; i!= m_numAtoms; ++i){
	if (!visited[i]){
	  if( m_atoms[i]=="cc" || 
	      m_atoms[i]=="ce" ||
	      m_atoms[i]=="cp" || 
	      m_atoms[i]=="nc" || 
	      m_atoms[i]=="ne" || 
	      m_atoms[i]=="pc" || 
	      m_atoms[i]=="pe"){
	    FOR_BONDS_OF_ATOM(bond,mol.GetAtom(i+1)) {
	      a = bond->GetBeginAtom();
	      b = bond->GetEndAtom();
	      ia = a->GetIdx() - 1;
	      ib = b->GetIdx() - 1;
	      if (m_atoms[ia]=="cc"||m_atoms[ia]=="cd"&&m_atoms[ib]=="cc"||m_atoms[ib]=="cd"){
		BO = bond->GetBondOrder();
		if ( (BO > 1  && m_atoms[ia]==m_atoms[ib]) || (BO==1 && m_atoms[ia]!=m_atoms[ib]) ){
		  if (!visited[ia]){
		    m_atoms[ia]="cd";
		  }
		  else if (!visited[ib]){
		    m_atoms[ib]=("cd");
		  }
		}
		visited.SetBitOn(ia);
		visited.SetBitOn(ib);
	      }
	      if (m_atoms[ia]=="ce"||m_atoms[ia]=="cf"&&(m_atoms[ib]=="ce"||m_atoms[ib]=="cf")){
		BO = bond->GetBondOrder();
		if ( (BO > 1  && m_atoms[ia]==m_atoms[ib]) || (BO==1 && m_atoms[ia]!=m_atoms[ib])){
		  if (!visited[ia]){
		    m_atoms[ia]=("cf");
		  }
		  else if (!visited[ib]){
		    m_atoms[ib]=("cf");
		  }
		}
		visited.SetBitOn(ia);
		visited.SetBitOn(ib);
	      }
	      if ((m_atoms[ia]=="cp"||m_atoms[ia]=="cq") && (m_atoms[ib]=="cp"||m_atoms[ib]=="cq") ){
		BO = bond->GetBondOrder();
		if ( (BO > 1  && m_atoms[ia]==m_atoms[ib]) || (BO==1 && m_atoms[ia]!=m_atoms[ib])){
		  if (!visited[ia]){
		    m_atoms[ia]=("cq");
		  }
		  else if (!visited[ib]){
		    m_atoms[ib]=("cq");
		  }
		}
		visited.SetBitOn(ia);
		visited.SetBitOn(ib);
	      }
	      if ((m_atoms[ia]=="nc"||m_atoms[ia]=="nd") && (m_atoms[ib]=="nc"||m_atoms[ib]=="nd") ) {
		BO = bond->GetBondOrder();
		if ( (BO > 1  && m_atoms[ia]==m_atoms[ib]) || (BO==1 && m_atoms[ia]!=m_atoms[ib]) ){
		  if (!visited[ia]){
		    m_atoms[ia]=("nd");
		  }
		  else if (!visited[ib]){
		    m_atoms[ib]=("nd");
		  }
		}
		visited.SetBitOn(ia);
		visited.SetBitOn(ib);
	      }
	      if ((m_atoms[ia]=="ne"|| m_atoms[ia]=="nf" ) && (m_atoms[ib] == "ne" || m_atoms[ib]=="nf")) {
		BO = bond->GetBondOrder();
		if ( (BO > 1  && m_atoms[ia]==m_atoms[ib]) || (BO==1 && m_atoms[ia]!=m_atoms[ib]) ){
		  if (!visited[ia]){
		    m_atoms[ia]=("nf");
		  }
		  else if (!visited[ib]){
		    m_atoms[ib]=("nf");
		  }
		}
		visited.SetBitOn(ia);
		visited.SetBitOn(ib);
	      }
	      if ((m_atoms[ia]=="pc" ||m_atoms[ia]=="pd")&&(m_atoms[ib]=="pc"||m_atoms[ib]=="pd")){
		BO = bond->GetBondOrder();
		if ( (BO > 1  && m_atoms[ia]==m_atoms[ib]) || (BO==1 && m_atoms[ia]!=m_atoms[ib])){
		  if (!visited[ia]){
		    m_atoms[ia]=("pd");
		  }
		  else if (!visited[ib]){
		    m_atoms[ib]=("pd");
		  }
		}
		visited.SetBitOn(ia);
		visited.SetBitOn(ib);
	      }
	      if ((m_atoms[ia]=="pe"||m_atoms[ia]=="pf")&&(m_atoms[ib]=="pe"||m_atoms[ib]=="pf")){
		BO = bond->GetBondOrder();
		if ( (BO > 1  && m_atoms[ia]==m_atoms[ib]) || (BO==1 && m_atoms[ia]!=m_atoms[ib])){
		  if (!visited[ia]){
		    m_atoms[ia]=("pf");
		  }
		  else if (!visited[ib]){
		    m_atoms[ib]=("pf");
		  }
		}
		visited.SetBitOn(ia);
		visited.SetBitOn(ib);
	      }
	    }
	  }
	  visited.SetBitOn(i);
	}
      }

      /*
      m_bonds.clear();
      m_bonds.reserve(mol.NumBonds());
      OBFFType::BondIdentifier bondID;
      FOR_BONDS_OF_MOL(bond,const_cast<OBMol&>(mol)){
	bondID.iA = bond->GetBeginAtom()->GetIdx() - 1;
	bondID.iB = bond->GetEndAtom()->GetIdx() - 1;
	bondID.name=MakeBondName(m_atoms[bondID.iA],m_atoms[bondID.iB]);
	m_bonds.push_back(bondID);
      }

      size_t ic, id;
      m_angles.clear();
      OBFFType::AngleIdentifier angleID;
      FOR_ANGLES_OF_MOL(angle,const_cast<OBMol&>(mol)){
	angleID.iB = (*angle)[0];
	angleID.iA = (*angle)[1];
	angleID.iC = (*angle)[2];
	angleID.name=MakeAngleName(m_atoms[angleID.iA], m_atoms[angleID.iB], m_atoms[angleID.iC]);
	m_angles.push_back(angleID);
      }

      m_torsions.clear();
      OBFFType::TorsionIdentifier torsionID;
      FOR_TORSIONS_OF_MOL(t,const_cast<OBMol&>(mol)) {
	torsionID.iA = (*t)[0];
	torsionID.iB = (*t)[1];
	torsionID.iC = (*t)[2];
	torsionID.iD = (*t)[3];
	torsionID.name=MakeTorsionName(m_atoms[torsionID.iA], m_atoms[torsionID.iB], m_atoms[torsionID.iC], m_atoms[torsionID.iD]);
	m_torsions.push_back(torsionID);
      }

      m_oops.clear();
      OBFFType::OOPIdentifier oopID;
      FOR_ATOMS_OF_MOL(atom, const_cast<OBMol&>(mol)) {
	b = (OBAtom*) &*atom;
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
	      oopID.name=MakeOOPName(m_atoms[oopID.iA], m_atoms[oopID.iB], m_atoms[oopID.iC], m_atoms[oopID.iD]);
	      m_oops.push_back(oopID);
	    }
	  }
	}
      }
      m_oops.clear();
      */
     
     /* 
      OBAtom *c, *d;
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
      */

      return true;
    }

    bool GAFFType::ValidateTypes(GAFFParameterDB * pdatabase)
    {
      vector<OBParameterDBTable::Query> query;
      std::vector<OBVariant> row;
      bool valid(true), valid_tmp;
      string name;
      vector<string> vs, names;
      map<string,string> aliases;
      map<string,string>::iterator itr2;

      //check if atom type parameters in database
      //check if bond-types are in database
      BondIdentifier bond;
      vector<BondIdentifier> bonds_cleaned;
      bonds_cleaned.reserve(m_bonds.size());
      OBParameterDBTable * pTable = (pdatabase->GetTable("Bond Harmonic"));
      aliases.clear();
      for(vector<OBFFType::BondIdentifier>::const_iterator itr=m_bonds.begin();itr!=m_bonds.end();++itr){
	bond=*itr;
	valid_tmp=true;
	itr2=aliases.find(bond.name);
	if (itr2 != aliases.end()){
	  bond.name=itr2->second;
	  if (bond.name != "-") bonds_cleaned.push_back(bond);
	}
	else {
	  name = bond.name;
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(name)));
	  row = pTable->FindRow(query);
	  if (row.size() == 0){
	    obErrorLog.ThrowError(__FUNCTION__, "bond type: " + bond.name + " not present in database", obInfo);
	    tokenize(vs, bond.name, "-", 2);
	    names=MakeAlternativeBondNames(vs[0],vs[1]);
	    valid_tmp=false;
	    for(unsigned int i=0; (i!=names.size() && !valid_tmp); ++i) {
	      query.clear();
	      query.push_back( OBParameterDBTable::Query(0, names[i]));
	      row = pTable->FindRow(query);
	      if (row.size() != 0) {
		valid_tmp=true;
		bond.name = names[i];
	      }
	    }
	  }
	  if (valid_tmp){
	    obErrorLog.ThrowError(__FUNCTION__, "Using bond type: " + bond.name + " instead", obInfo);
	    bonds_cleaned.push_back(bond);
	    }
	  else {
	    valid = false;
	    obErrorLog.ThrowError(__FUNCTION__, "Please supply parameters for bond type " + bond.name, obInfo);
	    bond.name = "-";
	  }
	}
	aliases.insert(pair<string, string>(name, bond.name));
      }
      m_bonds=bonds_cleaned;

      AngleIdentifier angle;
      vector<AngleIdentifier> angles_cleaned;
      angles_cleaned.reserve(m_angles.size());
      pTable = (pdatabase->GetTable("Angle Harmonic"));
      aliases.clear();
      for(vector<OBFFType::AngleIdentifier>::const_iterator itr=m_angles.begin();itr!=m_angles.end();++itr){
	angle = *itr;
	valid_tmp=true;
	itr2=aliases.find(angle.name);
	if (itr2 != aliases.end()){
	  angle.name=itr2->second;
	  if (angle.name != "--") angles_cleaned.push_back(angle);
	}
	else {
	  name = angle.name;
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(name)));
	  row = pTable->FindRow(query);
	  if (row.size() == 0){
	    obErrorLog.ThrowError(__FUNCTION__, "angle type: " + angle.name + " not present in database", obInfo);
	    tokenize(vs, angle.name, "-", 3);
	    names=MakeAlternativeAngleNames(vs[0],vs[1],vs[2]);
	    valid_tmp=false;
	    for(unsigned int i=0; (i!=names.size() && !valid_tmp); ++i) {
	      query.clear();
	      query.push_back( OBParameterDBTable::Query(0, names[i]));
	      row = pTable->FindRow(query);
	      if (row.size() != 0) {
		valid_tmp=true;
		angle.name = names[i];
	      }
	    }
	  }
	  if (valid_tmp){
	    obErrorLog.ThrowError(__FUNCTION__, "Using angle type: " + angle.name + " instead", obInfo);
	    angles_cleaned.push_back(angle);
	  }
	  else {
	    valid = false;
	    obErrorLog.ThrowError(__FUNCTION__, "Please supply parameters for angle type " + angle.name, obInfo);
	    angle.name = "--";
	  }
	}
	aliases.insert(pair<string, string>(name, angle.name));
      }
      m_angles=angles_cleaned;

      TorsionIdentifier torsion;
      vector<TorsionIdentifier> torsions_cleaned;
      torsions_cleaned.reserve(m_torsions.size());
      pTable = (pdatabase->GetTable("Torsion Harmonic"));
      aliases.clear();
      for(vector<OBFFType::TorsionIdentifier>::const_iterator itr=m_torsions.begin();itr!=m_torsions.end();++itr){
	torsion = *itr;
	valid_tmp=true;
	itr2=aliases.find(torsion.name);
	if (itr2 != aliases.end()){
	  torsion.name=itr2->second;
	  if (torsion.name != "--") torsions_cleaned.push_back(torsion);
	}
	else {
	  name = torsion.name;
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(name)));
	  row = pTable->FindRow(query);
	  if (row.size() == 0){
	    obErrorLog.ThrowError(__FUNCTION__, "torsion type: " + torsion.name + " not present in database", obInfo);
	    tokenize(vs, torsion.name, "-",4);
	    names=MakeAlternativeTorsionNames(vs[0],vs[1],vs[2],vs[3]);
	    valid_tmp=false;
	    for(unsigned int i=0; (i!=names.size() && !valid_tmp); ++i) {
	      query.clear();
	      query.push_back( OBParameterDBTable::Query(0, names[i]));
	      row = pTable->FindRow(query);
	      if (row.size() != 0) {
		valid_tmp=true;
		torsion.name = names[i];
	      }
	    }
	  }
	  if (valid_tmp){
	    obErrorLog.ThrowError(__FUNCTION__, "Using torsion type: " + torsion.name + " instead", obInfo);
	    torsions_cleaned.push_back(torsion);
	  }
	  else {
	    valid = false;
	    obErrorLog.ThrowError(__FUNCTION__, "Please supply parameters for torsion type " + torsion.name, obInfo);
	    torsion.name="---";
	  }
	}
	  aliases.insert(pair<string, string>(name, torsion.name));
      }
      m_torsions=torsions_cleaned;

      OOPIdentifier oop;
      vector<OOPIdentifier> oops_cleaned;
      oops_cleaned.reserve(m_oops.size());
      pTable = (pdatabase->GetTable("Torsion Harmonic OOP"));
      aliases.clear();
      for(vector<OBFFType::OOPIdentifier>::iterator itr=m_oops.begin();itr!=m_oops.end();++itr){
	oop = *itr;
	valid_tmp=true;
	itr2=aliases.find(oop.name);
	if (itr2 != aliases.end()){
	  oop.name=itr2->second;
	  if (oop.name != "---") {
	    oops_cleaned.push_back(oop);
	  }
	}
	else {
	  name = oop.name;
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(name)));
	  row = pTable->FindRow(query);
	  if (row.size() == 0){
	    obErrorLog.ThrowError(__FUNCTION__, "oop type: " + oop.name + " not present in database", obInfo);
	    tokenize(vs, oop.name, "-", 4);
	    names=MakeAlternativeOOPNames(vs[0],vs[1],vs[2],vs[3]);
	    valid_tmp=false;
	    for(unsigned int i=0; (i!=names.size() && !valid_tmp); ++i) {
	      query.clear();
	      query.push_back( OBParameterDBTable::Query(0, names[i]));
	      row = pTable->FindRow(query);
	      if (row.size() != 0) {
		valid_tmp=true;
		oop.name = names[i];
	      }
	    }
	  }
	  if (valid_tmp){
	    obErrorLog.ThrowError(__FUNCTION__, "Using oop type: " + oop.name + " instead", obInfo);
	    oops_cleaned.push_back(oop);
	  }
	  else {
	    obErrorLog.ThrowError(__FUNCTION__, "Assuming no OOP potential for: " + oop.name, obInfo);
	    oop.name="---";
	  }
	}
	aliases.insert(pair<string, string>(name, oop.name));
      }
      m_oops=oops_cleaned;

      vector<AtomIdentifier> atoms_cleaned;
      AtomIdentifier atom;
      atoms_cleaned.reserve(m_atoms.size());
      pTable = (pdatabase->GetTable("LJ6_12"));
      aliases.clear();
      for(vector<OBFFType::AtomIdentifier>::const_iterator itr=m_atoms.begin();itr!=m_atoms.end();++itr){
	atom=*itr;
	valid_tmp=true;
	itr2=aliases.find(atom);
	if (itr2 != aliases.end()){
	  atom=itr2->second;
	  if (atom != "-") atoms_cleaned.push_back(atom);
	}
	else {
	  name = atom;
	  query.clear();
	  query.push_back( OBParameterDBTable::Query(0, OBVariant(name)));
	  row = pTable->FindRow(query);
	  if (row.size() == 0){
	    obErrorLog.ThrowError(__FUNCTION__, "atom type: " + atom + " not present in database", obInfo);
	    names=MakeAlternativeAtomNames(atom);
	    valid_tmp=false;
	    for(unsigned int i=0; (i!=names.size() && !valid_tmp); ++i) {
	      query.clear();
	      query.push_back( OBParameterDBTable::Query(0, names[i]));
	      row = pTable->FindRow(query);
	      if (row.size() != 0) {
		valid_tmp=true;
		atom = names[i];
	      }
	    }
	  }
	  if (valid_tmp){
	    obErrorLog.ThrowError(__FUNCTION__, "Using atom type: " + atom + " instead", obInfo);
	    atoms_cleaned.push_back(atom);
	  }
	  else {
	    valid = false;
	    obErrorLog.ThrowError(__FUNCTION__, "Please supply parameters for atom type " + atom , obInfo);
	    atom = "-";
	  }
	  aliases.insert(pair<string, string>(name, atom ));
	}
      }
      m_atoms=atoms_cleaned;

      return valid;
    }

    /*
    bool GAFFType::IsConnected(const size_t &idxA, const size_t &idxB) const
    {
      return (m_Connected.find((idxA-1)+m_numAtoms*(idxB-1))!=m_Connected.end());
    }

    bool GAFFType::IsOneThree(const size_t &idxA, const size_t &idxB) const
    {
      return (m_OneThree.find((idxA-1)+m_numAtoms*(idxB-1))!=m_OneThree.end());
    }

    bool GAFFType::IsOneFour(const size_t &idxA, const size_t &idxB) const
    {
      return (m_OneFour.find((idxA-1)+m_numAtoms*(idxB-1))!=m_OneFour.end());
    }


    const string & GAFFType::GetAtomType(const size_t & idx) const
    {
      return m_atoms[idx-1];
    }

    const vector<OBFFType::AtomIdentifier> & GAFFType::GetAtoms() const
    {
      return m_atoms;
    }

    const vector<OBFFType::BondIdentifier> & GAFFType::GetBonds() const
    {
      return m_bonds;
    }

    const vector<OBFFType::AngleIdentifier> & GAFFType::GetAngles() const
    {
      return m_angles;
    }

    const vector<OBFFType::TorsionIdentifier> & GAFFType::GetTorsions() const
    {
      return m_torsions;
    }

    const vector<OBFFType::OOPIdentifier> & GAFFType::GetOOPs() const
    {
      return m_oops;
    }
    */

    vector<string> GAFFType::MakeAlternativeAtomNames(string aName)
    {
      vector<string> names;
      names.push_back("X");
      return names;
    }

    std::string GAFFType::MakeBondName(const OBMol &mol, unsigned int iA, unsigned int iB)
    {
      vector<string> names;
      names.reserve(2);
      names.push_back(m_atoms.at(iA));
      names.push_back(m_atoms.at(iB));
      sort(names.begin(), names.end());
      return names[0]+ "-" + names[1];
    }
    
    std::string GAFFType::MakeAngleName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC)
    {
      vector<string> names;
      names.reserve(2);
      names.push_back(m_atoms.at(iA));
      names.push_back(m_atoms.at(iC));
      sort(names.begin(), names.end());
      return names[0]+ "-" + m_atoms.at(iC) + "-" + names[1];
    }
      
    std::string GAFFType::MakeTorsionName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD)
    {
      if (m_atoms.at(iA) < m_atoms.at(iD) || (m_atoms.at(iA) == m_atoms.at(iD) && m_atoms.at(iB) < m_atoms.at(iC)))
	return m_atoms.at(iA) + "-" + m_atoms.at(iB) + "-" + m_atoms.at(iC) + "-" + m_atoms.at(iD);
      else
	return m_atoms.at(iD) + "-" + m_atoms.at(iC) + "-" + m_atoms.at(iB) + "-" + m_atoms.at(iA);
 
    }

    
    std::string GAFFType::MakeOOPName(const OBMol &mol, unsigned int iA, unsigned int iB, unsigned int iC, unsigned int iD)
    {
      vector<string> names;
      names.reserve(3);
      names.push_back(m_atoms.at(iA));
      names.push_back(m_atoms.at(iC));
      names.push_back(m_atoms.at(iD));
      sort(names.begin(), names.end());
      return names[0]+ "-" + m_atoms.at(iB) + "-" + names[1] + "-" + names[2];
    }


    string GAFFType::MakeBondName(string aName, string bName)
    {
      vector<string> names;
      names.reserve(2);
      names.push_back(aName);
      names.push_back(bName);
      sort(names.begin(), names.end());
      return names[0]+ "-" + names[1];
    }

    vector<string> GAFFType::MakeAlternativeBondNames(string aName, string bName)
    {
      vector<string> names;
      names.push_back(MakeBondName(aName,"X"));
      names.push_back(MakeBondName(bName,"X"));
      names.push_back(MakeBondName("X","X"));
      return names;
    }

    string GAFFType::MakeAngleName(string  aName, string  bName, string  cName)
    {
      vector<string> names;
      names.reserve(2);
      names.push_back(aName);
      names.push_back(cName);
      sort(names.begin(), names.end());
      return names[0]+ "-" + bName + "-" + names[1];
    }

    vector<string> GAFFType::MakeAlternativeAngleNames(string  aName, string  bName, string  cName)
    {
      vector<string> names;
      names.push_back(MakeAngleName(aName,bName,"X"));
      names.push_back(MakeAngleName("X",bName,cName));
      names.push_back(MakeAngleName(aName,"X",cName));
      names.push_back(MakeAngleName(aName,"X","X"));
      names.push_back(MakeAngleName("X",bName,"X"));
      names.push_back(MakeAngleName("X","X",cName));
      names.push_back(MakeAngleName("X","X","X"));
      return names;
    }

    string GAFFType::MakeTorsionName(string  aName, string  bName, string  cName, string  dName)
    {
      if (aName < dName|| (aName == dName && bName < cName))
	return aName + "-" + bName + "-" + cName + "-" + dName;
      else
	return dName + "-" + cName + "-" + bName + "-" + aName;
    }
    
    vector<string> GAFFType::MakeAlternativeTorsionNames(string  aName, string  bName, string  cName, string  dName)
    {
      vector<string> names;
      names.push_back(MakeTorsionName(aName,bName,cName,"X"));
      names.push_back(MakeTorsionName(aName,bName,"X",dName));
      names.push_back(MakeTorsionName(aName,"X",cName,dName));
      names.push_back(MakeTorsionName("X",bName,cName,dName));
      names.push_back(MakeTorsionName("X",bName,cName,"X"));
      names.push_back(MakeTorsionName(aName,"X",cName,"X"));
      names.push_back(MakeTorsionName(aName,bName,"X","X"));
      names.push_back(MakeTorsionName("X",bName,"X",dName));
      names.push_back(MakeTorsionName(aName,"X","X",dName));
      names.push_back(MakeTorsionName("X","X",cName,dName));
      names.push_back(MakeTorsionName(aName,"X","X","X"));
      names.push_back(MakeTorsionName("X",bName,"X","X"));
      names.push_back(MakeTorsionName("X","X",cName,"X"));
      names.push_back(MakeTorsionName("X","X","X",dName));
      names.push_back(MakeTorsionName("X","X","X","X"));
      return names;
    }

    string GAFFType::MakeOOPName(string  aName, string  bName, string  cName, string  dName)
    {
      vector<string> names;
      names.reserve(3);
      names.push_back(aName);
      names.push_back(cName);
      names.push_back(dName);
      sort(names.begin(), names.end());
      return names[0]+ "-" + bName + "-" + names[1] + "-" + names[2];
    }

    vector<string> GAFFType::MakeAlternativeOOPNames(string  aName, string  bName, string  cName, string  dName)
    {
      vector<string> names;
      names.push_back(MakeOOPName(aName,bName,cName,"X"));
      names.push_back(MakeOOPName(aName,bName,"X",dName));
      names.push_back(MakeOOPName(aName,"X",cName,dName));
      names.push_back(MakeOOPName("X",bName,cName,dName));
      names.push_back(MakeOOPName("X",bName,cName,"X"));
      names.push_back(MakeOOPName(aName,"X",cName,"X"));
      names.push_back(MakeOOPName(aName,bName,"X","X"));
      names.push_back(MakeOOPName("X",bName,"X",dName));
      names.push_back(MakeOOPName(aName,"X","X",dName));
      names.push_back(MakeOOPName("X","X",cName,dName));
      names.push_back(MakeOOPName(aName,"X","X","X"));
      names.push_back(MakeOOPName("X",bName,"X","X"));
      names.push_back(MakeOOPName("X","X",cName,"X"));
      names.push_back(MakeOOPName("X","X","X",dName));
      names.push_back(MakeOOPName("X","X","X","X"));
      return names;
    }
  } 
}
