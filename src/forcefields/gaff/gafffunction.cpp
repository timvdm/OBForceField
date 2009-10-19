/**********************************************************************
obfunction.h - Base class for force fields and scroring functions which
               depend on atom positions.
 
Copyright (C) 2008,2009 by Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the sGNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <OBForceField>
#include <OBLogFile>
#include <GAFF>

#include <openbabel/mol.h>

namespace OpenBabel {
  namespace OBFFs {

    class GAFFFunction : public OBForceField 
    {
    public:
      GAFFFunction();
      ~GAFFFunction();
      
      std::string GetName() const
      {
        return "GAFF";
      }
      
      bool Setup(/*const*/ OBMol &mol);
      void Compute(Computation computation = Value);
      double GetValue() const;
      
      std::string GetUnit() const { return "kcal/mol"; }
      
      bool HasAnalyticalGradients() const 
      { 
        return true; 
      }
    protected:
      void ProcessOptions(std::vector<Option> &options);
      std::string GetDefaultOptions() const;
      
      GAFFParameterDB *p_database;
      GAFFTypeRules *p_gaffTypeRules;
      GAFFType *p_gaffType;
      OBChargeMethod *p_charge;
      bool m_HaveCreatedDB, m_HaveCreatedTypeRules, m_HaveCreatedType, m_HaveCreatedCharge; 
    };

    GAFFFunction::GAFFFunction() 
      : m_HaveCreatedDB(false), m_HaveCreatedType(false), m_HaveCreatedCharge(false)
    {
      AddTerm(new BondHarmonic(this));
      AddTerm(new AngleHarmonic(this));
      AddTerm(new TorsionHarmonic(this));
      AddTerm(new TorsionHarmonic(this,"Torsion Harmonic OOP"));
      AddTerm(new LJ6_12(this, 0.5, LJ6_12::geometric));
      AddTerm(new Coulomb(this, 0.8333));
      //     AddTerm(new Coulomn(this, m_common));
    }

    GAFFFunction::~GAFFFunction(){
      if (m_HaveCreatedDB)
	delete p_database;
      if (m_HaveCreatedTypeRules){
	delete p_gaffTypeRules;
      }
      if (m_HaveCreatedType){
	delete p_gaffType;
      }
      if (m_HaveCreatedCharge)
	delete p_charge;
    }

    bool GAFFFunction::Setup(/*const*/ OBMol &mol)
    {
      p_gaffType = (GAFFType *) GetOBFFType();
      if (p_gaffType==NULL){
	p_gaffTypeRules = new GAFFTypeRules("../data/gaff.prm"); //here should be the default value to gaff.prm 
	if (p_gaffTypeRules==NULL)
	  return false;
	else {
	  m_HaveCreatedTypeRules = true;
	  p_gaffType = new GAFFType(p_gaffTypeRules);
	  if (p_gaffType==NULL)
	    return false;
	  else {
	    m_HaveCreatedType = true;
	    SetOBFFType(p_gaffType);
	  }
	}
      }
      else if (! (p_gaffType->IsInitialized()) ) {
	p_gaffTypeRules = new GAFFTypeRules("../data/gaff.prm"); //here should be the default value to gaff.prm 
	if (p_gaffTypeRules==NULL)
	  return false;
	else {
	  m_HaveCreatedTypeRules = true;
	  p_gaffType->SetGAFFTypeRules(p_gaffTypeRules);
	}
      }	
      p_database = (GAFFParameterDB *) GetParameterDB();
      if (p_database==NULL){
	p_database = new GAFFParameterDB("../data/gaff.dat");;
	if (p_database==NULL)
	  return false;
	else {
	  m_HaveCreatedDB = true;
	  SetParameterDB(p_database);
	}
      }
      p_charge = (OBChargeMethod *) GetOBChargeMethod();
      if (p_charge==NULL){
	p_charge = new OBGasteiger;
	if (p_charge==NULL)
	  return false;
	else {
	  m_HaveCreatedCharge = true;
	  SetOBChargeMethod(p_charge);
	}
      }

      p_gaffType->SetTypes(mol);
      p_gaffType->ValidateTypes(p_database);
      p_charge->ComputeCharges(mol);

      OBFunction::Setup(mol);
    }

    void GAFFFunction::Compute(Computation computation)
    {
      if (computation == OBFunction::Gradients)
	for (unsigned int idx = 0; idx < m_gradients.size(); ++idx)
	  m_gradients[idx] = Eigen::Vector3d::Zero();

      std::vector<OBFunctionTerm*>::iterator term;
      for (term = m_terms.begin(); term != m_terms.end(); ++term)
	(*term)->Compute(computation);
    }
  
    double GAFFFunction::GetValue() const
    {
      double energy = 0.0;
   
      std::vector<OBFunctionTerm*> terms = GetTerms();
      std::vector<OBFunctionTerm*>::iterator term;
      for (term = terms.begin(); term != terms.end(); ++term)
	energy += (*term)->GetValue();

      return energy;
    }
      
    std::string GAFFFunction::GetDefaultOptions() const
    {
      std::stringstream ss;
      ss << "# parameters = gaff" << std::endl;
      ss << std::endl;
      ss << "################" << std::endl;
      ss << "# Bonded Terms #" << std::endl;
      ss << "################" << std::endl;
      ss << std::endl;
      ss << "# By default, all bonded terms are enabled." << std::endl;
      ss << "# bonded = [bond] [angle] [torsion] [oop] | [none]" << std::endl;
      ss << "bonded = bond angle torsion oop" << std::endl;
      ss << std::endl;
      ss << "######################" << std::endl;
      ss << "# Van der Waals Term #" << std::endl;
      ss << "######################" << std::endl;
      ss << std::endl;
      ss << "# vdwterm = allpair | none" << std::endl;
      ss << "vdwterm = allpair" << std::endl;
      ss << std::endl;
      return ss.str();
    }
     
    void GAFFFunction::ProcessOptions(std::vector<Option> &options)
    {
      enum BondedTerm {
	BondedBond    = (1<<0),
	BondedAngle   = (1<<1),
	BondedTorsion = (1<<2),
	BondedOOP     = (1<<3),
      };
      int bondedterm = 0;
      bool isBondFound = false;

      enum VdWTerm {
	VdWNone,
	VdWAllPair
      };
      int vdwterm = VdWAllPair;

      enum ElectroTerm {
	ElectroNone,
	ElectroAllPair
      };
      int electroterm = ElectroAllPair;

      OBLogFile *logFile = GetLogFile();
      logFile->Write("Processing GAFF options...\n");
 
      std::vector<Option>::iterator option;
      for (option = options.begin(); option != options.end(); ++option) {

	if ((*option).name == "bonded") {
	  isBondFound=true;
	  if ((*option).value.find("non") != std::string::npos) {
	    logFile->Write("  Disabling all bonded interactions ...\n");
	    bondedterm = 0;        
	  }
	  if ((*option).value.find("all") != std::string::npos) {
	    bondedterm = BondedBond | BondedAngle | BondedTorsion | BondedOOP;
	  }
	  if ((*option).value.find("bond") != std::string::npos) {
	    bondedterm |= BondedBond;        
	  }
	  if ((*option).value.find("angle") != std::string::npos) {
	    bondedterm |= BondedAngle; 
	  }
	  if ((*option).value.find("torsion") != std::string::npos) {
	    bondedterm |= BondedTorsion; 
	  }
	  if ((*option).value.find("oop") != std::string::npos) {
	    bondedterm |= BondedOOP; 
	  }
 	}
	
	if ((*option).name == "vdwterm") {
	  if ((*option).value == "allpair") {
	    vdwterm = VdWAllPair;
	  } else if ((*option).value == "none") {
	    vdwterm = VdWNone;
	  } else {
	    std::stringstream ss;
	    ss << "Invalid value for option: " << (*option).name << " = " << (*option).value << std::endl;
	    logFile->Write(ss.str());
	  }
	}

	if ((*option).name == "electroterm") {
	  if ((*option).value == "allpair")
	    electroterm = ElectroAllPair;
	  if ((*option).value == "none") {
	    electroterm = ElectroNone;
	  } else {
	    std::stringstream ss;
	    ss << "Invalid value for option: " << (*option).name << " = " << (*option).value << std::endl;
	    logFile->Write(ss.str());
	  }
	}
      }
      // use default if option for bonded interaction is not supplied
      isBondFound ? : bondedterm = BondedBond | BondedAngle | BondedTorsion | BondedOOP;


      // remove previous terms
      RemoveAllTerms();
      // add new bonded terms
      if (bondedterm & BondedBond){
	AddTerm(new BondHarmonic(this));
	logFile->Write("  Enabling bond stretching term...\n");
      }
      if (bondedterm & BondedAngle){
	AddTerm(new AngleHarmonic(this));
	logFile->Write("  Enabling angle bending term...\n");
      }
      if (bondedterm & BondedTorsion) {
	AddTerm(new TorsionHarmonic(this));
	logFile->Write("  Enabling torsion term...\n");
      }
      if (bondedterm & BondedOOP) {
	AddTerm(new TorsionHarmonic(this,"Torsion Harmonic OOP"));
	logFile->Write("  Enabling out of plane term...\n");
      }
      // van der waals term
      switch (vdwterm) {
      case VdWNone:
	logFile->Write("  Disabling Van der Waals term\n");
	break;
      case VdWAllPair:
      default:
	AddTerm(new LJ6_12(this, 0.5, LJ6_12::geometric));
	logFile->Write("  Using all-pairs Van der Waals term\n");
	break;
      }
      // electrostatic term
      switch (electroterm) {
      case ElectroNone:
	logFile->Write("  Disabling Van der electrostatic term\n");
	break;
      case ElectroAllPair:
      default:
	logFile->Write("  Using all-pairs electrostatic term\n");
	AddTerm(new Coulomb(this, 0.8333));
	break;
      }
    }
 
    class GAFFFunctionFactory : public OBFunctionFactory
    {
    public:
      GAFFFunctionFactory() : OBFunctionFactory()
      {
      }
      std::string GetName() const
      {
        return "GAFF";
      }
      virtual OBFunction* NewInstance()
      {
        return new GAFFFunction;
      }
    };

    GAFFFunctionFactory theGAFFFunctionFactory;

  } // OBFFs
} // OpenBabel

