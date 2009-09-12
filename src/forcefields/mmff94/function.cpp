/**********************************************************************
obfunction.h - Base class for force fields and scroring functions which
               depend on atom positions.
 
Copyright (C) 2008,2009 by Tim Vandermeersch
 
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

#include <OBForceField>
#include <OBLogFile>
#include <OBParameterDB>

#include <openbabel/mol.h>

#include "common.h"
#include "parameter.h"

// terms
#include "bond.h"
#include "angle.h"
#include "strbnd.h"
#include "torsion.h"
#include "oop.h"
#include "vdw.h"
#include "electro.h"
#include "vdw_opencl.h"
#include "electro_opencl.h"


namespace OpenBabel {
namespace OBFFs {

  class MMFF94Function : public OBForceField 
  {
    public:
      MMFF94Function();
      //virtual ~MMFF94Function();
      
      std::string GetName() const
      {
        return "MMFF94";
      }
      
      bool Setup(/*const*/ OBMol &mol);
      void Compute(Computation computation = Value);
      double GetValue() const;
      
      std::string GetUnit() const { return "kcal/mol"; }
      
      std::vector<std::string> GetAtomTypes() const;
      std::vector<double> GetFormalCharges() const;
      std::vector<double> GetPartialCharges() const;
 
      bool HasAnalyticalGradients() const 
      { 
        return true; 
      }
    protected:
      void ProcessOptions(std::vector<Option> &options);
      std::string GetDefaultOptions() const;

      MMFF94Common *m_common;
  };

  MMFF94Function::MMFF94Function() : m_common(new MMFF94Common)
  {
    SetParameterDB(m_common->GetParameterDB());
    AddTerm(new MMFF94BondTerm(this, m_common));
    AddTerm(new MMFF94AngleTerm(this, m_common));
    AddTerm(new MMFF94StrBndTerm(this, m_common));
    AddTerm(new MMFF94TorsionTerm(this, m_common));
    AddTerm(new MMFF94OutOfPlaneTerm(this, m_common));
    AddTerm(new MMFF94VDWTerm(this, m_common));
    AddTerm(new MMFF94ElectroTerm(this, m_common));
  }

  bool MMFF94Function::Setup(/*const*/ OBMol &mol)
  {
    if (!m_common->SetTypes(mol))
      return false;

    PrintAtomTypes();

    m_common->SetFormalCharges(mol);
    PrintFormalCharges();
    m_common->SetPartialCharges(mol);
    PrintPartialCharges();

    // call setup for all terms
    OBFunction::Setup(mol);
  }

  void MMFF94Function::Compute(Computation computation)
  {
    if (computation == OBFunction::Gradients)
      for (unsigned int idx = 0; idx < m_gradients.size(); ++idx)
        m_gradients[idx] = Eigen::Vector3d::Zero();

    std::vector<OBFunctionTerm*>::iterator term;
    for (term = m_terms.begin(); term != m_terms.end(); ++term)
      (*term)->Compute(computation);
  }
  
  double MMFF94Function::GetValue() const
  {
    double energy = 0.0;
   
    std::vector<OBFunctionTerm*> terms = GetTerms();
    std::vector<OBFunctionTerm*>::iterator term;
    for (term = terms.begin(); term != terms.end(); ++term)
      energy += (*term)->GetValue();

    return energy;
  }
      
  std::string MMFF94Function::GetDefaultOptions() const
  {
    std::stringstream ss;
    ss << "# parameters = mmff94 | mmff94s" << std::endl;
    ss << "parameters = mmff94" << std::endl;
    ss << std::endl;
    ss << "################" << std::endl;
    ss << "# Bonded Terms #" << std::endl;
    ss << "################" << std::endl;
    ss << std::endl;
    ss << "# By default, all bonded terms are enabled." << std::endl;
    ss << "# bonded = [bond] [angle] [strbnd] [torsion] [oop] | none" << std::endl;
    ss << "bonded = bond angle strbnd torsion oop" << std::endl;
    ss << std::endl;
    ss << "######################" << std::endl;
    ss << "# Van der Waals Term #" << std::endl;
    ss << "######################" << std::endl;
    ss << std::endl;
    ss << "# vdwterm = allpair | rvdw | opencl | none" << std::endl;
    ss << "vdwterm = allpair" << std::endl;
    ss << std::endl;
    ss << "# rvdw = <double>" << std::endl;
    ss << "rvdw = 8.0" << std::endl;
    ss << std::endl;
    return ss.str();
  }
     
  void MMFF94Function::ProcessOptions(std::vector<Option> &options)
  {
    enum BondedTerm {
      BondedBond    = (1<<0),
      BondedAngle   = (1<<1),
      BondedStrBnd  = (1<<2),
      BondedTorsion = (1<<3),
      BondedOOP     = (1<<4),
    };
    int bondedterm = BondedBond | BondedAngle | BondedStrBnd | BondedTorsion | BondedOOP;

    enum VdwTerm {
      VdwAllPair,
      VdwCutOff,
      VdwOpenCL,
      VdwNone
    };
    int vdwterm = VdwAllPair;

    enum ElectroTerm {
      ElectroAllPair,
      ElectroCutOff,
      ElectroOpenCL,
      ElectroNone
    };
    int electroterm = ElectroAllPair;

    OBLogFile *logFile = GetLogFile();
    logFile->Write("Processing MMFF94 options...\n");
 
    std::vector<Option>::iterator option;
    for (option = options.begin(); option != options.end(); ++option) {

      if ((*option).name == "bonded") {
        bondedterm = 0;
        if ((*option).value.find("bond") != std::string::npos) {
          logFile->Write("  Enabling bond stretching term...\n");
          bondedterm |= BondedBond;        
        }
        if ((*option).value.find("angle") != std::string::npos) {
          logFile->Write("  Enabling angle bending term...\n");
          bondedterm |= BondedAngle; 
        }
        if ((*option).value.find("strbnd") != std::string::npos) {
          logFile->Write("  Enabling stretch bending term...\n");
          bondedterm |= BondedStrBnd; 
        }
        if ((*option).value.find("torsion") != std::string::npos) {
          logFile->Write("  Enabling torsion term...\n");
          bondedterm |= BondedTorsion; 
        }
        if ((*option).value.find("oop") != std::string::npos) {
          logFile->Write("  Enabling out of plane term...\n");
          bondedterm |= BondedOOP; 
        }
 
      }

      if ((*option).name == "vdwterm") {
        if ((*option).value == "allpair") {
          vdwterm = VdwAllPair;
          logFile->Write("  Using all-pairs Van der Waals term\n");
        } else if ((*option).value == "rvdw") {
          vdwterm = VdwCutOff;
          logFile->Write("  Using cut-off Van der Waals term\n");
        } else if ((*option).value == "opencl") {
          vdwterm = VdwOpenCL;
          logFile->Write("  Using OpenCL Van der Waals term\n");
        } else if ((*option).value == "none") {
          vdwterm = VdwNone;
          logFile->Write("  Van der Waals term set to none\n");
        } else {
          std::stringstream ss;
          ss << "Invalid value for option: " << (*option).name << " = " << (*option).value << std::endl;
          logFile->Write(ss.str());
        }
      }

      if ((*option).name == "electroterm") {
        if ((*option).value == "allpair")
          electroterm = ElectroAllPair;
        if ((*option).value == "rele")
          electroterm = ElectroCutOff;
        if ((*option).value == "opencl")
          electroterm = ElectroOpenCL;
        if ((*option).value == "none")
          electroterm = ElectroNone;        
      }
 
    }

    // remove previous terms
    RemoveAllTerms();
    // add new bonded terms
    if (bondedterm & BondedBond)
      AddTerm(new MMFF94BondTerm(this, m_common));
    if (bondedterm & BondedAngle)
      AddTerm(new MMFF94AngleTerm(this, m_common));
    if (bondedterm & BondedStrBnd)
      AddTerm(new MMFF94StrBndTerm(this, m_common));
    if (bondedterm & BondedTorsion)
      AddTerm(new MMFF94TorsionTerm(this, m_common));
    if (bondedterm & BondedOOP)
      AddTerm(new MMFF94OutOfPlaneTerm(this, m_common));
    // van der waals term
    switch (vdwterm) {
      case VdwNone:
        break;
      case VdwOpenCL:
        AddTerm(new MMFF94VDWTermOpenCL(this, m_common));
        break;
      case VdwAllPair:
      default:
        AddTerm(new MMFF94VDWTerm(this, m_common));
        break;
    }
    // electrostatic term
    switch (electroterm) {
      case VdwNone:
        break;
      case ElectroOpenCL:
        AddTerm(new MMFF94ElectroTermOpenCL(this, m_common));
        break;
      case VdwAllPair:
      default:
        AddTerm(new MMFF94ElectroTerm(this, m_common));
        break;
    }
  }
 
  std::vector<std::string> MMFF94Function::GetAtomTypes() const
  {
    std::vector<std::string> types;
    for (unsigned int i = 0; i < m_common->m_types.size(); ++i) {
      std::stringstream ss;
      ss << m_common->m_types.at(i);
      types.push_back(ss.str());
    }
    return types;
  }
   
  std::vector<double> MMFF94Function::GetFormalCharges() const
  {
    return m_common->m_fCharges;
  }
  
  std::vector<double> MMFF94Function::GetPartialCharges() const
  {
    return m_common->m_pCharges;
  }

  class MMFF94FunctionFactory : public OBFunctionFactory
  {
    public:
      MMFF94FunctionFactory() : OBFunctionFactory()
      {
      }
      std::string GetName() const
      {
        return "MMFF94";
      }
      virtual OBFunction* NewInstance()
      {
        return new MMFF94Function;
      }
  };

  MMFF94FunctionFactory theMMFF94FunctionFactory;

} // OBFFs
} // OpenBabel

