/**********************************************************************
obfunction.cpp - Base class for force fields and scroring functions which
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

#include <OBFunction>
#include <OBFunctionTerm>
#include <OBLogFile>

#include <openbabel/mol.h>
#include <iostream>
#include <iterator>
using namespace std;

namespace OpenBabel {
namespace OBFFs {

  OBFunction::OBFunction() : m_logfile(new OBLogFile), m_parameterDB(0)
  {
  }

  OBFunction::~OBFunction()
  {
    std::vector<OBFunctionTerm*>::iterator term;
    for (term = m_terms.begin(); term != m_terms.end(); ++term)
      delete &*term;
    delete m_logfile;
  }

  bool OBFunction::Setup(/*const*/ OBMol &mol)
  {
    m_positions.resize(mol.NumAtoms());
    FOR_ATOMS_OF_MOL (atom, mol)
      m_positions[atom->GetIdx()-1] = Eigen::Vector3d(atom->GetVector().AsArray());

    m_gradients.resize(mol.NumAtoms(), Eigen::Vector3d::Zero());

    std::vector<OBFunctionTerm*>::iterator term;
    for (term = m_terms.begin(); term != m_terms.end(); ++term)
      (*term)->Setup(mol);
  }

  void OBFunction::SetParameterDB(OBParameterDB *db) 
  { 
    m_parameterDB = db; 
  }

  void OBFunction::AddTerm(OBFunctionTerm *term)
  {
    if (term)
      m_terms.push_back(term);
  }
      
  const std::vector<OBFunctionTerm*>& OBFunction::GetTerms() const
  {
    return m_terms;
  }

  std::string OBFunction::GetOptions() const
  {
    if (m_options.empty())
      return GetDefaultOptions();
    return m_options;  
  }
  
  void OBFunction::SetOptions(const std::string &options)
  {
    int lineCount = 0;
    std::vector<Option> voptions;
    std::string line;
    std::stringstream ss(options);
    while (std::getline(ss, line, '\n')) {
      lineCount++;
      std::istringstream iss(line);
      std::vector<std::string> tokens;
      copy(std::istream_iterator<string>(iss), 
         std::istream_iterator<std::string>(), 
         back_inserter<std::vector<std::string> >(tokens));

      if (tokens.size() >= 1) {
        if (tokens.at(0)[0] == '#')
          continue;
      }

      if (tokens.size() < 3) {
        std::stringstream msg;
        msg << "line " << lineCount << ": All options should have at least 3 tokens (e.g. 'rvdw = 9.0', ...)" << endl;
        m_logfile->Write(msg.str());
        continue;
      }

      if (tokens.at(1) != "=") {
        std::stringstream msg;
        msg << "line " << lineCount << ": No '=' character found. All options should be of the form 'name = value' (e.g. 'rvdw = 9.0', ...)" << endl;
        m_logfile->Write(msg.str());
        continue;
      }

      voptions.push_back( Option(lineCount, tokens.at(0), tokens.at(2)) );
    }

    ProcessOptions(voptions);

    m_options = options;  
  }
 


  OBFunctionFactory::OBFunctionFactory()
  {
    OBFunctionFactory::GetFactoriesRef().push_back(this);
  }
  
  std::vector<OBFunctionFactory*> OBFunctionFactory::GetFactories()
  {
    return OBFunctionFactory::GetFactoriesRef();
  }

  std::vector<OBFunctionFactory*>& OBFunctionFactory::GetFactoriesRef()
  {
    static std::vector<OBFunctionFactory*> factories = std::vector<OBFunctionFactory*>();
    return factories;
  }

  OBFunctionFactory* OBFunctionFactory::GetFactory(const std::string &name)
  {
    std::vector<OBFunctionFactory*> factories(GetFactoriesRef());
    for (unsigned int i = 0; i < factories.size(); ++i)
      if (factories.at(i)->GetName() == name)
        return factories.at(i);
    return 0;
  }

 
}
} // end namespace OpenBabel


//! @file obfunction.cpp
//! @brief Handle OBFunction class
