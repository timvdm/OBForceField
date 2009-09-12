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
  
  void OBFunction::RemoveAllTerms(bool deleteTerms)
  {
    std::vector<OBFunctionTerm*>::iterator term;
    for (term = m_terms.begin(); term != m_terms.end(); ++term)
      delete *term;
    m_terms.clear();
  }
      
  const std::vector<OBFunctionTerm*>& OBFunction::GetTerms() const
  {
    return m_terms;
  }

  //  
  //         f(1) - f(0)
  // f'(0) = -----------      f(1) = f(0+h)
  //              h
  //
  Eigen::Vector3d OBFunction::NumericalDerivative(unsigned int index)
  {
    double e_orig, e_plus_delta, delta, dx, dy, dz;
    delta = 1.0e-5;

    const Eigen::Vector3d va = m_positions.at(index);
    Compute(OBFunction::Value);
    e_orig = GetValue();
    
    // X direction
    m_positions[index].x() += delta;
    Compute(OBFunction::Value);
    e_plus_delta = GetValue();
    dx = (e_plus_delta - e_orig) / delta;
    
    // Y direction
    m_positions[index].x() = va.x();
    m_positions[index].y() += delta;
    Compute(OBFunction::Value);
    e_plus_delta = GetValue();
    dy = (e_plus_delta - e_orig) / delta;
    
    // Z direction
    m_positions[index].y() = va.y();
    m_positions[index].z() += delta;
    Compute(OBFunction::Value);
    e_plus_delta = GetValue();
    dz = (e_plus_delta - e_orig) / delta;

    // reset coordinates to original
    m_positions[index].z() = va.z();

    return Eigen::Vector3d(-dx, -dy, -dz);
  }

  //  
  //         f(2) - 2f(1) + f(0)
  // f'(0) = -------------------      f(1) = f(0+h)
  //                 h^2              f(1) = f(0+2h)
  //
  Eigen::Vector3d OBFunction::NumericalSecondDerivative(unsigned int index)
  {
    double e_0, e_1, e_2, delta, dx, dy, dz;
    delta = 1.0e-5;

    const Eigen::Vector3d va = m_positions.at(index);

    // calculate f(0)
    Compute(OBFunction::Value);
    e_0 = GetValue();
    
    // 
    // X direction
    //
    
    // calculate f(1)
    m_positions[index].x() += delta;
    Compute(OBFunction::Value);
    e_1 = GetValue();

    // calculate f(2)
    m_positions[index].x() += delta;
    Compute(OBFunction::Value);
    e_2 = GetValue();
    
    dx = (e_2 - 2 * e_1 + e_0) / (delta * delta);
    m_positions[index].x() = va.x();
    
    // 
    // Y direction
    //
    
    // calculate f(1)
    m_positions[index].y() += delta;
    Compute(OBFunction::Value);
    e_1 = GetValue();

    // calculate f(2)
    m_positions[index].y() += delta;
    Compute(OBFunction::Value);
    e_2 = GetValue();

    dy = (e_2 - 2 * e_1 + e_0) / (delta * delta);
    m_positions[index].y() = va.y();

    // 
    // Z direction
    //
    
    // calculate f(1)
    m_positions[index].z() += delta;
    Compute(OBFunction::Value);
    e_1 = GetValue();

    // calculate f(2)
    m_positions[index].z() += delta;
    Compute(OBFunction::Value);
    e_2 = GetValue();

    dz = (e_2 - 2 * e_1 + e_0) / (delta * delta);
    m_positions[index].z() = va.z();


    return Eigen::Vector3d(-dx, -dy, -dz);
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
