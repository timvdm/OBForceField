/*********************************************************************
forcefieldmmff94.cpp - MMFF94 force field

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

#include "mmffparameter.h"

#include <openbabel/locale.h>
#include <openbabel/tokenst.h>
#include <openbabel/oberror.h>
#include <openbabel/mol.h>
#include <fstream>


using namespace std;

namespace OpenBabel {
namespace OBFFs {

  MMFF94ParameterDB::MMFF94ParameterDB(const std::string &filename) : 
    OBFFParameterDB("MMFF94"), m_initialized(false), m_filename(filename)
  {
    ParseParamFile();
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Parse parameter files
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
 
  bool MMFF94ParameterDB::ParseParamFile()
  {
    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    vector<string> vs;
    char buffer[80];
   
    // open data/_parFile
    ifstream ifs;
    if (OpenDatafile(ifs, m_filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open parameter file", obError);
      cout << "Cannot open parameter file: " << m_filename << endl;
      return false;
    }
 
    std::string::size_type p = m_filename.rfind('/');
    std::string path = m_filename.substr(0, p+1);
        
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "#", 1)) continue;
      
      tokenize(vs, buffer);
      if (vs.size() < 2)
        continue;

      if (vs[0] == "prop")
        ParseParamProp(path + vs[1]);
      if (vs[0] == "def")
        ParseParamDef(path + vs[1]);
      if (vs[0] == "bond")
        ParseParamBond(path + vs[1]);
      if (vs[0] == "ang")
        ParseParamAngle(path + vs[1]);
      if (vs[0] == "bndk")
        ParseParamBndk(path + vs[1]);
      if (vs[0] == "chg")
        ParseParamCharge(path + vs[1]);
      if (vs[0] == "dfsb")
        ParseParamDfsb(path + vs[1]);
      if (vs[0] == "oop")
        ParseParamOOP(path + vs[1]);
      if (vs[0] == "pbci")
        ParseParamPbci(path + vs[1]);
      if (vs[0] == "stbn")
        ParseParamStrBnd(path + vs[1]);
      if (vs[0] == "tor")
        ParseParamTorsion(path + vs[1]);
      if (vs[0] == "vdw")
        ParseParamVDW(path + vs[1]);
    }
	
    if (ifs)
      ifs.close();
  
    // return the locale to the original one
    obLocale.RestoreLocale();
    return true;
  }
  
  bool MMFF94ParameterDB::ParseParamBond(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
   
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("name");
    header.push_back("class");
    header.push_back("type1");
    header.push_back("type2");
    header.push_back("kb");
    header.push_back("r0");

    OBFFTable *table = AddTable("Bond Parameters", header);
    
    // open data/mmffbond.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffbond.par", obError);
      cout << "Cannot open parameter file: " << filename << endl;
      return false;
    }

    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      std::stringstream ss;
      ss << vs[0] << ":" << vs[1] << "-" << vs[2];
      cout << ss.str() << endl;
      parameter.push_back(OBVariant(ss.str(), "name"));
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "class"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type2"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "kb"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "r0"));
      table->AddRow(parameter);
    }
	
    if (ifs)
      ifs.close();
  
    return 0;
  }

  
  bool MMFF94ParameterDB::ParseParamBndk(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("species1");
    header.push_back("species2");
    header.push_back("r0-ref");
    header.push_back("kb-ref");
    
    OBFFTable *table = AddTable("Empirical Bond Parameters", header);
 
    // open data/mmffbndk.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffbndk.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "species1"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "species2"));
      parameter.push_back(OBVariant(strtod(vs[2].c_str(), 0), "r0-ref"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "kb-ref"));
      table->AddRow(parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94ParameterDB::ParseParamAngle(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("class");
    header.push_back("type1");
    header.push_back("type2");
    header.push_back("type3");
    header.push_back("ka");
    header.push_back("theta0");
    
    OBFFTable *table = AddTable("Angle Parameters", header);
 
    // open data/mmffang.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffang.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "class"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type2"));
      parameter.push_back(OBVariant(atoi(vs[3].c_str()), "type3"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "ka"));
      parameter.push_back(OBVariant(strtod(vs[5].c_str(), 0), "theta0"));
      table->AddRow(parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
   
  bool MMFF94ParameterDB::ParseParamStrBnd(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("class");
    header.push_back("type1");
    header.push_back("type2");
    header.push_back("type3");
    header.push_back("kbaIJK");
    header.push_back("kbaKJI");
    
    OBFFTable *table = AddTable("Stretch-Bend Parameters", header);
 
    // open data/mmffstbn.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffstbn.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "class"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type2"));
      parameter.push_back(OBVariant(atoi(vs[3].c_str()), "type3"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "kbaIJK"));
      parameter.push_back(OBVariant(strtod(vs[5].c_str(), 0), "kbaKJI"));
      table->AddRow(parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94ParameterDB::ParseParamDfsb(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("periodic-table-row1");
    header.push_back("periodic-table-row2");
    header.push_back("periodic-table-row3");
    header.push_back("Fijk");
    header.push_back("Fkji");
    
    OBFFTable *table = AddTable("Empirical Stretch-Bend Parameters", header);
 
    // open data/mmffdfsb.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffdfsb.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "periodic-table-row1"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "periodic-table-row2"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "periodic-table-row3"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "Fijk"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "Fkji"));
      table->AddRow(parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94ParameterDB::ParseParamOOP(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("type1");
    header.push_back("type2");
    header.push_back("type3");
    header.push_back("type4");
    header.push_back("koop");
    
    OBFFTable *table = AddTable("Out-Of-Plane Parameters", header);
 
    // open data/mmffoop.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffoop.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type2"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type3"));
      parameter.push_back(OBVariant(atoi(vs[3].c_str()), "type4"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "koop"));
      table->AddRow(parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94ParameterDB::ParseParamTorsion(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("class");
    header.push_back("type1");
    header.push_back("type2");
    header.push_back("type3");
    header.push_back("type4");
    header.push_back("V1");
    header.push_back("V2");
    header.push_back("V3");
    
    OBFFTable *table = AddTable("Torsion Parameters", header);
 
    // open data/mmfftor.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmfftor.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "class"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type2"));
      parameter.push_back(OBVariant(atoi(vs[3].c_str()), "type3"));
      parameter.push_back(OBVariant(atoi(vs[4].c_str()), "type4"));
      parameter.push_back(OBVariant(strtod(vs[5].c_str(), 0), "V1"));
      parameter.push_back(OBVariant(strtod(vs[6].c_str(), 0), "V2"));
      parameter.push_back(OBVariant(strtod(vs[7].c_str(), 0), "V3"));
      table->AddRow(parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
   
  bool MMFF94ParameterDB::ParseParamVDW(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("type");
    header.push_back("alpha-i");
    header.push_back("N-i");
    header.push_back("A-i");
    header.push_back("G-i");
    header.push_back("HBD/HBA");
    
    OBFFTable *table = AddTable("Van der Waals Parameters", header);
 
    // open data/mmffvdw.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffvdw.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "type"));
      parameter.push_back(OBVariant(strtod(vs[1].c_str(), 0), "alpha-i"));
      parameter.push_back(OBVariant(strtod(vs[2].c_str(), 0), "N-i"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "A-i"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "G-i"));
      int HBDorHBA;
      if (EQn(vs[5].c_str(), "-", 1))
        HBDorHBA = 0;
      else if (EQn(vs[5].c_str(), "D", 1))
        HBDorHBA = 1; // hydrogen bond donor
      else if (EQn(vs[5].c_str(), "A", 1))
        HBDorHBA = 2; // hydrogen bond acceptor
      parameter.push_back(OBVariant(HBDorHBA, "HBD/HBA"));
      table->AddRow(parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
   
  bool MMFF94ParameterDB::ParseParamCharge(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("class");
    header.push_back("type1");
    header.push_back("type2");
    header.push_back("bci");
           
    OBFFTable *table = AddTable("Charge Parameters", header);
    
    // open data/mmffchg.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffchg.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "class"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type2"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "bci"));
      table->AddRow(parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
 
  bool MMFF94ParameterDB::ParseParamPbci(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("type");
    header.push_back("pbci");
    header.push_back("fcadj");
 
    OBFFTable *table = AddTable("Partial Bond Charge Increments", header);
    
    // open data/mmffpbci.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffpbci", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type"));
      parameter.push_back(OBVariant(strtod(vs[2].c_str(), 0), "pbci"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "fcadj"));
      table->AddRow(parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
 
  bool MMFF94ParameterDB::ParseParamProp(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("atype");
    header.push_back("aspec");
    header.push_back("crd");
    header.push_back("val");
    header.push_back("pilp");
    header.push_back("mltb");
    header.push_back("arom");
    header.push_back("lin");
    header.push_back("sbmb");
         
    OBFFTable *table = AddTable("Atom Properties", header);
    
    // open data/mmffprop.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffprop.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);

      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "type"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "at.no"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "crd"));
      parameter.push_back(OBVariant(atoi(vs[3].c_str()), "val"));
      parameter.push_back(OBVariant(vs[4] == "1", "pilp"));
      parameter.push_back(OBVariant(atoi(vs[5].c_str()), "mltb"));
      parameter.push_back(OBVariant(vs[6] == "1", "arom"));
      parameter.push_back(OBVariant(vs[7] == "1", "linh"));
      parameter.push_back(OBVariant(vs[8] == "1", "sbmb"));
      table->AddRow(parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94ParameterDB::ParseParamDef(const std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    std::vector<std::string> header;
    header.push_back("level1");
    header.push_back("level2");
    header.push_back("level3");
    header.push_back("level4");
    header.push_back("level5");
         
    OBFFTable *table = AddTable("Atom Type Levels", header);

    // open data/mmffdef.par
    ifstream ifs;
    if (OpenDatafile(ifs, filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open mmffdef.par", obError);
      return false;
    }
    
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "*", 1)) continue;
      if (EQn(buffer, "$", 1)) continue;
	
      tokenize(vs, buffer);
      
      parameter.clear();
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "level1"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "level2"));
      parameter.push_back(OBVariant(atoi(vs[3].c_str()), "level3"));
      parameter.push_back(OBVariant(atoi(vs[4].c_str()), "level4"));
      parameter.push_back(OBVariant(atoi(vs[5].c_str()), "level5"));
      table->AddRow(parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
} // end namespace OpenBabel
}
//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field

