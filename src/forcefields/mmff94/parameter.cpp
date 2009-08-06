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

#include <OBParameterDB>

#include <openbabel/locale.h>
#include <openbabel/tokenst.h>
#include <openbabel/oberror.h>
#include <openbabel/mol.h>
#include <fstream>

#include "parameter.h"

using namespace std;

namespace OpenBabel {
namespace OBFFs {

  const int numTables = 12;

  template<unsigned int idx> std::string TableToString();
  template<> std::string TableToString<0>() { return "MMFF94 Atom Properties"; }
  template<> std::string TableToString<1>() { return "MMFF94 Atom Type Levels"; }
  template<> std::string TableToString<2>() { return "MMFF94 Bond Stretching Parameters"; }
  template<> std::string TableToString<3>() { return "MMFF94 Angle Bending Parameters"; }
  template<> std::string TableToString<4>() { return "MMFF94 Stretch Bend Parameters"; }
  template<> std::string TableToString<5>() { return "MMFF94 Torsion Parameters"; }
  template<> std::string TableToString<6>() { return "MMFF94 Out-Of-Plane Parameters"; }
  template<> std::string TableToString<7>() { return "MMFF94 Van der Waals Parameters"; }
  template<> std::string TableToString<8>() { return "MMFF94 Charge Parameters"; }
  template<> std::string TableToString<9>() { return "MMFF94 Bond Stretching Emprirical Rules"; }
  template<> std::string TableToString<10>() { return "MMFF94 Stretch Bending Emprirical Rules"; }
  template<> std::string TableToString<11>() { return "MMFF94 Partial Bond Charge Increments"; }

  template<int N> unsigned int StringToTable(const std::string &table)
  {
    if (table == TableToString<N>()) 
      return N;
    return StringToTable<N-1>(table);
  }
  template<> unsigned int StringToTable<0>(const std::string &table)
  {
    if (table == TableToString<0>()) 
      return 0;
    return numTables; // not found, return numTables
  }


      MMFF94SimpleParameterDB::MMFF94SimpleParameterDB(const std::string &filename) : m_initialized(false), m_filename(filename)
      {
        ParseParamFile();
      }

      unsigned int MMFF94SimpleParameterDB::NumTables() const 
      { 
        return numTables; 
      }
      
      std::vector<std::string> MMFF94SimpleParameterDB::GetTables() const
      {
        std::vector<std::string> tables;
        tables.push_back(TableToString<0>());
        tables.push_back(TableToString<1>());
        tables.push_back(TableToString<2>());
        tables.push_back(TableToString<3>());
        tables.push_back(TableToString<4>());
        tables.push_back(TableToString<5>());
        tables.push_back(TableToString<6>());
        tables.push_back(TableToString<7>());
        tables.push_back(TableToString<8>());
        tables.push_back(TableToString<9>());
        tables.push_back(TableToString<10>());
        tables.push_back(TableToString<11>());
        return tables;            
      }

      unsigned int MMFF94SimpleParameterDB::GetTable(const std::string &table) const
      {
        return StringToTable<numTables-1>(table);
      }

      unsigned int MMFF94SimpleParameterDB::NumRows(unsigned int table) const
      {
        if (table >= m_tables.size())
          return 0;
        return m_tables.at(table)->size();      
      }
      
      unsigned int MMFF94SimpleParameterDB::NumColumns(unsigned int table) const
      {
        switch (table) {
          case AtomProperties:
            return 9;
          case AtomTypeLevels:
            return 5;
          case BondParameters:
            return 5;
          case AngleParameters:
            return 6;
          case StretchBendParameters:
            return 6;
          case TorsionParameters:
            return 8;
          case OutOfPlaneParameters:
            return 5;
          case VanDerWaalsParameters:
            return 6;
          case ChargeParameters:
            return 4;
          case EmpiricalBondRules:
            return 4;
          case EmpiricalStretchBendRules: 
            return 5;
          case PartialBondChargeIncrements:
            return 3;
          default: 
            return 0;
        }
      }
      
      std::vector<std::string> MMFF94SimpleParameterDB::GetHeader(unsigned int table) const
      {
        std::vector<std::string> header;
        switch (table) {
          case AtomProperties:
            header.push_back("atype");
            header.push_back("aspec");
            header.push_back("crd");
            header.push_back("val");
            header.push_back("pilp");
            header.push_back("mltb");
            header.push_back("arom");
            header.push_back("lin");
            header.push_back("sbmb");
            break;
          case AtomTypeLevels:
            header.push_back("level1");
            header.push_back("level2");
            header.push_back("level3");
            header.push_back("level4");
            header.push_back("level5");
            break;
          case BondParameters:
            header.push_back("class");
            header.push_back("type1");
            header.push_back("type2");
            header.push_back("kb");
            header.push_back("r0");
            break;
          case AngleParameters:
            header.push_back("class");
            header.push_back("type1");
            header.push_back("type2");
            header.push_back("type3");
            header.push_back("ka");
            header.push_back("theta0");
            break;
          case StretchBendParameters:
            header.push_back("class");
            header.push_back("type1");
            header.push_back("type2");
            header.push_back("type3");
            header.push_back("kbaIJK");
            header.push_back("kbaKJI");
            break;
          case TorsionParameters:
            header.push_back("class");
            header.push_back("type1");
            header.push_back("type2");
            header.push_back("type3");
            header.push_back("type4");
            header.push_back("V1");
            header.push_back("V2");
            header.push_back("V3");
            break;
          case OutOfPlaneParameters:
            header.push_back("type1");
            header.push_back("type2");
            header.push_back("type3");
            header.push_back("type4");
            header.push_back("koop");
            break;
          case VanDerWaalsParameters:
            header.push_back("type");
            header.push_back("alpha-i");
            header.push_back("N-i");
            header.push_back("A-i");
            header.push_back("G-i");
            header.push_back("HBD/HBA");
            break;
          case ChargeParameters:
            header.push_back("class");
            header.push_back("type1");
            header.push_back("type2");
            header.push_back("bci");
            break;
          case EmpiricalBondRules:
            header.push_back("species1");
            header.push_back("species2");
            header.push_back("r0-ref");
            header.push_back("kb-ref");
            break;
          case EmpiricalStretchBendRules: 
            header.push_back("periodic-table-row1");
            header.push_back("periodic-table-row2");
            header.push_back("periodic-table-row3");
            header.push_back("Fijk");
            header.push_back("Fkji");
            break;
          case PartialBondChargeIncrements:
            header.push_back("type");
            header.push_back("pbci");
            header.push_back("fcadj");
            break;
          default: 
            break;
        }

        return header;      
      }

      std::vector<OBVariant::Type> MMFF94SimpleParameterDB::GetTypes(unsigned int table) const
      {
        std::vector<OBVariant::Type> header;
        switch (table) {
          case AtomProperties:
            header.push_back(OBVariant::Int); // atype
            header.push_back(OBVariant::Int); // aspec
            header.push_back(OBVariant::Int); // crd
            header.push_back(OBVariant::Int); // val
            header.push_back(OBVariant::Bool); // pilp
            header.push_back(OBVariant::Int); // mltb
            header.push_back(OBVariant::Bool); // arom
            header.push_back(OBVariant::Bool); // lin
            header.push_back(OBVariant::Bool); // sbmb
            break;
          case AtomTypeLevels:
            header.push_back(OBVariant::Int); // level1
            header.push_back(OBVariant::Int); // level2
            header.push_back(OBVariant::Int); // level3
            header.push_back(OBVariant::Int); // level4
            header.push_back(OBVariant::Int); // level5
            break;
          case BondParameters:
            header.push_back(OBVariant::Int); // class
            header.push_back(OBVariant::Int); // type1
            header.push_back(OBVariant::Int); // type2
            header.push_back(OBVariant::Double); // kb
            header.push_back(OBVariant::Double); // r0
            break;
          case AngleParameters:
            header.push_back(OBVariant::Int); // class
            header.push_back(OBVariant::Int); // type1
            header.push_back(OBVariant::Int); // type2
            header.push_back(OBVariant::Int); // type3
            header.push_back(OBVariant::Double); // ka
            header.push_back(OBVariant::Double); // theta0
            break;
          case StretchBendParameters:
            header.push_back(OBVariant::Int); // class
            header.push_back(OBVariant::Int); // type1
            header.push_back(OBVariant::Int); // type2
            header.push_back(OBVariant::Int); // type3
            header.push_back(OBVariant::Double); // kbaIJK
            header.push_back(OBVariant::Double); // bbaKJI
            break;
          case TorsionParameters:
            header.push_back(OBVariant::Int); // class
            header.push_back(OBVariant::Int); // type1
            header.push_back(OBVariant::Int); // type2
            header.push_back(OBVariant::Int); // type3
            header.push_back(OBVariant::Int); // type4
            header.push_back(OBVariant::Double); // V1
            header.push_back(OBVariant::Double); // V2
            header.push_back(OBVariant::Double); // V3
            break;
          case OutOfPlaneParameters:
            header.push_back(OBVariant::Int); // type1
            header.push_back(OBVariant::Int); // type2
            header.push_back(OBVariant::Int); // type3
            header.push_back(OBVariant::Int); // type4
            header.push_back(OBVariant::Double); // koop
            break;
          case VanDerWaalsParameters:
            header.push_back(OBVariant::Int); // type
            header.push_back(OBVariant::Double); // alpha-i 
            header.push_back(OBVariant::Double); // N-i
            header.push_back(OBVariant::Double); // A-i
            header.push_back(OBVariant::Double); // G-i
            header.push_back(OBVariant::Int); // HBD/HBA
            break;
          case ChargeParameters:
            header.push_back(OBVariant::Int); // class
            header.push_back(OBVariant::Int); // type1
            header.push_back(OBVariant::Int); // type2
            header.push_back(OBVariant::Double); // bci
            break;
          case EmpiricalBondRules:
            header.push_back(OBVariant::Int); // species1
            header.push_back(OBVariant::Int); // species2
            header.push_back(OBVariant::Double); // r0-ref
            header.push_back(OBVariant::Double); // kb-ref
            break;
          case EmpiricalStretchBendRules: 
            header.push_back(OBVariant::Int); // periodic-table-row1
            header.push_back(OBVariant::Int); // periodic-table-row2
            header.push_back(OBVariant::Int); // periodic-table-row3
            header.push_back(OBVariant::Double); // Fijk
            header.push_back(OBVariant::Double); // Fkji
            break;
          case PartialBondChargeIncrements:
            header.push_back(OBVariant::Int); // type1
            header.push_back(OBVariant::Double); // Fijk
            header.push_back(OBVariant::Double); // fcadj
            break;
          default: 
            break;
        }

        return header;      
      }
      
      bool MMFF94SimpleParameterDB::VerifyTypes(unsigned int table, const std::vector<OBVariant> &values) 
      {
        std::vector<OBVariant::Type> types = GetTypes(table);
        if (types.size() != values.size())
          return false;

        for (unsigned int i = 0; i < types.size(); ++i)
          if (types.at(i) != values.at(i).GetType())
            return false;

        return true;
      }

      const std::vector<OBVariant>& MMFF94SimpleParameterDB::GetRow(unsigned int table, unsigned int row) const
      {
        if (table >= m_tables.size())
          return emptyRow;
        if (row >= NumRows(table))
          return emptyRow;
        return m_tables.at(table)->at(row);
      }

      bool MMFF94SimpleParameterDB::AddRow(unsigned int table, const std::vector<OBVariant> &values)
      {
        if (table >= m_tables.size())
          return false;
        if (values.size() != NumColumns(table))
          return false;
        if (!VerifyTypes(table, values))
          return false;

        m_tables.at(table)->push_back(values);      
        return true;
      }
 
      const std::vector<OBVariant>& MMFF94SimpleParameterDB::FindRow(unsigned int table, const std::vector<Query> &query, bool *swapped) const
      {
        unsigned int numRows = NumRows(table);
        unsigned int numColumns = NumColumns(table);

        unsigned int swapCount = 0;
        // make sure the qury contains valid columns
        for (unsigned int i = 0; i < query.size(); ++i) {
          if (query.at(i).column >= numColumns)
            return emptyRow;
          if (query.at(i).swap)
            swapCount++;
        }
      
        if (!swapCount) {
          for (unsigned int i = 0; i < numRows; ++i) {
            bool match = true;
            // check each query entry
            for (unsigned int j = 0; j < query.size(); ++j) {
              if (m_tables.at(table)->at(i).at(query.at(j).column) != query.at(j).value) {
                match = false;
                break;
              }
            }

            if (match)
              return m_tables.at(table)->at(i);
          }
        } else {
          // construct swapped_query
          std::vector<Query> swapped_query = query;
          for (unsigned int i = 0; i < query.size(); ++i) {
            if (query.at(i).swap) {
              if (swapCount == 4) {
                swapped_query[i  ].column = query[i+3].column;
                swapped_query[i+1].column = query[i+2].column;
                swapped_query[i+2].column = query[i+1].column;
                swapped_query[i+3].column = query[i  ].column;
              } else {
                swapped_query[i].column = query[i+swapCount-1].column;
                swapped_query[i+swapCount-1].column = query[i].column;
              }
              break;
            }
          }

          for (unsigned int i = 0; i < numRows; ++i) {
            bool match = true;
            // check query
            for (unsigned int j = 0; j < query.size(); ++j) {
              if (m_tables.at(table)->at(i).at(query.at(j).column) != query.at(j).value) {
                match = false;
                break;
              }
            }

            if (match) {
              if (swapped)
                *swapped = false;
              return m_tables.at(table)->at(i);
            }
 
            match = true;
            // check swapped query
            for (unsigned int j = 0; j < swapped_query.size(); ++j) {
              if (m_tables.at(table)->at(i).at(swapped_query.at(j).column) != swapped_query.at(j).value) {
                match = false;
                break;
              }
            }
 
            if (match) {
              if (swapped)
                *swapped = true;
              return m_tables.at(table)->at(i);
            }
          }
        }
     
        return emptyRow;
      }
      
  const std::vector<std::vector<OBVariant> >& MMFF94SimpleParameterDB::GetAllRows(unsigned int table) const
  {
    return *m_tables.at(table);
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Parse parameter files
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
 
  bool MMFF94SimpleParameterDB::ParseParamFile()
  {
    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    vector<string> vs;
    char buffer[80];
   
    m_tables.clear();
    for (unsigned int i = 0; i < 12; ++i)
      m_tables.push_back(new std::vector< std::vector<OBVariant> >);

    // open data/_parFile
    ifstream ifs;
    if (OpenDatafile(ifs, m_filename).length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open parameter file", obError);
      cout << "Cannot open parameter file: " << m_filename << endl;
      return false;
    }
        
    while (ifs.getline(buffer, 80)) {
      if (EQn(buffer, "#", 1)) continue;
      
      tokenize(vs, buffer);
      if (vs.size() < 2)
        continue;

      if (vs[0] == "prop")
        ParseParamProp(vs[1]);
      if (vs[0] == "def")
        ParseParamDef(vs[1]);
      if (vs[0] == "bond")
        ParseParamBond(vs[1]);
      if (vs[0] == "ang")
        ParseParamAngle(vs[1]);
      if (vs[0] == "bndk")
        ParseParamBndk(vs[1]);
      if (vs[0] == "chg")
        ParseParamCharge(vs[1]);
      if (vs[0] == "dfsb")
        ParseParamDfsb(vs[1]);
      if (vs[0] == "oop")
        ParseParamOOP(vs[1]);
      if (vs[0] == "pbci")
        ParseParamPbci(vs[1]);
      if (vs[0] == "stbn")
        ParseParamStrBnd(vs[1]);
      if (vs[0] == "tor")
        ParseParamTorsion(vs[1]);
      if (vs[0] == "vdw")
        ParseParamVDW(vs[1]);
    }
	
    if (ifs)
      ifs.close();
  
    // return the locale to the original one
    obLocale.RestoreLocale();
    return true;
  }
  
  bool MMFF94SimpleParameterDB::ParseParamBond(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
   
    std::vector<OBVariant> parameter;
    
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
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "class"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type2"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "kb"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "r0"));
      AddRow(BondParameters, parameter);
    }
	
    if (ifs)
      ifs.close();
  
    return 0;
  }

  
  bool MMFF94SimpleParameterDB::ParseParamBndk(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type2"));
      parameter.push_back(OBVariant(strtod(vs[2].c_str(), 0), "r0-ref"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "kb-ref"));
      AddRow(EmpiricalBondRules, parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94SimpleParameterDB::ParseParamAngle(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(AngleParameters, parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
   
  bool MMFF94SimpleParameterDB::ParseParamStrBnd(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(StretchBendParameters, parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94SimpleParameterDB::ParseParamDfsb(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      parameter.push_back(OBVariant(atoi(vs[0].c_str()), "type1"));
      parameter.push_back(OBVariant(atoi(vs[1].c_str()), "type2"));
      parameter.push_back(OBVariant(atoi(vs[2].c_str()), "type3"));
      parameter.push_back(OBVariant(strtod(vs[3].c_str(), 0), "kbaIJK"));
      parameter.push_back(OBVariant(strtod(vs[4].c_str(), 0), "kbaKJI"));
      AddRow(EmpiricalStretchBendRules, parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94SimpleParameterDB::ParseParamOOP(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(OutOfPlaneParameters, parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94SimpleParameterDB::ParseParamTorsion(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(TorsionParameters, parameter);
    }
	
    if (ifs)
      ifs.close();
 
    return 0;
  }
   
  bool MMFF94SimpleParameterDB::ParseParamVDW(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(VanDerWaalsParameters, parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
   
  bool MMFF94SimpleParameterDB::ParseParamCharge(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(ChargeParameters, parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
 
  bool MMFF94SimpleParameterDB::ParseParamPbci(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(PartialBondChargeIncrements, parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
 
  bool MMFF94SimpleParameterDB::ParseParamProp(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(AtomProperties, parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
  
  bool MMFF94SimpleParameterDB::ParseParamDef(std::string &filename)
  {
    vector<string> vs;
    char buffer[80];
    
    std::vector<OBVariant> parameter;
    
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
      AddRow(AtomTypeLevels, parameter);
    }

    if (ifs)
      ifs.close();
 
    return 0;
  }
} // end namespace OpenBabel
}
//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field

