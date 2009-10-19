/*********************************************************************
General Forcefield database class

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

namespace OpenBabel {
  namespace OBFFs {

    class OBFFParameterDB;

    class OBFFTable : public OBParameterDBTable
    {
    public:
      OBFFTable(const std::string &tableName, const std::vector<std::string> &header);
      OBFFTable(const std::string &tableName);
      unsigned int NumRows() const;
      unsigned int NumColumns() const;
      std::vector<std::string> GetHeader() const;
      std::vector<OBVariant::Type> GetTypes() const;
      bool VerifyTypes(const std::vector<OBVariant> &values);
      const std::vector<OBVariant>& GetRow(unsigned int row) const;
      bool AddRow(const std::vector<OBVariant> &values);
      const std::vector<OBVariant>& FindRow(const std::vector<Query> &query, bool *swapped = 0) const;
      std::vector< std::vector<OBVariant> > FindRows(const std::vector<Query> &query, bool *swapped = 0) const;
      const std::vector<std::vector<OBVariant> >& GetAllRows() const;
    private:
      std::string _name; 
      std::vector<std::string> _header;
      std::vector<OBVariant::Type> _types;
      std::vector< std::vector<OBVariant> > _rows;
      unsigned int _numColumns;
      const std::vector<OBVariant> _emptyRow;
      friend class OBFFParameterDB;
    };
    
    class OBFFParameterDB : public OBParameterDB
    {
    public:
      OBFFParameterDB(const std::string &databaseName = "");
      ~OBFFParameterDB();
      OBFFParameterDB(const OBFFParameterDB &rhs) {}; //no copying of databases allowed
      unsigned int NumTables() const;
      std::vector<std::string> GetTables() const;
      OBFFTable *  GetTable(const std::string &tableName) const;
      OBFFTable *  AddTable(const std::string &tableName, const std::vector<std::string> &header);
      OBFFTable *  AddTable(const std::string &tableName);
    private:
      std::string _name;
      std::vector<OBFFTable *> _tables;
    };
    
  }
}
