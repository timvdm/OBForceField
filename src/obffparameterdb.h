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
      /**
       * Constructor specifying the @p tableName and @p header.
       */
      OBFFTable(const std::string &tableName, const std::vector<std::string> &header);
      /**
       * Constructor specifying the @p tableName.
       */
      OBFFTable(const std::string &tableName);
      /**
       * @return The number of rows in this table.
       */
      unsigned int NumRows() const;
      unsigned int NumColumns() const;
      /**
       * @return std::vector containing the column headers as std:sstring
       */
      std::vector<std::string> GetHeader() const;
      /**
       * @return std::vector containing the column types as OBVariant::Type.
       */
      std::vector<OBVariant::Type> GetTypes() const;
      /**
       * @return True if all columns types in @p values match this tables column types.
       * @sa GetTypes()
       */
      bool VerifyTypes(const std::vector<OBVariant> &values);
      /**
       * Get the row with index @p row.
       */
      const std::vector<OBVariant>& GetRow(unsigned int row) const;
      /**
       * Add a row to this table containing @p values.
       */
      bool AddRow(const std::vector<OBVariant> &values);
      /**
       * Find a row that matches the specified @p query. If the query was matched in
       * reverse order (e.g. 1-2-3 vs. 3-2-1), the swapped flag will be set when non-zero.
       * @return A constant reference to the row.
       */
      const std::vector<OBVariant>& FindRow(const std::vector<Query> &query, bool *swapped = 0) const;
      /**
       * Find a row that matches the specified @p query. If the query was matched in
       * reverse order (e.g. 1-2-3 vs. 3-2-1), the swapped flag will be set when non-zero.
       * @return A copy of the row.
       */
      std::vector< std::vector<OBVariant> > FindRows(const std::vector<Query> &query, bool *swapped = 0) const;
      /**
       * Get a constant reference to all rows in this table. This function can 
       * be used in graphical user interfaces to quickly access all data for display.
       */
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
      /**
       * Constructor specifying the @p databaseName.
       */
      OBFFParameterDB(const std::string &databaseName = "");
      /**
       * Destructor.
       */
      ~OBFFParameterDB();

      OBFFParameterDB(const OBFFParameterDB &rhs) {}; //no copying of databases allowed
      /**
       * @return The number of tables in this database.
       */
      unsigned int NumTables() const;
      /**
       * @return std::vector containing the table names as std::string.
       */
      std::vector<std::string> GetTables() const;
      /**
       * Get the table with the specified @p tableName.
       */
      OBFFTable* GetTable(const std::string &tableName) const;
      /**
       * Add a table named @p tableName and @p header.
       * @return A pointer to the newly added table.
       */
      OBFFTable* AddTable(const std::string &tableName, const std::vector<std::string> &header);
      /**
       * Add a table named @p tableName.
       * @return A pointer to the newly added table.
       */
      OBFFTable* AddTable(const std::string &tableName);
    private:
      std::string _name;
      std::vector<OBFFTable *> _tables;
    };
    
  }
}
