#ifndef OBFFS_PARAMETERDB_H
#define OBFFS_PARAMETERDB_H

#include <OBVariant>

#include <string>
#include <vector>

namespace OpenBabel {
  namespace OBFFs {
    
    /**
     * @class OBParameterDB
     * @brief Abstract class for simple databases
     *
     */
    class OBParameterDBTable
    {
    public:
      
      struct Query {
	Query(int _column, const OBVariant &_value, bool _swap = false) : column(_column), value(_value), swap(_swap)
        {
        }
        int column;
        OBVariant value;
        bool swap;
      };

      static Query MakeQuery(int column, const OBVariant &value);
      
      virtual unsigned int NumRows() const = 0;
      /**
       * Get the number of rows for the table with index @p table.
       */
      virtual unsigned int NumColumns() const = 0;
      /**
       * Get the header for each column in the table with index @p table.
       */
      virtual std::vector<std::string> GetHeader() const = 0;
      /**
       * Get the type for @p table.
       */
      virtual std::vector<OBVariant::Type> GetTypes() const = 0;
      /**
       * Chech if the types in @p values are correct for @p tabl.
       */
      virtual bool VerifyTypes(const std::vector<OBVariant> &values) = 0;
      /**
       * Get a row from a table.
       */
      virtual const std::vector<OBVariant>& GetRow(unsigned int row) const = 0;
      virtual bool AddRow(const std::vector<OBVariant> &values) = 0;
      
      virtual const std::vector<OBVariant>& FindRow(const std::vector<Query> &query, bool *swapped = 0) const = 0;

      virtual std::vector< std::vector<OBVariant> > FindRows(const std::vector<Query> &query, bool *swapped = 0) const = 0; //return all matches 
      
      virtual const std::vector<std::vector<OBVariant> >& GetAllRows() const = 0;
    };
    
    class OBParameterDB
    {
    public:
      virtual unsigned int NumTables() const = 0;
      /**
       * Get the names for the tables in this database.
       */
      virtual std::vector<std::string> GetTables() const = 0;
      /**
       * Get the unsigned int for the table named @table or NumTables() if there is no such table.
       */
      virtual OBParameterDBTable * GetTable(const std::string &tableName) const = 0;
      virtual OBParameterDBTable * AddTable(const std::string &tableName) = 0;
    };

  } // namespace
}

#endif
