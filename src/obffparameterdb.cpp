/*********************************************************************
OBFF parameter database

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

#include "obffparameterdb.h"

#include <string>
#include <vector>

using namespace std;

namespace OpenBabel {
  namespace OBFFs {

    OBFFTable::OBFFTable(const string &tableName, const vector<string> &header)
      :_name(tableName), _header(header), _numColumns(header.size()) {}

    OBFFTable::OBFFTable(const string &tableName)
      :_name(tableName), _numColumns(0) {}
      
    unsigned int OBFFTable::NumRows() const
    {
      return _rows.size();
    }

    unsigned int OBFFTable::NumColumns() const
    {
      return _numColumns;
    }
          
    vector<string> OBFFTable::GetHeader() const
    {
      return _header;      
    }
    
    vector<OBVariant::Type> OBFFTable::GetTypes() const
    {
      if (_rows.size() != 0)
	return _types;
      vector<OBVariant::Type> types;
      return types;
    }

    bool OBFFTable::VerifyTypes(const vector<OBVariant> &values) 
    {
      if (_types.size() != values.size())
	return false;

      for (unsigned int i = 0; i != _types.size(); ++i)
	if (values[i].GetType() != _types[i])
	  return false;

      return true;
    }

    const vector<OBVariant>& OBFFTable::GetRow(unsigned int i) const
    {
      if (i < _rows.size())
	return _rows[i];
      else
	return _emptyRow;
    }

    bool OBFFTable::AddRow(const vector<OBVariant> &values)
    {
      if (_rows.size()!=0) {
	if (values.size()!=_numColumns)
	  return false;
	for (unsigned int i = 0; i != _types.size(); ++i)
	  if (values[i].GetType() != _types[i])
	    return false;
      }
      else {
	for (unsigned int i = 0; i != values.size(); ++i)
	  _types.push_back(values[i].GetType());
	_numColumns = values.size();
      }
      _rows.push_back(values);
      return true;
    }

    const vector<OBVariant>& OBFFTable::FindRow(const vector<Query> &query, bool *swapped) const
    {
      unsigned int numRows = NumRows();

      unsigned int swapCount = 0;
      // make sure the query contains valid columns
      for (unsigned int i = 0; i < query.size(); ++i) {
	if (query.at(i).column >= _numColumns)
	  return _emptyRow;
	if (query.at(i).swap)
	  swapCount++;
      }
      
      if (!swapCount) {
	for (unsigned int i = 0; i < numRows; ++i) {
	  bool match = true;
	  // check each query entry
	  for (unsigned int j = 0; j < query.size(); ++j) {
	    if (_rows[i].at(query.at(j).column) != query.at(j).value) {
	      match = false;
	      break;
	    }
	  }

	  if (match)
	    return _rows[i];
	}
      } else {
	// construct swapped_query
	vector<Query> swapped_query = query;
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
	    if (_rows[i].at(query.at(j).column) != query.at(j).value) {
	      match = false;
	      break;
	    }
	  }

	  if (match) {
	    if (swapped)
	      *swapped = false;
	    return _rows[i];
	  }
 
	  match = true;
	  // check swapped query
	  for (unsigned int j = 0; j < swapped_query.size(); ++j) {
	    if (_rows[i].at(swapped_query.at(j).column) != swapped_query.at(j).value) {
	      match = false;
	      break;
	    }
	  }
 
	  if (match) {
	    if (swapped)
	      *swapped = true;
	    return _rows[i];
	  }
	}
      }
     
      return _emptyRow;
    }

    vector< vector<OBVariant> > OBFFTable::FindRows(const vector<Query> &query, bool *swapped) const
    {
      unsigned int numRows = NumRows();
      vector< vector<OBVariant> > rows;

      unsigned int swapCount = 0;
      // make sure the query contains valid columns
      for (unsigned int i = 0; i < query.size(); ++i) {
	if (query.at(i).column >= _numColumns)
	  return rows;
	if (query.at(i).swap)
	  swapCount++;
      }
      
      if (!swapCount) {
	for (unsigned int i = 0; i < numRows; ++i) {
	  bool match = true;
	  // check each query entry
	  for (unsigned int j = 0; j < query.size(); ++j) {
	    if (_rows[i].at(query.at(j).column) != query.at(j).value) {
	      match = false;
	      break;
	    }
	  }

	  if (match)
	    rows.push_back(_rows[i]);
	}
      } else {
	// construct swapped_query
	vector<Query> swapped_query = query;
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
	    if (_rows[i].at(query.at(j).column) != query.at(j).value) {
	      match = false;
	      break;
	    }
	  }

	  if (match) {
	    if (swapped)
	      *swapped = false;
	    rows.push_back(_rows[i]);
	  }
 
	  match = true;
	  // check swapped query
	  for (unsigned int j = 0; j < swapped_query.size(); ++j) {
	    if (_rows[i].at(swapped_query.at(j).column) != swapped_query.at(j).value) {
	      match = false;
	      break;
	    }
	  }
 
	  if (match) {
	    if (swapped)
	      *swapped = true;
	    rows.push_back(_rows[i]);
	  }
	}
      }
     
      return rows;
    }
      
    const vector<vector<OBVariant> >& OBFFTable::GetAllRows() const
    {
      return _rows;
    }
    
    OBFFParameterDB::OBFFParameterDB(const string &databaseName)
      :_name(databaseName) {}
    
    OBFFParameterDB::~OBFFParameterDB(){
      vector<OBFFTable *>::iterator itr;
      for(itr=_tables.begin(); itr!=_tables.end(); ++itr)
	delete *itr;
    }

    unsigned int OBFFParameterDB::NumTables() const 
    { 
      return _tables.size();
    }
    
    vector<string> OBFFParameterDB::GetTables() const
    {
      vector<string> tableNames;
      vector<OBFFTable *>::const_iterator itr;
      for(itr = _tables.begin(); itr != _tables.end(); ++itr)
	tableNames.push_back((*itr)->_name);
      return tableNames;
    }
    
    OBFFTable * OBFFParameterDB::GetTable(const string &tableName) const
    {
      vector<OBFFTable *>::const_iterator itr;
      for(itr = _tables.begin(); itr !=_tables.end(); ++itr)
	if ((*itr)->_name == tableName) 
	  return *itr;
      return NULL;
    }

    OBFFTable * OBFFParameterDB::AddTable(const string &tableName)
    {
      OBFFTable * p = new OBFFTable(tableName);
      _tables.push_back(p);
      return p;
    }

    OBFFTable * OBFFParameterDB::AddTable(const string &tableName, const vector<string> &header)
    {
      OBFFTable * p = new OBFFTable(tableName, header);
      _tables.push_back(p);
      return p;
    }


  } // end namespace OpenBabel
}
//! \brief OBFF parameter database
