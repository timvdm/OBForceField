#include <OBVariant>
#include "obtest.h"

#include "../src/forcefields/mmff94/parameter.h"

using namespace OpenBabel::OBFFs;

using namespace std;

void testTables(OBParameterDB *database)
{
  OB_ASSERT( database->NumTables() == 12 );
      
  std::vector<std::string> tables = database->GetTables();
  OB_ASSERT( tables.size() == 12 );

  for (unsigned int i = 0; i < 12; ++i)
    OB_ASSERT( database->GetTable(tables[i]) == i);
}

void testAtomTypeLevels(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Atom Type Levels");
  OB_ASSERT( database->NumRows(table) == 95 );
  OB_ASSERT( database->NumColumns(table) == 5 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Int );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsInt() == 0 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(20)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 20 );
  OB_ASSERT( row.at(1).AsInt() == 20 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsInt() == 0 );
}

void testAtomProperties(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Atom Properties");
  OB_ASSERT( database->NumRows(table) == 95 );
  OB_ASSERT( database->NumColumns(table) == 9 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 9 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 9 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Bool );
  OB_ASSERT( types.at(5) == OBVariant::Int );
  OB_ASSERT( types.at(6) == OBVariant::Bool );
  OB_ASSERT( types.at(7) == OBVariant::Bool );
  OB_ASSERT( types.at(8) == OBVariant::Bool );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 9 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 6 );
  OB_ASSERT( row.at(2).AsInt() == 4 );
  OB_ASSERT( row.at(3).AsInt() == 4 );
  OB_ASSERT( row.at(4).AsInt() == 0 );
  OB_ASSERT( row.at(5).AsInt() == 0 );
  OB_ASSERT( row.at(6).AsInt() == 0 );
  OB_ASSERT( row.at(7).AsInt() == 0 );
  OB_ASSERT( row.at(8).AsInt() == 0 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(9)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 9 );
  OB_ASSERT( row.at(0).AsInt() == 9 );
  OB_ASSERT( row.at(1).AsInt() == 7 );
  OB_ASSERT( row.at(2).AsInt() == 2 );
  OB_ASSERT( row.at(3).AsInt() == 3 );
  OB_ASSERT( row.at(4).AsInt() == 0 );
  OB_ASSERT( row.at(5).AsInt() == 2 );
  OB_ASSERT( row.at(6).AsInt() == 0 );
  OB_ASSERT( row.at(7).AsInt() == 0 );
  OB_ASSERT( row.at(8).AsInt() == 1 );
}

void testBondParameters(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Bond Stretching Parameters");
  OB_ASSERT( database->NumRows(table) == 493 );
  OB_ASSERT( database->NumColumns(table) == 5 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsDouble() == 4.258 );
  OB_ASSERT( row.at(4).AsDouble() == 1.508 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(1)) );
  query.push_back( OBParameterDB::Query(2, OBVariant(12)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 12 );
  OB_ASSERT( row.at(3).AsDouble() == 2.974 );
  OB_ASSERT( row.at(4).AsDouble() == 1.773 );

  // try swapped query
  query.clear();
  query.push_back( OBParameterDB::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(12), true) );
  query.push_back( OBParameterDB::Query(2, OBVariant(1), true) );
  bool swapped;
  row = database->FindRow(table, query, &swapped);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( swapped == true );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 12 );
  OB_ASSERT( row.at(3).AsDouble() == 2.974 );
  OB_ASSERT( row.at(4).AsDouble() == 1.773 );






}

void testAngleParameters(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Angle Bending Parameters");
  OB_ASSERT( database->NumRows(table) == 2342 );
  OB_ASSERT( database->NumColumns(table) == 6 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 6 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 6 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  OB_ASSERT( types.at(5) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 0 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 0 );
  OB_ASSERT( row.at(4).AsDouble() == 0.000 );
  OB_ASSERT( row.at(5).AsDouble() == 108.900 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(5)) );
  query.push_back( OBParameterDB::Query(2, OBVariant(2)) );
  query.push_back( OBParameterDB::Query(3, OBVariant(6)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 5 );
  OB_ASSERT( row.at(2).AsInt() == 2 );
  OB_ASSERT( row.at(3).AsInt() == 6 );
  OB_ASSERT( row.at(4).AsDouble() == 0.589 );
  OB_ASSERT( row.at(5).AsDouble() == 108.757 );
}

void testStretchBendParameters(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Stretch Bend Parameters");
  OB_ASSERT( database->NumRows(table) == 282 );
  OB_ASSERT( database->NumColumns(table) == 6 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 6 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 6 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  OB_ASSERT( types.at(5) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsDouble() == 0.206 );
  OB_ASSERT( row.at(5).AsDouble() == 0.206 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(1)) );
  query.push_back( OBParameterDB::Query(2, OBVariant(20)) );
  query.push_back( OBParameterDB::Query(3, OBVariant(5)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 20 );
  OB_ASSERT( row.at(3).AsInt() == 5 );
  OB_ASSERT( row.at(4).AsDouble() == 0.290 );
  OB_ASSERT( row.at(5).AsDouble() == 0.098 );
}

void testTorsionParameters(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Torsion Parameters");
  OB_ASSERT( database->NumRows(table) == 926 );
  OB_ASSERT( database->NumColumns(table) == 8 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 8 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 8 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Int );
  OB_ASSERT( types.at(5) == OBVariant::Double );
  OB_ASSERT( types.at(6) == OBVariant::Double );
  OB_ASSERT( types.at(7) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 8 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 0 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsInt() == 0 );
  OB_ASSERT( row.at(5).AsDouble() == 0.000 );
  OB_ASSERT( row.at(6).AsDouble() == 0.000 );
  OB_ASSERT( row.at(7).AsDouble() == 0.300 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(5)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(1)) );
  query.push_back( OBParameterDB::Query(2, OBVariant(1)) );
  query.push_back( OBParameterDB::Query(3, OBVariant(8)) );
  query.push_back( OBParameterDB::Query(4, OBVariant(1)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 8 );
  OB_ASSERT( row.at(0).AsInt() == 5 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 8 );
  OB_ASSERT( row.at(4).AsInt() == 1 );
  OB_ASSERT( row.at(5).AsDouble() == 0.115 );
  OB_ASSERT( row.at(6).AsDouble() == -0.390 );
  OB_ASSERT( row.at(7).AsDouble() == 0.658 );


  // try swapped query
  query.clear();
  query.push_back( OBParameterDB::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(5), true) );
  query.push_back( OBParameterDB::Query(2, OBVariant(1), true) );
  query.push_back( OBParameterDB::Query(3, OBVariant(1), true) );
  query.push_back( OBParameterDB::Query(4, OBVariant(5), true) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 8 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 5 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsInt() == 5 );
  OB_ASSERT( row.at(5).AsDouble() == 0.284 );
  OB_ASSERT( row.at(6).AsDouble() == -1.386 );
  OB_ASSERT( row.at(7).AsDouble() == 0.314 );

}

void testOutOfPlaneParameters(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Out-Of-Plane Parameters");
  OB_ASSERT( database->NumRows(table) == 117 );
  OB_ASSERT( database->NumColumns(table) == 5 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 2 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 2 );
  OB_ASSERT( row.at(4).AsDouble() == 0.030 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(5)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(30)) );
  query.push_back( OBParameterDB::Query(2, OBVariant(20)) );
  query.push_back( OBParameterDB::Query(3, OBVariant(30)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 5 );
  OB_ASSERT( row.at(1).AsInt() == 30 );
  OB_ASSERT( row.at(2).AsInt() == 20 );
  OB_ASSERT( row.at(3).AsInt() == 30 );
  OB_ASSERT( row.at(4).AsDouble() == 0.008 );
}

void testVanDerWaalsParameters(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Van der Waals Parameters");
  OB_ASSERT( database->NumRows(table) == 95 );
  OB_ASSERT( database->NumColumns(table) == 6 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 6 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 6 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Double );
  OB_ASSERT( types.at(2) == OBVariant::Double );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  OB_ASSERT( types.at(5) == OBVariant::Int );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsDouble() == 1.050 );
  OB_ASSERT( row.at(2).AsDouble() == 2.490 );
  OB_ASSERT( row.at(3).AsDouble() == 3.890 );
  OB_ASSERT( row.at(4).AsDouble() == 1.282 );
  OB_ASSERT( row.at(5).AsInt() == 0 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(12)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 12 );
  OB_ASSERT( row.at(1).AsDouble() == 2.300 );
  OB_ASSERT( row.at(2).AsDouble() == 5.100 );
  OB_ASSERT( row.at(3).AsDouble() == 3.320 );
  OB_ASSERT( row.at(4).AsDouble() == 1.345 );
  OB_ASSERT( row.at(5).AsInt() == 2 );
}

void testChargeParameters(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Charge Parameters");
  OB_ASSERT( database->NumRows(table) == 498 );
  OB_ASSERT( database->NumColumns(table) == 4 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 4 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 4 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsDouble() == 0.0000 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(1)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(9)) );
  query.push_back( OBParameterDB::Query(2, OBVariant(78)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 9 );
  OB_ASSERT( row.at(2).AsInt() == 78 );
  OB_ASSERT( row.at(3).AsDouble() == 0.1380 );
}

void testBondEmpiricalRules(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Bond Stretching Emprirical Rules");
  OB_ASSERT( database->NumRows(table) == 58 );
  OB_ASSERT( database->NumColumns(table) == 4 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 4 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 4 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Double );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 6 );
  OB_ASSERT( row.at(2).AsDouble() == 1.084 );
  OB_ASSERT( row.at(3).AsDouble() == 5.15 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(6)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(9)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 6 );
  OB_ASSERT( row.at(1).AsInt() == 9 );
  OB_ASSERT( row.at(2).AsDouble() == 1.353 );
  OB_ASSERT( row.at(3).AsDouble() == 6.20 );
}

void testStrBndEmpiricalRules(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Stretch Bending Emprirical Rules");
  OB_ASSERT( database->NumRows(table) == 30 );
  OB_ASSERT( database->NumColumns(table) == 5 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 0 );
  OB_ASSERT( row.at(3).AsDouble() == 0.15 );
  OB_ASSERT( row.at(4).AsDouble() == 0.15 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(1)) );
  query.push_back( OBParameterDB::Query(1, OBVariant(1)) );
  query.push_back( OBParameterDB::Query(2, OBVariant(3)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 3 );
  OB_ASSERT( row.at(3).AsDouble() == 0.30 );
  OB_ASSERT( row.at(4).AsDouble() == 0.50 );
}

void testPartialBondChargeIncrements(OBParameterDB *database)
{
  unsigned int table = database->GetTable("MMFF94 Partial Bond Charge Increments");
  OB_ASSERT( database->NumRows(table) == 98 );
  OB_ASSERT( database->NumColumns(table) == 3 );
  
  // check the header
  std::vector<std::string> header = database->GetHeader(table);
  OB_ASSERT( header.size() == 3 );

  // check the row types
  std::vector<OBVariant::Type> types = database->GetTypes(table);
  OB_ASSERT( types.size() == 3 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Double );
  OB_ASSERT( types.at(2) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = database->GetRow(table, 0);
  OB_ASSERT( row.size() == 3 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsDouble() == 0.000 );
  OB_ASSERT( row.at(2).AsDouble() == 0.000 );

  // try query
  std::vector<OBParameterDB::Query> query;
  query.push_back( OBParameterDB::Query(0, OBVariant(35)) );
  row = database->FindRow(table, query);
  OB_ASSERT( row.size() == 3 );
  OB_ASSERT( row.at(0).AsInt() == 35 );
  OB_ASSERT( row.at(1).AsDouble() == -0.456 );
  OB_ASSERT( row.at(2).AsDouble() == 0.500 );
}

int main()
{
  MMFF94SimpleParameterDB *database = new MMFF94SimpleParameterDB("/home/timvdm/OBForceField/laptop/data/mmff94.ff");
  testTables(database);
  testAtomProperties(database);
  testAtomTypeLevels(database);
  testBondParameters(database);
  testAngleParameters(database);
  testStretchBendParameters(database);
  testTorsionParameters(database);
  testOutOfPlaneParameters(database);
  testVanDerWaalsParameters(database);
  testChargeParameters(database);
  testBondEmpiricalRules(database);
  testStrBndEmpiricalRules(database);
  testPartialBondChargeIncrements(database);





}
