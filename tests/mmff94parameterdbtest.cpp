#include <OBVariant>
#include "obtest.h"

#include "../src/forcefields/mmff94/mmffparameter.h"

using namespace OpenBabel::OBFFs;

using namespace std;

void testTables(OBParameterDB *database)
{
  OB_ASSERT( database->NumTables() == 12 );
      
  std::vector<std::string> tables = database->GetTables();
  OB_ASSERT( tables.size() == 12 );

  OB_ASSERT( database->GetTable("Bond Parameters") );
  OB_ASSERT( database->GetTable("Angle Parameters") );
  OB_ASSERT( database->GetTable("Stretch-Bend Parameters") );
  OB_ASSERT( database->GetTable("Torsion Parameters") );
  OB_ASSERT( database->GetTable("Out-Of-Plane Parameters") );
  OB_ASSERT( database->GetTable("Charge Parameters") );
  OB_ASSERT( database->GetTable("Partial Bond Charge Increments") );
  OB_ASSERT( database->GetTable("Van der Waals Parameters") );
  OB_ASSERT( database->GetTable("Atom Properties") );
  OB_ASSERT( database->GetTable("Atom Type Levels") );
  OB_ASSERT( database->GetTable("Empirical Bond Parameters") );

}

void testAtomTypeLevels(OBParameterDB *database)
{
  OBParameterDBTable *table = database->GetTable("Atom Type Levels");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 95 );
  OB_ASSERT( table->NumColumns() == 5 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Int );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsInt() == 0 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(20)) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 20 );
  OB_ASSERT( row.at(1).AsInt() == 20 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsInt() == 0 );
}

void testAtomProperties(OBParameterDB *database)
{
  OBParameterDBTable *table = database->GetTable("Atom Properties");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 95 );
  OB_ASSERT( table->NumColumns() == 9 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 9 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
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
  std::vector<OBVariant> row = table->GetRow(0);
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
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(9)) );
  row = table->FindRow(query);
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
  OBParameterDBTable *table = database->GetTable("Bond Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 493 );
  OB_ASSERT( table->NumColumns() == 6 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 6 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 6 );
  OB_ASSERT( types.at(0) == OBVariant::String );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  OB_ASSERT( types.at(5) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(1).AsInt() == 0 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsDouble() == 4.258 );
  OB_ASSERT( row.at(5).AsDouble() == 1.508 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant("0:1-12")) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(1).AsInt() == 0 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 12 );
  OB_ASSERT( row.at(4).AsDouble() == 2.974 );
  OB_ASSERT( row.at(5).AsDouble() == 1.773 );

  // try swapped query
  /*
  query.clear();
  query.push_back( OBParameterDBTable::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(12), true) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(1), true) );
  bool swapped;
  row = table->FindRow(query, &swapped);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( swapped == true );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 12 );
  OB_ASSERT( row.at(3).AsDouble() == 2.974 );
  OB_ASSERT( row.at(4).AsDouble() == 1.773 );
  */

}

void testAngleParameters(OBParameterDB *database)
{
  OBParameterDBTable *table = database->GetTable("Angle Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 2342 );
  OB_ASSERT( table->NumColumns() == 6 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 6 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 6 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  OB_ASSERT( types.at(5) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 0 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 0 );
  OB_ASSERT( row.at(4).AsDouble() == 0.000 );
  OB_ASSERT( row.at(5).AsDouble() == 108.900 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(5)) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(2)) );
  query.push_back( OBParameterDBTable::Query(3, OBVariant(6)) );
  row = table->FindRow(query);
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
  OBParameterDBTable *table = database->GetTable("Stretch-Bend Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 282 );
  OB_ASSERT( table->NumColumns() == 6 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 6 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 6 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  OB_ASSERT( types.at(5) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 1 );
  OB_ASSERT( row.at(4).AsDouble() == 0.206 );
  OB_ASSERT( row.at(5).AsDouble() == 0.206 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(1)) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(20)) );
  query.push_back( OBParameterDBTable::Query(3, OBVariant(5)) );
  row = table->FindRow(query);
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
  OBParameterDBTable *table = database->GetTable("Torsion Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 926 );
  OB_ASSERT( table->NumColumns() == 8 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 8 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
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
  std::vector<OBVariant> row = table->GetRow(0);
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
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(5)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(1)) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(1)) );
  query.push_back( OBParameterDBTable::Query(3, OBVariant(8)) );
  query.push_back( OBParameterDBTable::Query(4, OBVariant(1)) );
  row = table->FindRow(query);
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
  query.push_back( OBParameterDBTable::Query(0, OBVariant(0)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(5), true) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(1), true) );
  query.push_back( OBParameterDBTable::Query(3, OBVariant(1), true) );
  query.push_back( OBParameterDBTable::Query(4, OBVariant(5), true) );
  row = table->FindRow(query);
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
  OBParameterDBTable *table = database->GetTable("Out-Of-Plane Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 117 );
  OB_ASSERT( table->NumColumns() == 5 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Int );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 2 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsInt() == 2 );
  OB_ASSERT( row.at(4).AsDouble() == 0.030 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(5)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(30)) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(20)) );
  query.push_back( OBParameterDBTable::Query(3, OBVariant(30)) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 5 );
  OB_ASSERT( row.at(1).AsInt() == 30 );
  OB_ASSERT( row.at(2).AsInt() == 20 );
  OB_ASSERT( row.at(3).AsInt() == 30 );
  OB_ASSERT( row.at(4).AsDouble() == 0.008 );
}

void testVanDerWaalsParameters(OBParameterDB *database)
{
  OBParameterDBTable *table = database->GetTable("Van der Waals Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 95 );
  OB_ASSERT( table->NumColumns() == 6 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 6 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 6 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Double );
  OB_ASSERT( types.at(2) == OBVariant::Double );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  OB_ASSERT( types.at(5) == OBVariant::Int );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 6 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsDouble() == 1.050 );
  OB_ASSERT( row.at(2).AsDouble() == 2.490 );
  OB_ASSERT( row.at(3).AsDouble() == 3.890 );
  OB_ASSERT( row.at(4).AsDouble() == 1.282 );
  OB_ASSERT( row.at(5).AsInt() == 0 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(12)) );
  row = table->FindRow(query);
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
  OBParameterDBTable *table = database->GetTable("Charge Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 498 );
  OB_ASSERT( table->NumColumns() == 4 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 4 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 4 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 1 );
  OB_ASSERT( row.at(3).AsDouble() == 0.0000 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(1)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(9)) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(78)) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 9 );
  OB_ASSERT( row.at(2).AsInt() == 78 );
  OB_ASSERT( row.at(3).AsDouble() == 0.1380 );
}

void testBondEmpiricalRules(OBParameterDB *database)
{
  OBParameterDBTable *table = database->GetTable("Empirical Bond Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 58 );
  OB_ASSERT( table->NumColumns() == 4 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 4 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 4 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Double );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 6 );
  OB_ASSERT( row.at(2).AsDouble() == 1.084 );
  OB_ASSERT( row.at(3).AsDouble() == 5.15 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(6)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(9)) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 4 );
  OB_ASSERT( row.at(0).AsInt() == 6 );
  OB_ASSERT( row.at(1).AsInt() == 9 );
  OB_ASSERT( row.at(2).AsDouble() == 1.353 );
  OB_ASSERT( row.at(3).AsDouble() == 6.20 );
}

void testStrBndEmpiricalRules(OBParameterDB *database)
{
  OBParameterDBTable *table = database->GetTable("Empirical Stretch-Bend Parameters");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 30 );
  OB_ASSERT( table->NumColumns() == 5 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Int );
  OB_ASSERT( types.at(2) == OBVariant::Int );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 0 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 0 );
  OB_ASSERT( row.at(3).AsDouble() == 0.15 );
  OB_ASSERT( row.at(4).AsDouble() == 0.15 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(1)) );
  query.push_back( OBParameterDBTable::Query(1, OBVariant(1)) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant(3)) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsInt() == 1 );
  OB_ASSERT( row.at(2).AsInt() == 3 );
  OB_ASSERT( row.at(3).AsDouble() == 0.30 );
  OB_ASSERT( row.at(4).AsDouble() == 0.50 );
}

void testPartialBondChargeIncrements(OBParameterDB *database)
{
  OBParameterDBTable *table = database->GetTable("Partial Bond Charge Increments");
  OB_REQUIRE( table );
  OB_ASSERT( table->NumRows() == 98 );
  OB_ASSERT( table->NumColumns() == 3 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 3 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 3 );
  OB_ASSERT( types.at(0) == OBVariant::Int );
  OB_ASSERT( types.at(1) == OBVariant::Double );
  OB_ASSERT( types.at(2) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 3 );
  OB_ASSERT( row.at(0).AsInt() == 1 );
  OB_ASSERT( row.at(1).AsDouble() == 0.000 );
  OB_ASSERT( row.at(2).AsDouble() == 0.000 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant(35)) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 3 );
  OB_ASSERT( row.at(0).AsInt() == 35 );
  OB_ASSERT( row.at(1).AsDouble() == -0.456 );
  OB_ASSERT( row.at(2).AsDouble() == 0.500 );
}

int main()
{
  std::cout << string(TESTDATADIR) + string("../data/mmff94.ff") << std::endl;
  MMFF94ParameterDB *database = new MMFF94ParameterDB(string(TESTDATADIR) + string("../data/mmff94.ff"));
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
