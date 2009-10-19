#include <OBVariant>
#include "obtest.h"

#include "../src/forcefields/gaff/gaffparameter.h"

using namespace OpenBabel::OBFFs;

using namespace std;

 void testTables(OBParameterDB *database)
{
  OB_ASSERT( database->NumTables() == 7 );
      
  std::vector<std::string> tables = database->GetTables();
  OB_ASSERT( tables.size() == 7 );
}

void testAtomProperties(OBParameterDB *database)
{
  OBParameterDBTable * table = database->GetTable("Atom Properties");
  OB_ASSERT( table->NumRows() == 70 );
  OB_ASSERT( table->NumColumns() == 3 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 3 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 3 );
  OB_ASSERT( types.at(0) == OBVariant::String );
  OB_ASSERT( types.at(1) == OBVariant::Double );
  OB_ASSERT( types.at(2) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 3 );
  OB_ASSERT( row.at(0).AsString() == "c" );
  OB_ASSERT( row.at(1).AsDouble() == 12.01 );
  OB_ASSERT( row.at(2).AsDouble() == 0.616 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back(OBParameterDBTable::Query(0, OBVariant("h4")) );
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 3 );
  OB_ASSERT( row.at(0).AsString() == "h4" );
  OB_ASSERT( row.at(1).AsDouble() == 1.008 );
  OB_ASSERT( row.at(2).AsDouble() == 0.135 );
}

void testBondParameters(OBParameterDB *database)
{
  OBParameterDBTable * table = database->GetTable("Bond Harmonic");
  OB_ASSERT( table->NumRows() == 748 );
  OB_ASSERT( table->NumColumns() == 5 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::String );
  OB_ASSERT( types.at(1) == OBVariant::String );
  OB_ASSERT( types.at(2) == OBVariant::String );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsString() == "hw-ow" );
  OB_ASSERT( row.at(1).AsString() == "ow" );
  OB_ASSERT( row.at(2).AsString() == "hw" );
  OB_ASSERT( row.at(3).AsDouble() == 553.0 );
  OB_ASSERT( row.at(4).AsDouble() == 0.9572 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant("cx-oh")) );
  row = table->FindRow(query);
  OB_ASSERT( row.at(0).AsString() == "cx-oh" );
  OB_ASSERT( row.at(1).AsString() == "cx" );
  OB_ASSERT( row.at(2).AsString() == "oh" );
  OB_ASSERT( row.at(3).AsDouble() == 387.4 );
  OB_ASSERT( row.at(4).AsDouble() == 1.361 );

  // try query
  query.clear();
  query.push_back( OBParameterDBTable::Query(1, OBVariant("cx")) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant("oh")) );
  row = table->FindRow(query);
  OB_ASSERT( row.at(0).AsString() == "cx-oh" );
  OB_ASSERT( row.at(1).AsString() == "cx" );
  OB_ASSERT( row.at(2).AsString() == "oh" );
  OB_ASSERT( row.at(3).AsDouble() == 387.4 );
  OB_ASSERT( row.at(4).AsDouble() == 1.361 );

  // try swapped query
  query.clear();
  query.push_back( OBParameterDBTable::Query(1, OBVariant("oh"), true) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant("cx"), true) );
  query.push_back( OBParameterDBTable::Query(3, OBVariant(387.4)) );
  bool swapped;
  row = table->FindRow(query, &swapped);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( swapped == true );
  OB_ASSERT( row.at(0).AsString() == "cx-oh" );
  OB_ASSERT( row.at(1).AsString() == "cx" );
  OB_ASSERT( row.at(2).AsString() == "oh" );
  OB_ASSERT( row.at(3).AsDouble() == 387.4 );
  OB_ASSERT( row.at(4).AsDouble() == 1.361 );
}

void testAngleParameters(OBParameterDB *database)
{
  OBParameterDBTable * table = database->GetTable("Angle Bending");
  OB_ASSERT( table->NumRows() == 3513 );
  OB_ASSERT( table->NumColumns() == 5 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 5 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 5 );
  OB_ASSERT( types.at(0) == OBVariant::String );
  OB_ASSERT( types.at(1) == OBVariant::String );
  OB_ASSERT( types.at(2) == OBVariant::String );
  OB_ASSERT( types.at(3) == OBVariant::Double );
  OB_ASSERT( types.at(4) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(0).AsString() == "hw" );
  OB_ASSERT( row.at(1).AsString() == "ow" );
  OB_ASSERT( row.at(2).AsString() == "hw" );
  OB_ASSERT( row.at(3).AsDouble() == 100 );
  OB_ASSERT( row.at(4).AsDouble() == 104.52 );

  // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(0, OBVariant("s"), true ));
  query.push_back( OBParameterDBTable::Query(1, OBVariant("c")) );
  query.push_back( OBParameterDBTable::Query(2, OBVariant("c2"), true ));
  row = table->FindRow(query);
  OB_ASSERT( row.size() == 5 );
  OB_ASSERT( row.at(3).AsDouble() == 64.7 );
  OB_ASSERT( row.at(4).AsDouble() == 119.16 );
}

void testTorsionParameters(OBParameterDB *database)
{
  OBParameterDBTable * table = database->GetTable("Torsion Harmonic");
  OB_ASSERT( table->NumRows() == 698 );
  OB_ASSERT( table->NumColumns() == 8 );
  
  // check the header
  std::vector<std::string> header = table->GetHeader();
  OB_ASSERT( header.size() == 8 );

  // check the row types
  std::vector<OBVariant::Type> types = table->GetTypes();
  OB_ASSERT( types.size() == 8 );
  OB_ASSERT( types.at(0) == OBVariant::String );
  OB_ASSERT( types.at(1) == OBVariant::String );
  OB_ASSERT( types.at(2) == OBVariant::String );
  OB_ASSERT( types.at(3) == OBVariant::String );
  OB_ASSERT( types.at(4) == OBVariant::String );
  OB_ASSERT( types.at(5) == OBVariant::Double );
  OB_ASSERT( types.at(6) == OBVariant::Double );
  OB_ASSERT( types.at(7) == OBVariant::Double );
  
  // try first entry
  std::vector<OBVariant> row = table->GetRow(0);
  OB_ASSERT( row.size() == 8 );
  OB_ASSERT( row.at(0).AsString() == "X-c-c-X" );
  OB_ASSERT( row.at(1).AsString() == "X" );
  OB_ASSERT( row.at(2).AsString() == "c" );
  OB_ASSERT( row.at(3).AsString() == "c" );
  OB_ASSERT( row.at(4).AsString() == "X" );
  OB_ASSERT( row.at(5).AsDouble() == 0.3 );
  OB_ASSERT( row.at(6).AsDouble() == -1.0 );
  OB_ASSERT( row.at(7).AsDouble() == 2.0 );

   // try query
  std::vector<OBParameterDBTable::Query> query;
  query.push_back( OBParameterDBTable::Query(2, OBVariant("ca")) );
  row = table->FindRow(query);
  OB_ASSERT( row.at(0).AsString() == "X-ca-ca-X" );
  OB_ASSERT( row.at(1).AsString() == "X" );
  OB_ASSERT( row.at(2).AsString() == "ca" );
  OB_ASSERT( row.at(3).AsString() == "ca" );
  OB_ASSERT( row.at(4).AsString() == "X" );
  OB_ASSERT( row.at(5).AsDouble() == 3.625 );
  OB_ASSERT( row.at(6).AsDouble() == -1.0 );
  OB_ASSERT( row.at(7).AsDouble() == 2 );

//   OB_ASSERT( row.size() == 8 );
//   OB_ASSERT( row.at(0).AsInt() == 5 );
//   OB_ASSERT( row.at(1).AsInt() == 1 );
//   OB_ASSERT( row.at(2).AsInt() == 1 );
//   OB_ASSERT( row.at(3).AsInt() == 8 );
//   OB_ASSERT( row.at(4).AsInt() == 1 );
//   OB_ASSERT( row.at(5).AsDouble() == 0.115 );
//   OB_ASSERT( row.at(6).AsDouble() == -0.390 );
//   OB_ASSERT( row.at(7).AsDouble() == 0.658 );


  query.clear();
  query.push_back( OBParameterDBTable::Query(3, OBVariant("ca")) );
  std::vector< std::vector<OBVariant> > rows;
  std::vector< std::vector<OBVariant> >::const_iterator itr;
  rows = table->FindRows(query);
  for(itr=rows.begin();itr!=rows.end();++itr){
   cout << itr->at(0).AsString() << " " << itr->at(1).AsString() << " " << itr->at(2).AsString() << " " << itr->at(3).AsString() << " " << itr->at(4).AsString();
   cout << " " << itr->at(5).AsDouble() << " " << itr->at(6).AsDouble() << " " << itr->at(7).AsDouble() << endl;
  }
  cout << endl;

  query.clear();
  query.push_back( OBParameterDBTable::Query(0, OBVariant("X-n2-n2-X")) );
  rows = table->FindRows(query);
  for(itr=rows.begin();itr!=rows.end();++itr){
   cout << itr->at(0).AsString() << " " << itr->at(1).AsString() << " " << itr->at(2).AsString() << " " << itr->at(3).AsString() << " " << itr->at(4).AsString();
   cout << " " << itr->at(5).AsDouble() << " " << itr->at(6).AsDouble() << " " << itr->at(7).AsDouble() << endl;
  }


  


//   // try swapped query
//   query.clear();
//   query.push_back( OBParameterDBTable::Query(0, OBVariant(0)) );
//   query.push_back( OBParameterDBTable::Query(1, OBVariant(5), true) );
//   query.push_back( OBParameterDBTable::Query(2, OBVariant(1), true) );
//   query.push_back( OBParameterDBTable::Query(3, OBVariant(1), true) );
//   query.push_back( OBParameterDBTable::Query(4, OBVariant(5), true) );
//   row = table->FindRow(query);
//   OB_ASSERT( row.size() == 8 );
//   OB_ASSERT( row.at(0).AsInt() == 0 );
//   OB_ASSERT( row.at(1).AsInt() == 5 );
//   OB_ASSERT( row.at(2).AsInt() == 1 );
//   OB_ASSERT( row.at(3).AsInt() == 1 );
//   OB_ASSERT( row.at(4).AsInt() == 5 );
//   OB_ASSERT( row.at(5).AsDouble() == 0.284 );
//   OB_ASSERT( row.at(6).AsDouble() == -1.386 );
//   OB_ASSERT( row.at(7).AsDouble() == 0.314 );

}

// void testOutOfPlaneParameters(OBParameterDB *database)
// {
//   unsigned int table = database->GetTable("GAFF Out-Of-Plane Parameters");
//   OB_ASSERT( table->NumRows() == 117 );
//   OB_ASSERT( table->NumColumns() == 5 );
  
//   // check the header
//   std::vector<std::string> header = table->GetHeader();
//   OB_ASSERT( header.size() == 5 );

//   // check the row types
//   std::vector<OBVariant::Type> types = table->GetTypes();
//   OB_ASSERT( types.size() == 5 );
//   OB_ASSERT( types.at(0) == OBVariant::Int );
//   OB_ASSERT( types.at(1) == OBVariant::Int );
//   OB_ASSERT( types.at(2) == OBVariant::Int );
//   OB_ASSERT( types.at(3) == OBVariant::Int );
//   OB_ASSERT( types.at(4) == OBVariant::Double );
  
//   // try first entry
//   std::vector<OBVariant> row = table->GetRow(0);
//   OB_ASSERT( row.size() == 5 );
//   OB_ASSERT( row.at(0).AsInt() == 1 );
//   OB_ASSERT( row.at(1).AsInt() == 2 );
//   OB_ASSERT( row.at(2).AsInt() == 1 );
//   OB_ASSERT( row.at(3).AsInt() == 2 );
//   OB_ASSERT( row.at(4).AsDouble() == 0.030 );

//   // try query
//   std::vector<OBParameterDBTable::Query> query;
//   query.push_back( OBParameterDBTable::Query(0, OBVariant(5)) );
//   query.push_back( OBParameterDBTable::Query(1, OBVariant(30)) );
//   query.push_back( OBParameterDBTable::Query(2, OBVariant(20)) );
//   query.push_back( OBParameterDBTable::Query(3, OBVariant(30)) );
//   row = table->FindRow(query);
//   OB_ASSERT( row.size() == 5 );
//   OB_ASSERT( row.at(0).AsInt() == 5 );
//   OB_ASSERT( row.at(1).AsInt() == 30 );
//   OB_ASSERT( row.at(2).AsInt() == 20 );
//   OB_ASSERT( row.at(3).AsInt() == 30 );
//   OB_ASSERT( row.at(4).AsDouble() == 0.008 );
// }

// void testVanDerWaalsParameters(OBParameterDB *database)
// {
//   unsigned int table = database->GetTable("GAFF Van der Waals Parameters");
//   OB_ASSERT( table->NumRows() == 95 );
//   OB_ASSERT( table->NumColumns() == 6 );
  
//   // check the header
//   std::vector<std::string> header = table->GetHeader();
//   OB_ASSERT( header.size() == 6 );

//   // check the row types
//   std::vector<OBVariant::Type> types = table->GetTypes();
//   OB_ASSERT( types.size() == 6 );
//   OB_ASSERT( types.at(0) == OBVariant::Int );
//   OB_ASSERT( types.at(1) == OBVariant::Double );
//   OB_ASSERT( types.at(2) == OBVariant::Double );
//   OB_ASSERT( types.at(3) == OBVariant::Double );
//   OB_ASSERT( types.at(4) == OBVariant::Double );
//   OB_ASSERT( types.at(5) == OBVariant::Int );
  
//   // try first entry
//   std::vector<OBVariant> row = table->GetRow(0);
//   OB_ASSERT( row.size() == 6 );
//   OB_ASSERT( row.at(0).AsInt() == 1 );
//   OB_ASSERT( row.at(1).AsDouble() == 1.050 );
//   OB_ASSERT( row.at(2).AsDouble() == 2.490 );
//   OB_ASSERT( row.at(3).AsDouble() == 3.890 );
//   OB_ASSERT( row.at(4).AsDouble() == 1.282 );
//   OB_ASSERT( row.at(5).AsInt() == 0 );

//   // try query
//   std::vector<OBParameterDBTable::Query> query;
//   query.push_back( OBParameterDBTable::Query(0, OBVariant(12)) );
//   row = table->FindRow(query);
//   OB_ASSERT( row.size() == 6 );
//   OB_ASSERT( row.at(0).AsInt() == 12 );
//   OB_ASSERT( row.at(1).AsDouble() == 2.300 );
//   OB_ASSERT( row.at(2).AsDouble() == 5.100 );
//   OB_ASSERT( row.at(3).AsDouble() == 3.320 );
//   OB_ASSERT( row.at(4).AsDouble() == 1.345 );
//   OB_ASSERT( row.at(5).AsInt() == 2 );
// }

// void testChargeParameters(OBParameterDB *database)
// {
//   unsigned int table = database->GetTable("GAFF Charge Parameters");
//   OB_ASSERT( table->NumRows() == 498 );
//   OB_ASSERT( table->NumColumns() == 4 );
  
//   // check the header
//   std::vector<std::string> header = table->GetHeader();
//   OB_ASSERT( header.size() == 4 );

//   // check the row types
//   std::vector<OBVariant::Type> types = table->GetTypes();
//   OB_ASSERT( types.size() == 4 );
//   OB_ASSERT( types.at(0) == OBVariant::Int );
//   OB_ASSERT( types.at(1) == OBVariant::Int );
//   OB_ASSERT( types.at(2) == OBVariant::Int );
//   OB_ASSERT( types.at(3) == OBVariant::Double );
  
//   // try first entry
//   std::vector<OBVariant> row = table->GetRow(0);
//   OB_ASSERT( row.size() == 4 );
//   OB_ASSERT( row.at(0).AsInt() == 0 );
//   OB_ASSERT( row.at(1).AsInt() == 1 );
//   OB_ASSERT( row.at(2).AsInt() == 1 );
//   OB_ASSERT( row.at(3).AsDouble() == 0.0000 );

//   // try query
//   std::vector<OBParameterDBTable::Query> query;
//   query.push_back( OBParameterDBTable::Query(0, OBVariant(1)) );
//   query.push_back( OBParameterDBTable::Query(1, OBVariant(9)) );
//   query.push_back( OBParameterDBTable::Query(2, OBVariant(78)) );
//   row = table->FindRow(query);
//   OB_ASSERT( row.size() == 4 );
//   OB_ASSERT( row.at(0).AsInt() == 1 );
//   OB_ASSERT( row.at(1).AsInt() == 9 );
//   OB_ASSERT( row.at(2).AsInt() == 78 );
//   OB_ASSERT( row.at(3).AsDouble() == 0.1380 );
// }

// void testBondEmpiricalRules(OBParameterDB *database)
// {
//   unsigned int table = database->GetTable("GAFF Bond Stretching Emprirical Rules");
//   OB_ASSERT( table->NumRows() == 58 );
//   OB_ASSERT( table->NumColumns() == 4 );
  
//   // check the header
//   std::vector<std::string> header = table->GetHeader();
//   OB_ASSERT( header.size() == 4 );

//   // check the row types
//   std::vector<OBVariant::Type> types = table->GetTypes();
//   OB_ASSERT( types.size() == 4 );
//   OB_ASSERT( types.at(0) == OBVariant::Int );
//   OB_ASSERT( types.at(1) == OBVariant::Int );
//   OB_ASSERT( types.at(2) == OBVariant::Double );
//   OB_ASSERT( types.at(3) == OBVariant::Double );
  
//   // try first entry
//   std::vector<OBVariant> row = table->GetRow(0);
//   OB_ASSERT( row.size() == 4 );
//   OB_ASSERT( row.at(0).AsInt() == 1 );
//   OB_ASSERT( row.at(1).AsInt() == 6 );
//   OB_ASSERT( row.at(2).AsDouble() == 1.084 );
//   OB_ASSERT( row.at(3).AsDouble() == 5.15 );

//   // try query
//   std::vector<OBParameterDBTable::Query> query;
//   query.push_back( OBParameterDBTable::Query(0, OBVariant(6)) );
//   query.push_back( OBParameterDBTable::Query(1, OBVariant(9)) );
//   row = table->FindRow(query);
//   OB_ASSERT( row.size() == 4 );
//   OB_ASSERT( row.at(0).AsInt() == 6 );
//   OB_ASSERT( row.at(1).AsInt() == 9 );
//   OB_ASSERT( row.at(2).AsDouble() == 1.353 );
//   OB_ASSERT( row.at(3).AsDouble() == 6.20 );
// }

// void testStrBndEmpiricalRules(OBParameterDB *database)
// {
//   unsigned int table = database->GetTable("GAFF Stretch Bending Emprirical Rules");
//   OB_ASSERT( table->NumRows() == 30 );
//   OB_ASSERT( table->NumColumns() == 5 );
  
//   // check the header
//   std::vector<std::string> header = table->GetHeader();
//   OB_ASSERT( header.size() == 5 );

//   // check the row types
//   std::vector<OBVariant::Type> types = table->GetTypes();
//   OB_ASSERT( types.size() == 5 );
//   OB_ASSERT( types.at(0) == OBVariant::Int );
//   OB_ASSERT( types.at(1) == OBVariant::Int );
//   OB_ASSERT( types.at(2) == OBVariant::Int );
//   OB_ASSERT( types.at(3) == OBVariant::Double );
//   OB_ASSERT( types.at(4) == OBVariant::Double );
  
//   // try first entry
//   std::vector<OBVariant> row = table->GetRow(0);
//   OB_ASSERT( row.size() == 5 );
//   OB_ASSERT( row.at(0).AsInt() == 0 );
//   OB_ASSERT( row.at(1).AsInt() == 1 );
//   OB_ASSERT( row.at(2).AsInt() == 0 );
//   OB_ASSERT( row.at(3).AsDouble() == 0.15 );
//   OB_ASSERT( row.at(4).AsDouble() == 0.15 );

//   // try query
//   std::vector<OBParameterDBTable::Query> query;
//   query.push_back( OBParameterDBTable::Query(0, OBVariant(1)) );
//   query.push_back( OBParameterDBTable::Query(1, OBVariant(1)) );
//   query.push_back( OBParameterDBTable::Query(2, OBVariant(3)) );
//   row = table->FindRow(query);
//   OB_ASSERT( row.size() == 5 );
//   OB_ASSERT( row.at(0).AsInt() == 1 );
//   OB_ASSERT( row.at(1).AsInt() == 1 );
//   OB_ASSERT( row.at(2).AsInt() == 3 );
//   OB_ASSERT( row.at(3).AsDouble() == 0.30 );
//   OB_ASSERT( row.at(4).AsDouble() == 0.50 );
// }

// void testPartialBondChargeIncrements(OBParameterDB *database)
// {
//   unsigned int table = database->GetTable("GAFF Partial Bond Charge Increments");
//   OB_ASSERT( table->NumRows() == 98 );
//   OB_ASSERT( table->NumColumns() == 3 );
  
//   // check the header
//   std::vector<std::string> header = table->GetHeader();
//   OB_ASSERT( header.size() == 3 );

//   // check the row types
//   std::vector<OBVariant::Type> types = table->GetTypes();
//   OB_ASSERT( types.size() == 3 );
//   OB_ASSERT( types.at(0) == OBVariant::Int );
//   OB_ASSERT( types.at(1) == OBVariant::Double );
//   OB_ASSERT( types.at(2) == OBVariant::Double );
  
//   // try first entry
//   std::vector<OBVariant> row = table->GetRow(0);
//   OB_ASSERT( row.size() == 3 );
//   OB_ASSERT( row.at(0).AsInt() == 1 );
//   OB_ASSERT( row.at(1).AsDouble() == 0.000 );
//   OB_ASSERT( row.at(2).AsDouble() == 0.000 );

//   // try query
//   std::vector<OBParameterDBTable::Query> query;
//   query.push_back( OBParameterDBTable::Query(0, OBVariant(35)) );
//   row = table->FindRow(query);
//   OB_ASSERT( row.size() == 3 );
//   OB_ASSERT( row.at(0).AsInt() == 35 );
//   OB_ASSERT( row.at(1).AsDouble() == -0.456 );
//   OB_ASSERT( row.at(2).AsDouble() == 0.500 );
// }

int main()
{
  std::cout << string(TESTDATADIR) + string("../data/gaff.dat") << std::endl;
  GAFFParameterDB *database = new GAFFParameterDB(string(TESTDATADIR) + string("../data/gaff.dat"));
  testTables(database);
  testAtomProperties(database);
  testBondParameters(database);
  // testAngleParameters(database);
  // testStretchBendParameters(database);
  testTorsionParameters(database);
  // testOutOfPlaneParameters(database);
  // testVanDerWaalsParameters(database);
  // testChargeParameters(database);
  // testBondEmpiricalRules(database);
  // testStrBndEmpiricalRules(database);
  // testPartialBondChargeIncrements(database);
}
