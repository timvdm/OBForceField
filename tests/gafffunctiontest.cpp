#include <OBFunction>
#include <OBLogFile>
#include "obtest.h"
#include <GAFF>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using OpenBabel::OBMol;
using OpenBabel::OBConversion;

using namespace OpenBabel::OBFFs;

using namespace std;


int main()
{
  OBFunctionFactory *gaff_factory = OBFunctionFactory::GetFactory("GAFF");
  OB_ASSERT( gaff_factory != 0);

  OBFunction *gaff_function = gaff_factory->NewInstance();
  OB_ASSERT( gaff_function != 0);

  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("mol");

  std::ifstream ifs;
  //ifs.open("hexane.xyz");
  ifs.open("acetone.mol");
  conv.Read(&mol, &ifs);
  ifs.close();

  cout << "num atoms = " << mol.NumAtoms() << endl;

  gaff_function->GetLogFile()->SetOutputStream(&std::cout);

  GAFFParameterDB gaff_parameterDB("../data/gaff.dat");
  GAFFTypeRules gaff_typerules("../data/gaff.prm");
  GAFFType gaff_type(& gaff_typerules);
  OBGasteiger chargeMethod;

  gaff_function->SetParameterDB(& gaff_parameterDB);
  gaff_function->SetOBFFType(& gaff_type);
  gaff_function->SetOBChargeMethod(& chargeMethod);

  gaff_function->Setup(mol);
  gaff_function->Compute();

  cout << "E = " << gaff_function->GetValue() << endl;
  
  cout << "Options:" << endl;
  cout << gaff_function->GetOptions() << endl;
 
  std::stringstream ss;
  ss << "vdwterm = allpair";
  gaff_function->SetOptions(ss.str());


}
