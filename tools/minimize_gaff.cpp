#include <OBFunction>
#include <OBLogFile>
#include <OBMinimize>
//#include <GAFF>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/obconversion.h>

using OpenBabel::OBMol;
using OpenBabel::OBAtom;
using OpenBabel::OBConversion;
using OpenBabel::OBFormat;

using namespace OpenBabel::OBFFs;
using namespace std;


int main(int argc, char **argv)
{
  OBFunctionFactory *factory = OBFunctionFactory::GetFactory("GAFF");
  OBFunction *function = factory->NewInstance();
  if (!function) {
    cout << "ERROR: could not find GAFF function" << endl;
    return -1;
  }

  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <filename in> <filename out> <filename options>" << endl;
    return -1;    
  }

  OBMol mol;
  OBConversion conv;
  OBFormat *format_in = conv.FormatFromExt(argv[1]);
  if (!format_in || !conv.SetInFormat(format_in)) {
    cout << "ERROR: could not find format for file " << argv[1] << endl;
    return -1;
  }

  std::ifstream ifs;
  ifs.open(argv[1]);
  conv.Read(&mol, &ifs);
  ifs.close();

  cout << "# atoms = " << mol.NumAtoms() << endl;
  function->GetLogFile()->SetOutputStream(&std::cout);
  function->GetLogFile()->SetLogLevel(OBLogFile::Low);

  bool toCout;
  OBFormat *format_out;
  if (argc == 2) {
    format_out=format_in;
    toCout=true;
  } else {
    format_out = conv.FormatFromExt(argv[2]);
    toCout=false;
  }

  if (!format_out || !conv.SetOutFormat(format_out)) {
    cout << "ERROR: could not find format for file " << argv[2] << endl;
    return -1;
  }
  
  // read options file
  if (argc == 4) {
    std::ifstream cifs;
    cifs.open(argv[3]);
    std::stringstream options;
    std::string line;
    while (std::getline(cifs, line))
      options << line << std::endl;
    function->SetOptions(options.str());
  }

//   GAFFParameterDB gaff_parameterDB("../data/gaff.dat");
//   GAFFTypeRules gaff_typerules("../data/gaff.prm");
//   GAFFType gaff_type(& gaff_typerules);
//   OBGasteiger chargeMethod;
//   function->SetParameterDB(& gaff_parameterDB);
//   function->SetOBFFType(& gaff_type);
//   function->SetOBChargeMethod(& chargeMethod);

  function->Setup(mol);

  OBMinimize minimize(function);

  minimize.ConjugateGradients(50);

  function->CopyPositionsToMol(mol);

  if (toCout){
    cout << endl;
    conv.Write(&mol, &cout);
  }
  else {
    std::ofstream ofs;
    ofs.open(argv[2]);
    conv.Write(&mol, &ofs);
    ofs.close();
  }

  return 1;
}
