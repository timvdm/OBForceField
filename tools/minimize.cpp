#include <OBFunction>
#include <OBLogFile>
#include <OBMinimize>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using OpenBabel::OBMol;
using OpenBabel::OBConversion;
using OpenBabel::OBFormat;

using namespace OpenBabel::OBFFs;
using namespace std;


int main(int argc, char **argv)
{
  OBFunctionFactory *factory = OBFunctionFactory::GetFactory("MMFF94");
  OBFunction *function = factory->NewInstance();
  if (!function) {
    cout << "ERROR: could not find MMFF94 function" << endl;
    return -1;
  }

  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <filename>" << endl;
    return -1;    
  }

  OBMol mol;
  OBConversion conv;
  OBFormat *format = conv.FormatFromExt(argv[1]);
  if (!format || !conv.SetInFormat(format)) {
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

  // read options file
  if (argc == 3) {
    std::ifstream cifs;
    cifs.open(argv[2]);
    std::stringstream options;
    std::string line;
    while (std::getline(cifs, line))
      options << line << std::endl;

    function->SetOptions(options.str());
  }


  function->Setup(mol);

  OBMinimize minimize(function);

  minimize.SteepestDescent(50);


}
