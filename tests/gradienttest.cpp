#include <OBFunction>
#include <OBLogFile>
#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

using OpenBabel::OBMol;
using OpenBabel::OBConversion;

using namespace OpenBabel::OBFFs;

using namespace std;

Eigen::Vector3d ValidateGradientError(const Eigen::Vector3d &numgrad, const Eigen::Vector3d &anagrad)
{
  double errx, erry, errz;

  if (fabs(numgrad.x()) < 1.0)
    errx = numgrad.x() * fabs(numgrad.x() - anagrad.x()) * 100;
  else
    errx = fabs((numgrad.x() - anagrad.x()) / numgrad.x()) * 100;

  if (fabs(numgrad.y()) < 1.0)
    erry = numgrad.y() * fabs(numgrad.y() - anagrad.y()) * 100;
  else
    erry = fabs((numgrad.y() - anagrad.y()) / numgrad.y()) * 100;

  if (fabs(numgrad.z()) < 1.0)
    errz = numgrad.z() * fabs(numgrad.z() - anagrad.z()) * 100;
  else
    errz = fabs((numgrad.z() - anagrad.z()) / numgrad.z()) * 100;

  errx = fabs(errx);
  erry = fabs(erry);
  errz = fabs(errz);

  return Eigen::Vector3d(errx, erry, errz);
}

bool ValidateGradients(OBFunction *function)
{
  Eigen::Vector3d numgrad, anagrad, err;

  bool passed = true; // set to false if any component fails

  unsigned int numAtoms = function->NumParticles();

  cout << "----------------------------------------------------------------------------------------" << endl;
  cout << "                                                                                        " << endl;
  cout << "  VALIDATE GRADIENTS                                                                    " << endl;
  cout << "                                                                                        " << endl;
  cout << "ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERROR (%)    " << endl;
  cout << "----------------------------------------------------------------------------------------" << endl;
  //     "XX       (000.000, 000.000, 000.000)  (000.000, 000.000, 000.000)  (00.00, 00.00, 00.00)"

  for (unsigned int i = 0; i < numAtoms; ++i) {

    // OBFF_ENERGY
    numgrad = function->NumericalDerivative(i);
    function->Compute(OBFunction::Gradients);
    anagrad = function->GetGradients()[i];
    err = ValidateGradientError(numgrad, anagrad);

    char _logbuf[1001];
    snprintf(_logbuf, 1000, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", i, numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    cout << _logbuf;

    /*
    // OBFF_EBOND
    numgrad = NumericalDerivative(&*a, OBFF_EBOND);
    ClearGradients();
    E_Bond(); // compute
    anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
    err = ValidateGradientError(numgrad, anagrad);

    snprintf(_logbuf, BUFF_SIZE, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    OBFFLog(_logbuf);
    if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
      passed = false;

    // OBFF_EANGLE
    numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
    ClearGradients();
    E_Angle(); // compute
    anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
    err = ValidateGradientError(numgrad, anagrad);

    snprintf(_logbuf, BUFF_SIZE, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    OBFFLog(_logbuf);
    if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
      passed = false;

    // OBFF_ESTRBND
    numgrad = NumericalDerivative(&*a, OBFF_ESTRBND);
    ClearGradients();
    E_StrBnd(); // compute
    anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
    err = ValidateGradientError(numgrad, anagrad);

    snprintf(_logbuf, BUFF_SIZE, "    strbnd  (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    OBFFLog(_logbuf);
    if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
      passed = false;

    // OBFF_ETORSION
    numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
    ClearGradients();
    E_Torsion(); // compute
    anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
    err = ValidateGradientError(numgrad, anagrad);

    snprintf(_logbuf, BUFF_SIZE, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    OBFFLog(_logbuf);
    if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
      passed = false;

    // OBFF_EOOP
    numgrad = NumericalDerivative(&*a, OBFF_EOOP);
    ClearGradients();
    E_OOP(); // compute
    anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
    err = ValidateGradientError(numgrad, anagrad);

    snprintf(_logbuf, BUFF_SIZE, "    oop     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    OBFFLog(_logbuf);
    // disable OOP gradient validation for now -- some small errors, but nothing major
    //      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
    //        passed = false;

    // OBFF_EVDW
    numgrad = NumericalDerivative(&*a, OBFF_EVDW);
    ClearGradients();
    E_VDW(); // compute
    anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
    err = ValidateGradientError(numgrad, anagrad);

    snprintf(_logbuf, BUFF_SIZE, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    OBFFLog(_logbuf);
    if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
      passed = false;

    // OBFF_EELECTROSTATIC
    numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
    ClearGradients();
    E_Electrostatic(); // compute
    anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
    err = ValidateGradientError(numgrad, anagrad);

    snprintf(_logbuf, BUFF_SIZE, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
        anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
    OBFFLog(_logbuf);
    if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
      passed = false;
    */
  }

  return passed; // did we pass every single component?
}


int main()
{
  OBFunctionFactory *mmff94_factory = OBFunctionFactory::GetFactory("MMFF94");
  OB_ASSERT( mmff94_factory != 0);

  OBFunction *mmff94_function = mmff94_factory->NewInstance();
  OB_ASSERT( mmff94_function != 0);

  OBMol mol;
  OBConversion conv;
  conv.SetInFormat("sdf");

  std::ifstream ifs;
  //ifs.open("hexane.xyz");
  ifs.open("aceton.sdf");
  conv.Read(&mol, &ifs);
  ifs.close();

  cout << "num atoms = " << mol.NumAtoms() << endl;

  mmff94_function->GetLogFile()->SetOutputStream(&std::cout);
  mmff94_function->GetLogFile()->SetLogLevel(OBLogFile::Low);



  std::stringstream options;
  options << "bonded = none" << std::endl;
  //options << "vdwterm = none" << std::endl;
  options << "vdwterm = opencl" << std::endl;
  options << "electroterm = none" << std::endl;
  //options << "electroterm = opencl" << std::endl;

  mmff94_function->SetOptions(options.str());
  mmff94_function->Setup(mol);



  ValidateGradients(mmff94_function);
/*
  

  cout << "Options:" << endl;
  cout << mmff94_function->GetOptions() << endl;
 
  std::stringstream ss;
  ss << "mmff_vdw = allpair";
  mmff94_function->SetOptions(ss.str());

*/
}
