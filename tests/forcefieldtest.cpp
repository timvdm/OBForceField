#include <OBFunctionTerm>
#include <OBLogFile>


#include <openbabel/atom.h>
#include <openbabel/mol.h>

using namespace OpenBabel;
using namespace OpenBabel::OBFFs;

class MMFF94
{

  template<typename T>
  static OBFunctionTerm::Parameter FindParameter1(T value, const OBFunctionTerm::Parameters &parameter, int column = 0)
  {
    OBFunctionTerm::Parameter par;

    for (unsigned int row = 0; row < parameter.size(); row++)
      if (value == parameter.at(row).at(column)) {
        return parameter.at(row);
      }

    return 0;
  }
 
 /* 
  OBFFParameter* OBForceFieldMMFF94::GetParameter2Atom(int a, int b, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b)) || 
          ((a == parameter[idx].b) && (b == parameter[idx].a))) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }
  
  OBFFParameter* OBForceFieldMMFF94::GetParameter3Atom(int a, int b, int c, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;

    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c)) || 
          ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a))) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetTypedParameter2Atom(int ffclass, int a, int b, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
      
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (ffclass == parameter[idx]._ipar[0])) || 
          ((a == parameter[idx].b) && (b == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0]))) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetTypedParameter3Atom(int ffclass, int a, int b, int c, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && (ffclass == parameter[idx]._ipar[0])) || 
          ((a == parameter[idx].c) && (b == parameter[idx].b) && (c == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0])) ) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

  OBFFParameter* OBForceFieldMMFF94::GetTypedParameter4Atom(int ffclass, int a, int b, int c, int d, std::vector<OBFFParameter> &parameter)
  {
    OBFFParameter *par;
    
    for (unsigned int idx=0; idx < parameter.size(); idx++)
      if (((a == parameter[idx].a) && (b == parameter[idx].b) && (c == parameter[idx].c) && 
           (d == parameter[idx].d) && (ffclass == parameter[idx]._ipar[0])) 
          // || ((a == parameter[idx].d) && (b == parameter[idx].c) && (c == parameter[idx].b) && 
          //   (d == parameter[idx].a) && (ffclass == parameter[idx]._ipar[0])) 
          ) 
        {
          par = &parameter[idx];
          return par;
        }

    return NULL;
  }

*/

  //
  // MMFF part V - page 620
  //
  // BTij is 1 when:
  // a) single bond between atoms i and j, both i and j are not aromatic and both types have sbmb set in mmffprop.par, or
  // b) bewtween two aromatic atoms, but the bond is not aromatic (e.g. connecting bond in biphenyl)
  //
  int OBForceFieldMMFF94::GetBondType(OBAtom* a, OBAtom* b)
  {
    /*
    OBMol *mol = GetMolecule();
    
    if (!mol->GetBond(a,b)->IsSingle())
      return 0;
    
    if (!mol->GetBond(a,b)->IsAromatic())
        */
    OBBond *bond = a->GetBond(b);

    if (bond->IsSingle())
      return 0;

    if (bond->IsAromatic())
      if (HasAromSet(atoi(a->GetType())) && HasAromSet(atoi(b->GetType())))
        return 1;
      
    if (HasSbmbSet(atoi(a->GetType())) && HasSbmbSet(atoi(b->GetType())))
      return 1;
    
    return 0;
  }
 
};

class MMFF94BondTerm : public OBFunctionTerm
{
  public:
    struct Calculation
    {
      unsigned int idx1, idx2;
      double kb, r0;
      char bondType;
    };

    MMFF94BondTerm(OBFunction *func) : OBFunctionTerm(func)
    {
    }
    
    std::string GetName() const 
    { 
      return "MMFF94 bond stretching"; 
    }
    
    int NumParameterRows() const 
    { 
      return m_parameters.size(); 
    }
    
    int NumParameterColumns() const 
    { 
      return 5; 
    }
    
    std::vector<OBVariant::Type> GetParameterTypes() const
    {
      std::vector<OBVariant::Type> types;
      types.push_back(OBVariant::Int); // bond type
      types.push_back(OBVariant::Int); // atom type 1
      types.push_back(OBVariant::Int); // atom type 2
      types.push_back(OBVariant::Double); // kb
      types.push_back(OBVariant::Double); // r0
      return types;      
    }
    
    std::vector<std::string> GetParameterHeaders() const
    {
      std::vector<std::string> headers;
      headers.push_back("bond type");
      headers.push_back("atom type 1");
      headers.push_back("atom type 2");
      headers.push_back("kb");
      headers.push_back("r0");
      return headers;      
    }
    
    double Compute(OBFunction::Computation computation = OBFunction::Value) const
    {
   
    }

    bool SetMolecule(OBMol *mol)
    {
      OBLogFile *logfile = m_function->GetLogFile();
      OBAtom *a, *b;
    
      if (logfile->IsLow())
        logfile->Write("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");
    
      // 
      // Bond Calculations
      //
      // no "step-down" procedure
      // MMFF part V - page 625 (empirical rule)
      //
      if (logfile->IsLow())
        logfile->Write("SETTING UP BOND CALCULATIONS...\n");
 
      Calculation *bondcalc;
      int bondtype;

      m_calculations.clear();
    
      FOR_BONDS_OF_MOL(bond, mol) {
        a = bond->GetBeginAtom();
        b = bond->GetEndAtom();
     
        bondtype = GetBondType(a, b);

        parameter = GetTypedParameter2Atom(bondtype, atoi(a->GetType()), atoi(b->GetType()), _ffbondparams); // from mmffbond.par
        if (parameter == NULL) {
          parameter = GetParameter2Atom(a->GetAtomicNum(), b->GetAtomicNum(), _ffbndkparams); // from mmffbndk.par - emperical rules
          if (parameter == NULL) { 
            if (logfile->IsLow()) {
              // This should never happen
              snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR BOND %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
              logfile->Write(_logbuf);
            }
            return false;
          } else {
            if (logfile->IsLow()) {
              snprintf(_logbuf, BUFF_SIZE, "   USING EMPIRICAL RULE FOR BOND STRETCHING %d-%d (IDX)...\n", a->GetIdx(), b->GetIdx());
              logfile->Write(_logbuf);
            }

            double rr, rr2, rr4, rr6;
            bondcalc.a = a;
            bondcalc.b = b;
            bondcalc.r0 = GetRuleBondLength(a, b); 

            rr = parameter->_dpar[0] / bondcalc.r0; // parameter->_dpar[0]  = r0-ref
            rr2 = rr * rr;
            rr4 = rr2 * rr2;
            rr6 = rr4 * rr2;

            bondcalc.kb = parameter->_dpar[1] * rr6; // parameter->_dpar[1]  = kb-ref
            bondcalc.bt = bondtype;

            _bondcalculations.push_back(bondcalc);
          }
        } else {
          bondcalc.a = a;
          bondcalc.b = b;
          bondcalc.kb = parameter->_dpar[0];
          bondcalc.r0 = parameter->_dpar[1];
          bondcalc.bt = bondtype;

          _bondcalculations.push_back(bondcalc);
        }
      }

    }



};


int main()
{



}
