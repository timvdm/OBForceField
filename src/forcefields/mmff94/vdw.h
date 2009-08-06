#include <OBFunctionTerm>

namespace OpenBabel {
namespace OBFFs {

  class MMFF94Common;

  class MMFF94VDWTerm : public OBFunctionTerm
  {
    public:
      struct Calculation
      {
        unsigned int idx1, idx2;
        double epsilon, R_AB, R_AB7;
      };

      MMFF94VDWTerm(OBFunction *function, MMFF94Common *common);
      std::string GetName() const { return "MMFF94 Van der Waals term"; }
      bool Setup(/*const*/ OBMol &molecule);
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      double GetValue() const { return m_value; }
    protected:
      MMFF94Common *m_common;
      std::vector<Calculation> m_calcs;
      double m_value;
  };

} // OBFFs
} // OpenBabel
