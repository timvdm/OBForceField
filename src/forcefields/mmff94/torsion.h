#include <OBFunctionTerm>

namespace OpenBabel {
namespace OBFFs {

  class MMFF94Common;

  class MMFF94TorsionTerm : public OBFunctionTerm
  {
    public:
      struct Calculation
      {
        unsigned int idx1, idx2, idx3, idx4;
        double v1, v2, v3;
        int torsionType;
      };

      MMFF94TorsionTerm(OBFunction *function, MMFF94Common *common);
      std::string GetName() const { return "MMFF94 torsion term"; }
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
