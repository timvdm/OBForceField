#include <OBFunctionTerm>

namespace OpenBabel {
namespace OBFFs {

  class MMFF94Common;

  class MMFF94StrBndTerm : public OBFunctionTerm
  {
    public:
      struct Calculation
      {
        unsigned int idx1, idx2, idx3;
        double kba123, kba321, theta0, rab0, rbc0;
        char strbndType;
      };

      MMFF94StrBndTerm(OBFunction *function, MMFF94Common *common);
      std::string GetName() const { return "MMFF94 stretch bending term"; }
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
