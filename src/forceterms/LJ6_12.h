#include <OBFunction>
#include <OBFunctionTerm>

namespace OpenBabel {
  namespace OBFFs {

    class LJ6_12 : public OBFunctionTerm
    {
    public:
      enum MixingRule {geometric, arithmetic, sixthpower};

      struct Index
      {
	unsigned int iA, iB;
      };
      struct Parameter
      {
	double epsilon, sigma;
      };
      LJ6_12(OBFunction *function, const double factorOneFour = 0.5, const LJ6_12::MixingRule rule = geometric, const std::string tableName="LJ6_12");
      ~LJ6_12();
      std::string GetName() const { return m_name; }
      bool Setup();
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      double GetValue() const { return m_value; }
    protected:
      template <MixingRule rule>
      static void Mix(double & sigma, double & epsilon, const double & sigma_1,  const double & epsilon_1,  const double & sigma_2,  const double & epsilon_2);
    private:
      static const std::string m_name;
      const std::string m_tableName;
      unsigned int m_numPairs;
      Parameter *  m_calcs;
      Index * m_i;
      double m_value;
      void (*m_Mix)(double &, double &, const double &,  const double &,  const double &,  const double &);
      const double m_factorOneFour;
    };

  } // OBFFs
} // OpenBabel
