#include <OBFunction>
#include <OBFunctionTerm>

namespace OpenBabel {
  namespace OBFFs {

    class Coulomb : public OBFunctionTerm
    {
    public:
      struct Index
      {
	unsigned int iA, iB;
      };
      struct Parameter
      {
	double qq;
      };
      Coulomb(OBFunction *function, const double factorOneFour = 0.8333, const double relativePermittivity = 1.0);
      ~Coulomb();
      std::string GetName() const { return m_name; }
      bool Setup();
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      double GetValue() const { return m_value; }
    private:
      static const std::string m_name;
      unsigned int m_numPairs;
      Parameter *  m_calcs;
      Index * m_i;
      double m_value;
      const double m_relativePermittivity;
      const double m_factorOneFour;
    };

  } // OBFFs
} // OpenBabel
