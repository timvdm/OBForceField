#include <OBFunctionTerm>

namespace OpenBabel {
  namespace OBFFs {

    class TorsionHarmonic : public OBFunctionTerm
    {
    public:
      struct Index
      {
	unsigned int iA, iB, iC, iD;
      };
      struct Parameter
      {
	double K, d, n;
      };
      TorsionHarmonic(OBFunction *function, std::string tableName="Torsion Harmonic");
      ~TorsionHarmonic();
      std::string GetName() const { return m_name;}
      bool Setup();
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      double GetValue() const { return m_value;}
    private:
      static const std::string m_name;
      const std::string m_tableName;
      unsigned int m_numTorsions;
      Parameter *  m_calcs;
      Index * m_i;
      double m_value;
    };
    
  } // OBFFs
} // OpenBabel
