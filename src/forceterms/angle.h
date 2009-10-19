#include <OBFunctionTerm>

namespace OpenBabel {
  namespace OBFFs {

    class AngleHarmonic : public OBFunctionTerm
    {
    public:
      struct Index
      {
	unsigned int iA, iB, iC;
      };
      struct Parameter
      {
	double K, theta0;
      };
      AngleHarmonic(OBFunction *function, std::string tableName="Angle Harmonic");
      ~AngleHarmonic();
      std::string GetName() const { return m_name;}
      bool Setup();
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      double GetValue() const { return m_value;}
    private:
      static const std::string m_name;
      const std::string m_tableName;
      unsigned int m_numAngles;
      Parameter *  m_calcs;
      Index * m_i;
      double m_value;
    };
    
  } // OBFFs
} // OpenBabel
