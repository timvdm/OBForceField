#include <OBFunctionTerm>

namespace OpenBabel {
  namespace OBFFs {

    class BondHarmonic : public OBFunctionTerm
    {
    public:
      struct Index
      {
	unsigned int iA, iB;
      };
      struct Parameter
      {
	double K, r0;
      };
      BondHarmonic(OBFunction *function, std::string tableName="Bond Harmonic");
      ~BondHarmonic();
      std::string GetName() const { return m_name; }
      bool Setup();
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      double GetValue() const { return m_value; }
    private:
      static const std::string m_name;
      const std::string m_tableName;
      unsigned int m_numBonds;
      Parameter *  m_calcs;
      Index * m_i;
      double m_value;
    };

    class BondClass2 : public OBFunctionTerm
    {
    public:
      struct Index
      {
	unsigned int iA, iB;
      };
      struct Parameter
      {
	double K2, K3, K4, r0;
      };
      BondClass2(OBFunction *function, std::string tableName="Bond Class 2");
      ~BondClass2();
      std::string GetName() const { return m_name; }
      bool Setup();
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      double GetValue() const { return m_value; }
    private:
      static const std::string m_name;
      const std::string m_tableName;
      unsigned int m_numBonds;
      Parameter *  m_calcs;
      Index * m_i;
      double m_value;
    };
    
  } // OBFFs
} // OpenBabel
