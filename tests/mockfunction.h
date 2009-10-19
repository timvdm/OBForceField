#include <OBFunction>

namespace OpenBabel {
  namespace OBFFs {

    class MockFunction : public OBFunction
    {
      public:
        MockFunction(unsigned int numParticles) : OBFunction()
        {
          m_positions.resize(numParticles, Eigen::Vector3d::Zero());
          m_gradients.resize(numParticles, Eigen::Vector3d::Zero());
        }
        std::string GetName() const
        {
          return "MockFunction";
        }
        void Compute(Computation computation = Value)
        {
        }
        double GetValue() const
        {
          return 0.0;
        }
        std::string GetUnit() const
        {
          return "kJ/mol";
        }
        void ProcessOptions(std::vector<Option> &options)
        {
        }
        std::string GetDefaultOptions() const
        {
          return "";
        }
    };


  }
}
