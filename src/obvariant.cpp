#include <OBVariant>
#include <sstream>

namespace OpenBabel {
namespace OBFFs {

  template<typename T>
  T OBVariant::AsT() const
  {
    switch (m_type) {
      case Int:
        return m_int;
      case Double:
        return m_double;
      case Bool:
        return m_bool;
    }
  }

  int OBVariant::AsInt() const
  {
    return AsT<int>();
  }

  double OBVariant::AsDouble() const
  {
    return AsT<double>();
  }

  bool OBVariant::AsBool() const
  {
    return AsT<bool>();
  }

  std::string OBVariant::AsString() const
  {
    std::stringstream ss;
    switch (m_type) {
      case Int:
        ss << m_int;
        break;
      case Double:
        ss << m_double;
        break;
      case Bool:
        if (m_bool == true)
          ss << "True";
        else 
          ss << "False";
        break;
    }

    return ss.str();
  }
      
  bool OBVariant::operator==(const OBVariant &other)
  {
    if (m_type != other.m_type)
      return false;
    switch (m_type) {
      case Int:
        return (m_int == other.m_int);
      case Double:
        return (m_double == other.m_double);
      case Bool:
        return (m_bool == other.m_bool);
      default:
        return false;  
    }
  }
  
  bool OBVariant::operator!=(const OBVariant &other)
  {
    return !(*this == other);
  }

} // OBFFs
} // OpenBabel

