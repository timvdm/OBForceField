#include <OBVariant>
#include <sstream>
#include <stdlib.h>

namespace OpenBabel {
namespace OBFFs {

  template<typename T>
  T OBVariant::AsT() const
  {
    switch (m_type) {
      case Int:
        return ( static_cast< holder<int> *>(p_value) -> m_value) ;
      case Double:
        return ( static_cast< holder<double> * >(p_value) -> m_value) ;
      case Bool:
        return ( static_cast< holder<bool> * >(p_value) -> m_value) ;
    }
  }

  int OBVariant::AsInt() const
  {
    if (m_type==String)
      return  atoi((static_cast< holder<std::string> * >(p_value) -> m_value).c_str()) ;
    else
      return AsT<int>();
  }

  double OBVariant::AsDouble() const
  {
    if (m_type==String)
      return  strtod((static_cast< holder<std::string> * >(p_value) -> m_value).c_str(),0) ;
    else
      return AsT<double>();
  }

  bool OBVariant::AsBool() const
  {
    if (m_type==String){
      std::string s(static_cast< holder<std::string> * >(p_value) -> m_value);
      return (s=="true"||s=="True"||s=="1"||s=="t") ;
    }
    else
      return AsT<bool>();
  }

  std::string OBVariant::AsString() const
  {
    std::stringstream ss;
    switch (m_type) {
      case Int:
        ss << ( static_cast< holder<int> * >(p_value) -> m_value);
        break;
      case Double:
        ss << ( static_cast< holder<double> * >(p_value) -> m_value);
        break;
      case Bool:
        if ( static_cast< holder<bool> * >(p_value) -> m_value == true)
          ss << "True";
        else 
          ss << "False";
        break;
      case String:
          ss << ( static_cast< holder<std::string> * >(p_value) -> m_value);
        break;
    }
    return ss.str();
  }
      
  bool OBVariant::operator==(const OBVariant &other) const
  {
    if (m_type != other.m_type)
      return false;
    switch (m_type) {
      case Int:
        return (static_cast< holder<int> * >(p_value) -> m_value ==  static_cast< holder<int> * >(other.p_value) -> m_value);
      case Double:
        return (static_cast< holder<double> * >(p_value) -> m_value ==  static_cast< holder<double> * >(other.p_value) -> m_value);
      case Bool:
        return (static_cast< holder<bool> * >(p_value) -> m_value ==  static_cast< holder<bool> * >(other.p_value) -> m_value);
      case String:
        return (static_cast< holder<std::string> * >(p_value) -> m_value ==  static_cast< holder<std::string> * >(other.p_value) -> m_value);
    }
  }
  
  bool OBVariant::operator!=(const OBVariant &other) const
  {
    return !(*this == other);
  }

} // OBFFs
} // OpenBabel

