/**********************************************************************
  OBVariant - Class to store variables (i.e. bool, int, double)

  Copyright (C) 2009 Tim Vandermeersch
  Some contributions copyright (C) 2009 Frank Peters 

  This file is part of the OpenBabel project.
  For more information, see <http://openbabel.openmolecules.net/>

  Avogadro is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Avogadro is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/

#ifndef OBFFS_VARIANT_H
#define OBFFS_VARIANT_H

#include <string>
#include <typeinfo>
#include <iostream>

namespace OpenBabel {
  namespace OBFFs {
    
    class OBVariant
    {
    public:
      enum Type {Int, Double, Bool, String};
    
      OBVariant(const int &value, const std::string &name = std::string()) 
	: m_name(name), p_value(new holder<int>(value)), m_type(Int) {}
      OBVariant(const double &value, const std::string &name = std::string()) 
	: m_name(name), p_value(new holder<double>(value)), m_type(Double) {}
      OBVariant(const bool &value, const std::string &name = std::string()) 
	: m_name(name), p_value(new holder<bool>(value)), m_type(Bool) {}
      OBVariant(const std::string &value, const std::string &name = std::string()) 
	: m_name(name), p_value(new holder<std::string>(value)), m_type(String) {}
      OBVariant(const char * value, const std::string &name = std::string()) 
	: m_name(name), p_value(new holder<std::string>(value)), m_type(String) {}

      OBVariant(const OBVariant& rhs)
	: m_name(rhs.m_name), m_type(rhs.m_type), p_value(rhs.p_value) {
	++(p_value->references);
      }

      OBVariant & operator= (const OBVariant& rhs){
	if (this == &rhs)
	  return *this;
	else {
	  if (--(p_value->references)==0) { 
	    delete p_value;
	  }
	  m_name = rhs.m_name; 
	  m_type = rhs.m_type;
	  p_value = rhs.p_value;
	  ++(p_value->references);
	  return *this;
	}
      }
      
      ~OBVariant() 
      {
	if (--(p_value->references)==0) { 
	  delete p_value;
	}
      }

      const std::string& GetName() const { return m_name;}
      Type GetType() const { return m_type;}
    
      int AsInt() const;
      double AsDouble() const;
      bool AsBool() const;
      std::string AsString() const;
    
      bool operator==(const OBVariant &other) const;
      bool operator!=(const OBVariant &other) const;
    
      class placeholder 
      { 
      public:
	placeholder() : references(1) {}
	virtual ~placeholder() {} 
	unsigned int references;
      };

      template<typename ValueType>
      class holder : public placeholder
      {
      public:
	holder(const ValueType & value)
	  : m_value(value) {}
	ValueType m_value;
      };

    private:
      template<typename T>
      T AsT() const;
      std::string m_name;
      Type m_type;
      placeholder * p_value;
    };
    
  } // OBFFs
} // OpenBabel

#endif
