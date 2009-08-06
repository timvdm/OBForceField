/**********************************************************************
  OBVariant - Class to store variables (i.e. bool, int, double)

  Copyright (C) 2009 Tim Vandermeersch

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

namespace OpenBabel {
namespace OBFFs {

  class OBVariant
  {
    public:
      enum Type {Int, Double, Bool};

      OBVariant(int value, const std::string &name = std::string()) 
        : m_name(name), m_int(value), m_type(Int) {}
      OBVariant(double value, const std::string &name = std::string()) 
        : m_name(name), m_double(value), m_type(Double) {}
      OBVariant(bool value, const std::string &name = std::string()) 
        : m_name(name), m_bool(value), m_type(Bool) {}
      
      const std::string& GetName() const { return m_name; }
      const Type GetType() const { return m_type; }
      
      int AsInt() const;
      double AsDouble() const;
      bool AsBool() const;
      std::string AsString() const;

      bool operator==(const OBVariant &other);
      bool operator!=(const OBVariant &other);

    private:
      template<typename T> T AsT() const;
 
      std::string m_name;
      Type m_type;
      union {
        int    m_int;
        double m_double;
        bool   m_bool;
      };
  };

} // OBFFs
} // OpenBabel

#endif
