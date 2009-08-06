/**********************************************************************
obfunctionterm.cpp - 
 
Copyright (C) 2009 by Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

***********************************************************************/

#include <OBFunctionTerm>
#include <OBFunction>

namespace OpenBabel {
namespace OBFFs {
 
  OBFunctionTerm::OBFunctionTerm(OBFunction *func) : m_function(func)
  {
  }

  OBFunctionTerm::~OBFunctionTerm()
  {}

}
} // end namespace OpenBabel


//! @file obfunctionterm.cpp
//! @brief 
