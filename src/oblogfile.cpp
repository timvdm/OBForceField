/**********************************************************************
oblogbase.cpp - Base class for classes that generate loggable output.

Copyright (C) 2008 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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

#include <OBLogFile>
#include <iostream>

namespace OpenBabel {
namespace OBFFs {
  

  OBLogFile::OBLogFile() : m_loglvl(High), m_logos(&std::cout)
  {
  }
  
  OBLogFile::~OBLogFile()
  {
  }

  bool OBLogFile::SetOutputStream(std::ostream* pos)
  {
    if(pos)
      m_logos = pos;
    else
      m_logos = &std::cout;
    
    return true;
  }
    
  void OBLogFile::Write(const std::string &msg)
  {
    if (!m_logos)
      return;
     
    *m_logos << msg;
  }

  void OBLogFile::Write(LogLevel lvl, const std::string &msg)
  {
    if (!m_logos || (m_loglvl < lvl))
      return;
     
    *m_logos << msg;
  }


}
} // end namespace OpenBabel


//! @file oblogbase.cpp
//! @brief 
