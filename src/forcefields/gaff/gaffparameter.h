/*********************************************************************
GAFF force field

Copyright (C) 2006-2009 by Tim Vandermeersch
 
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

#include <OBFFParameterDB>

namespace OpenBabel {
namespace OBFFs {
    
  class GAFFParameterDB : public OBFFParameterDB
  {
    public:
      GAFFParameterDB(const std::string &filename);      
      bool IsInitialized() const { return _initialized; }
      void EnsureInit() { if (!_initialized) ParseParamFile(); }
    private:
      bool ParseParamFile();
      bool _initialized;      
      std::string _filename;
  };
  
} // namespace OBFFs
} // namespace OpenBabel
