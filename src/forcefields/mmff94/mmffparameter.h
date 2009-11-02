/*********************************************************************
forcefieldmmff94.cpp - MMFF94 force field

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

  class MMFF94ParameterDB : public OBFFParameterDB
  {
    public:
      MMFF94ParameterDB(const std::string &filename);
      bool IsInitialized() const { return m_initialized; }
      void EnsureInit() { if (!m_initialized) ParseParamFile(); }

    private:
      bool ParseParamFile();
      bool ParseParamProp(const std::string &filename);
      bool ParseParamDef(const std::string &filename);
      bool ParseParamBond(const std::string &filename);
      bool ParseParamBndk(const std::string &filename);
      bool ParseParamAngle(const std::string &filename);
      bool ParseParamStrBnd(const std::string &filename);
      bool ParseParamDfsb(const std::string &filename);
      bool ParseParamOOP(const std::string &filename);
      bool ParseParamTorsion(const std::string &filename);
      bool ParseParamVDW(const std::string &filename);
      bool ParseParamCharge(const std::string &filename);
      bool ParseParamPbci(const std::string &filename);
      
      bool m_initialized;      
      std::string m_filename;
  };
 
}
}
