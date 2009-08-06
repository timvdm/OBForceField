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

#include <OBParameterDB>

namespace OpenBabel {
namespace OBFFs {

  class MMFF94SimpleParameterDB : public OBParameterDB
  {
    public:
      enum Tables {
        AtomProperties = 0,
        AtomTypeLevels,
        BondParameters,
        AngleParameters,
        StretchBendParameters,
        TorsionParameters,
        OutOfPlaneParameters,
        VanDerWaalsParameters,
        ChargeParameters,
        EmpiricalBondRules,
        EmpiricalStretchBendRules,
        PartialBondChargeIncrements,
        LastEmptyTable
      };

      MMFF94SimpleParameterDB(const std::string &filename);
      unsigned int NumTables() const;
      std::vector<std::string> GetTables() const;
      unsigned int GetTable(const std::string &table) const;
      unsigned int NumRows(unsigned int table) const;
      unsigned int NumColumns(unsigned int table) const;
      std::vector<std::string> GetHeader(unsigned int table) const;
      std::vector<OBVariant::Type> GetTypes(unsigned int table) const;
      bool VerifyTypes(unsigned int table, const std::vector<OBVariant> &values);
      const std::vector<OBVariant>& GetRow(unsigned int table, unsigned int row) const;
      bool AddRow(unsigned int table, const std::vector<OBVariant> &values);
      const std::vector<OBVariant>& FindRow(unsigned int table, const std::vector<Query> &query, bool *swapped = 0) const;
      const std::vector<std::vector<OBVariant> >& GetAllRows(unsigned int table) const;

    private:
      bool IsInitialized() const { return m_initialized; }
      void EnsureInit() { if (!m_initialized) ParseParamFile(); }

      bool ParseParamFile();
      bool ParseParamProp(std::string &filename);
      bool ParseParamDef(std::string &filename);
      bool ParseParamBond(std::string &filename);
      bool ParseParamBndk(std::string &filename);
      bool ParseParamAngle(std::string &filename);
      bool ParseParamStrBnd(std::string &filename);
      bool ParseParamDfsb(std::string &filename);
      bool ParseParamOOP(std::string &filename);
      bool ParseParamTorsion(std::string &filename);
      bool ParseParamVDW(std::string &filename);
      bool ParseParamCharge(std::string &filename);
      bool ParseParamPbci(std::string &filename);
 
      
      bool m_initialized;      
      std::string m_filename;
      std::vector< std::vector<std::vector<OBVariant> >* > m_tables;
      const std::vector<OBVariant> emptyRow;
      
  };
 
}
}
