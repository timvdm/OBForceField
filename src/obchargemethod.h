/**********************************************************************
 
Copyright (C) 2009 by Frank Peters
 
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
#include <vector>

namespace OpenBabel {

  class OBMol;

  namespace OBFFs {
  
    class OBChargeMethod
    {
    public:
      void CopyFromMol(OBMol & mol);
      bool CopyToMol(OBMol & mol) const;
      const std::vector<double> & GetFormalCharges() const;
      const std::vector<double> & GetPartialCharges() const;
      virtual bool ComputeCharges(OBMol & mol);
    protected:
      std::vector<double> m_partialCharges;
      std::vector<double> m_formalCharges;
    }; 
  }
}// namespace OpenBabel


