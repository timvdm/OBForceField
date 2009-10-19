/**********************************************************************

 
Copyright (C) 2006 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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
#include <OBChargeMethod>

namespace OpenBabel {

  class OBMol;
  class OBChargeMethod;

  namespace OBFFs {
  
    class OBGasteiger : public OBChargeMethod
    {
    public:
      bool ComputeCharges(OBMol & mol);
    }; 
  }
}// namespace OpenBabel


