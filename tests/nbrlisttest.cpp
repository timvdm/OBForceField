/**********************************************************************
  Molecule - NeighborListTest class provides unit testing for the Molecule class

  Copyright (C) 2008 Marcus D. Hanwell

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.openmolecules.net/>

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

#include <OBNbrList>
#include <OBFunction>

#include "obtest.h"
#include "mockfunction.h"

using namespace OpenBabel::OBFFs;

unsigned int test(OBFunction *function, int n, double r)
{
  OBNbrList *nbrList = new OBNbrList(function, r, false, n);
  
  unsigned int count = 0;
  for (unsigned int i = 0; i < function->NumParticles(); ++i) {
    std::vector<unsigned int> nbrs = nbrList->GetNbrs(i);
    count += nbrs.size();
  }

  return count;
}



int main()
{
  MockFunction *function = new MockFunction(1000);

  // create a 10x10x10 regular grid with atoms
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k) {
        function->GetPositions()[i * 100 + j * 10 + k] = Eigen::Vector3d((double)i, (double)j, (double)k);
      }

  // compute the correct number of pairs
  unsigned int correct5 = 0;
  unsigned int correct10 = 0;
  for (unsigned int i = 0; i < 1000; ++i) {
    for (unsigned int j = 0; j < 1000; ++j) {
      if (i >= j)
        continue;
      
      double r2 = ( function->GetPositions()[i] - function->GetPositions()[j] ).squaredNorm();

      if (r2 <= 25.0)
        correct5++;
      if (r2 <= 100.0)
        correct10++;
    }
  }
  
  unsigned int count;

  count = test(function, 1, 5.);
  OB_ASSERT(correct5 == count);

  count = test(function, 2, 5.);
  OB_ASSERT(correct5 == count);
  
  count = test(function, 3, 5.);
  OB_ASSERT(correct5 == count);
  
  count = test(function, 1, 10.);
  OB_ASSERT(correct10 == count);
  
  count = test(function, 2, 10.);
  OB_ASSERT(correct10 == count);
  
  count = test(function, 3, 10.);
  OB_ASSERT(correct10 == count);

  delete function;
}


