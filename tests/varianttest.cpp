#include <OBVariant>
#include "obtest.h"

using namespace OpenBabel::OBFFs;

using namespace std;

int main()
{
  OBVariant var1(3, "myvar");
  OB_ASSERT( var1.GetName() == "myvar" );
  OB_ASSERT( var1.AsInt() == 3 );
  OB_ASSERT( var1.AsDouble() == 3.0 );
  OB_ASSERT( var1.AsBool() == true );
  OB_ASSERT( var1.AsString() == "3" );

  OBVariant var2(3.1, "kb");
  OB_ASSERT( var2.GetName() == "kb" );
  OB_ASSERT( var2.AsInt() == 3 );
  OB_ASSERT( var2.AsDouble() == 3.1 );
  OB_ASSERT( var2.AsBool() == true );
  OB_ASSERT( var2.AsString() == "3.1" );

  OBVariant var3(true, "mybool");
  OB_ASSERT( var3.GetName() == "mybool" );
  OB_ASSERT( var3.AsInt() == 1 );
  OB_ASSERT( var3.AsBool() == true );
  OB_ASSERT( var3.AsDouble() == 1.0 );
  OB_ASSERT( var3.AsString() == "True" );

  OBVariant var4(0, "myvar");
  OB_ASSERT( var4.AsBool() == false );

  OBVariant var5("3.678", "myvar");
  OB_ASSERT( var5.GetName() == "myvar" );
  OB_ASSERT( var5.AsInt() == 3 );
  OB_ASSERT( var5.AsDouble() == 3.678 );
  OB_ASSERT( var5.AsBool() == false );
  OB_ASSERT( var5.AsString() == "3.678" );

}
