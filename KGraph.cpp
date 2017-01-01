//===-- KGraph.hpp -- graph representation supplement ---------------------===//
//
// This file is distributed under the GNU GPL v3 License. 
// See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "KGraph.hpp"

namespace KGR {

const char *recode (int color)
{
  switch (color)
    {
    default:
      return "black";
    case 1:
      return "red";
    case 2:
      return "blue";
    case 3:
      return "green";
    }
}

}

