//===-- KGFormats.cpp -- supplement to read/write graphs ------------------===//
//
// This file is distributed under the GNU GPL v3 License. 
// See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "KGFormats.hpp"

void 
trim (string &str)
{
  // trim leading spaces
  size_t startpos = str.find_first_not_of(" \t");
  if( string::npos != startpos )
    str = str.substr( startpos );

  // trim trailing spaces
  size_t endpos = str.find_last_not_of(" \t");
  if( string::npos != endpos )
    str = str.substr( 0, endpos+1 );
}

int 
update_vertices (string s, map<string, int>& vertices, int &vidx)
{
  if (vertices.find(s) != vertices.end())
    return vertices[s];
  vertices[s] = vidx;
  return vidx++;
}

