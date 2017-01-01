//===-- KGFormats.hpp -- read/write graphs to/from different formats ------===//
//
// This file is distributed under the GNU GPL v3 License. 
// See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// This file contains:
//
// out_dot_to_stream -- outputs G in dot format
//
// out_mps_to_stream -- outputs G in mps format for vertex cover LPVC approx
//
// read_graph_from_stream -- reads G from file in simplest form (vertex pairs)
//
//===----------------------------------------------------------------------===//

#ifndef GRAPH_KFMTS_GUARD__
#define GRAPH_KFMTS_GUARD__

#include "KGInc.hpp"

//===----------------------------------------------------------------------===//
//
// Writing graphs in given formats
//
//===----------------------------------------------------------------------===//

// dot format: https://en.wikipedia.org/wiki/DOT_(graph_description_language)
// useful for visualizations
template <typename G> void
out_dot_to_stream (ostream& stream, G& g)
{
  using VD = typename G::VertexDescriptor;
  int n = 0;
  map<VD, int> indexes;
  stream << "graph " << g.name() << "{" << endl;
  for (auto vd : g) {
    stream << "v" << n << "[" << vd->load << "]" << ";" << endl;
    indexes[vd] = n;
    n += 1;
  }

  // zero or one vertices corner case
  if (n < 2) {
    stream << "}" << endl;
    return;
  }

  for (auto vd : g) 
    for (auto ed = vd->arcs; ed != g.last_edge(); ed = ed->next)
      {
        int fst = indexes[vd];
        int snd = indexes[ed->tip];
        if (fst > snd) continue; // links are always symmetric
        stream << "v" << fst << " -- " << "v" << snd 
               << "[" << ed->load << "]"
               << endl;
      }

  stream << "}" << endl;
}

// mps format: https://en.wikipedia.org/wiki/MPS_(format)
// useful for LP approximations
template <typename G> void
out_mps_to_stream (ostream& stream, G& g)
{
  using VD = typename G::VertexDescriptor;
  stream << std::setw(14) << std::left << "NAME" << 
                           std::setw(0) << "BIPART" << endl;
  stream << "ROWS" << endl;
  stream << std::setw(4) << std::left << " N" << std::setw(0) << "COST" << endl;

  int n = 0;
  map<VD, int> indexes;
  set <pair<int, int>> proper_edges;

  for (auto vd : g)
    indexes[vd] = n++;

  for (auto vd : g) 
    for (auto ed = vd->arcs; ed != g.last_edge(); ed = ed->next)
      if (indexes[vd] < indexes[ed->tip])
        proper_edges.insert(make_pair(indexes[vd], indexes[ed->tip]));
  
  for (auto pe : proper_edges)
    stream << std::setw(4) << std::left << " G" << 
              std::setw(0) << "V" << pe.first <<
              std::setw(0) << "V" << pe.second << endl;

  stream << "COLUMNS" << endl;
  for (auto vd : g) 
    {
      int vidx = indexes[vd];
      string vstr = string("V") + to_string(vidx);
      stream << std::setw(4) << " " 
             << std::setw(10) << std::left << vstr
             << std::setw(20) << std::left << "COST"
             << std::setw(0) << "1" << endl;

      for (auto ed = vd->arcs; ed != g.last_edge(); ed = ed->next)
        {
          int iless = indexes[ed->tip];
          int ibigger = vidx;
          if (iless > ibigger) 
            std::swap (iless, ibigger);
          string estr = string("V") + to_string(iless) +
                        string("V") + to_string (ibigger);
          stream << std::setw(4) << " " 
                 << std::setw(10) << std::left << vstr
                 << std::setw(20) << std::left << estr
                 << std::setw(0) << "1" << endl;
        }
    }
  
  stream << "RHS" << endl;
  for (auto pe : proper_edges)
    {
      string estr = string("V") + to_string(pe.first) +
                    string("V") + to_string(pe.second);
      stream << std::setw(4) << " " 
             << std::setw(10) << std::left << "RHS1"
             << std::setw(20) << std::left << estr
             << std::setw(0) << "1" << endl;
    }

  stream << "BOUNDS" << endl;
  for (auto vd : g)
    {
      string vstr = string("V") + to_string(indexes[vd]);

      stream << std::setw(4) << " LO" 
             << std::setw(10) << std::left << "BND1"
             << std::setw(20) << std::left << vstr
             << std::setw(0) << "0" << endl;
    }
  stream << "ENDATA" << endl;
}

//===----------------------------------------------------------------------===//
//
// Reading graphs in given formats
//
//===----------------------------------------------------------------------===//

void trim (string &str);
int update_vertices (string s, map<string, int>& vertices, int &vidx);

// simplest format
// AnyVertexId AnyOtherId
template <typename G> void
read_graph_from_stream (istream& stream, G& g)
{
  g.cleanup();
  string line;
  int vidx = 0;
  map<string, int> vertices;
  map<int, set<int>> edges;

  while (getline (stream, line))
    {
      size_t found = line.find(" "); 
      assert (found != string::npos && 
              "You must separate vertices with space(s)");
      string lhs = line.substr(0, found);
      string rhs = line.substr(found + 1, string::npos);
      int lnum = update_vertices (lhs, vertices, vidx);
      int rnum = update_vertices (rhs, vertices, vidx);
      if (lnum > rnum) std::swap (lnum, rnum);
      if (edges.find (lnum) == edges.end())
        {
          set<int> s {rnum};
          edges[lnum] = s;          
        }
      else
        edges[lnum].insert(rnum);        
    }
  
  g.add_isolated(vidx);
  for (auto v1 : edges)
    for (auto v2 : v1.second)
      g.add_link (v1.first, v2);  
}

#endif
