//===-- KGraph_tests.cpp -- basic tests for different graph algorithms ----===//
//
// This file is distributed under the GNU GPL v3 License.
// See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#include "KGraph.hpp"
#include "KGAlg.hpp"

using KGR::noload;
using KGR::colorload;
using KGR::GraphBuilder;

int test_simple(void) {
  GraphBuilder<noload, noload> GN;
  ofstream ofs;

  GN.add_path(5);
  ofs.open("path.dot", ofstream::out | ofstream::trunc);
  ofs << GN << endl;
  ofs.close();
  GN.cleanup();

  GN.add_cycle(5);
  ofs.open("cycle.dot", ofstream::out | ofstream::trunc);
  ofs << GN << endl;
  ofs.close();
  GN.cleanup();

  GN.add_clique(5);
  ofs.open("clique.dot", ofstream::out | ofstream::trunc);
  ofs << GN << endl;
  ofs.close();
  GN.cleanup();

  GN.add_path(5);
  GN.add_cycle(5);
  GN.add_clique(5);
  ofs.open("three-component.dot", ofstream::out | ofstream::trunc);
  ofs << GN << endl;
  ofs.close();
  GN.cleanup();

  return 0;
}

int test_bipart(void) {
  bool is_ok;

  GraphBuilder<colorload, noload> GN;

  GN.add_path(5);
  is_ok = color_bipartite(GN);
  assert(is_ok);
  GN.add_cycle(6);
  is_ok = color_bipartite(GN);
  assert(is_ok);
  GN.add_cycle(5);
  is_ok = color_bipartite(GN);
  assert(!is_ok);
  GN.cleanup();

  GN.add_clique(5);
  is_ok = color_bipartite(GN);
  assert(!is_ok);
  GN.cleanup();

  GN.add_full_bipart(3, 5);
  is_ok = color_bipartite(GN);
  assert(is_ok);

  ofstream ofs;
  ofs.open("bipart.dot", ofstream::out | ofstream::trunc);
  ofs << GN << endl;
  ofs.close();
  GN.cleanup();

  GraphBuilder<colorload, colorload> GNC;
  GNC.add_full_bipart(3, 5);
  is_ok = color_bipartite(GNC);
  assert(is_ok);
  hopcroft_karp(GNC);
  ofs.open("hopcroft.dot", ofstream::out | ofstream::trunc);
  ofs << GNC << endl;
  ofs.close();

  matching_to_cover(GNC);
  ofs.open("hopcroft_cover.dot", ofstream::out | ofstream::trunc);
  ofs << GNC << endl;
  ofs.close();

  GNC.cleanup();
  return 0;
}

int test_vc(void) {
  GraphBuilder<colorload, colorload> GNC;
  ifstream ifs;
  ifs.open("us.inp", ifstream::in);
  read_graph_from_stream(ifs, GNC);
  ifs.close();

  vertex_2approx(GNC);

  ofstream ofs;
  ofs.open("us.dot", ofstream::out | ofstream::trunc);
  ofs << GNC << endl;
  ofs.close();

  ofs.open("us.mps", ofstream::out | ofstream::trunc);
  out_mps_to_stream(ofs, GNC);
  ofs.close();

  GNC.cleanup();

  GraphBuilder<colorload, colorload> GNFB;

  GNFB.add_full_bipart(3, 3);
  ofs.open("bip.mps", ofstream::out | ofstream::trunc);
  out_mps_to_stream(ofs, GNFB);
  ofs.close();
  GNFB.cleanup();

  GNFB.add_clique(3);
  ofs.open("triangle.mps", ofstream::out | ofstream::trunc);
  out_mps_to_stream(ofs, GNFB);
  ofs.close();

  ofs.open("triangle.dot", ofstream::out | ofstream::trunc);
  ofs << GNFB << endl;
  ofs.close();

  GNFB.duplicate_to_bipart();

  ofs.open("bipart_triangle.dot", ofstream::out | ofstream::trunc);
  ofs << GNFB << endl;
  ofs.close();

  GNFB.cleanup();

  return 0;
}

int main(void) {
  test_simple();
  test_bipart();
  test_vc();
}
