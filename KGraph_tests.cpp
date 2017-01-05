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

int vc_routine(string gname) {
  GraphBuilder<colorload, colorload> GNC;
  using VD = typename GraphBuilder<colorload, colorload>::VertexDescriptor;
  ifstream ifs;
  ofstream ofs;
  string ins = gname + string(".inp");
  string outs = gname + string("_joined.dot");
  ifs.open(ins, ifstream::in);
  read_graph_from_stream(ifs, GNC);
  ifs.close();
  GNC.duplicate_to_bipart([](VD vdst) { vdst->load.color = 1; });
  hopcroft_karp(GNC);
  matching_to_cover(GNC);
  GNC.join_from_bipart(
      [](VD vdst, VD vsrc) { vdst->load.color += vsrc->load.color; });
  ofs.open(outs, ofstream::out | ofstream::trunc);
  ofs << GNC << endl;
  ofs.close();
  GNC.cleanup();
}

int test_vc(void) {
  bool res;
  GraphBuilder<colorload, colorload> GNC;
  using VD = typename GraphBuilder<colorload, colorload>::VertexDescriptor;

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

  GNC.add_full_bipart(3, 3);
  ofs.open("bip.mps", ofstream::out | ofstream::trunc);
  out_mps_to_stream(ofs, GNC);
  ofs.close();
  GNC.cleanup();

  GNC.add_clique(3);
  int nisol = GNC.add_default_vertex();
  GNC.add_link(0, nisol);

  ofs.open("mod_triangle.mps", ofstream::out | ofstream::trunc);
  out_mps_to_stream(ofs, GNC);
  ofs.close();

  GNC.duplicate_to_bipart([](VD vdst) { vdst->load.color = 1; });
  hopcroft_karp(GNC);
  matching_to_cover(GNC);

  GNC.join_from_bipart(
      [](VD vdst, VD vsrc) { vdst->load.color += vsrc->load.color; });

  ofs.open("mod_triangle_joined.dot", ofstream::out | ofstream::trunc);
  ofs << GNC << endl;
  ofs.close();

  GNC.cleanup();

  vc_routine("petersen");
  vc_routine("chvatal");
  vc_routine("us");

  // brute-force test for petersen
  ifs.open("petersen.inp", ifstream::in);
  read_graph_from_stream(ifs, GNC);
  ifs.close();
  res = vertex_cover_brute(GNC, 5, [](VD vsrc) { return -1; });
  assert(!res); // no 5-cover
  res = vertex_cover_brute(GNC, 6, [](VD vsrc) { return -1; });
  assert(res); // existing 6-cover
  ofs.open("petersen_bruted_6.dot", ofstream::out | ofstream::trunc);
  ofs << GNC << endl;
  ofs.close();
  GNC.cleanup();

  return 0;
}

using VD = typename GraphBuilder<colorload, colorload>::VertexDescriptor;

int standart_cbf(VD vsrc) {
  int c = vsrc->load.color;
  return (c == 0) ? 0 : (c == 2) ? 1 : -1;
}

void standart_cmf(VD vdst, int c) { 
  vdst->load.color = (c > 0) ? 2 : 0; 
}

int test_bst(void) {
  int n;
  ofstream ofs;
  GraphBuilder<colorload, colorload> GNC;
  GNC.add_path(5);
  GNC.add_path(6);
  GNC.add_cycle(5);
  GNC.add_cycle(6);

  // mark all as kernel
  for (auto vd : GNC)
    vd->load.color = 1;

  n = vertex_cover_trivial(GNC, standart_cbf, standart_cmf);
  assert(n == 11);

  ofs.open("trivial_solved.dot", ofstream::out | ofstream::trunc);
  ofs << GNC << endl;
  ofs.close();
  GNC.cleanup();

  return 0;
}

int main(void) {
  test_simple();
  test_bipart();
  test_vc();
  test_bst();
}
