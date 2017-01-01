//===-- KGAlg.hpp -- some common algorithms, related to graphs ------------===//
//
// This file is distributed under the GNU GPL v3 License.
// See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// This file contains:
//
// color_bipartite -- properly {0,1}-color graph if it is bipartite
//                    or return false otherwise
//
// hopcroft_karp -- find maximum cardinality matching
//                  in bipartite {0,1}-colored graph
//
// matching_to_cover -- maximal cardinality matching to minimum vertex cover
//                      for bipartite graph with {0,1}-colored matching
//
// vertex_2approx -- find 2-approximation for vertex cover in general graph
//
// TODO: here shall also go crown decomposition and LPVC approximation
//
//===----------------------------------------------------------------------===//

#ifndef GRAPH_KALG_GUARD__
#define GRAPH_KALG_GUARD__

#include "KGInc.hpp"

// DFS-like coloring with additional stack, like Knuth alg7-B
template <typename G> bool color_bipartite(G &g) {
  using VD = typename G::VertexDescriptor;
  for (auto vd : g)
    vd->load.color = -1;

  // DFS-like approach.
  // Simpler one is possible, but this one allows to output loop on error
  forward_list<VD> stack;
  for (auto vd : g) {
    int &wc = vd->load.color;
    if (wc >= 0)
      continue;
    wc = 0;
    stack.push_front(vd);
    while (!stack.empty()) {
      VD u = stack.front();
      stack.pop_front();
      int &uc = u->load.color;
      assert(uc >= 0);

      for (auto a = u->arcs; a != g.last_edge(); a = a->next) {
        auto vd = a->tip;
        if (vd->load.color == -1) {
          vd->load.color = 1 - uc;
          stack.push_front(vd);
        }
        if (vd->load.color == uc) {
          // TODO: here we may also output odd-length loop
          //       it is already on stack
          return false;
        }
      }
    }
  }

  return true;
}

template <typename G, typename VD>
bool hk_bfs(G &g, vector<VD> &U, map<VD, VD> &PairU, map<VD, VD> &PairV,
            map<VD, int> &Dist);

template <typename G, typename VD>
bool hk_dfs(G &g, vector<VD> &U, vector<VD> &V, map<VD, VD> &PairU,
            map<VD, VD> &PairV, map<VD, int> &Dist, VD u);

// input is 0-1 colored bipartite graph
// with colorable edges
template <typename G> int hopcroft_karp(G &g) {
  using VD = typename G::VertexDescriptor;

  int matching = 0;
  auto nil = g.last_vertex();
  vector<VD> U, V;
  map<VD, VD> PairU, PairV;
  map<VD, int> Dist;

  for (auto vd : g)
    if (vd->load.color == 0) {
      U.push_back(vd);
      PairU[vd] = nil;
    } else {
      V.push_back(vd);
      PairV[vd] = nil;
    }

  while (hk_bfs(g, U, PairU, PairV, Dist))
    for (auto ud : U)
      if (PairU[ud] == nil)
        if (hk_dfs(g, U, V, PairU, PairV, Dist, ud))
          matching = matching + 1;

  return matching;
}

template <typename G, typename VD>
bool hk_bfs(G &g, vector<VD> &U, map<VD, VD> &PairU, map<VD, VD> &PairV,
            map<VD, int> &Dist) {
  auto nil = g.last_vertex();
  auto enil = g.last_edge();
  int inf = std::numeric_limits<int>::max();
  list<VD> Q;

  for (auto ud : U) {
    if (PairU[ud] == nil) {
      Dist[ud] = 0;
      Q.push_back(ud);
    } else
      Dist[ud] = inf;
  }

  Dist[nil] = inf;

  while (!Q.empty()) {
    auto u = Q.front();
    Q.pop_front();
    if (Dist[u] < Dist[nil])
      for (auto e = u->arcs; e != enil; e = e->next) {
        VD v = e->tip;
        if (Dist[PairV[v]] == inf) {
          Dist[PairV[v]] = Dist[u] + 1;
          Q.push_back(PairV[v]);
        }
      }
  }
  return (Dist[nil] != inf);
}

template <typename G, typename VD>
bool hk_dfs(G &g, vector<VD> &U, vector<VD> &V, map<VD, VD> &PairU,
            map<VD, VD> &PairV, map<VD, int> &Dist, VD u) {
  auto nil = g.last_vertex();
  auto enil = g.last_edge();
  int inf = std::numeric_limits<int>::max();
  if (u == nil)
    return true;
  for (auto e = u->arcs; e != enil; e = e->next) {
    VD v = e->tip;
    if (Dist[PairV[v]] == Dist[u] + 1)
      if (hk_dfs(g, U, V, PairU, PairV, Dist, PairV[v])) {
        PairV[v] = u;
        PairU[u] = v;
        e->load.color = 1;
        auto ev = g.get_sibling(e, u);
        assert(ev != enil);
        ev->load.color = 1;
        return true;
      }
  }
  Dist[u] = inf;
  return false;
}

template <typename G, typename VD> bool vertex_unmatched(G &g, VD u);

template <typename G, typename VD> void remove_matching(G &g, VD u, int newc);

// input is 0-1 edge colored bipartite graph
// with colorable vertices
template <typename G> int matching_to_cover(G &g) {
  using VD = typename G::VertexDescriptor;
  auto enil = g.last_edge();
  bool has_unmatched = true;
  int vcsz = 0;

  for (auto vd : g)
    vd->load.color = -1;

  while (has_unmatched) {
    has_unmatched = false;
    for (auto vd : g)
      if ((vd->load.color == -1) && vertex_unmatched(g, vd)) {
        has_unmatched = true;

        // color it 0
        vd->load.color = 0;

        // color all matched siblings 1
        for (auto e = vd->arcs; e != enil; e = e->next)
          if ((e->tip->load.color == -1) && !vertex_unmatched(g, e->tip)) {
            e->tip->load.color = 1;
            vcsz += 1;
            remove_matching(g, e->tip, 0);
          }
      }
  }

  for (auto vd : g)
    if (vd->load.color == -1) {
      assert(!vertex_unmatched(g, vd));
      vd->load.color = 1;
      vcsz += 1;
      remove_matching(g, vd, 0);
    }

  return vcsz;
}

template <typename G, typename VD> bool vertex_unmatched(G &g, VD u) {
  auto enil = g.last_edge();
  for (auto e = u->arcs; e != enil; e = e->next)
    if (e->load.color == 1)
      return false;
  return true;
}

template <typename G, typename VD> void remove_matching(G &g, VD u, int newc) {
  auto enil = g.last_edge();
  for (auto e = u->arcs; e != enil; e = e->next)
    if (e->load.color == 1) {
      e->load.color = 0;
      assert(e->tip->load.color == -1);
      e->tip->load.color = newc;
      return;
    }
  assert(0 && "we should not be here: vertex supposed to be matched");
}

// 2-approximation for vertex cover
template <typename G> void vertex_2approx(G &g) {
  auto enil = g.last_edge();

  for (auto vd : g)
    vd->load.color = 0;

  for (auto vd : g) {
    if (vd->load.color == 1)
      continue;

    for (auto e = vd->arcs; e != enil; e = e->next)
      if (e->tip->load.color == 0) {
        vd->load.color = 1;
        e->tip->load.color = 1;
        e->load.color = 2;
        break;
      }
  }
}

#endif
