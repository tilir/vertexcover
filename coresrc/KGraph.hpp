//===-- KGraph.hpp -- graph representation experiments --------------------===//
//
// This file is distributed under the GNU GPL v3 License.
// See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// Graph may be mutable or immutable (TODO: memory-efficient immutable graph)
// Vertices and edges might have some load (like color, weight, etc) or not
// TODO: specialization for noload without load at all
//
// Author knows about BGL
// TODO: maybe some performance comparisons with BGL here
//
//===----------------------------------------------------------------------===//

#ifndef GRAPH_KG_GUARD__
#define GRAPH_KG_GUARD__

#include "KGFormats.hpp"

namespace KGR {

//------------------------------------------------------------------------------
//
//  Load types
//
//------------------------------------------------------------------------------

struct noload {
  friend ostream &operator<<(ostream &stream, const noload &) { return stream; }
};

const char *recode(int color);

struct colorload {
  int color;
  friend ostream &operator<<(ostream &stream, const colorload &l) {
    stream << "color=\"" << recode(l.color) << "\"";
    return stream;
  }
};

// TODO: think about immutable graph
// idea is: fixed vertex array, fixed edge array, everything on stack, etc.

//------------------------------------------------------------------------------
//
//  Mutable graph
//
//------------------------------------------------------------------------------

// base for vertex for mutable graph
template <typename VL, typename ET> struct IVertex {
  VL load{};
  ET *arcs = nullptr;

  // protected dtor to not make it virtual
  // now attempt to do from outside:
  //   IVertex *base = new Vertex;
  //   delete base; /* possible memory leak here */
  // will fail because we can not delete with protected dtor
protected:
  ~IVertex() {}
};

// edge for mutable graph
template <typename EL, typename VT> struct Edge final {
  EL load;
  VT *tip = nullptr;
  Edge *next = nullptr;
  Edge(EL l, VT *t) : load(l), tip(t) {}
};

template <typename VT, typename EL>
static inline void add_link_to(VT *v1, VT *v2, EL l) {
  using ET = Edge<EL, VT>;
  ET *e12 = new ET(l, v2);
  v1->link_to(v2, e12);
}

template <typename VT, typename EL>
static inline void link(VT *v1, VT *v2, EL l) {
  assert(v1 && v2 && "Linking to null vertex is bad idea");
  // TODO: may be allocate edges in 2-blocks as in Knuth 4A
  //       but this will complicate deletion
  using ET = Edge<EL, VT>;
  add_link_to(v1, v2, l);
  add_link_to(v2, v1, l);
}

// mutable graph, heap only
template <typename VL, typename EL> class GraphBuilder final {
  // vertex for this graph
  struct Vertex : public IVertex<VL, Edge<EL, Vertex>> {
    using ET = Edge<EL, Vertex>;
    void link_to(Vertex *v, ET *edge) {
      assert(edge->tip == v);
      // without this-> we have unqualified lookup!
      edge->next = this->arcs;
      this->arcs = edge;
    }
  };
  vector<Vertex *> vertices_;

public:
  GraphBuilder() = default;
  GraphBuilder(const GraphBuilder &) = delete;
  GraphBuilder &operator=(const GraphBuilder &) = delete;
  ~GraphBuilder() { cleanup(); }

  // general interface
public:
  using VertexDescriptor = Vertex *;
  using VertexIterator = typename vector<Vertex *>::iterator;
  using VT = Vertex;
  using ET = typename Vertex::ET;
  using EdgeDescriptor = ET *;
  const char *name() const { return "G"; }
  int nvertices() { return vertices_.size(); }
  VT *front() { return vertices_.front(); }
  VT *back() { return vertices_.back(); }
  VertexIterator begin() { return vertices_.begin(); }
  VertexIterator end() { return vertices_.end(); }
  VertexDescriptor last_vertex() { return nullptr; }
  EdgeDescriptor last_edge() { return nullptr; }
  ET *get_edge(VT *u, VT *v) {
    assert(u && v && "Edge for null is bad idea");
    for (auto eu = u->arcs; eu != nullptr; eu = eu->next)
      if (eu->tip == v)
        return eu;
    return nullptr;
  }
  ET *get_sibling(ET *e, VT *u) {
    assert(e->tip != u);
    assert(e && u && "Sibling for null is bad idea too");
    // TODO: 2-blocks for edges will make this O(1) instead O(n)
    //       like: return ((e % 4) == 0) ? (e + 4) : (e - 4)
    //       but this effectiveness cost some uglyness and protability
    for (auto ev = e->tip->arcs; ev != nullptr; ev = ev->next)
      if (ev->tip == u)
        return ev;
    return nullptr;
  }
  int degree(VT *u) {
    assert(u != nullptr);
    int deg = 0;
    for (auto eu = u->arcs; eu != nullptr; eu = eu->next)
      deg += 1;
    return deg;
  }

  // modifiable specifics
public:
  int add_default_vertex(void) {
    VT *vert = new Vertex();
    vertices_.push_back(vert);
    return vertices_.size() - 1;
  }

  // link with default load
  void add_link(int i, int j) {
    assert(i >= 0 && i < (int)vertices_.size());
    assert(j >= 0 && j < (int)vertices_.size());
    link(vertices_[i], vertices_[j], EL{});
  }

  void partial_cleanup(int nstart, int nend) {
    assert(nstart < nend);
    assert(nstart >= 0);
    assert(nend <= vertices_.size());
    for (auto vi = nstart; vi != nend; ++vi) {
      for (ET *e = vertices_[vi]->arcs; e != nullptr;) {
        ET *tmp = e->next;
        // TODO:
        // ET *back = get_sibling (e, vertices_[vi]);
        // delete back;
        delete e;
        e = tmp;
      }
      delete vertices_[vi];
    }
    vertices_.erase(vertices_.begin() + nstart, vertices_.begin() + nend);
  }

  void cleanup() {
    for (auto v : vertices_) {
      for (ET *e = v->arcs; e != nullptr;) {
        ET *tmp = e->next;
        delete e;
        e = tmp;
      }
      delete v;
    }
    vertices_.clear();
  }

  // add n isolated vertices
  void add_isolated(int n) {
    for (int vcount = 0; vcount < n; ++vcount)
      add_default_vertex();
  }

  void add_path(int n) {
    VT *vcurr = nullptr;
    for (int vcount = 0; vcount < n; ++vcount) {
      int inext = add_default_vertex();
      VT *vnext = vertices_[inext];
      if (vcurr)
        link(vcurr, vnext, EL{});
      vcurr = vnext;
    }
  }

  void add_cycle(int n) {
    int start = vertices_.size();
    add_path(n);
    assert(n > 2);
    add_link(start, start + n - 1);
  }

  void add_clique(int n) {
    assert(n > 2);
    int start = vertices_.size();
    add_path(n);
    int process = n / 2;
    for (int i = start; i != start + process; ++i)
      for (int j = i + 2; j != start + n; ++j)
        add_link(i, j);
  }

  void add_full_bipart(int n, int m) {
    int start = vertices_.size();
    add_isolated(n + m);
    for (int i = start; i != start + n; ++i)
      for (int j = start + n; j != start + n + m; ++j)
        add_link(i, j);
  }

  // duplicates current graph to create bipartite for LPVC
  template <typename C> void duplicate_to_bipart(C colors_callback) {
    int start = vertices_.size();
    assert(start > 0 && "Not good idea doing this on empty graph");
    add_isolated(start);
    int n = 0;
    map<VertexDescriptor, int> indexes;
    for (auto vd : vertices_)
      indexes[vd] = n++;

    for (int i = 0; i != start; ++i)
      for (auto ed = vertices_[i]->arcs; ed != nullptr; ed = ed->next) {
        int nold = indexes[ed->tip];
        int nnew = nold + start;
        ed->tip = vertices_[nnew];
        add_link_to(vertices_[nnew], vertices_[i], EL{});
      }

    for (int i = 0; i != start; ++i)
      colors_callback(vertices_[i]);
  }

  // brings {0,1}-colored bipartite back to {0,1,2}-colored graph
  // color 1 is for 1/2 vertices of core task
  template <typename C> void join_from_bipart(C colors_callback) {
    int n = 0;
    map<VertexDescriptor, int> indexes;
    for (auto vd : vertices_)
      indexes[vd] = n++;

    int nall = vertices_.size();
    assert((nall % 2) == 0);
    int nhalf = nall / 2;
    for (int idx = 0; idx != nhalf; ++idx) {
      for (auto ed = vertices_[idx]->arcs; ed != nullptr; ed = ed->next) {
        int tipold = indexes[ed->tip];
        assert(tipold > nhalf - 1);
        int tipnew = tipold - nhalf;
        ed->tip = vertices_[tipnew];
      }
      colors_callback(vertices_[idx], vertices_[idx + nhalf]);
    }

    partial_cleanup(nhalf, nall);
    assert(vertices_.size() == nhalf);
  }

  friend ostream &operator<<(ostream &stream, GraphBuilder &g) {
    out_dot_to_stream(stream, g);
    return stream;
  }

  // TODO: dump to immutable
};
}

#endif
