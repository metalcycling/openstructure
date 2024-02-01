#include <boost/thread.hpp>

#include <ost/img/image_factory.hh>
#include <ost/mol/atom_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/img/image_handle.hh>
#include <ost/mol/residue_view.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/mol/alg/domain_find.hh>
#include <ost/mol/alg/adjacency_bitmap.hh>

using namespace ost::seq;
using namespace ost::mol;
using namespace ost::img;
using namespace ost;

namespace ost { namespace mol { namespace alg {

int clear_component(const AdjacencyBitmap& bmap, std::vector<int>& components,
                     int x, int comp_id, int offset) {
  int count=1;
  components[offset+x]=-1;
  for (size_t y=x+1; y<bmap.size(); ++y) {
    if (components[offset+y]==comp_id && bmap.get(x, y) && bmap.defined(x, y)) {
      count+=clear_component(bmap, components, y, comp_id, offset);
    }
  }
  return count;
}

int label_component(const AdjacencyBitmap& bmap, std::vector<int>& components,
                    int x, int comp_id, int offset) {
  int count=1;
  components[offset+x]=comp_id;
  for (size_t y=x+1; y<bmap.size(); ++y) {
    if (components[offset+y]==-2 && bmap.get(x, y) && bmap.defined(x, y)) {
      count+=label_component(bmap, components, y, comp_id, offset);
    }
  }
  return count;
}

// use template to force generation of two loops, one without image update
// and one with image update. we for sure don't want to have the image update 
// when running DomainFind in fast mode
//
// this is the single-threaded variant
template <bool update_img>
void update_edge_weights_st(AdjacencyBitmap* adj1p, 
                            AdjacencyBitmap* adj2p,
                            img::ImageHandle& im,
                            int n, Real tolerance,
                            int offset) {
  for (size_t i=0; i<adj1p->size(); ++i) {
    for (size_t j=i+1; j<adj1p->size(); ++j) {
      OverlapResult ov = adj1p->Overlap(i, j);
      if (ov.defined>10) {
        // fixme: avoid promotion to float. totally unneccessary.
        Real set = (static_cast<Real>(ov.nominator)/
                    static_cast<Real>(ov.denominator ? ov.denominator : 1));
        adj2p->set(i, j, set>=tolerance);
        adj2p->set(j, i, set>=tolerance);
        if (update_img) {
          im.SetReal(img::Point(i+offset, j+offset, n), set);
          im.SetReal(img::Point(j+offset, i+offset, n), set);
        }
      } else {
        adj2p->set(i, j, false, false);
        adj2p->set(j, i, false, false);
        if (update_img) {
          im.SetReal(img::Point(i+offset, j+offset, n), 0);
          im.SetReal(img::Point(j+offset, i+offset, n), 0);
        }
      }
    }
  }
}

template <bool update_img>
struct update_edge_weights_mt {

  AdjacencyBitmap* adj1p;
  AdjacencyBitmap* adj2p;
  Real tolerance;
  size_t start;
  size_t end;

  img::ImageHandle im;
  int n;
  int offset;
  update_edge_weights_mt(AdjacencyBitmap* a1p, AdjacencyBitmap* a2p,
                         Real tol, size_t s, size_t e,
                         img::ImageHandle i, int nth,
                         int o):
    adj1p(a1p), adj2p(a2p), tolerance(tol), start(s), 
    end(std::min(a1p->size(), e)), im(i), n(nth),
    offset(o) {}

  void operator()() {
    for (size_t i=start; i<end; ++i) {
      for (size_t j=0; j<adj1p->size(); ++j) {
        if (i==j) continue;
        OverlapResult ov = adj1p->Overlap(i, j);
        if (ov.defined>10) {
          // fixme: avoid promotion to float. totally unneccessary.
          Real set = (static_cast<Real>(ov.nominator)/
                      static_cast<Real>(ov.denominator ? ov.denominator : 1));
          adj2p->set(i, j, set>=tolerance);
          if (update_img) {
            im.SetReal(img::Point(i+offset, j+offset, n), set);
          }
        } else {
          adj2p->set(i, j, false, false);
          if (update_img) {
            im.SetReal(img::Point(i+offset, j+offset, n), 0.0);
          }
        }
      }
    }
  }
};

AdjacencyBitmap optimize_domains(AdjacencyBitmap adj1, Real tolerance,
                                 img::ImageHandle& im, 
                                 bool update_img,
                                 int num_threads,
                                 int offset,
                                 int full_size) {

  const static int MAX_ITERATIONS=10;
  AdjacencyBitmap adj2(adj1.size());

  AdjacencyBitmap* adj1p=&adj1;
  AdjacencyBitmap* adj2p=&adj2;
  int n=MAX_ITERATIONS;
  if (update_img) {
    im = img::CreateImage(img::Size(full_size, full_size, n+1));
    for (size_t i=0;i<adj1.size();++i) {
      for (size_t j=0;j<adj1.size();++j) {
        if (i==j) {
          im.SetReal(img::Point(i+offset, j+offset, n), 0.0);
        } else {
          im.SetReal(img::Point(i+offset, j+offset, n), 
                     adj1.get(i, j));
        }
      }
    }
  }

  if (num_threads>1) {
    while (n>0) {
      --n;
      boost::thread_group tg;
      int rows_per_thread= ((adj1p->size()+1)/num_threads);
      int adj_rows_per_thread = rows_per_thread;
      if (int modulo = rows_per_thread % 64) {
        adj_rows_per_thread = adj_rows_per_thread+64-modulo;
      }
      for (int i=0;i<num_threads;++i) {
        if (size_t(i*adj_rows_per_thread) > adj1p->size()) {
          continue;
        }
        if (update_img) {
          tg.create_thread(update_edge_weights_mt<true>(adj1p, adj2p, tolerance,
                                                  (i+0)*adj_rows_per_thread, 
                                                  (i+1)*adj_rows_per_thread,
                                                  im, n, offset));
        } else {
          tg.create_thread(update_edge_weights_mt<false>(adj1p, adj2p, tolerance,
                                                  (i+0)*adj_rows_per_thread, 
                                                  (i+1)*adj_rows_per_thread,
                                                  im, n, offset));
        }
      }
      tg.join_all();
      std::swap(adj1p, adj2p);
    }
    return *adj1p;
  }
  while (n>0) {
    --n;
    if (update_img) {
      update_edge_weights_st<true>(adj1p, adj2p, im, n, tolerance, offset);
    } else {
      update_edge_weights_st<false>(adj1p, adj2p, im, n, tolerance, offset);
    }
    std::swap(adj1p, adj2p);
  }
  return *adj1p;
}

void connected_components(const AdjacencyBitmap& adj_bitmap,
                         std::vector<int>& components, int offset) {
  int comp_id=0;
  for (int i=0; i<offset; ++i) {
    components[i] = -1;
  }
  for (size_t i=0; i<adj_bitmap.size(); ++i) {
    if (components[offset+i]!=-2)  {
      continue;
    }
    int count = label_component(adj_bitmap, components, i, comp_id, offset);
    if (count>16) {
      comp_id+=1;
    } else {
      clear_component(adj_bitmap, components, i, comp_id, offset);
    }
  }
  for (size_t i=offset+adj_bitmap.size(); i<components.size(); ++i) {
    components[i] = -1;
  }
}

DomainsPtr Domains::FromAln(ost::seq::AlignmentHandle aln,
                            Real tolerance, Real radius, Real threshold,
                            int idx_a, int idx_b) {
  SimilarityMatrix sim = ContactOverlapMap(aln, tolerance, OM_ABSOLUTE,
                                           radius);
  return Domains::FromCOM(sim, aln, threshold, idx_a, idx_b);
}

DomainsPtr Domains::FromCOM(const SimilarityMatrix& sim, AlignmentHandle aln,
                            Real tolerance, int idx_a, int idx_b)
{
  bool update_img = false;
  AdjacencyBitmap adj1(sim, 0.5);
  img::ImageHandle im;

  AdjacencyBitmap final_adj(optimize_domains(adj1, tolerance, im, 
                                             update_img, 1, 0,
                                             adj1.size()));

  // identify connected components in the graph
  std::vector<int> components(adj1.size(), -2);
  connected_components(final_adj, components, 0);

  DomainsMap domains_map;
  if (idx_a>=0 && idx_b>=0) {
    for (int i=0; i<aln.GetLength(); ++i) {
      if (components[i]==-1) {
        continue;
      }
      DomainsMap::iterator k = domains_map.find(components[i]);
      if (k == domains_map.end()) {
        std::pair<DomainsMap::iterator, bool> k2;
        k2 = domains_map.insert(std::make_pair(components[i], 
                                               Domain(EntityViewList())));
        k = k2.first;
        k->second.views.push_back(aln.GetSequence(idx_a).GetAttachedView().CreateEmptyView());
        k->second.views.push_back(aln.GetSequence(idx_b).GetAttachedView().CreateEmptyView());
      }
      if (ResidueView r=aln.GetResidue(idx_a, i)) {
        k->second.views[0].AddResidue(r.GetHandle(), ViewAddFlag::INCLUDE_ALL);
      }
      if (ResidueView r=aln.GetResidue(idx_b, i)) {
        k->second.views[1].AddResidue(r.GetHandle(), ViewAddFlag::INCLUDE_ALL);
      }
    }
  }

  int num_free = 0;

  Real free_d = 0, free_s = 0;

  for (size_t i=0; i<adj1.size(); ++i) {
    if (components[i]==-1) {
      num_free+=1;
    }
    for (size_t j=i+1; j<adj1.size(); ++j) {
      if (!final_adj.defined(i, j)) {
        continue;
      }
      if (components[i] == -1 || components[j] == -1) {
       free_d+=1;
       if (sim.Get(i, j)>=0)
          free_s+=sim.Get(i, j);
        continue;
      }
      Real s = sim.Get(i, j);
      if (s<0) {
        continue;
      }
      if (components[i] == components[j]) {
        domains_map.find(components[i])->second.intra_s+= adj1.get(i, j);
        domains_map.find(components[i])->second.intra_d+= 1;
      } else {
        domains_map.find(components[i])->second.inter_s+=adj1.get(i, j);
        domains_map.find(components[i])->second.inter_d+= 1;
        domains_map.find(components[j])->second.inter_s+= adj1.get(i, j);
        domains_map.find(components[j])->second.inter_d+= 1;
      }
    }
  }

  DomainsPtr domains(new Domains);
  for (DomainsMap::iterator i = domains_map.begin(), 
       e = domains_map.end(); i!=e; ++i) {
    domains->domains_.push_back(i->second);
  }
  domains->num_free_ = num_free;
  domains->components_ = components;
  domains->free_d_ = free_d;
  domains->free_s_ = free_s;
  if (im.IsValid())
     domains->adj_map_ = im;
  return domains;
}


std::pair<int, int> get_defined_region(const DistanceMatrix& dmat_a,
                                       const DistanceMatrix& dmat_b) {
  int start = 0;
  for (int i=0; i<dmat_a.GetSize(); ++i) {
    if (dmat_a.Get(i, i)>=0.0 && dmat_b.Get(i, i)>=0.0) {
      break;
    }
    start++;
  }
  if (start == dmat_a.GetSize()) {
    // the two don't have an overlap
    return std::make_pair(start, start);
  }
  int end = dmat_a.GetSize();
  for (int i=dmat_a.GetSize()-1; i>=0; --i) {
    if (dmat_a.Get(i, i)>=0.0 && dmat_b.Get(i, i)>=0.0) {
      break;
    }
    end--;
  }
  return std::make_pair(start, end);
}

DomainsPtr Domains::FromDistMats(const DistanceMatrix& dmat_a, 
                                 const DistanceMatrix& dmat_b,
                                 Real tolerance, Real radius, 
                                 Real threshold, int num_threads,
                                 bool defined_only,
                                 bool update_img)
{
  // for distance matrices derived from target-template alignments, C-
  // and N-terminal regions of the distance map are undefined, since
  // they are not covered in the alignment. To avoid doing a full-blown
  // O(n^3) update on the full distance matrix, we first identify regions
  // of the distance map which are defined and perform the update only on
  // that part of the map.
  std::pair<int, int> region(0, dmat_a.GetSize());
  if (defined_only)
    region = get_defined_region(dmat_a, dmat_b);

  DomainsPtr domains(new Domains);

  if (region.first == region.second) {
    // there is no overlap whatsoever between the two distance matrices.
    // all present residues will thus be assigned to the unassigned
    // component (designated as -1).
    domains->components_.resize(dmat_a.GetSize(), -1);
    return domains;
  }

  // derive initial connectivity bitmap from the two distance matrices
  AdjacencyBitmap adj1(dmat_a, dmat_b, tolerance, radius, region.first,
                       region.second);

  img::ImageHandle im;

  // optimize domains
  AdjacencyBitmap final_adj(optimize_domains(adj1, threshold, im, 
                                             update_img, num_threads,
                                             region.first,
                                             dmat_a.GetSize()));
  // identify connected components in the graph
  std::vector<int> components(dmat_a.GetSize(), -2);
  connected_components(final_adj, components, region.first);
  domains->components_ = components;
  if (im.IsValid())
     domains->adj_map_ = im;
  return domains;

}

}}}

