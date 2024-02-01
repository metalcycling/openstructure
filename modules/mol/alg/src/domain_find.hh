#ifndef OST_MOL_ALG_DOMAIN_FIND_HH
#define OST_MOL_ALG_DOMAIN_FIND_HH

#include <ost/seq/alignment_handle.hh>

#include "module_config.hh"
#include "contact_overlap.hh"

namespace ost { namespace mol { namespace alg {

struct Domain {
  Domain(const ost::mol::EntityViewList& v): 
    inter_s(0), inter_d(0), intra_s(0), intra_d(0), views(v) {}
  Domain(): inter_s(0), inter_d(0), intra_s(0), intra_d(0) {}

  Real             inter_s;
  Real             inter_d;
  Real             intra_s;
  Real             intra_d;

  bool operator==(const Domain& rhs) const {
    return views==rhs.views && inter_s == rhs.inter_s &&
           inter_d == rhs.inter_d && intra_s == rhs.intra_s &&
           intra_d == rhs.intra_d;
  }
  ost::mol::EntityViewList   views;
};


class Domains;

typedef boost::shared_ptr<Domains> DomainsPtr;

class DLLEXPORT_OST_MOL_ALG Domains {
public:

  static DomainsPtr FromAln(ost::seq::AlignmentHandle aln,
                            Real tolerance, Real radius, Real threshold,
                            int idx_a=-1, int idx_b=-1);

  static DomainsPtr FromCOM(const SimilarityMatrix& sim, 
                            ost::seq::AlignmentHandle aln,
                            Real tolerance, int idx_a=-1, int idx_b=-1);

  // version of domain find which operates on two distances matrices
  // we on purpose only return the connected components, since the
  // distance matrices may not need to correpsond to an actual
  // entity view.
  static DomainsPtr FromDistMats(const DistanceMatrix& dmat_a,
                                 const DistanceMatrix& dmat_b,
                                 Real tolerance,
                                 Real radius,
                                 Real threshold,
                                 int num_threads,
                                 bool defined_only,
                                 bool adj_map);
  size_t GetNumFreeResidues() const { return num_free_; }

  Real GetFreeD() const { return free_d_; }
  Real GetFreeS() const { return free_s_; }

  const std::vector<Domain>& GetDomains() const { return domains_; }
  const std::vector<int>& GetComponents() const { return components_; }

  ost::img::ImageHandle GetAdjMap() const { return adj_map_; }
private:
  Domains(): num_free_(0), free_s_(0), free_d_(0) {}
  // holds the steps of the algorithm. mainly for debuggin purposes
  ost::img::ImageHandle     adj_map_;
  int                       num_free_;
  Real                      free_s_;
  Real                      free_d_;
  std::vector<Domain>       domains_;
  std::vector<int>          components_;
};

typedef std::map<int, Domain> DomainsMap;

}}}
#endif

