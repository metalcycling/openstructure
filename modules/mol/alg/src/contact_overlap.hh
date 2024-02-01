#ifndef OST_MOL_ALG_CONTACT_HH
#define OST_MOL_ALG_CONTACT_HH

#include <boost/shared_ptr.hpp>

#include <ost/mol/entity_view.hh>
#include <ost/img/image.hh>
#include <ost/seq/alignment_handle.hh>


#include "module_config.hh"
#include "similarity_matrix.hh"

namespace ost { namespace mol { namespace alg {

struct Contact {
  Contact(): atom_a(), atom_b() {}
  Contact(ost::mol::AtomHandle a, ost::mol::AtomHandle b):
    atom_a(a), atom_b(b)
  { }
  ost::mol::AtomHandle atom_a;
  ost::mol::AtomHandle atom_b;

  bool operator==(const Contact& rhs) const
  {
    return atom_a==rhs.atom_a && atom_b==rhs.atom_b;
  }

  bool operator!=(const Contact& rhs) const
  {
    return !this->operator==(rhs);
  }
};

typedef std::vector<Contact> ContactList;

/// \brief returns a a list of contacts between C-alpha atoms
ContactList DLLEXPORT_OST_MOL_ALG
Contacts(const ost::mol::EntityView& ent, Real min_dist,
         Real max_dist);

typedef enum {
 OM_RELATIVE,
 OM_ABSOLUTE
} OverlapMode;
/// \brief calculates and assigns the contact node degree as the int
///   property 'degree' for each atom
void DLLEXPORT_OST_MOL_ALG CalcContactNodeDegree(const ContactList& contacts);

struct LDT {
  Real ha;
  Real ma;
  Real la;
};

class DistanceMatrix;

typedef boost::shared_ptr<DistanceMatrix> DistanceMatrixPtr;

class DLLEXPORT_OST_MOL_ALG DistanceMatrix: public ost::TriMatrix<Real> {
public:
  static DistanceMatrixPtr Create(const ost::seq::ConstSequenceHandle& s);
  static DistanceMatrixPtr FromAln(const ost::seq::AlignmentHandle& a,
                                   int index=1);
private:
  DistanceMatrix(const ost::seq::ConstSequenceHandle& sequence);
  DistanceMatrix(const ost::seq::AlignmentHandle& aln, int index);

};


Real DLLEXPORT_OST_MOL_ALG ContactOverlap(const DistanceMatrixPtr& lhs,
                              const DistanceMatrixPtr& rhs,
                              Real tolerance,
                              OverlapMode mode=OM_RELATIVE,
                              Real dist_limit=1e6);


Real DLLEXPORT_OST_MOL_ALG
LocalDistanceTest(const ost::mol::EntityView& mdl,
                  const ost::mol::EntityView& ref,
                  Real max_dist,
                  Real tolerance=0.02,
                  bool only_complete=false);

SimilarityMatrix DLLEXPORT_OST_MOL_ALG
ContactOverlapMap(const ost::seq::AlignmentHandle& target_tpl_aln,
                  Real tolerance, int idx_a, int idx_b);

SimilarityMatrix DLLEXPORT_OST_MOL_ALG
ContactOverlapMap(const ost::seq::AlignmentHandle& target_tpl_aln,
                  Real tolerance, OverlapMode mode=OM_RELATIVE,
                  Real max_dist=1e6);

SimilarityMatrix DLLEXPORT_OST_MOL_ALG
ContactOverlapMap(const DistanceMatrix& dmat1, const DistanceMatrix& dmat2,
                  Real tolerance, OverlapMode mode=OM_RELATIVE,
                  Real max_dist=1e6);


//ost::img::ImageHandle 
void DLLEXPORT_OST_MOL_ALG
DomainsFromCOM(const SimilarityMatrix& sim, ost::seq::AlignmentHandle aln,
               Real tolerance, int idx_a=-1, int idx_b=-1);

}}}

#endif
