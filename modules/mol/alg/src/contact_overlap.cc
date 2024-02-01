#include <ost/log.hh>
#include <ost/mol/mol.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/merge_pairwise_alignments.hh>
#include "contact_overlap.hh"

using namespace ost;
using namespace ost::mol;
using namespace ost::seq;

namespace ost { namespace mol { namespace alg {

ContactList Contacts(const ost::mol::EntityView& ent, Real min_dist,
                     Real max_dist)
{
  Real min_sqr=min_dist*min_dist;
  Real max_sqr=max_dist*max_dist;
  ContactList contacts;
  ResidueViewList residues=ent.GetResidueList();
  for (ResidueViewList::const_iterator
    i=residues.begin(), e=residues.end(); i!=e; ++i) {
    AtomView atom_a=i->FindAtom("CA");
    if (!atom_a.IsValid()) {
      continue;
    }
    for (ResidueViewList::const_iterator
         j=i+1, e2=residues.end(); j!=e2; ++j) {
      AtomView atom_b=j->FindAtom("CA");
      if (!atom_b.IsValid()) {
        continue;
      }
      Real d=geom::Length2(atom_a.GetPos()-atom_b.GetPos());
      if (min_sqr<d && max_sqr>d) {
        contacts.push_back(Contact(atom_a.GetHandle(), atom_b.GetHandle()));
      }
    }

  }
  return contacts;
}

namespace {

struct AtomHashCmp {
  bool operator()(const AtomHandle& a,
                  const AtomHandle& b) const {
    return a.GetHashCode()<b.GetHashCode();
  }
};

}

void CalcContactNodeDegree(const ContactList& contacts)
{
  std::map<AtomHandle, int, AtomHashCmp> degrees;
  for (ContactList::const_iterator
       i=contacts.begin(), e=contacts.end(); i!=e; ++i) {
    // yes, this is sane. upon first accessing a key, the value is initialised
    // to 0 and *then* incremented.
    degrees[i->atom_a]+=1;
    degrees[i->atom_b]+=1;
  }
  for (std::map<AtomHandle, int>::iterator
       i=degrees.begin(), e=degrees.end(); i!=e; ++i) {
    AtomHandle atom=i->first;
    atom.SetIntProp("degree", i->second);
  }
}

Real Reward(Real ref_dist, Real obs_dist, Real tolerance,
            OverlapMode mode=OM_RELATIVE) {
  assert(ref_dist>0.0);
  Real norm_diff=mode==OM_RELATIVE ?
                       std::abs((obs_dist-ref_dist)/sqrt(obs_dist+ref_dist)) :
                       std::abs(ref_dist-obs_dist);
  if (norm_diff<tolerance*0.5) {
    return 1.0;
  }
  Real one_over_tol=0.5/tolerance;
  return std::max(0.0, 1.0-(norm_diff-tolerance*0.5)*one_over_tol);
}

Real LocalDistanceTest(const EntityView& mdl, const EntityView& ref,
                       Real max_dist, Real tolerance, bool only_complete)
{
  Real max_sqr=pow(max_dist+max_dist*tolerance, 2);
  ContactList contacts;
  ResidueViewList residues=ref.GetResidueList();
  Real in_ref=0.0;
  Real in_mdl=0.0;
  for (ResidueViewList::const_iterator
    i=residues.begin(), e=residues.end(); i!=e; ++i) {
    AtomView atom_a=i->FindAtom("CA");
    if (!atom_a.IsValid()) {
      continue;
    }
    for (ResidueViewList::const_iterator
         j=i+1, e2=residues.end(); j!=e2; ++j) {
      AtomView atom_b=j->FindAtom("CA");
      if (!atom_b.IsValid()) {
        continue;
      }
      Real d=geom::Length2(atom_a.GetPos()-atom_b.GetPos());
      if (max_sqr<d) {
        continue;
      }
      if (!only_complete) {
        in_ref+=1.0;
      }

      AtomView atom_c=mdl.FindAtom(atom_a.GetResidue().GetChain().GetName(),
                                   atom_a.GetResidue().GetNumber(), "CA");
      if (!atom_c.IsValid()) {
        continue;
      }
      AtomView atom_d=mdl.FindAtom(atom_b.GetResidue().GetChain().GetName(),
                                   atom_b.GetResidue().GetNumber(), "CA");
      if (!atom_d.IsValid()) {
        continue;
      }
      if (only_complete) {
        in_ref+=1.0;
      }
      in_mdl+=Reward(sqrt(d), geom::Length(atom_c.GetPos()-atom_d.GetPos()),
                     tolerance);
    }
  }
  if (in_ref==0.0) {
    return 0.0;
  }
  return in_mdl/in_ref;
}

namespace {
void calculate_com(const ost::seq::AlignmentHandle& target_tpl_aln,
                   Real tolerance, int idx_a, int idx_b,
                   SimilarityMatrix& overlap, OverlapMode mode,
                   Real max_dist)
{
  ConstSequenceHandle s1=target_tpl_aln.GetSequence(idx_a);
  ConstSequenceHandle s2=target_tpl_aln.GetSequence(idx_b);
  if (!(s1.HasAttachedView() && s2.HasAttachedView())) {
    throw Error("Both sequences must have a view attached");
  }
  EntityView av=s1.GetAttachedView();
  for (int x=0; x<s1.GetLength(); ++x) {
    if (s1[x]=='-' || s2[x]=='-') {
      continue;
    }
    ResidueView res_a=s1.GetResidue(x);
    if (!res_a) {
      continue;
    }
    AtomView ca_a=res_a.FindAtom("CA");
    if (!ca_a) {
      continue;
    }
    for (int y=x+1; y<s1.GetLength(); ++y) {
      if (s1[y]=='-' || s2[y]=='-') {
        continue;
      }
      ResidueView res_b=s1.GetResidue(y);
      AtomView ca_b=res_b.FindAtom("CA");
      if (!ca_b.IsValid()) {
        continue;
      }
      ResidueView res_c=s2.GetResidue(x);
      AtomView ca_c=res_c.FindAtom("CA");
      if (!ca_c.IsValid()) {
        continue;
      }
      ResidueView res_d=s2.GetResidue(y);
      AtomView ca_d=res_d.FindAtom("CA");
      if (!ca_d.IsValid()) {
        continue;
      }
      Real da=geom::Length(ca_a.GetPos()-ca_b.GetPos());
      Real db=geom::Length(ca_c.GetPos()-ca_d.GetPos());
      if (da>max_dist && db>max_dist) {
        continue;
      }
      Real reward=Reward(da, db, tolerance, mode);
      if (overlap.Get(x, y)>=0) {
        reward=std::min(reward, overlap.Get(x, y));
      }
      overlap.Set(x, y, reward);
    }
  }
}

}
DistanceMatrixPtr DistanceMatrix::FromAln(const AlignmentHandle& aln,
                                          int index) {
  if (!aln.GetSequence(index).HasAttachedView()) {
    throw Error("Sequence must have a view attached");
  }
  return DistanceMatrixPtr(new DistanceMatrix(aln, index));

}

DistanceMatrixPtr DistanceMatrix::Create(const ConstSequenceHandle& s)
{
  if (!s.HasAttachedView()) {
    throw Error("Sequence must have a view attached");
  }
  return DistanceMatrixPtr(new DistanceMatrix(s));
}

SimilarityMatrix ContactOverlapMap(const DistanceMatrix& dmat1,
                                   const DistanceMatrix& dmat2,
                                   Real tolerance, OverlapMode mode,
                                   Real max_dist)
{
  if (dmat1.GetSize()!=dmat2.GetSize()) {
    throw Error("Distance matrices must have equal size");
  }
  SimilarityMatrix sim(dmat1.GetSize(), -1.0);
  for (int i=0; i<dmat1.GetSize(); ++i) {
    for (int j=i+1; j<dmat1.GetSize(); ++j) {
      Real a=dmat1.Get(i, j);
      Real b=dmat2.Get(i, j);
      if ((a<0 || b<0.0) || (a>max_dist && a>max_dist)) {
        continue;
      }
      sim(i, j)=Reward(a, b, tolerance, mode);
    }
  }
  return sim;
}

DistanceMatrix::DistanceMatrix(const AlignmentHandle& a, int index):
  ost::TriMatrix<Real>(a.GetSequence(0).GetGaplessString().size(), -1.0)
{
   EntityView av=a.GetSequence(index).GetAttachedView();
   ConstSequenceHandle trg_seq = a.GetSequence(0);
   ConstSequenceHandle tpl_seq = a.GetSequence(index);
   for (int x=0; x<a.GetLength(); ++x) {
     if (trg_seq[x]=='-' || tpl_seq[x]=='-') {
       continue;
     }
     ResidueView res_a=tpl_seq.GetResidue(x);
     if (!res_a) {
       continue;
     }
     AtomView ca_a=res_a.FindAtom("CA");
     if (!ca_a) {
       continue;
     }
     for (int y=x; y<a.GetLength(); ++y) {
       if (trg_seq[y]=='-' || tpl_seq[y] =='-') {
         continue;
       }
       ResidueView res_b=tpl_seq.GetResidue(y);
       if (!res_b.IsValid()) {
         continue;
       }
       AtomView ca_b=res_b.FindAtom("CA");
       if (!ca_b.IsValid()) {
         continue;
       }
       this->Set(trg_seq.GetResidueIndex(x), 
                 trg_seq.GetResidueIndex(y), 
                 geom::Length(ca_a.GetPos()-ca_b.GetPos()));
     }
   }
}

DistanceMatrix::DistanceMatrix(const ConstSequenceHandle& s):
  ost::TriMatrix<Real>(s.GetLength(), -1.0)
{
   EntityView av=s.GetAttachedView();
   for (int x=0; x<s.GetLength(); ++x) {
     if (s[x]=='-') {
       continue;
     }
     ResidueView res_a=s.GetResidue(x);
     if (!res_a) {
       continue;
     }
     AtomView ca_a=res_a.FindAtom("CA");
     if (!ca_a) {
       continue;
     }
     for (int y=x; y<s.GetLength(); ++y) {
       if (s[y]=='-') {
         continue;
       }
       ResidueView res_b=s.GetResidue(y);
       if (!res_b.IsValid()) {
         continue;
       }
       AtomView ca_b=res_b.FindAtom("CA");
       if (!ca_b.IsValid()) {
         continue;
       }
       this->Set(x, y, geom::Length(ca_a.GetPos()-ca_b.GetPos()));
     }
   }
}

SimilarityMatrix ContactOverlapMap(const AlignmentHandle& target_tpl_aln,
                                   Real tolerance, OverlapMode mode, Real max_dist)
{
  SimilarityMatrix overlap(target_tpl_aln.GetLength(), -1.0);
  for (int i=1; i<target_tpl_aln.GetCount(); ++i) {
    calculate_com(target_tpl_aln, tolerance, 0, i, overlap, mode, max_dist);
  }
  return overlap;
}

SimilarityMatrix ContactOverlapMap(const AlignmentHandle& target_tpl_aln,
                                   Real tolerance,
                                   int idx_a, int idx_b)
{
    SimilarityMatrix overlap(target_tpl_aln.GetLength(), -1.0);
    calculate_com(target_tpl_aln, tolerance, idx_a, idx_b, overlap,
                  OM_ABSOLUTE, 1e6);
    return overlap;
}


Real ContactOverlap(const DistanceMatrixPtr& lhs,
                    const DistanceMatrixPtr& rhs,
                    Real tolerance, OverlapMode mode,
                    Real dist_limit)
{
  Real nom=0.0, denom=0.0;
  for (int i=0; i<lhs->GetSize(); ++i) {
    for (int j=i+1; j<lhs->GetSize(); ++j) {
      Real a=lhs->Get(i, j);
      Real b=rhs->Get(i, j);
      if (a>dist_limit && b>dist_limit) {
        continue;
      }
      if (a>=0.0 && b>=0.0) {
        denom+=1.0;
        nom+=Reward(a, b, tolerance, mode);
      }
    }
  }
  return denom>0 ? nom/denom : 0.0;
}

int LabelComponents(const SimilarityMatrix& sim, std::vector<int>& components,
                     int x, int comp_id, Real tolerance)
{
  int count=1;
  components[x]=comp_id;
  for (int y=x+1; y<sim.GetSize(); ++y) {
    if (sim.Get(x, y)<tolerance || components[y]>-1) {
      continue;
    }
    count+=LabelComponents(sim, components, y, comp_id, tolerance);
  }
  return count;
}


//img::ImageHandle 
void DomainsFromCOM(const SimilarityMatrix& sim,
                                AlignmentHandle aln,
                                Real tolerance, int idx_a, int idx_b)
{
  const static int MAX_ITERATIONS=20;
  const static Real EPSILON=1e-6;
  SimilarityMatrix sim1=sim;
  SimilarityMatrix sim2(sim1.GetSize(), 0.0);
  SimilarityMatrix* sim1p=&sim1;
  SimilarityMatrix* sim2p=&sim2;
  img::ImageHandle im/*=img::CreateImage(img::Size(sim1.GetSize(),
                                                 sim1.GetSize(), MAX_ITERATIONS+1))*/;

  /*for (int i=0; i<sim.GetSize(); ++i) {
   for (int j=i+1; j<sim.GetSize(); ++j) {
     im.SetReal(img::Point(i, j, MAX_ITERATIONS), sim1.Get(i, j));
   }
  }*/
  // iteratively update similarity based on direct neighbours
  int n=MAX_ITERATIONS;
  std::vector<bool> fixed_residues(sim1.GetSize(), false);
  Real last_error=0.0;
  while (n>0) {
    --n;
    Real error=0.0;
    for (int i=0; i<sim.GetSize(); ++i) {
      /*if (fixed_residues[i]) {
        continue;
      }*/
      fixed_residues[i]=true;
      for (int j=i+1; j<sim.GetSize(); ++j) {
        int nom=0, denom=0;
        for (int k=0; k<sim.GetSize(); ++k) {
          Real ik=sim1p->Get(i, k);
          Real jk=sim1p->Get(j, k);
          if (jk<0 || ik <0) {
            continue;
          }
          if (ik>tolerance) {
            denom+=1;
            if (jk>tolerance) {
              nom+=1;
            }
          } else if (jk>tolerance) {
            denom+=1;
          }
        }
        if (denom) {
          Real newval=Real(nom)/denom;
          if (!(newval>0.95 || newval<0.05)) {
            fixed_residues[i]=false;
          }
          sim2p->Set(i, j, newval);
        } else {
          sim2p->Set(i, j, 0);
        }
        error+=std::abs(sim2p->Get(i, j)-sim1p->Get(i, j));
        //sim.SetReal(img::Point(i, j, n), sim2p->Get(i, j));
      }
    }
    error/=sim.GetSize()*(sim.GetSize()-1)/2;
    std::swap(sim1p, sim2p);
    if (std::abs(error-last_error)<EPSILON) {
      break;
    }
    last_error=error;
  }
  // run connected component algorithm on similarity matrix to identify
  // conserved domains.
  std::vector<int> components(sim.GetSize(), -1);
  int comp_id=0;
  for (int i=0; i<sim.GetSize(); ++i) {
    if (components[i]!=-1)  {
      continue;
    }
    if (LabelComponents(*sim2p, components, i, comp_id, tolerance)>1) {
      comp_id+=1;
    } {
      components[i]=-1;
    }
  }

  for (int i=0; i<aln.GetLength(); ++i) {
    if (idx_a>=0 && idx_b>=0) {
      if (ResidueView r=aln.GetResidue(idx_a, i)) {
        r.SetIntProp("comp_id", components[i]);
      }
      if (ResidueView r=aln.GetResidue(idx_b, i)) {
        r.SetIntProp("comp_id", components[i]);
      }
    } else {
      for (int j=0; j<aln.GetCount(); ++j) {
        if (ResidueView r=aln.GetResidue(j, i)) {
          r.SetIntProp("comp_id", components[i]);
        }
      }
    }

  }
}

}}}
