#include <cmath>

#include <ost/message.hh>
#include <ost/mol/mol.hh>
#include <ost/log.hh>
#include "construct_cbeta.hh"

namespace ost { namespace mol { namespace alg {

namespace detail {



class CBetaInserter:public EntityVisitor {
public:
  CBetaInserter(mol::EntityHandle ent, bool include_gly);

  virtual bool VisitResidue(const mol::ResidueHandle& res);
private:
  XCSEditor edi_;
  bool      include_gly_;
};
CBetaInserter::CBetaInserter(mol::EntityHandle ent, bool include_gly) {
  edi_=ent.EditXCS(BUFFERED_EDIT);
  include_gly_=include_gly;
}

bool CBetaInserter::VisitResidue(const mol::ResidueHandle& res)
{
  if (!res.IsPeptideLinking()) {
    return false;
  }
  try { 
    if(res.GetOneLetterCode()=='G' && include_gly_==false) {
      return false;
    }
    if (res.FindAtom("CB").IsValid()) {
      return false;
    }
    geom::Vec3 cbeta_pos=CBetaPosition(res);
    edi_.InsertAtom(res, "CB", cbeta_pos, "C");
  }
  catch(...) {
    LOG_WARNING("residue " << res.GetQualifiedName() 
                << " doesn't have enough backbone atoms to reconstruct CBeta position");
    return false;
  }
  return false;
}

} // namespace


// Todo: Applying on an entity does not work since the handle is used and the constructed
// C-betas are there not present in the resulting entity-reference -> use handle onyl or
// ensure that C-betas are added to the view as well
void ConstructCBetas(mol::EntityHandle& entity_handle, bool include_gly)
{
  detail::CBetaInserter cbi(entity_handle, include_gly);
  entity_handle.Apply(cbi);
}


geom::Vec3 CBetaPosition(const geom::Vec3& n_pos, const geom::Vec3& ca_pos, 
                         const geom::Vec3& c_pos, Real bond_length) {

   geom::Vec3 v1=geom::Normalize(ca_pos-n_pos);
   geom::Vec3 v2=geom::Normalize(ca_pos-c_pos);
   geom::Vec3 in_plane_v=geom::Normalize(v1+v2);
   geom::Plane p(ca_pos ,n_pos, c_pos);
   // rotate around vector perpendicular  to p and in_plane_v
   geom::Vec3 axis=geom::Normalize(geom::Cross(p.GetNormal(), in_plane_v));
   geom::Mat3 rot_mat=geom::AxisRotation(axis, (-54/180.0)*M_PI);
   return ca_pos+rot_mat*in_plane_v*bond_length;                           
}


geom::Vec3 CBetaPosition(const ResidueHandle& residue, Real bond_length)
{
  AtomHandle ca=residue.FindAtom("CA");
  AtomHandle n=residue.FindAtom("N");
  AtomHandle c=residue.FindAtom("C");
  if (!(ca.IsValid() && c.IsValid() && n.IsValid())) {
  std::stringstream ss;
  ss << "residue " << residue.GetQualifiedName() 
     << " doesn't have enough backbone atoms to reconstruct CBeta position";
    throw ost::Error(ss.str());
  }
  return CBetaPosition(n.GetPos(), ca.GetPos(), c.GetPos(), bond_length);
}


}}} // ns


