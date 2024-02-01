//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/mol/entity_view.hh>
#include <ost/mol/impl/entity_impl.hh>
#include <ost/mol/impl/chain_impl.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/atom_impl.hh>

#include <ost/mol/query.hh>
#include <ost/mol/mol.hh>
#include "bounds.hh"

using namespace ost;
using namespace ost::mol;

#include <ost/export_helper/generic_property_def.hh>

/* Including NumPy headers produces compiler warnings. The ones about "Using
   deprecated NumPy API..." we can not get rid of. The highest NumPy version we
   support is 1.6 while the non-deprecated API starts with version 1.7.
   Also see the comment in modules/gfx/pymod/export_primlist.cc for further
   information.
*/
#if OST_NUMPY_SUPPORT_ENABLED
#include <numpy/numpyconfig.h>
#define NPY_NO_DEPRECATED_API NPY_1_6_API_VERSION
#include <numpy/arrayobject.h>
#endif

namespace {
EntityHandle create1() {  return CreateEntity(); }

typedef EntityView (EntityHandle::*StringMethod)(const String&, uint) const;
typedef BondHandle (EntityHandle::*Connect1)(const AtomHandle&,
                                             const AtomHandle&);
typedef BondHandle (EntityHandle::*Connect2)(const AtomHandle&,
                                             const AtomHandle&, Real,
                                             Real, Real);

typedef EntityView (EntityHandle::*QueryMethod)(const Query&, uint) const;

Real (EntityHandle::*get_angle1)(const AtomHandle&, const AtomHandle&, const AtomHandle&) const = &EntityHandle::GetAngle;
Real (EntityHandle::*get_angle2)(const AtomView&, const AtomView&, const AtomView&) const = &EntityHandle::GetAngle;
                                      
XCSEditor depr_request_xcs_editor(EntityHandle e, EditMode m)
{
  WARN_DEPRECATED("EntityHandle::RequestXCSEditor is deprecated. Use "
                  "EntityHandle::EditXCS instead");
  return e.EditXCS(m);
}

ICSEditor depr_request_ics_editor(EntityHandle e, EditMode m)
{
  WARN_DEPRECATED("EntityHandle::RequestICSEditor is deprecated. Use "
                  "EntityHandle::EditICS instead");
  return e.EditICS(m);
}


#if OST_NUMPY_SUPPORT_ENABLED

bool less_index(const mol::AtomHandle& a1, const mol::AtomHandle& a2)
{
  return a1.GetIndex()<a2.GetIndex();
}
PyObject* get_pos2(EntityHandle& entity, bool id_sorted)
{
  npy_intp dims[]={entity.GetAtomCount(),3};
  PyObject* na = PyArray_SimpleNew(2,dims,NPY_FLOAT);
  npy_float* nad = reinterpret_cast<npy_float*>(PyArray_DATA(na));
  if(id_sorted) {
    AtomHandleList alist = entity.GetAtomList();
    std::sort(alist.begin(),alist.end(),less_index);
    for(AtomHandleList::const_iterator it=alist.begin();it!=alist.end();++it,nad+=3) {
      geom::Vec3 pos=(*it).GetPos();
      nad[0]=static_cast<npy_float>(pos[0]);
      nad[1]=static_cast<npy_float>(pos[1]);
      nad[2]=static_cast<npy_float>(pos[2]);
    }
  } else {
    impl::EntityImplPtr ei=entity.Impl();
    for(impl::ChainImplList::iterator cit=ei->GetChainList().begin();
        cit!=ei->GetChainList().end();++cit) {
      for (impl::ResidueImplList::iterator rit = (*cit)->GetResidueList().begin(),
          ret = (*cit)->GetResidueList().end(); rit != ret; ++rit) {
            
        for (impl::AtomImplList::iterator ait = (*rit)->GetAtomList().begin(), 
            aet = (*rit)->GetAtomList().end(); ait != aet; ++ait, nad+=3) {

          geom::Vec3 pos=(*ait)->TransformedPos();
          nad[0]=static_cast<npy_float>(pos[0]);
          nad[1]=static_cast<npy_float>(pos[1]);
          nad[2]=static_cast<npy_float>(pos[2]);
    }}}
  }
  return na;
}

PyObject* get_pos1(EntityHandle& entity)
{
  return get_pos2(entity,true);
}

#endif

geom::Mat4 depr_get_transformation_matrix(const EntityHandle& eh)
{
  return eh.GetTransformationMatrix();
}

bool depr_is_transformation_identity(const EntityHandle& eh)
{
  return eh.IsTransformationIdentity();
}

} // ns

void export_Entity()
{
#if OST_NUMPY_SUPPORT_ENABLED
  // The following define enforces no return value when calling import_array
  #undef NUMPY_IMPORT_ARRAY_RETVAL
  #define NUMPY_IMPORT_ARRAY_RETVAL
  import_array();
#endif

  class_<EntityBase> ent_base("EntityBase", no_init);
  ent_base
    .def(self_ns::str(self))
    .def("GetName", &EntityBase::GetName,
        return_value_policy<copy_const_reference>())
    .def("SetName", &EntityBase::SetName)
    .def("IsValid", &EntityBase::IsValid)
    .add_property("valid", &EntityBase::IsValid)
  ;
  generic_prop_def<EntityBase>(ent_base);

  EntityView (EntityHandle::*select1)(const Query&) const = &EntityHandle::Select;
  EntityView (EntityHandle::*select2)(const Query&, QueryFlags) const = &EntityHandle::Select;
  EntityView (EntityHandle::*select3)(const String&) const = &EntityHandle::Select;
  EntityView (EntityHandle::*select4)(const String&, QueryFlags) const = &EntityHandle::Select;
  
  class_<EntityHandle, bases<EntityBase> >("EntityHandle", init<>())
    .def("Select",select1)
    .def("Select",select2)
    .def("Select",select3)
    .def("Select",select4)
    .def("SetDefaultQueryFlags",&EntityHandle::SetDefaultQueryFlags)
    .def("GetDefaultQueryFlags",&EntityHandle::GetDefaultQueryFlags)
    .add_property("default_query_flags",
                  &EntityHandle::GetDefaultQueryFlags,
                  &EntityHandle::SetDefaultQueryFlags)
    .def("FindChain", &EntityHandle::FindChain)
    .def("FindResidue", &EntityHandle::FindResidue)
    .def("FindAtom", &EntityHandle::FindAtom)
    .def("GetChainList", &EntityHandle::GetChainList)
    .def("CreateFullView",&EntityHandle::CreateFullView)
    .def("CreateEmptyView",&EntityHandle::CreateEmptyView)
    .def("GetMass", &EntityHandle::GetMass)
    .def("GetCenterOfMass", &EntityHandle::GetCenterOfMass)
    .def("GetCenterOfAtoms", &EntityHandle::GetCenterOfAtoms)
    .def("GetAtomPosList", &EntityHandle::GetAtomPosList, arg("ordered_by_index")=false)
    .def("GetGeometricCenter", geom_center<EntityHandle>)
    .add_property("geometric_center", geom_center<EntityHandle>)

    .add_property("geometric_end", geom_end<EntityHandle>)
    .add_property("geometric_start", geom_start<EntityHandle>)
    .def("GetGeometricStart", geom_start<EntityHandle>)
    .def("GetGeometricEnd", geom_end<EntityHandle>)
    .def("GetBoundarySize", geom_size<EntityHandle>)
    .def("GetResidueCount", &EntityHandle::GetResidueCount)
    .def("GetAtomCount", &EntityHandle::GetAtomCount)
    .def("GetBondCount", &EntityHandle::GetBondCount)
    .def("GetChainCount", &EntityHandle::GetChainCount)
    .def("GetHandle", &EntityHandle::GetHandle)
    .add_property("handle", &EntityHandle::GetHandle)
    .add_property("chain_count", &EntityHandle::GetChainCount)
    .add_property("residue_count", &EntityHandle::GetResidueCount)
    .add_property("atom_count", &EntityHandle::GetAtomCount)
    .add_property("mass", &EntityHandle::GetMass)
    .add_property("center_of_mass", &EntityHandle::GetCenterOfMass)
    .add_property("center_of_atoms", &EntityHandle::GetCenterOfAtoms)
    .def("FindWithin", &EntityHandle::FindWithin)
    .def("GetAngle", get_angle1)
    .def("GetAngle", get_angle2)
    .def("FindTorsion", &EntityHandle::FindTorsion)
    .def("Copy", &EntityHandle::Copy)
    .def("GetBondList", &EntityHandle::GetBondList)
    .def("GetAtomList", &EntityHandle::GetAtomList)
    .def("GetResidueList", &EntityHandle::GetResidueList)
    .add_property("residues", &EntityHandle::GetResidueList)
    .add_property("atoms", &EntityHandle::GetAtomList)
    .add_property("chains", &EntityHandle::GetChainList)
    .add_property("bonds", &EntityHandle::GetBondList)
    .def("GetBounds", &EntityHandle::GetBounds)
    .add_property("bounds", &EntityHandle::GetBounds)
    .def("GetTransformationMatrix", depr_get_transformation_matrix)
    .def("IsTransformationIdentity",depr_is_transformation_identity)
    .def("GetTransform",&EntityHandle::GetTransform)
    .def("SetTransform",&EntityHandle::SetTransform)
    .add_property("transform",&EntityHandle::GetTransform,&EntityHandle::SetTransform)
    .def("HasTransform",&EntityHandle::HasTransform)
    .def("ClearTransform",&EntityHandle::ClearTransform)
    .def("FixTransform",&EntityHandle::FixTransform)
    .def("EditXCS", &EntityHandle::EditXCS, arg("mode")=UNBUFFERED_EDIT)
    .def("EditICS", &EntityHandle::EditICS, arg("mode")=UNBUFFERED_EDIT)
    .def("RequestXCSEditor", &depr_request_xcs_editor, arg("mode")=UNBUFFERED_EDIT)
    .def("RequestICSEditor", &depr_request_ics_editor, arg("mode")=UNBUFFERED_EDIT)
    .def("GetHashCode", &EntityHandle::GetHashCode) 
    .def(self==self)
    .def(self!=self)
#if OST_NUMPY_SUPPORT_ENABLED
    .def("GetPositions",get_pos1)
    .def("GetPositions",get_pos2)
    .add_property("positions",get_pos1)
#endif
  ;

  def("CreateEntity",create1);

  class_<EntityHandleList>("EntityHandleList", init<>())
    .def(vector_indexing_suite<EntityHandleList>())
  ;
}

