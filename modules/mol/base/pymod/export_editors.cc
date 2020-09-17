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

#include <ost/mol/mol.hh>
#include <ost/mol/impl/entity_impl.hh>
#include <ost/mol/impl/chain_impl.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/atom_impl.hh>
using namespace ost;
using namespace ost::mol;

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

BondHandle (EditorBase::*connect_a)(const AtomHandle&, 
                                    const AtomHandle&)=&EditorBase::Connect;
BondHandle (EditorBase::*connect_c)(const AtomHandle&, 
                                    const AtomHandle&,
                                    unsigned char)=&EditorBase::Connect;
BondHandle (EditorBase::*connect_b)(const AtomHandle&, const AtomHandle&, 
                                    Real, Real, Real)=&EditorBase::Connect;
BondHandle (EditorBase::*connect_d)(const AtomHandle&, const AtomHandle&, 
                                    Real, Real, Real,
                                    unsigned char)=&EditorBase::Connect;

ChainHandle (EditorBase::*insert_chain_a)(const String& chain_name)=&EditorBase::InsertChain;
ChainHandle (EditorBase::*insert_chain_b)(const String& chain_name,
                                          ChainHandle chain, bool deep)=&EditorBase::InsertChain;
                                    
ResidueHandle (EditorBase::*append_a)(ChainHandle ch, 
                                      const ResidueKey&)=&EditorBase::AppendResidue;
ResidueHandle (EditorBase::*append_b)(ChainHandle ch, const ResidueKey&, 
                                      const ResNum&)=&EditorBase::AppendResidue;
ResidueHandle (EditorBase::*append_c)(ChainHandle ch, ResidueHandle residue,
                                      bool deep)=&EditorBase::AppendResidue;

AtomHandle (EditorBase::*insert_atom_a)(ResidueHandle residue, const String& name,
                                        const geom::Vec3& pos, const String& ele,
                                        Real occupancy, Real b_factor,
                                        bool is_hetatm)=&EditorBase::InsertAtom;
AtomHandle (EditorBase::*insert_atom_b)(ResidueHandle residue, AtomHandle atom)=&EditorBase::InsertAtom;


AtomHandle (EditorBase::*insert_alt_atom_a)(ResidueHandle residue, const String& name,
                                            const String& alt_group, const geom::Vec3& pos,
                                            const String& ele, Real occ,
                                            Real b_factor)=&EditorBase::InsertAltAtom;
AtomHandle (EditorBase::*insert_alt_atom_b)(ResidueHandle residue, AtomHandle atom,
                                            const String& alt_group)=&EditorBase::InsertAltAtom;



void (ICSEditor::*set_torsion_a)(TorsionHandle, Real, bool)=&ICSEditor::SetTorsionAngle;
void (ICSEditor::*set_torsion_b)(const AtomHandle&, const AtomHandle&,
                                 const AtomHandle&, const AtomHandle&,
                                 Real, bool)=&ICSEditor::SetTorsionAngle;

void (ICSEditor::*rotate_torsion_a)(TorsionHandle, Real, bool)=&ICSEditor::RotateTorsionAngle;
void (ICSEditor::*rotate_torsion_b)(const AtomHandle&, const AtomHandle&,
                                    const AtomHandle&, const AtomHandle&,
                                    Real, bool)=&ICSEditor::RotateTorsionAngle;

void (EditorBase::*renumber_chain_a)(ChainHandle,const ResNumList&)=&EditorBase::RenumberChain;
void (EditorBase::*renumber_chain_b)(const ChainHandle&,int, bool)=&EditorBase::RenumberChain;
#if OST_NUMPY_SUPPORT_ENABLED
template<typename T, bool O>
void set_pos2_nc_t(XCSEditor& e, const AtomHandleList& alist, PyArrayObject* na)
{
  size_t count=0;
  for(AtomHandleList::const_iterator ait=alist.begin();ait!=alist.end();++ait,++count) {
    if(O) {
      e.SetAtomOriginalPos(*ait,geom::Vec3(static_cast<Real>(*reinterpret_cast<T*>(PyArray_GETPTR2(na,count,0))),
                                           static_cast<Real>(*reinterpret_cast<T*>(PyArray_GETPTR2(na,count,1))),
                                           static_cast<Real>(*reinterpret_cast<T*>(PyArray_GETPTR2(na,count,2)))));
    } else {
      e.SetAtomTransformedPos(*ait,geom::Vec3(static_cast<Real>(*reinterpret_cast<T*>(PyArray_GETPTR2(na,count,0))),
                                              static_cast<Real>(*reinterpret_cast<T*>(PyArray_GETPTR2(na,count,1))),
                                              static_cast<Real>(*reinterpret_cast<T*>(PyArray_GETPTR2(na,count,2)))));
    }
  }
}

template<bool O>
void set_pos2_t(XCSEditor& e, const AtomHandleList& alist, object pyobj)
{
  size_t acount = alist.size();
  
  if(!PyArray_Check(pyobj.ptr())) {
    throw Error("expected a numpy array");
    return;
  }
  PyArrayObject* na=reinterpret_cast<PyArrayObject*>(pyobj.ptr());
  
  if(PyArray_NDIM(na)!=2 || PyArray_DIM(na,0)!=int(acount) || PyArray_DIM(na,1)!=3) {
    throw Error("expected a numpy array of shape (NAtoms, 3)");
    return;
  }
  
  if(PyArray_ISCONTIGUOUS(na)) {
    if(PyArray_TYPE(na)==NPY_FLOAT) {
      if(O) {
        e.SetAtomOriginalPos(alist,reinterpret_cast<float*>(PyArray_DATA(na)));
      } else {
        e.SetAtomTransformedPos(alist,reinterpret_cast<float*>(PyArray_DATA(na)));
      }
    } else if(PyArray_TYPE(na)==NPY_DOUBLE) {
      if(O) {
        e.SetAtomOriginalPos(alist,reinterpret_cast<double*>(PyArray_DATA(na)));
      } else {
        e.SetAtomTransformedPos(alist,reinterpret_cast<double*>(PyArray_DATA(na)));
      }
    } else {
      throw Error("expected a numpy array of type float or double");
      return;
    }
  } else {
    // non-contiguous
#if 0
    throw Error("expected contiguous numpy array");
#else
    if(PyArray_TYPE(na)==NPY_FLOAT) {
      set_pos2_nc_t<float,O>(e,alist,na);
    } else if(PyArray_TYPE(na)==NPY_DOUBLE) {
      set_pos2_nc_t<double,O>(e,alist,na);
    } else {
      throw Error("expected a numpy array of type float or double");
      return;
    }
#endif
  }
}
#endif

void set_pos(XCSEditor& e, object o1, object o2, bool trans)
{
  extract<AtomHandle> eah(o1);
  extract<geom::Vec3> ev3(o2);
  if(eah.check() && ev3.check()) {
    if(trans) {
      e.SetAtomTransformedPos(eah(),ev3());
    } else {
      e.SetAtomOriginalPos(eah(),ev3());
    }
    return;
  }

#if OST_NUMPY_SUPPORT_ENABLED

  extract<AtomHandleList> eal(o1);
  if(eal.check()) {
    if(trans) {
      set_pos2_t<false>(e,eal(),o2);
    } else {
      set_pos2_t<true>(e,eal(),o2);
    }
    return;
  }

  std::map<unsigned long,AtomHandle> amap;
  impl::EntityImplPtr ei=e.GetEntity().Impl();
  for(impl::ChainImplList::iterator cit=ei->GetChainList().begin();
      cit!=ei->GetChainList().end();++cit) {
    for (impl::ResidueImplList::iterator rit = (*cit)->GetResidueList().begin(),
         ret = (*cit)->GetResidueList().end(); rit != ret; ++rit) {
           
      for (impl::AtomImplList::iterator ait = (*rit)->GetAtomList().begin(), 
           aet = (*rit)->GetAtomList().end(); ait != aet; ++ait) {

        amap[(*ait)->GetIndex()]=*ait;
      }
    }
  }

  AtomHandleList alist;
  for(int i=0;i<len(o1);++i) {
    int gid = extract<int>(o1[i]);
    std::map<unsigned long,AtomHandle>::iterator ait=amap.find(static_cast<unsigned long>(gid));
    alist.push_back(ait==amap.end() ? AtomHandle() : ait->second);
  }

  if(trans) {
    set_pos2_t<false>(e,alist,o2);
  } else {
    set_pos2_t<true>(e,alist,o2);
  }

#else
  throw Error("SetAtom*Pos(...,ndarray) not available, because numpy support not compiled in");
#endif
}

void set_o_pos(XCSEditor& e, object o1, object o2)
{
  set_pos(e,o1,o2,false);
}

void set_t_pos(XCSEditor& e, object o1, object o2)
{
  set_pos(e,o1,o2,true);
}

}

void export_Editors()
{
#if OST_NUMPY_SUPPORT_ENABLED
  // The following define enforces no return value when calling import_array
  #undef NUMPY_IMPORT_ARRAY_RETVAL
  #define NUMPY_IMPORT_ARRAY_RETVAL
  import_array();
#endif


  class_<EditorBase>("EditorBase", no_init)
    .def("InsertChain", insert_chain_a)
    .def("InsertChain", insert_chain_b,(arg("chain_name"),arg("chain"), arg("deep")=false))
    .def("InsertAtom", insert_atom_a,
         (arg("residue"), arg("name"), arg("pos"), arg("element")="", 
          arg("occupancy")=1.0, arg("b_factor")=0.0, arg("is_hetatm")=false))
    .def("InsertAtom", insert_atom_b)
    .def("InsertAltAtom", insert_alt_atom_a)
    .def("InsertAltAtom", insert_alt_atom_b)
    .def("DeleteResidue", &EditorBase::DeleteResidue)
    .def("DeleteChain", &EditorBase::DeleteChain)
    .def("DeleteAtom", &EditorBase::DeleteAtom)
    .def("DeleteBond", &EditorBase::DeleteBond)
    .def("DeleteBonds", &EditorBase::DeleteBonds)
    .def("InsertResidueBefore", &EditorBase::InsertResidueBefore)
    .def("InsertResidueAfter", &EditorBase::InsertResidueAfter)
    .def("AppendResidue", append_a)
    .def("AppendResidue", append_b)
    .def("AppendResidue", append_c, (arg("chain"), arg("residue"),
                                     arg("deep")=false))
    .def("Connect", connect_a)
    .def("Connect", connect_b)    
    .def("Connect", connect_c)
    .def("Connect", connect_d)    
    .def("RenameChain", &EditorBase::RenameChain)
    .def("SetChainType", &EditorBase::SetChainType)
    .def("SetChainDescription", &EditorBase::SetChainDescription)
    .def("RenameResidue", &EditorBase::RenameResidue)
    .def("SetResidueNumber", &EditorBase::SetResidueNumber)
    .def("RenameAtom", &EditorBase::RenameAtom)
    .def("AddTorsion", &EditorBase::AddTorsion)
    .def("ReorderResidues",&EditorBase::ReorderResidues)
    .def("ReorderAllResidues",&EditorBase::ReorderAllResidues)
    .def("RenumberAllResidues",&EditorBase::RenumberAllResidues)
    .def("RenumberChain",renumber_chain_a)
    .def("RenumberChain",renumber_chain_b)
    .def("Prune", &EditorBase::Prune)
  ;
  
  void (XCSEditor::*apply_transform1)(const geom::Mat4&) = &XCSEditor::ApplyTransform;
  void (XCSEditor::*apply_transform2)(const geom::Transform&) = &XCSEditor::ApplyTransform;
  void (XCSEditor::*set_transform1)(const geom::Mat4&) = &XCSEditor::SetTransform;
  void (XCSEditor::*set_transform2)(const geom::Transform&) = &XCSEditor::SetTransform;

  class_<XCSEditor, bases<EditorBase> >("XCSEditor", no_init)
    .def("SetAtomPos", set_t_pos)
    .def("SetAtomTransformedPos", set_t_pos)
    .def("SetAtomOriginalPos", set_o_pos)
    .def("ApplyTransform", apply_transform1)
    .def("ApplyTransform", apply_transform2)
    .def("SetTransform", set_transform1)
    .def("SetTransform", set_transform2)
    .def("UpdateICS", &XCSEditor::UpdateICS)
    .def("ForceUpdate", &XCSEditor::ForceUpdate)
    .def("__exit__", &XCSEditor::ForceUpdate)    
    .def("__del__", &XCSEditor::ForceUpdate)    
  ;
  
  class_<ICSEditor, bases<EditorBase> >("ICSEditor", no_init)
    .def("SetAngle", &ICSEditor::SetAngle)
    .def("SetBondLength", &ICSEditor::SetBondLength)
    .def("SetTorsionAngle", set_torsion_a,
          (arg("torsion"),arg("angle"), arg("update_others")=true))
    .def("SetTorsionAngle", set_torsion_b,
         (arg("atom_a"), arg("atom_b"), arg("atom_c"), arg("atom_d"),
          arg("angle"), arg("update_others")))
    .def("RotateTorsionAngle", rotate_torsion_a,
          (arg("torsion"),arg("angle"), arg("update_others")=true))
    .def("RotateTorsionAngle", rotate_torsion_b,
         (arg("atom_a"), arg("atom_b"), arg("atom_c"), arg("atom_d"),
          arg("angle"), arg("update_others")))
    .def("UpdateXCS", &ICSEditor::UpdateXCS)
    .def("__exit__", &ICSEditor::UpdateXCS)
  ;  
}
