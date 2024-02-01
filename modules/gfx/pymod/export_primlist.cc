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
using namespace boost::python;

#include <boost/shared_ptr.hpp>

#include <ost/message.hh>
#include <ost/gfx/prim_list.hh>
using namespace ost;
using namespace ost::gfx;

/* The following #include triggers a deprecation warning upon compilation:

In file included from /unibas/lcs-software/software/Python/2.7.11-goolf-1.7.20/lib/python2.7/site-packages/numpy/core/include/numpy/ndarraytypes.h:1781:0,
                 from /unibas/lcs-software/software/Python/2.7.11-goolf-1.7.20/lib/python2.7/site-packages/numpy/core/include/numpy/ndarrayobject.h:18,
                 from /unibas/lcs-software/software/Python/2.7.11-goolf-1.7.20/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h:4,
                 from /home/zohixe92/build/ost-develop/modules/gfx/pymod/export_primlist.cc:30:
/unibas/lcs-software/software/Python/2.7.11-goolf-1.7.20/lib/python2.7/site-packages/numpy/core/include/numpy/npy_1_7_deprecated_api.h:15:2: warning: #warning "Using deprecated NumPy API, disable it by " "#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION" [-Wcpp]
 #warning "Using deprecated NumPy API, disable it by " \
  ^

It looks like just including <numpy/arrayobject.h>, without actually using anything from it, already triggers the
warning. This could be a side effect of other includes in export_primlist.cc though, such as boost.

The Numpy C API documentation <https://docs.scipy.org/doc/numpy-1.15.1/reference/c-api.deprecations.html> is very
unclear about what is deprecated exactly (the whole API or only direct data access?) The Array API doc
<https://docs.scipy.org/doc/numpy-1.16.0/reference/c-api.array.html> doesn't clearly state what to import instead,
and while some functions are marked as deprecated I couldn't see that we are using any of them. Last, the migration sed
script <https://github.com/numpy/numpy/blob/master/tools/replace_old_macros.sed> didn't trigger any changes.

The warnings appear to be safe though, behavior has been checked in <https://jira.biozentrum.unibas.ch/browse/SCHWED-3149>
and Numpy guarantees to maintain them until the next major release (2.0). So it is safe to ignore them for now.
They can be silenced with the following #define, however it was decided to keep them for now as functionality *will*
disappear in the future.

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

In the future it might be desirable to move to the newer API in Boost::Python::NumPy
which uses modern C++ instead of old-style C macros.
<https://www.boost.org/doc/libs/1_64_0/libs/python/doc/html/numpy/reference/ndarray.html>
*/
#if OST_NUMPY_SUPPORT_ENABLED
#include <numpy/arrayobject.h>
#endif

namespace {
  void add_mesh(PrimList& p, object ova, object ona, object oca, object oia)
  {
#if OST_NUMPY_SUPPORT_ENABLED
    if(!PyArray_Check(ova.ptr())) {
      throw Error("ova is not a numpy array");
    }
    PyArrayObject* va=reinterpret_cast<PyArrayObject*>(ova.ptr());
    if(!PyArray_ISCONTIGUOUS(va)) {
      throw Error("expected vertex array to be contiguous");
    }
    if(PyArray_TYPE(va)!=NPY_FLOAT) {
      throw Error("expected vertex array to be of dtype=float32");
    }
    size_t v_size=PyArray_SIZE(va);
    if(v_size%3!=0) {
      throw Error("expected vertex array size to be divisible by 3");
    }
    size_t v_count=v_size/3;
    float* vp=reinterpret_cast<float*>(PyArray_DATA(va));
    float* np=0;
    float* cp=0;
    if(!ona.is_none()) {
      if(!PyArray_Check(ona.ptr())) {
        throw Error("ona is not a numpy array");
      }
      PyArrayObject* na=reinterpret_cast<PyArrayObject*>(ona.ptr());
      if(!PyArray_ISCONTIGUOUS(na)) {
        throw Error("expected normal array to be contiguous");
      }
      if(PyArray_TYPE(na)!=NPY_FLOAT) {
        throw Error("expected normal array to be of dtype=float32");
      }
      if((size_t)PyArray_SIZE(na)!=v_size) {
        throw Error("expected normal array size to match vertex array size");
      }
      np=reinterpret_cast<float*>(PyArray_DATA(na));
    }
    if(!oca.is_none()) {
      if(!PyArray_Check(oca.ptr())) {
        throw Error("oca is not a numpy array");
      }
      PyArrayObject* ca=reinterpret_cast<PyArrayObject*>(oca.ptr());
      if(!PyArray_ISCONTIGUOUS(ca)) {
        throw Error("expected color array to be contiguous");
      }
      if(PyArray_TYPE(ca)!=NPY_FLOAT) {
        throw Error("expected color array to be of dtype=float32");
      }
      if((size_t)PyArray_SIZE(ca)!=v_count*4) {
        throw Error("expected color array size to equal vertex-count x 4");
      }
      cp=reinterpret_cast<float*>(PyArray_DATA(ca));
    }
    if(!PyArray_Check(oia.ptr())) {
      throw Error("oia is not a numpy array");
    }
    PyArrayObject* ia=reinterpret_cast<PyArrayObject*>(oia.ptr());
    if(!PyArray_ISCONTIGUOUS(ia)) {
      throw Error("expected vertex array to be contiguous");
    }
    if(PyArray_TYPE(ia)!=NPY_UINT) {
      throw Error("expected vertex array to be of dtype=uint32");
    }
    size_t i_size=PyArray_SIZE(ia);
    unsigned int* ip=reinterpret_cast<unsigned int*>(PyArray_DATA(ia));

    p.AddMesh(vp,np,cp,v_count,ip,i_size/3);
#else
    throw Error("AddMesh requires compiled-in numpy support");
#endif
  }
}

void export_primlist()
{
#if OST_NUMPY_SUPPORT_ENABLED
  // The following define enforces no return value when calling import_array
  #define NUMPY_IMPORT_ARRAY_RETVAL
  import_array(); // magic handshake for numpy module
#endif

  class_<PrimList, bases<GfxObj>, boost::shared_ptr<PrimList>, boost::noncopyable>("PrimList", init<const String& >())
    .def("Clear",&PrimList::Clear)
    .def("_add_line",&PrimList::AddLine)
    .def("_add_point",&PrimList::AddPoint)
    .def("_add_sphere",&PrimList::AddSphere)
    .def("_add_cyl",&PrimList::AddCyl)
    .def("_add_text",&PrimList::AddText)
    .def("AddMesh",add_mesh)
    .def("SetColor",&PrimList::SetColor)
    .def("SetDiameter",&PrimList::SetDiameter)
    .def("SetRadius",&PrimList::SetRadius)
  ;

}
