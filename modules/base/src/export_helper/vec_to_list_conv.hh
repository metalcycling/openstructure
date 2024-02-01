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


#ifndef VEC_LIST_CONV_HH
#define VEC_LIST_CONV_HH

#include <boost/python.hpp>
#include <vector>
namespace ost{

/// \brief helper to convert between python list tuple and std::vecot
/// 
/// Usage:
///
/// boost::python::to_python_converter<std::vector<T>,
  ///                                  VectorToListConverter<T> >()
                        

template<typename T>
struct VectorToListConverter {
  static PyObject* convert(const std::vector<T>& vec) {
    boost::python::list l;
    for(typename std::vector<T>::const_iterator it=vec.begin();it!=vec.end();++it){
      l.append(*it);
    }
    return boost::python::incref(l.ptr());
  }
};

}

#endif
