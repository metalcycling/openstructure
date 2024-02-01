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

namespace ost {

/// \brief helper to convert between python tuple and std::pair
/// 
/// Usage:
/// 
/// boost::python::to_python_converter<std::pair<T1, T2>, 
///                                    PairToTupleConverter<T1, T2> >()
template<class T1, class T2>
struct PairToTupleConverter {
  static PyObject* convert(const std::pair<T1, T2>& pair) {
    boost::python::tuple t=boost::python::make_tuple<T1,T2>(pair.first,
                                                            pair.second);
    return boost::python::incref(t.ptr());
  }
};
}

