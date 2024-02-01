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
#include <ost/mol/mm/index.hh>

using namespace boost::python;

namespace{

template<typename T>
int Get(T index, uint i){
  return index[i];
}

template<typename T>
void Set(T& index, uint i, uint v){
  index[i] = v;
}
}

void export_Index()
{
  class_<ost::mol::mm::Index<1> >("Index1", init<uint>())
    .def("__getitem__",Get<ost::mol::mm::Index<1> >)
    .def("__setitem__",Set<ost::mol::mm::Index<1> >)
  ;

  class_<ost::mol::mm::Index<2> >("Index2", init<uint, uint>())
    .def("__getitem__",Get<ost::mol::mm::Index<2> >)
    .def("__setitem__",Set<ost::mol::mm::Index<2> >)
  ;

  class_<ost::mol::mm::Index<3> >("Index3", init<uint, uint, uint>())
    .def("__getitem__",Get<ost::mol::mm::Index<3> >)
    .def("__setitem__",Set<ost::mol::mm::Index<3> >)
  ;

  class_<ost::mol::mm::Index<4> >("Index4", init<uint, uint, uint,uint>())
    .def("__getitem__",Get<ost::mol::mm::Index<4> >)
    .def("__setitem__",Set<ost::mol::mm::Index<4> >)
  ;

  class_<ost::mol::mm::Index<5> >("Index5", init<uint, uint, uint, uint, uint>())
    .def("__getitem__",Get<ost::mol::mm::Index<5> >)
    .def("__setitem__",Set<ost::mol::mm::Index<5> >)
  ;

  class_<ost::mol::mm::Index<6> >("Index6", init<uint, uint, uint, uint, uint, uint>())
    .def("__getitem__",Get<ost::mol::mm::Index<6> >)
    .def("__setitem__",Set<ost::mol::mm::Index<6> >)
  ;

  class_<ost::mol::mm::Index<7> >("Index7", init<uint, uint, uint, uint, uint, uint, uint>())
    .def("__getitem__",Get<ost::mol::mm::Index<7> >)
    .def("__setitem__",Set<ost::mol::mm::Index<7> >)
  ;
}
