//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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

/*
  Author: Andreas Schenk
*/

#include <boost/python.hpp>

using namespace boost::python;

#include <ost/img/phase.hh>


void export_Phase()
{
  using namespace ost::img;
  using namespace ost;
  class_<Phase> ("Phase" ,init<>())
    .def(init<Complex>())//init<Complex> has to be before init<Real> because of implicit Real->Complex conversion
    .def(init<Real>())
    .def(init<Phase>())
    .def(self+self)
    .def(self-self)
    .def(self+=self)
    .def(self-=self)
    .def(self*=Real())
    .def(self/=Real())
    .def(self*Real())
    .def(self/Real())
    .def(float_(self)) 
    .def(self_ns::str(self))
    ;
}

