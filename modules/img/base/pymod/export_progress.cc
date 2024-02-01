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
#include <ost/img/progress.hh>

using namespace boost::python;

//wrapper functions to facilitate passing of void ptr
void pregister1(ost::img::Progress&  p,PyObject* x, unsigned long steps, unsigned long delta) { p.Register(x,steps,delta); }
void pregister2(ost::img::Progress& p,PyObject* x, unsigned long steps) { p.Register(x,steps,1); }
void pderegister(ost::img::Progress& p,PyObject* x) { p.DeRegister(x); }
void padvance1(ost::img::Progress& p,PyObject* x,unsigned long step) { p.AdvanceProgress(x,step); }
void padvance2(ost::img::Progress& p,PyObject* x) { p.AdvanceProgress(x,1); }



void export_Progress()
{
  using namespace ost::img;


  class_<Progress, boost::noncopyable>("ProgressInst", no_init)
    .def("Register",pregister1)
    .def("Register",pregister2)
    .def("DeRegister",pderegister)
    .def("AdvanceProgress",padvance1)
    .def("AdvanceProgress",padvance2)
    ;

  // global functions
  def("Progress",Progress::Instance,
      return_value_policy<reference_existing_object>() );

}
