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
#include <boost/shared_ptr.hpp>
#include <ost/mol/mm/observer.hh>

using namespace boost::python;


void export_Observers()
{

  class_<ost::mol::mm::Observer, boost::noncopyable>("Observer", no_init);

  class_<ost::mol::mm::TrajObserver, bases<ost::mol::mm::Observer> >("TrajObserver", init<int>())
    .def("GetTraj", &ost::mol::mm::TrajObserver::GetTraj)
  ;

  //noncopyable, because TrajWriter contains a unique ofstream
  class_<ost::mol::mm::TrajWriter, boost::noncopyable, bases<ost::mol::mm::Observer> >("TrajWriter", init<int, String, String>())
    .def("Finalize", &ost::mol::mm::TrajWriter::Finalize)
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::ObserverPtr>();
  boost::python::register_ptr_to_python<ost::mol::mm::TrajObserverPtr>();
  boost::python::register_ptr_to_python<ost::mol::mm::TrajWriterPtr>();
}