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
#include <boost/python/exception_translator.hpp>
using namespace boost::python;

#include <ost/mol/query.hh>
#include <ost/mol/mol.hh>
using namespace ost;
using namespace ost::mol;

void translate_QueryError(const QueryError& e) {
  PyErr_SetString(PyExc_RuntimeError, e.GetFormattedMessage().c_str());
}

String error_msg(Query& q)
{
  return QueryError(q.GetQueryString(), 
                    q.GetErrorDescription()).GetFormattedMessage();
}

void export_Query()
{
  
  register_exception_translator<QueryError>(translate_QueryError);
  
  class_<QueryErrorDesc>("QueryErrorDesc", no_init)
    .def_readonly("msg", &QueryErrorDesc::msg)
    .def_readonly("range", &QueryErrorDesc::range)
  ;
  
  enum_<QueryFlag::Flag>("QueryFlag")
    .value("NO_BONDS", QueryFlag::NO_BONDS)
    .value("EXCLUSIVE_BONDS", QueryFlag::EXCLUSIVE_BONDS)
    .value("MATCH_RESIDUES", QueryFlag::MATCH_RESIDUES)
    .export_values()
  ;
  
  class_<Query>("Query", init<String>())
    .def("GetQueryString",&Query::GetQueryString,
         return_value_policy<copy_const_reference>())
    .def("IsValid", &Query::IsValid)
    .def("IsAtomSelected", &Query::IsAtomSelected)    
    .def("IsResidueSelected", &Query::IsResidueSelected)
    .def("IsChainSelected", &Query::IsChainSelected)    
    .add_property("valid", &Query::IsValid)
    .add_property("error", &error_msg)
    .add_property("string",  
                  make_function(&Query::GetQueryString,
                                return_value_policy<copy_const_reference>()))
    .def("GetErrorDescription", &Query::GetErrorDescription, 
         return_value_policy<copy_const_reference>())
  ;
  
  def("QueryQuoteName", &QueryQuoteName);
}
