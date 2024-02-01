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
#include <vector>
#include <ost/geom/export_helper/vector.hh>
#include <ost/base.hh>
#include <ost/string_ref.hh>
#include <ost/platform.hh>
#include <ost/message.hh>
#include <ost/version.hh>

using namespace boost::python;

void export_Logger();
void export_GenericProp();
void export_Range();
void export_Units();

// exception translator
void translator(const ost::Error& x) {
  PyErr_SetString(PyExc_Exception, x.what() );
}


struct stringref_to_python_string
{
  static PyObject* convert(ost::StringRef const& s)
  {
    return boost::python::incref(boost::python::object(s.str()).ptr());
  }
};

struct stringref_from_python_string
{
  stringref_from_python_string()
  {
    boost::python::converter::registry::push_back(&convertible,
                                                  &construct,
                                      boost::python::type_id<ost::StringRef>());
  }
  
  static void* convertible(PyObject* obj_ptr)
  {
    if (!PyUnicode_Check(obj_ptr)) return 0;
    return obj_ptr;
  }
  
  static void construct(PyObject* obj_ptr,
                 boost::python::converter::rvalue_from_python_stage1_data* data)
  {
    PyObject* temp_bytes = PyUnicode_AsEncodedString(obj_ptr, "UTF-8", "strict");
    char* value = NULL;
    if (temp_bytes != NULL) {
      value = PyBytes_AS_STRING(temp_bytes);
      value = strdup(value);
      Py_DECREF(temp_bytes);
    }
    if (value == NULL) boost::python::throw_error_already_set();
    void* storage = (
         (boost::python::converter::rvalue_from_python_storage<ost::StringRef>*)
                     data)->storage.bytes;
    new (storage) ost::StringRef(value, strlen(value));
    data->convertible = storage;
  }
};


BOOST_PYTHON_MODULE(_ost_base)
{
  boost::python::to_python_converter<ost::StringRef,
                                     stringref_to_python_string>();
  stringref_from_python_string();
  register_exception_translator<ost::Error>(&translator);

  def("SetPrefixPath", &ost::SetPrefixPath);
  def("GetPrefixPath", &ost::GetPrefixPath);
  def("GetSharedDataPath", &ost::GetSharedDataPath);  
  
  scope().attr("VERSION")=OST_VERSION_STRING;
  scope().attr("VERSION_MAJOR")=OST_VERSION_MAJOR;
  scope().attr("VERSION_MINOR")=OST_VERSION_MINOR;  
  scope().attr("VERSION_PATCH")=OST_VERSION_PATCH;
  scope().attr("WITH_NUMPY")= OST_NUMPY_SUPPORT_ENABLED;
  export_Logger();
  export_Range();
  export_Units();
  typedef std::vector<Real> FloatList;
  class_<FloatList>("FloatList", init<>())
    .def(vector_indexing_suite<FloatList>())
    .def(geom::VectorAdditions<FloatList>())
  ;
  
  class_<std::vector<String> >("StringList", init<>())
    .def(vector_indexing_suite<std::vector<String>, true>())
  ;

  class_<std::vector<char> >("CharList", init<>())
    .def(vector_indexing_suite<std::vector<char>, true>())
  ;
  
  typedef std::vector<int> IntList;
  class_<std::vector<int> >("IntList", init<>())
    .def(vector_indexing_suite<std::vector<int> >())
    .def(geom::VectorAdditions<IntList>())
  ;

  typedef std::vector< FloatList > FloatMatrix;
  class_<FloatMatrix>("FloatMatrix", init<>())
  .def(vector_indexing_suite<FloatMatrix>())
  ;

  typedef std::vector< FloatMatrix > FloatMatrix3;
  class_<FloatMatrix3>("FloatMatrix3", init<>())
  .def(vector_indexing_suite<FloatMatrix3>())
  ;

  typedef std::vector< FloatMatrix3 > FloatMatrix4;
  class_<FloatMatrix4>("FloatMatrix4", init<>())
  .def(vector_indexing_suite<FloatMatrix4>())
  ;

  ;
}
