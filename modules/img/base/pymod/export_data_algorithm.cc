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
  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>
using namespace boost::python;

#include <ost/img/algorithm.hh>
#include <ost/message.hh>

namespace ost { namespace img {

// wrapper for const algorithm
class PyNonModAlgorithm: public NonModAlgorithm {
public:
  PyNonModAlgorithm(const PyNonModAlgorithm& a):
    NonModAlgorithm(a) {}

  PyNonModAlgorithm(const String& s):
    NonModAlgorithm(s) {}

  virtual void Visit(const ConstImageHandle& ih) {VisitImage(ih);}
  virtual void VisitImage(const ConstImageHandle&) {
    throw(Error("NonModAlgorithm Base: VisitImage not overridden"));
  }
};

class WrapPyNonModAlgorithm: public PyNonModAlgorithm {
public:
  WrapPyNonModAlgorithm(PyObject* self, const PyNonModAlgorithm& a):
    PyNonModAlgorithm(a),
    self_(self) {}

  virtual void VisitImage(const ConstImageHandle& ih) {
    call_method<void>(self_, "VisitImage", ih); 
  }

private:
  PyObject* const self_;
};

// wrapper for mod ip algorithm
class PyModIPAlgorithm: public ModIPAlgorithm {
public:
  PyModIPAlgorithm(const PyModIPAlgorithm& a):
    ModIPAlgorithm(a) {}

  PyModIPAlgorithm(const String& s):
    ModIPAlgorithm(s) {}

  virtual void Visit(ImageHandle& ih) {VisitImage(ih);}
  virtual void VisitImage(ImageHandle&) {
    throw(Error("ModIPAlgorithm Base: VisitImage not overridden"));
  }
};

class WrapPyModIPAlgorithm: public PyModIPAlgorithm {
public:
  WrapPyModIPAlgorithm(PyObject* self, const PyModIPAlgorithm& a):
    PyModIPAlgorithm(a),
    self_(self) {}

  virtual void VisitImage(ImageHandle& ih) {
    call_method<void>(self_, "VisitImage", ih); 
  }

private:
  PyObject* const self_;
};

// wrapper for mod ip const algorithm
class PyConstModIPAlgorithm: public ConstModIPAlgorithm {
public:
  PyConstModIPAlgorithm(const PyConstModIPAlgorithm& a):
    ConstModIPAlgorithm(a) {}

  PyConstModIPAlgorithm(const String& s):
    ConstModIPAlgorithm(s) {}

  virtual void Visit(ImageHandle& ih) const {VisitImage(ih);}
  virtual void VisitImage(ImageHandle&) const {
    throw(Error("ConstModIPAlgorithm Base: VisitImage not overridden"));
  }
};

class WrapPyConstModIPAlgorithm: public PyConstModIPAlgorithm {
public:
  WrapPyConstModIPAlgorithm(PyObject* self, const PyConstModIPAlgorithm& a):
    PyConstModIPAlgorithm(a),
    self_(self) {}

  virtual void VisitImage(ImageHandle& ih) const {
    call_method<void>(self_, "VisitImage", ih); 
  }

private:
  PyObject* const self_;
};


// wrapper for mod op algorithm
class PyModOPAlgorithm: public ModOPAlgorithm {
public:
  PyModOPAlgorithm(const PyModOPAlgorithm& a):
    ModOPAlgorithm(a) {}

  PyModOPAlgorithm(const String& s):
    ModOPAlgorithm(s) {}

  virtual ImageHandle Visit(const ConstImageHandle& ih) {return VisitImage(ih);}
  virtual ImageHandle VisitImage(const ConstImageHandle& ih) {
    throw(Error("ModOPAlgorithm Base: VisitImage not overridden"));
  }
};

class WrapPyModOPAlgorithm: public PyModOPAlgorithm {
public:
  WrapPyModOPAlgorithm(PyObject* self, const PyModOPAlgorithm& a):
    PyModOPAlgorithm(a),
    self_(self) {}

  virtual ImageHandle VisitImage(const ConstImageHandle& ih) {
    object o = call_method<object>(self_, "VisitImage", ih); 
    return extract<ImageHandle>(o);
  }

private:
  PyObject* const self_;
};

class PyConstModOPAlgorithm: public ConstModOPAlgorithm {
public:
  PyConstModOPAlgorithm(const PyConstModOPAlgorithm& a):
    ConstModOPAlgorithm(a) {}

  PyConstModOPAlgorithm(const String& s):
    ConstModOPAlgorithm(s) {}

  virtual ImageHandle Visit(const ConstImageHandle& ih) const {return VisitImage(ih);}
  virtual ImageHandle VisitImage(const ConstImageHandle& ih) const {
    throw(Error("ConstModOPAlgorithm Base: VisitImage not overridden"));
  }
};

class WrapPyConstModOPAlgorithm: public PyConstModOPAlgorithm {
public:
  WrapPyConstModOPAlgorithm(PyObject* self, const PyConstModOPAlgorithm& a):
    PyConstModOPAlgorithm(a),
    self_(self) {}

  virtual ImageHandle VisitImage(const ConstImageHandle& ih) const {
    object o = call_method<object>(self_, "VisitImage", ih); 
    return extract<ImageHandle>(o);
  }

private:
  PyObject* const self_;
};

}} // namespace

void export_DataAlgorithm()
{

  using namespace ost::img;
  // algorithm wrapper
  class_<AlgorithmBase, boost::noncopyable>("AlgorithmBase", no_init)
    .def("GetName",&AlgorithmBase::GetName,
	 return_value_policy<return_by_value>())
    ;

  class_<NonModAlgorithm, boost::noncopyable>("NonModAlgorithmBase", no_init);

  class_<PyNonModAlgorithm, WrapPyNonModAlgorithm, bases<NonModAlgorithm>, boost::noncopyable>("NonModAlgorithm", init<const String&>() );


  class_<ModIPAlgorithm, boost::noncopyable>("ModIPAlgorithmBase", no_init);

  class_<PyModIPAlgorithm, WrapPyModIPAlgorithm, bases<ModIPAlgorithm>, boost::noncopyable>("ModIPAlgorithm", init<const String&>() );

  class_<ConstModIPAlgorithm, boost::noncopyable>("ConstModIPAlgorithmBase", no_init);

  class_<PyConstModIPAlgorithm, WrapPyConstModIPAlgorithm, bases<ConstModIPAlgorithm>, boost::noncopyable>("ConstModIPAlgorithm", init<const String&>() );

  class_<ModOPAlgorithm, boost::noncopyable>("ModOPAlgorithmBase", no_init);

  class_<PyModOPAlgorithm, WrapPyModOPAlgorithm, bases<ModOPAlgorithm>, boost::noncopyable>("ModOPAlgorithm", init<const String&>() );

  class_<ConstModOPAlgorithm, boost::noncopyable>("ConstModOPAlgorithmBase", no_init);

  class_<PyConstModOPAlgorithm, WrapPyConstModOPAlgorithm, bases<ConstModOPAlgorithm>, boost::noncopyable>("ConstModOPAlgorithm", init<const String&>() );
}
