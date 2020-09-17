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
#ifndef OST_GUI_SIP_HANDLER_HH
#define OST_GUI_SIP_HANDLER_HH

#include <boost/python.hpp>
#include <iostream>
#include <QWidget>

#include <ost/message.hh>

using namespace boost::python;

/*
   Author: Ansgar Philippsen, Marco Biasini
 */

namespace ost { namespace gui {


template <class O> object get_py_qobject(O* cpp_object)
{
  if (cpp_object != NULL){
    static object pyqt5_module=import("PyQt5.QtCore");
    static object sip_module=import("sip");
    QObject* qobject = qobject_cast<QObject*>(cpp_object);
    unsigned long addr = reinterpret_cast<unsigned long>(qobject);
    object py_qobject = pyqt5_module.attr("QObject");
    object object = sip_module.attr("wrapinstance")(addr, py_qobject);
    return object;
  }
  return object();
};


template <class O> O* get_cpp_qobject(object py_object)
{
  if(py_object.ptr() != Py_None){
    if(PyObject_HasAttrString(py_object.ptr(), "qobject")){
      py_object = py_object.attr("qobject");
    }
    static object pyqt5_module=import("PyQt5.QtCore");
    static object sip_module=import("sip");
    unsigned long addr = extract<unsigned long>(sip_module.attr("unwrapinstance")(py_object));
    if(addr){
      return reinterpret_cast<O*>(addr);
    }
  }
  return NULL;
};

class SipHandlerBase {
public:
  virtual ~SipHandlerBase() {}  
  virtual unsigned long GetSipHandle() const { return 0; }
};

template <typename T>
class SipHandler : public SipHandlerBase {
public:
  explicit SipHandler(T* instance):
    instance_(instance)
  { }

  virtual unsigned long GetSipHandle() const 
  {
    // cast to QWidget first, before doing reinterpret cast to avoid problems
    // with multiple inheritance.
    return reinterpret_cast<unsigned long>((QWidget*)instance_);
  }
  T* Me() { 
    if(!instance_) throw Error("SipHandler is no longer valid");
    return instance_; 
  }
  const T* Me() const { 
    if(!instance_) throw Error("SipHandler is no longer valid");
    return instance_; 
  }
  void Release(T* p){
    if(Me()==p){
      instance_=0;
    }
  }
private:
  T* instance_;
};
 
}}

#endif
