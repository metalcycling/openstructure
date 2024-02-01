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
#include <vector>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <ost/gui/tools/tool.hh>
#include <ost/gui/tools/tool_manager.hh>
#include <ost/gui/tools/tool_options_win.hh>
#include "sip_handler.hh"
#include "qptr.hh"

using namespace boost::python;
using namespace ost::gui;
using namespace ost;

  
struct WrappedTool : public Tool
{
    WrappedTool(PyObject *p, const String& name): 
      Tool(name), self(p) 
    { }

    virtual void Click(const MouseEvent& event) 
    {
      try {
        return call_method<void, MouseEvent>(self, "Click", event);
      } catch(error_already_set& ) {
        PyErr_Print();
      }
    }
    
    void MouseEventDefault(const MouseEvent&) 
    {
    }
    
    bool CanOperateOn(gfx::NodePtrList nodes)
    {
      try {
        return call_method< bool, gfx::NodePtrList >(self, "CanOperateOn", nodes);
      } catch(error_already_set& ) {
        PyErr_Print();
      }
      return false;
    }
    bool CanOperateOnDefault(gfx::NodePtrList) {return false; }

    QIcon GetIcon() const
    {
      try {
        String loc = call_method<String>(self, "GetIconPath");
        QIcon icon = QIcon(loc.c_str());
        return icon;
      } catch(error_already_set& ) {
        PyErr_Print();
      }
      QIcon icon = QIcon();
      return icon;
    }
 private:
    PyObject* self;
};
                                

  
void tm_add_tool(ToolManager& tm, QPtr<WrappedTool> tool)
{
  tm.AddTool(tool.get());
  tool.Release();
}

void to_add_option_a(ToolOptions& to, QPtr<ToolOption> opt)
{
  to.AddOption(opt.get());
  opt.Release();
}

void to_add_option_b(ToolOptions& to, const String& name,
                     QPtr<ToolOption> opt)
{
  to.AddOption(name, opt.get());
  opt.Release();
}


object qpoint_to_bp_object (QPoint& p)
{
  static object bp_module =import("PyQt5.QtCore");
  return bp_module.attr("QPoint")(p.x(),p.y());
}

object get_pos_wrapper(MouseEvent& me)
{
  QPoint pos=me.GetPos();
  return qpoint_to_bp_object(pos);
}

object get_last_pos_wrapper(MouseEvent& me)
{
  QPoint last_pos=me.GetPos();
  return qpoint_to_bp_object(last_pos);
}

object get_delta_wrapper(MouseEvent& me)
{
  QPoint delta=me.GetPos();
  return qpoint_to_bp_object(delta);
}


ToolOption* (ToolOptions::*get_option_a)(const String&, 
                                           const String&) const=&ToolOptions::GetOption;
ToolOption* (ToolOptions::*get_option_b)(const String&) const=&ToolOptions::GetOption;

void export_Tool()
{
  enum_<MouseEvent::Button>("MouseButton")
    .value("LeftButton", MouseEvent::LeftButton)
    .value("RightButton", MouseEvent::RightButton)
  ;

  class_<MouseEvent>("MouseEvent", no_init)
    .def("GetLastPos", &get_last_pos_wrapper)
    .def("GetPos", &get_pos_wrapper)
    .def("GetDelta", &get_delta_wrapper)
    .def("IsShiftPressed", &MouseEvent::IsShiftPressed)    
    .def("IsControlPressed", &MouseEvent::IsControlPressed)
    .def("IsMetaPressed", &MouseEvent::IsMetaPressed)        
    .add_property("delta", &get_delta_wrapper)
    .add_property("pos", &get_pos_wrapper)
    .add_property("last_pos", &get_last_pos_wrapper)
    .def("GetButtons", &MouseEvent::GetButtons)
    .add_property("buttons", &MouseEvent::GetButtons)
  ;
  
  enum_<ToolOption::Type>("ToolOptionType")
    .value("INT", ToolOption::INT)
    .value("FLOAT", ToolOption::FLOAT)
    .value("ENUM", ToolOption::ENUM)
  ;  

  class_<ToolOption, QPtr<ToolOption>,
         boost::noncopyable>("ToolOption", no_init)
    .def("GetVerboseName", &ToolOption::GetVerboseName,
         return_value_policy<copy_const_reference>())
    .def("GetKey", &ToolOption::GetKey,
         return_value_policy<copy_const_reference>())
    .def("GetType", &ToolOption::GetType)
    .add_property("type", &ToolOption::GetType)
    .add_property("key", 
                  make_function(&ToolOption::GetKey, 
                                return_value_policy<copy_const_reference>()))
    .add_property("verbose_name", 
                  make_function(&ToolOption::GetVerboseName, 
                                return_value_policy<copy_const_reference>()))
  ;
  //register_ptr_to_python<QPtr<ToolOption> >();
  class_<ToolOptionInt, QPtr<ToolOptionInt>, boost::noncopyable, 
         bases<ToolOption> >("ToolOptionInt", 
                              init<const String&, const String&, int>())
    .def(init<const String&, const String&, int, int, int>())
    .def("GetValue", &ToolOptionInt::GetValue)
    .def("SetValue", &ToolOptionInt::SetValue)
    .add_property("value", &ToolOptionInt::GetValue, 
                  &ToolOptionInt::SetValue)
  ;
  implicitly_convertible<QPtr<ToolOptionInt>, 
                         QPtr<ToolOption> >();
  class_<ToolOptionFloat, QPtr<ToolOptionFloat>, boost::noncopyable, 
         bases<ToolOption> >("ToolOptionFloat", 
                             init<const String&, const String&, float>())
    .def(init<const String&, const String&, float, float, float>())
    .def("GetValue", &ToolOptionFloat::GetValue)
    .def("SetValue", &ToolOptionFloat::SetValue)
    .add_property("value", &ToolOptionFloat::GetValue, 
                  &ToolOptionFloat::SetValue)
  ;  
  implicitly_convertible<QPtr<ToolOptionFloat>, 
                         QPtr<ToolOption> >();  
  class_<ToolOptionEnum, QPtr<ToolOptionEnum>, bases<ToolOption>, 
         boost::noncopyable>("ToolOptionEnum", 
                             init<const String&, const String&>())  
    .def("GetIndex", &ToolOptionEnum::GetIndex)
    .def("SetIndex", &ToolOptionEnum::SetIndex)
    .def("GetValue", &ToolOptionEnum::GetValue)    
    .def("Add", &ToolOptionEnum::Add)
    .add_property("value", &ToolOptionEnum::GetValue)
    .add_property("index", &ToolOptionEnum::GetIndex, &ToolOptionEnum::SetIndex)
  ;
  implicitly_convertible<QPtr<ToolOptionEnum>, 
                         QPtr<ToolOption> >();

  class_<ToolOptions, boost::noncopyable>("ToolOptions", no_init)
    .def("AddOption", to_add_option_a)
    .def("AddOption", to_add_option_b)
    .def("GetOption", get_option_a, 
         return_value_policy<reference_existing_object>())
    .def("GetOption", get_option_b,
         return_value_policy<reference_existing_object>())    
  ;
  
  register_ptr_to_python<ToolOptions*>();  
  class_<Tool, QPtr<WrappedTool>, WrappedTool, 
         boost::noncopyable>("Tool", init<const String&>())
    .def("Click", &WrappedTool::MouseEventDefault)
    .def("DoubleClick", &WrappedTool::MouseEventDefault)    
    .def("CanOperateOn", &WrappedTool::CanOperateOnDefault)
    .def("GetIcon", &WrappedTool::GetIcon)
    .def("GetToolOptions", &WrappedTool::GetToolOptions,
         return_value_policy<reference_existing_object>())
    .add_property("tool_options", 
                  make_function(&Tool::GetToolOptions,
                                return_value_policy<reference_existing_object>()))
    .def("GetName", &Tool::GetName, return_value_policy<copy_const_reference>())
  ;
  register_ptr_to_python<Tool*>();
  class_<ToolManager, boost::noncopyable>("ToolManager", no_init)
    .def("Instance", &ToolManager::Instance, 
         return_value_policy<reference_existing_object>()).staticmethod("Instance")
    .def("AddTool", &tm_add_tool)
    .def("FindToolByName", &ToolManager::FindToolByName,
         return_value_policy<reference_existing_object>())
  ;
  class_<ToolOptionsWin, boost::noncopyable>("ToolOptionsWin", no_init)
    .def("Show", &ToolOptionsWin::show)
    .def("Hide", &ToolOptionsWin::hide)
    .def("GetQObject",&get_py_qobject<ToolOptionsWin>)
    .add_property("qobject", &get_py_qobject<ToolOptionsWin>)
  ;
}

