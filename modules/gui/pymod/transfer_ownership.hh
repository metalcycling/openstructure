#ifndef OST_GUI_PYMOD_TRANSFER_OWNERSHIP_HH
#define OST_GUI_PYMOD_TRANSFER_OWNERSHIP_HH

/*
  Author: Marco Biasini
 */

#include <boost/python.hpp>


namespace ost { namespace gui {

inline void TransferOwnership(QWidget* obj)
{
  boost::python::object main_module=boost::python::import("__main__");  
  boost::python::object qt_gui=boost::python::import("PyQt5.QtGui");
  boost::python::object pobj(reinterpret_cast<unsigned long>(obj));
  boost::python::object sip=main_module.attr("sip");
  boost::python::object type=qt_gui.attr("QWidget");
  boost::python::object wrapped_obj=sip.attr("wrapinstance")(pobj, type);
  sip.attr("transferto")(wrapped_obj, boost::python::object());
}

}}

#endif
