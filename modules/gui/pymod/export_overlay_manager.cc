/////////////////////////////////////////////////////////
// IPLT - Image Processing Library & Toolkit
// Copyright (c) 2003-2007 University of Basel
//
// The usage terms and conditions of this software 
// are governed by the GNU General Public License, 
// as described in the accompanying text file LICENSE.
// Please review this information carefully and contact
// info@iplt.org should questions arise.
/////////////////////////////////////////////////////////

/*
  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>

#include <ost/gui/data_viewer/overlay_manager.hh>
#include <ost/gui/data_viewer/overlay_base.hh>


using namespace boost::python;
using namespace ost;
using namespace ost::img;
using namespace ost::img::gui;

void export_overlay_manager()
{
  void (OverlayManager::*ActivateOverlay1)(const String& name) = &OverlayManager::ActivateOverlay;
  void (OverlayManager::*ActivateOverlay2)(int id) = &OverlayManager::ActivateOverlay;
  void (OverlayManager::*SetOverlayVisibility1)(const String& name, bool visible) = &OverlayManager::SetOverlayVisibility;
  void (OverlayManager::*SetOverlayVisibility2)(int id, bool visible ) = &OverlayManager::SetOverlayVisibility;
  OverlayPtr (OverlayManager::*RetrieveOverlay1)(const String& name) = &OverlayManager::RetrieveOverlay;
  OverlayPtr (OverlayManager::*RetrieveOverlay2)(int id) = &OverlayManager::RetrieveOverlay;
  bool (OverlayManager::*IsVisible1)(const String& name) = &OverlayManager::IsVisible;
  bool (OverlayManager::*IsVisible2)(int id) = &OverlayManager::IsVisible;
  bool (OverlayManager::*IsActive1)(const String& name) = &OverlayManager::IsActive;
  bool (OverlayManager::*IsActive2)(int id) = &OverlayManager::IsActive;

  class_<OverlayManager, OverlayManagerPtr, boost::noncopyable>("OverlayManager", no_init)
    .def("AddOverlay",&OverlayManager::AddOverlay)
    .def("GetActiveOverlay",&OverlayManager::GetActiveOverlay)
    .def("ActivateOverlay",ActivateOverlay1)
    .def("ActivateOverlay",ActivateOverlay2)
    .def("SetOverlayVisibility",SetOverlayVisibility1)
    .def("SetOverlayVisibility",SetOverlayVisibility2)
    .def("RetrieveOverlay",RetrieveOverlay1)
    .def("RetrieveOverlay",RetrieveOverlay2)
    .def("IsVisible",IsVisible1)
    .def("IsVisible",IsVisible2)
    .def("IsActive",IsActive1)
    .def("IsActive",IsActive2)
    .def("GetOverlayName",&OverlayManager::GetOverlayName)
    .def("GetOverlayId",&OverlayManager::GetOverlayId)
    ;
}
