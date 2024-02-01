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
  Authors: Ansgar Philippsen
*/

#ifndef IMG_GUI_OVERLAY_MANAGER_OBSERVER_H
#define IMG_GUI_OVERLAY_MANAGER_OBSERVER_H

#include "overlay_manager_fw.hh"

namespace ost { namespace img { namespace gui {

//! Abstract observer base class for the OverlayManager
class OverlayManagerObserver {
public:
  virtual ~OverlayManagerObserver() {}

  //! called on overlay lock change
  virtual void OnLockChange(OverlayManager* m, int id, bool lock) {}
  //! called on overlay visibility change
  virtual void OnVisibilityChange(OverlayManager* m, int id, bool is_visible) {}
  //! called on overlay activation
  virtual void OnActivate(OverlayManager* m, int id) {}
  //! called on addition of an overlay
  virtual void OnAddOverlay(OverlayManager* m, int id) {}
  //! called on removal of an overlay
  virtual void OnRemoveOverlay(OverlayManager* m, int id) {}

  //! called if observed manager will go out of scope
  /*!
    do not attempt to call OverlayManager::Detach in this routine!
  */
  virtual void OnRelease(OverlayManager* m) {}
};

}}}  //ns

#endif
