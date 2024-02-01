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

/*
  Authors: Stefan Scheuber, Ansgar Philippsen
*/

#include "spnav_input.hh"

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <stdexcept>

#include <spnav.h>
#include <sys/un.h>
#include <sys/types.h>
#include <sys/socket.h>

namespace {

  // immitate spnav_open, but without the annoying perror littering 
  bool check_spnav() {
    int s=0;
    sockaddr_un addr;

    memset(&addr, 0, sizeof(sockaddr_un));
    addr.sun_family = AF_UNIX;
    strncpy(addr.sun_path, "/var/run/spnav.sock", sizeof(addr.sun_path));

    if((s = socket(PF_UNIX, SOCK_STREAM, 0)) == -1) {
      return false;
    }
    if(connect(s, (struct sockaddr*)&addr, sizeof(sockaddr_un)) == -1) {
      return false;
    }
    close(s);
    return true;
  }

}

namespace ost { namespace gui {

SpnavInput::SpnavInput(QObject* parent): 
  QThread(parent),
  valid_(false)
{
  qRegisterMetaType<geom::Mat4>("geom::Mat4");
  if(check_spnav() && spnav_open()!=-1) {
    valid_=true;
  } 
}

void SpnavInput::run()
{
  if(!valid_) return;
	spnav_event sev;

  while(spnav_wait_event(&sev)) {
    if(sev.type == SPNAV_EVENT_MOTION) {
      emit this->deviceTransformed(sev.motion.x,sev.motion.y,sev.motion.z, sev.motion.rx, sev.motion.ry, sev.motion.rz);
    } else if (sev.type == SPNAV_EVENT_BUTTON) {
      //printf("got button %s event b(%d)\n", sev.button.press ? "press" : "release", sev.button.bnum);
      if(sev.button.press)
        emit this->deviceButtonPressed(sev.button.bnum);
    }
  }
}

SpnavInput::~SpnavInput()
{
  if(valid_) {
    spnav_close();
    this->quit();
    this->terminate();
  }
}

SpnavInput* SpnavInput::Instance() {
  static SpnavInput inst;
  return &inst;
}

}} //ns
