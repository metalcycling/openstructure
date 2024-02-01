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

#include "admin.hh"

namespace ost { namespace gui {

#if defined(__APPLE__)

AdminRights::AdminRights(): auth_(NULL)
{
  
}
AdminRights::~AdminRights()
{
  this->Release();
}

bool AdminRights::Acquire()
{
  // Connect to Authorization Services.
  if (AuthorizationCreate(NULL, NULL, 0, &auth_)!=noErr) {
    return false;
  }
    
  static const AuthorizationFlags  kFlags =
                kAuthorizationFlagInteractionAllowed
              | kAuthorizationFlagExtendRights;
  AuthorizationItem   kActionRight = { "", 0, 0, 0 };
  AuthorizationRights kRights      = { 1, &kActionRight };

  // Request the application-specific right.
  return noErr==AuthorizationCopyRights(auth_, &kRights, NULL, kFlags, NULL);
}

void AdminRights::Release()
{
  if (auth_) {
    AuthorizationFree(auth_, kAuthorizationFlagDestroyRights);
    auth_=NULL;
  }
}

void AdminRights::CreateLink(const String& from, const String& to)
{
  static const char* minus_sf="-sf";
  // const_casts are gross, but at least that way we keep gcc happy
  char* const args[]={
    const_cast<char*>(minus_sf),
    const_cast<char*>(from.c_str()),
    const_cast<char*>(to.c_str()),
    NULL
  };
  AuthorizationExecuteWithPrivileges(auth_, "/bin/ln", 
                                     kAuthorizationFlagDefaults,
                                     args, NULL);
}

#else

AdminRights::AdminRights()
{ }

void AdminRights::CreateLink(const String& from, const String& to)
{
  
}
AdminRights::~AdminRights()
{ }

bool AdminRights::Acquire()
{ return false; }

void AdminRights::Release()
{ }
#endif


}}

