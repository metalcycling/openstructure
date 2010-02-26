//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
  messages and logs

  Author: Ansgar Philippsen
*/

#ifndef OST_MESSAGE_HH
#define OST_MESSAGE_HH

#include <exception>
#include <ost/module_config.hh>

namespace ost {

struct DLLEXPORT_OST_BASE Message: public std::exception {
  Message(const String& mesg,const String& prefix=""):
    _prefix(prefix), _mesg(mesg) {}
  virtual ~Message() throw() {}
  // exception interface
  virtual const char* what() const throw() {
    String msg = _prefix + ": " +_mesg;
    return msg.c_str();
  }

  String _prefix;
  String _mesg;
};

struct DLLEXPORT_OST_BASE Error: public Message {
  Error(const String& m): Message(m,"Error") {}
};

} // namespace




#endif // OST_MESSAGE_HH
