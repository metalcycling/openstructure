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
#ifndef OST_GFX_GL_INCLUDE_HH
#define OST_GFX_GL_INCLUDE_HH

#include <ost/config.hh>
/*
  wrapper for inclusion of OpenGL headers
  for various platforms and builds

  Author: Ansgar Philippsen
*/

#ifdef _MSC_VER
#include <windows.h>
# ifdef FindAtom
#  undef FindAtom
# endif
# ifdef AddAtom
#  undef AddAtom
# endif
#endif


# if defined(__APPLE__)
#  define GL_SILENCE_DEPRECATION
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  undef  GL_SILENCE_DEPRECATION
# else
#  include <GL/gl.h>
#  include <GL/glu.h>
# endif


#endif
