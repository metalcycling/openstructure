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
#ifndef PYTHON_SHELL_FW_HH
#define PYTHON_SHELL_FW_HH

/*
  Author: Andreas Schenk
 */

namespace ost { namespace gui {

// fw decl
class PythonShell;

enum BlockEditMode {
  EDITMODE_SINGLELINE,
  EDITMODE_MULTILINE_ACTIVE,
  EDITMODE_MULTILINE_INACTIVE
};
enum BlockType {
  BLOCKTYPE_OUTPUT=1,
  BLOCKTYPE_ERROR=2,
  BLOCKTYPE_CODE=4,
  BLOCKTYPE_ACTIVE=8,
  BLOCKTYPE_BLOCKEDIT=16,
  BLOCKTYPE_MULTILINE_SQ=32,
  BLOCKTYPE_MULTILINE_DQ=64
  };

struct GutterBlock{
  GutterBlock(int s, int e,BlockType bt):start(s),end(e),type(bt){}
  int start;
  int end;
  BlockType type;
};


}}//ns

#endif
