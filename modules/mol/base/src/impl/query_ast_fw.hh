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
#ifndef OST_QUERY_AST_FW_HH
#define OST_QUERY_AST_FW_HH

namespace ost { namespace mol { namespace impl {

      class WithinParam;
      class StringOrRegexParam;
      class Node;
      class LogicOpNode;
      class SelNode;

      typedef enum {
	LOP_OR, LOP_AND
      } LogicOP;

      typedef enum {
	COP_EQ, COP_NEQ, COP_GE, COP_LE, COP_LT, COP_GT
      } CompOP;
      
      
}}}

#endif
