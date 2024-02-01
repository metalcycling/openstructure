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
#ifndef SHELL_HISTORY_HH
#define SHELL_HISTORY_HH
/*
  Author: Andreas Schenk
 */

#include "python_shell_fw.hh"
#include <vector>
#include <QString>

namespace ost { namespace gui {


class ShellHistory
{
public:
 ShellHistory();
 void SetCurrentCommand(const QString& command, BlockEditMode block_edit_mode);
 void AppendCommand(const QString& command, BlockEditMode block_edit_mode);
 QString GetCommand();
 BlockEditMode GetCommandMode();
 bool AtEnd();
 void MoveToPreviousMatch();
 void MoveToNextMatch();
protected:
 std::vector<std::pair<QString,BlockEditMode> > history_list_;
 QString current_;
 BlockEditMode current_mode_;
 unsigned int index_;
};

}}//ns

#endif
