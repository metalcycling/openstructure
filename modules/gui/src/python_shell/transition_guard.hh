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

#ifndef TRANSITION_GUARD_HH
#define TRANSITION_GUARD_HH

#include "python_shell_fw.hh"
#include "python_interpreter.hh"
#include <QObject>

namespace ost { namespace gui {


//fw decl
  class ShellHistory;
  class PythonShellWidget;


class TransitionGuard: public QObject
{
Q_OBJECT

public:
  TransitionGuard();
  virtual bool check();
};



class HistoryGuard : public TransitionGuard
{
public:
  HistoryGuard(ShellHistory* history, BlockEditMode mode);
  virtual bool check();
protected:
   ShellHistory* history_;
   BlockEditMode mode_;
};

class EditPositionGuard : public TransitionGuard
{
public:
  enum FLAGS{
    ALWAYSTRUE=0,
    SMALLER=1,
    EQUAL=2,
    BIGGER=4,
    ANCHORSMALLER=8,
    ANCHOREQUAL=16,
    ANCHORBIGGER=32
  };
    EditPositionGuard(PythonShellWidget* shell,int flags1, int flags2=ALWAYSTRUE);
    virtual bool check();
protected:
    virtual bool check_flag(int flag);
    PythonShellWidget* shell_;
    int flags1_;
    int flags2_;
};

class BlockStatusGuard : public TransitionGuard
{
public:
    BlockStatusGuard(PythonShellWidget* shell, int status);
    virtual bool check();
protected:
    PythonShellWidget* shell_;
    int status_;
};



}} //ns


#endif // TRANSITION_GUARD_HH
