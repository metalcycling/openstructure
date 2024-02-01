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
  Author: Andreas Schenk
*/
#ifndef PROGRESS_HH_
#define PROGRESS_HH_

#include <vector>

#include <ost/img/module_config.hh>

namespace ost {

namespace img
{

class DLLEXPORT_OST_IMG_BASE Progress
{
public:
  static Progress& Instance();
  virtual ~Progress();
  void Register(const void* process, unsigned long steps, unsigned long delta=1);
  void DeRegister(const void* process);
  void AdvanceProgress(const void* process,unsigned long step=1);
  Real GetProgress();
  void AddCallback(void(*fpt)(void)){function_pointer_=fpt;}
  bool IsActive(){return ! processes_.empty();}
private:
  Progress();
  Real progress_;
  std::vector<const void *> processes_;
  std::vector<unsigned long> process_steps_;
  std::vector<unsigned long> process_deltas_;
  unsigned long delta_step_;
  void(*function_pointer_)(void);
};

}}//ns

#endif /*PROGRESS_HH_*/
