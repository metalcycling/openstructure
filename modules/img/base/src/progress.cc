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

#include <iostream>

#include <algorithm>       // STL algorithms class library
#include "progress.hh"

namespace ost { namespace img {

Progress& Progress::Instance() 
{
  static Progress instance;
  return instance;
}

Progress::Progress():
  progress_(0.0),
  processes_(),
  process_steps_(),
  process_deltas_(),
  delta_step_(0),
  function_pointer_()
{
}

Progress::~Progress()
{
}

void Progress::Register(const void* process, unsigned long steps, unsigned long delta)
{
  processes_.push_back(process);
  process_steps_.push_back(steps);
  process_deltas_.push_back(std::max<int>(1,delta));
  delta_step_=0;
  if(function_pointer_){
    function_pointer_();
  }
}

void Progress::DeRegister(const void* process)
{
  std::vector<const void *>::iterator pos=std::find(processes_.begin(),processes_.end(),process);
  if (pos !=processes_.end()){
    processes_.erase(pos,processes_.end());    
    process_steps_.erase(process_steps_.begin()+processes_.size(),process_steps_.end());
    process_deltas_.erase(process_deltas_.begin()+processes_.size(),process_deltas_.end());
  }
  if(processes_.size()==0){
    progress_=0.0;
  }
  if(function_pointer_){
    function_pointer_();
  }
  delta_step_=0;
}

void Progress::AdvanceProgress(const void* process,unsigned long step)
{
  //check if process forgot to deregister
  if(process!=processes_.back()){
    std::vector<const void *>::iterator pos=std::find(processes_.begin(),processes_.end(),process);
    if (pos+1 !=processes_.end()){
      processes_.erase(pos+1,processes_.end());    
      process_steps_.erase(process_steps_.begin()+processes_.size(),process_steps_.end());
    }
  }

  //calc progress
  Real prefac=1.0;

  for(std::vector<unsigned long>::iterator it=process_steps_.begin();it!=process_steps_.end()-1;++it){
    prefac*=1.0/static_cast<Real>(*it);
  }
  progress_+=prefac*static_cast<Real>(step)/static_cast<Real>(process_steps_.back());
  delta_step_+=1;
  if(delta_step_%process_deltas_.back()==0) {
    delta_step_=0;
    if(function_pointer_){
      function_pointer_();
    }
  }
}

Real Progress::GetProgress()
{
  return progress_;
}

}}//ns
