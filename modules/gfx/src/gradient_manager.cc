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
  Author: Stefan Scheuber
*/
#include "boost/filesystem/operations.hpp"

#include <ost/platform.hh>

#include <ost/info/info.hh>

#include "gradient_manager.hh"

namespace ost { namespace gfx {

String STORED_GRADIENT_PATH = "/scene/gradients.xml";

GradientManager& GradientManager::Instance()
{
  static GradientManager inst;
  return inst;
}

GradientManager::GradientManager()
{
  String path = GetSharedDataPath()+STORED_GRADIENT_PATH;
  if(boost::filesystem::exists(path)){
    info::InfoHandle info_handler = info::LoadInfo(path);
    if(info_handler.HasGroup("Gradients")){
      info::InfoGroup gradients = info_handler.GetGroup("Gradients");
      info::InfoGroupList list = gradients.GetGroups("Gradient");
      for(unsigned int i = 0; i< list.size();i++){
        info::InfoGroup info_gradient = list[i];
        String name = info_gradient.GetAttribute("Name");
        Gradient gradient = Gradient::GradientFromInfo(info_gradient);
        this->SetGradient(name,gradient);
      }
    }
  }
}

GradientManager::~GradientManager()
{

}

const std::vector<String> GradientManager::GetGradientNames() const
{
  std::vector<String> keys;
  for(std::map<String , Gradient>::const_iterator it = gradients_.begin(); it != gradients_.end(); ++it) {
    keys.push_back(it->first);
  }
  return keys;
}

const Gradient& GradientManager::GetGradient(String name) const
{
  std::map<String,Gradient>::const_iterator it=gradients_.find(name);
  if(gradients_.end() != it){
    return gradients_[name];
  }
  return empty_gradient_;
}

void GradientManager::SetGradient(String name, Gradient gradient)
{
  gradients_[name] = gradient;
}

void GradientManager::RemoveGradient(String name)
{
  std::map<String,Gradient>::const_iterator it=gradients_.find(name);
  if(gradients_.end() != it){
    gradients_.erase(name);
  }
}

}} //ns
