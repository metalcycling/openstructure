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
  Author: Ansgar Philippsen
*/

#include "data_algorithm.hh"

// algorithm base

namespace ost { namespace img {

AlgorithmBase::AlgorithmBase(const String& name):
  name_(name)
{}

AlgorithmBase::AlgorithmBase(const AlgorithmBase& b):
  name_(b.name_)
{}

AlgorithmBase::~AlgorithmBase()
{
}

AlgorithmBase& AlgorithmBase::operator=(const AlgorithmBase& b)
{
  name_=b.name_;
  return *this;
}

const String& AlgorithmBase::GetName() const 
{
  return name_;
}


// non mod algorithm

NonModAlgorithm::NonModAlgorithm(const String& name):
  AlgorithmBase(name) 
{}

NonModAlgorithm::NonModAlgorithm(const NonModAlgorithm& a):
  AlgorithmBase(a) 
{}

NonModAlgorithm& NonModAlgorithm::operator=(const NonModAlgorithm& a)
{
  AlgorithmBase::operator=(a);
  return *this;
}


// Mod-IP algorithm

ModIPAlgorithm::ModIPAlgorithm(const String& name):
  AlgorithmBase(name) 
{}

ModIPAlgorithm::ModIPAlgorithm(const ModIPAlgorithm& a):
  AlgorithmBase(a) 
{}

ModIPAlgorithm& ModIPAlgorithm::operator=(const ModIPAlgorithm& a)
{
  AlgorithmBase::operator=(a);
  return *this;
}

ConstModIPAlgorithm::ConstModIPAlgorithm(const String& name):
  AlgorithmBase(name) 
{}

ConstModIPAlgorithm::ConstModIPAlgorithm(const ConstModIPAlgorithm& a):
  AlgorithmBase(a) 
{}

ConstModIPAlgorithm& ConstModIPAlgorithm::operator=(const ConstModIPAlgorithm& a)
{
  AlgorithmBase::operator=(a);
  return *this;
}

// Mod-OP algorithm

ModOPAlgorithm::ModOPAlgorithm(const String& name):
  AlgorithmBase(name) 
{}

ModOPAlgorithm::ModOPAlgorithm(const ModOPAlgorithm& a):
  AlgorithmBase(a) 
{}

ModOPAlgorithm& ModOPAlgorithm::operator=(const ModOPAlgorithm& a)
{
  AlgorithmBase::operator=(a);
  return *this;
}


ConstModOPAlgorithm::ConstModOPAlgorithm(const String& name):
  AlgorithmBase(name) 
{}

ConstModOPAlgorithm::ConstModOPAlgorithm(const ConstModOPAlgorithm& a):
  AlgorithmBase(a) 
{}

ConstModOPAlgorithm& ConstModOPAlgorithm::operator=(const ConstModOPAlgorithm& a)
{
  AlgorithmBase::operator=(a);
  return *this;
}


}} // namespace
