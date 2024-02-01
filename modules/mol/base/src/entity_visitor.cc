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
#include "entity_visitor.hh"

namespace ost { namespace mol {

EntityVisitor::~EntityVisitor()
{}

void EntityVisitor::OnEntry()
{}

void EntityVisitor::OnExit()
{}

bool EntityVisitor::VisitChain(const ChainHandle& c) 
{
  return true;
}

bool EntityVisitor::VisitResidue(const ResidueHandle& r) 
{
  return true;
}

bool EntityVisitor::VisitAtom(const AtomHandle& a) 
{
  return true;
}

bool EntityVisitor::VisitBond(const BondHandle& b)
{
  return true;
}

bool EntityVisitor::VisitTorsion(const TorsionHandle& t)
{
  return true;
}

///// entity view visitor

EntityViewVisitor::~EntityViewVisitor()
{}

void EntityViewVisitor::OnEntry()
{}

void EntityViewVisitor::OnExit()
{}

bool EntityViewVisitor::VisitChain(const ChainView& c) 
{
  return true;
}

bool EntityViewVisitor::VisitResidue(const ResidueView& r) 
{
  return true;
}

bool EntityViewVisitor::VisitAtom(const AtomView& a) 
{
  return true;
}

bool EntityViewVisitor::VisitBond(const BondHandle& b)
{
  return true;
}


}} // ns

