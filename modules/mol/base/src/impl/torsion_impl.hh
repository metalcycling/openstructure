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
#ifndef TORSION_IMPL_HH
#define TORSION_IMPL_HH

#include <boost/enable_shared_from_this.hpp>

#include <ost/mol/module_config.hh>
#include <ost/geom/geom.hh>

#include "torsion_impl_fw.hh"
#include "connector_impl_fw.hh"
#include "atom_impl_fw.hh"
#include <ost/mol/entity_visitor_fw.hh>
#include "dihedral.hh"
namespace ost { namespace mol { namespace impl {

/// \internal
class TorsionImpl: public Dihedral, 
                   public boost::enable_shared_from_this<TorsionImpl> {
public:
  TorsionImpl(const String& name, 
              const AtomImplPtr& atom1, 
              const AtomImplPtr& atom2, 
              const AtomImplPtr& atom3,
              const AtomImplPtr& atom4);

  void Apply(EntityVisitor& h);
  
  const String& GetName() const;
  
  Real GetAngle() const;
  
  void SetAngle(Real a, bool up);

  void RotateAngle(Real angle, bool up);
private:
  String           name_;
  ConnectorImplW   conn2_;
};

}}} // ns

#endif
