//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2023 by the OpenStructure authors
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
#ifndef OST_MOL_ALG_BIOUNIT_HH
#define OST_MOL_ALG_BIOUNIT_HH

#include <ost/mol/entity_handle.hh>
#include <ost/io/mmcif_info.hh>

namespace ost { namespace mol { namespace alg {

struct BUInfo {

  BUInfo() { };

  BUInfo(const ost::io::MMCifInfoBioUnit& bu);

  void ToStream(std::ostream& stream) const;

  static BUInfo FromStream(std::istream& stream);

  String ToString() const;

  static BUInfo FromString(const String& s);

  const std::vector<std::vector<String> >& GetAUChains() const;

  const std::vector<std::vector<geom::Mat4> >& GetTransformations() const;

  void _InitTransforms() const;

  std::vector<String> au_chains;
  std::vector<int> chain_intvl;
  std::vector<std::vector<geom::Mat4> > operations;
  std::vector<int> op_intvl;

private:
  mutable std::vector<std::vector<String> > au_chains_;
  mutable std::vector<std::vector<geom::Mat4> > transforms_;
};

ost::mol::EntityHandle CreateBU(const ost::mol::EntityHandle& asu,
                                const ost::io::MMCifInfoBioUnit& bu);


ost::mol::EntityHandle CreateBU(const ost::mol::EntityHandle& asu,
                                const BUInfo& bu_info);

}}} // ns

#endif // OST_MOL_ALG_BIOUNIT_HH
