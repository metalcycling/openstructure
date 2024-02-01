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
#ifndef OST_MOL_ALG_CLASH_SCORE_HH
#define OST_MOL_ALG_CLASH_SCORE_HH

/*
  Author: Marco Biasini
 */

#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/alg/module_config.hh>

namespace ost { namespace mol { namespace alg {

/// \defgroup Clash Steric Clash Score Calculation
/// 
/// The used scoring function is a simple steric energy function with linear 
/// fallof that is capped for small distances to avoid explosion of energies.
/// This is identical to the clash function used in SCWRL 3

//@{ 
/// \brief calculate clash score between two entity views
/// 
/// For each atom of ent_a the interaction with atoms of ent_b calculated.
/// \return 0.0 if there are no clashes and a positive clash score otherwise.
/// \sa \ref the_hammer.py "The Hammer Example"
Real DLLEXPORT_OST_MOL_ALG ClashScore(const mol::EntityView& ent_a, const mol::EntityView& ent_b);

/// \brief calculate clash score between full entity and view
Real DLLEXPORT_OST_MOL_ALG ClashScore(const mol::EntityHandle& ent_a, 
                                 const mol::EntityView& ent_b);
//// \brief calculate clash score of one single atom
/// 
/// \return floating point between 0 and 10
/// \sa \ref the_hammer.py "The Hammer Example"
Real DLLEXPORT_OST_MOL_ALG ClashScore(const mol::AtomHandle& atom, 
                                 const mol::EntityView& ent_b);

/// \brief calculate steric energy of two atoms
Real DLLEXPORT_OST_MOL_ALG StericEnergy(const geom::Vec3& pos1, Real r1,
                                   const geom::Vec3& pos2, Real r2);
//@}

/// \example the_hammer.py
/// 
/// Dynamic recalculation of clash score for a moving object. The real-valued
/// clash score is then color-mapped onto the objects. 
}}}

#endif
