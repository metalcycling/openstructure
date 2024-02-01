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

/*
 * Author Juergen Haas
 */
#ifndef OST_SUPERPOSITION_HH
#define OST_SUPERPOSITION_HH

#include <ost/base.hh>
#include <ost/geom/geom.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/alg/module_config.hh>
#include <ost/mol/atom_view.hh>

namespace ost { namespace mol { 
  
class EntityHandle;

namespace alg {

/// \brief stores the number of cycles and the two final EntityViews
/// (in case IterativeSuperposition was applied),
/// the root mean square deviation of the two superposed objects and the 
/// transformation matrix (geom::Mat4) to superpose input1 onto input2.
struct DLLEXPORT_OST_MOL_ALG SuperpositionResult {
  SuperpositionResult()
    : ncycles(0), rmsd(0.0), transformation() {
  }
  int ncycles;
  Real rmsd;
  Real rmsd_superposed_atoms;
  Real fraction_superposed;
  geom::Mat4 transformation;
  mol::EntityView entity_view1;
  mol::EntityView entity_view2;
};


// the mean square minimizer is split into a private and a public class to avoid 
// exposing the lib Eigen classes to the outside world
class MeanSquareMinimizerImpl;


class DLLEXPORT_OST_MOL_ALG MeanSquareMinimizer {
public:
  static MeanSquareMinimizer FromAtomLists(const mol::AtomViewList& atoms,
                                           const mol::AtomViewList& atoms_ref);
  static MeanSquareMinimizer FromPointLists(const std::vector<geom::Vec3>& points,
                                            const std::vector<geom::Vec3>& points_ref);
  SuperpositionResult MinimizeOnce() const;

  SuperpositionResult IterativeMinimize(int ncycles=5, Real distance_threshold=2.0) const;

  ~MeanSquareMinimizer();
  

  MeanSquareMinimizer& operator=(const MeanSquareMinimizer& rhs);
  MeanSquareMinimizer(const MeanSquareMinimizer& rhs);
  void swap(MeanSquareMinimizer& rhs) {
    std::swap(rhs.impl_, impl_);
  }

protected:
  MeanSquareMinimizer(): impl_(NULL) {}
  MeanSquareMinimizerImpl* impl_;
};

/// \brief takes the corresponding atoms and superposes them
SuperpositionResult DLLEXPORT_OST_MOL_ALG SuperposeAtoms(const mol::AtomViewList& atoms1,
                                                    const mol::AtomViewList& atoms2,
                                                    bool apply_transform);
/// \brief superposes two entity views
SuperpositionResult DLLEXPORT_OST_MOL_ALG SuperposeSVD(const mol::EntityView& ev1,
                                                    const mol::EntityView& ev2,
                                                    bool apply_transform);

/// \brief superposes two pointlists
SuperpositionResult DLLEXPORT_OST_MOL_ALG SuperposeSVD(const geom::Vec3List& pl1,
                                                       const geom::Vec3List& pl2);

/// \brief iterative superposition
SuperpositionResult DLLEXPORT_OST_MOL_ALG IterativeSuperposeSVD(const mol::EntityView& ev1,
                                                                const mol::EntityView& ev2,
                                                                int max_cycles,
                                                                Real distance_threshold,
                                                                bool apply_transform);

/// \brief iterative superposition of two point lists
SuperpositionResult DLLEXPORT_OST_MOL_ALG IterativeSuperposeSVD(const geom::Vec3List& pl1,
                                                                const geom::Vec3List& pl2,
                                                                int max_cycles,
                                                                Real distance_threshold);



/// \brief calculates RMSD for two entity view 
Real DLLEXPORT_OST_MOL_ALG CalculateRMSD(const mol::EntityView& ev1,
                                         const mol::EntityView& ev2,
                                         const geom::Mat4& transformation);



}}}//ns
#endif
