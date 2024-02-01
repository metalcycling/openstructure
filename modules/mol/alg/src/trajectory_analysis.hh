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
 * Niklaus Johner
 */
#ifndef OST_TRAJECTORY_ANALYSIS_HH
#define OST_TRAJECTORY_ANALYSIS_HH

#include <ost/mol/alg/module_config.hh>

#include <ost/base.hh>
#include <ost/geom/geom.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/coord_group.hh>


namespace ost { namespace mol { namespace alg {
  
  geom::Vec3List DLLEXPORT_OST_MOL_ALG AnalyzeAtomPos(const CoordGroupHandle& traj, const AtomHandle& a1,unsigned int stride=1);
  geom::Vec3List DLLEXPORT_OST_MOL_ALG AnalyzeCenterOfMassPos(const CoordGroupHandle& traj, const EntityView& sele,unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeDistanceBetwAtoms(const CoordGroupHandle& traj, const AtomHandle& a1, const AtomHandle& a2,unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeAngle(const CoordGroupHandle& traj, const AtomHandle& a1, const AtomHandle& a2, const AtomHandle& a3,unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeDistanceBetwCenterOfMass(const CoordGroupHandle& traj, const EntityView& sele1, const EntityView& sele2,unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeDihedralAngle(const CoordGroupHandle& traj, const AtomHandle& a1, const AtomHandle& a2, const AtomHandle& a3, const AtomHandle& a4,unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeRMSD(const CoordGroupHandle& traj, const EntityView& reference_view, const EntityView& sele,unsigned int stride=1);
  Real DLLEXPORT_OST_MOL_ALG AnalyzeRMSF(const CoordGroupHandle& traj, const EntityView& selection,int from=0, int to=-1, unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeMinDistance(const CoordGroupHandle& traj, const EntityView& view1, const EntityView& view2,unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeMinDistanceBetwCenterOfMassAndView(const CoordGroupHandle& traj, const EntityView& view_cm, const EntityView& view_atoms,unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeAromaticRingInteraction(const CoordGroupHandle& traj, const EntityView& view_ring1, const EntityView& view_ring2,unsigned int stride=1);
  void  DLLEXPORT_OST_MOL_ALG AnalyzeAlphaHelixAxis(const CoordGroupHandle& traj, const EntityView& prot_seg, geom::Vec3List& directions, geom::Vec3List& centers, unsigned int stride=1);
  void  DLLEXPORT_OST_MOL_ALG AnalyzeBestFitLine(const CoordGroupHandle& traj, const EntityView& prot_seg, geom::Vec3List& directions, geom::Vec3List& centers, unsigned int stride=1);
  void  DLLEXPORT_OST_MOL_ALG AnalyzeBestFitPlane(const CoordGroupHandle& traj, const EntityView& prot_seg, geom::Vec3List& normals, geom::Vec3List& origins, unsigned int stride=1);
  EntityHandle  DLLEXPORT_OST_MOL_ALG CreateMeanStructure(const CoordGroupHandle& traj, const EntityView& selection, int from=0, int to=-1, unsigned int stride=1);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG AnalyzeHelicity(const CoordGroupHandle& traj, const EntityView& prot_seg, unsigned int stride=1);
}}}//ns
#endif
