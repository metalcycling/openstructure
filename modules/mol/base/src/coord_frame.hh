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
#ifndef OST_MOL_COORD_FRAME_HH
#define OST_MOL_COORD_FRAME_HH

/*
  Authors: Marco Biasini, Niklaus Johner, Ansgar Philippsen
 */
#include <ost/mol/module_config.hh>
#include <boost/shared_ptr.hpp>
#include <ost/geom/geom.hh>
#include <ost/mol/module_config.hh>
#include <ost/mol/entity_view.hh>
#include "atom_handle.hh"

namespace ost { namespace mol {

/*!
  A coordinate frame in a trajectory, containing atom positions and
  per-frame unit cell size and angles.
*/

/*
  TODO: 
  - move algorithmic code to separate functions in mol/alg
  - clean up mix of views and atom indices methods/functions
  - use existing UnitCell class
*/

class DLLEXPORT_OST_MOL CoordFrame : public geom::Vec3List{
private:
  geom::Vec3 ucell_size_;
  geom::Vec3 ucell_angles_;
public:
  typedef geom::Vec3List base_type;
  
  CoordFrame() : base_type() {}
  CoordFrame(size_t size, const geom::Vec3& value=geom::Vec3()) : base_type(size, value) {}
  CoordFrame(base_type::iterator b, base_type::iterator e): base_type(b, e) { }
  CoordFrame(const base_type& rhs) : base_type(rhs) { }
  CoordFrame(const std::vector<geom::Vec3>& rhs) : base_type(rhs) { }
  CoordFrame(const std::vector<geom::Vec3>& rhs,
             const geom::Vec3 box_size, 
             const geom::Vec3 box_angles) : base_type(rhs) {
    ucell_size_=box_size;
    ucell_angles_=box_angles;
  }

  void SetCellSize(const geom::Vec3& s) {
    ucell_size_=s;
  }

  geom::Vec3 GetCellSize() const {
    return this->ucell_size_;
  }

  void SetCellAngles(const geom::Vec3& a) {
    ucell_angles_=a;
  }

  geom::Vec3 GetCellAngles() const {
    return this->ucell_angles_;
  }
  
  geom::Vec3List GetCellVectors() const;
  
  geom::Vec3 GetAtomPos(const AtomHandle& atom) const;
  geom::Vec3 GetAtomPos(int atom_index) const;
  Real GetDistanceBetwAtoms(const AtomHandle& a1, const AtomHandle& a2) const;
  Real GetDistanceBetwAtoms(int atom1_index, int atom2_index) const;

  /*!
    Returns the angle between the three atoms in rad

    a2 is considered the central atom, i.e the result is the angle between 
    the vectors (a1.pos-a2.pos) and (a3.pos-a2.pos)
  */
  Real GetAngle(const AtomHandle& a1, const AtomHandle& a2, const AtomHandle& a3) const;

  /*!
    Returns the angle between the three atom indices in rad

    atom2 is considered the central atom, i.e the result is the angle between 
    the vectors (atom1.pos-atom2.pos) and (atom3.pos-atom2.pos)
  */
  Real GetAngle(int atom1_index, int atom2_index, int atom3_index) const;  

  //! Returns the Dihedral angle between the four atoms a1,a2,a3,a4
  Real GetDihedralAngle(const AtomHandle& a1, const AtomHandle& a2, const AtomHandle& a3, const AtomHandle& a4) const;

  //! Returns the Dihedral angle between the four atom indices a1,a2,a3,a4
  Real GetDihedralAngle(int a1_index, int a2_index, int a3_index, int a4_index) const; 

  //! Returns the position of the centor of mass of the atoms in the EntityView
  geom::Vec3 GetCenterOfMassPos(const mol::EntityView& sele) const;

  //! Calculates the center of mass given a list of indices and provided masses
  /*!
    the masses vector contains only the masses for the atoms indexed by indices, not for
    all atoms in the frame
  */
  geom::Vec3 GetCenterOfMassPos(std::vector<unsigned long>& indices, 
                                std::vector<Real>& masses) const;

  //! Returns the distance between the centers of mass of the two EntityViews
  Real GetDistanceBetwCenterOfMass(const mol::EntityView& sele1, const mol::EntityView& sele2) const;

  /*! Returns the distance betweem two custom center of masses
    
    Calculates the centers of mass of the two groups of atoms from given by index lists and
    custom masses. See comment regarding masses vector for GetCenterOfMassPos
  */
  Real GetDistanceBetwCenterOfMass(std::vector<unsigned long>& indices1, std::vector<Real>& masses1, 
                                   std::vector<unsigned long>& indices2, std::vector<Real>& masses2) const;

  /*!
    Returns the RMSD between the positions of the atoms whose indices are given in indices_sele 
    and the positions given in ref_pos
  */
  Real GetRMSD(const std::vector<geom::Vec3>& ref_pos, 
               const std::vector<unsigned long>& indices_sele) const;
  
  /*!
    Returns the rmsd between two EntityViews.

    The reference positions are taken directly from the reference_view,
    whereas they are taken from this CoordFrame for the sele_view
  */
  Real GetRMSD(const mol::EntityView& reference_view,
               const mol::EntityView& sele_view) const;

  //! Returns the minimal distance between two groups of atoms 
  /*!
    with indices given by index_list1 and index_list2
  */
  Real GetMinDistance(std::vector<unsigned long>& index_list1, 
                      std::vector<unsigned long>& index_list2) const;

  //! Returns the minimal distance between the atoms of two views (view1 and view2)
  Real GetMinDistance(const mol::EntityView& view1, const mol::EntityView& view2) const;

  /*!
    Returns the minimal distance between the center of mass of two groups of atoms
    
    atoms and their masses are defined by indices and masses for those atoms
  */
  Real GetMinDistBetwCenterOfMassAndView(std::vector<unsigned long>& indices_cm,
                                         std::vector<Real>& masses_cm,
                                         std::vector<unsigned long>& indices_atoms) const;
  //! Returns the minimal distance between the center of mass of two views
  Real GetMinDistBetwCenterOfMassAndView(const mol::EntityView& view_cm,
                                         const mol::EntityView& view_atoms) const;

  //! Returns the best fit line to atoms with indices in indices_ca
  geom::Line3 GetODRLine(std::vector<unsigned long>& indices_ca) const;

  //! Returns the best fit line to atoms in the EntityView view1
  geom::Plane GetODRPlane(std::vector<unsigned long>& indices_ca) const;

  //! Returns the normal to the best fit plane to atoms with indices in indices_ca
  geom::Line3 GetODRLine(const mol::EntityView& view1) const;
  geom::Plane GetODRPlane(const mol::EntityView& view1) const;

  /*!
    Returns a line which is the axis of a fitted cylinder to the atoms with indices given in indices_ca
    It is assumed that we fit an alpha-helix and that the CA atoms are oredered sequentially
    This is used for the initial guess of the helix axis
  */
  std::pair<geom::Line3, Real> FitCylinder(std::vector<unsigned long>& indices_ca) const;

  //! see FitCylinder(std::vector<unsigned long>&)
  std::pair<geom::Line3, Real> FitCylinder(const mol::EntityView& view1) const;

  /*!
    Returns the percentage of residues in the EntityView (segment) that are in an alpha-helix
    Each residue is assigned as being in an alpha-helix or not, based on and backbone torsion angles
    and the hydrogen bonding (O-N) (this term is used only if the segment contains at least 8 residues
    The entity view should be a single segment with no gaps and no missing N,CA,C,O backbone atoms
    the first and last residues will not be asessed for helicity
  */
  Real GetAlphaHelixContent(std::vector<unsigned long>& indices_ca, 
                            std::vector<unsigned long>& indices_c,
                            std::vector<unsigned long>& indices_o,
                            std::vector<unsigned long>& indices_n) const;

  //! see above
  Real GetAlphaHelixContent(const mol::EntityView& segment) const;
};

  typedef boost::shared_ptr<CoordFrame> CoordFramePtr;
  typedef std::vector<CoordFramePtr> CoordFrameList;

  // factory method
  // create a frame from a Vec3List containing the positions of the atoms, optionally with unit cell
  // TODO: why is this necessary? the ctor works fine 
  DLLEXPORT_OST_MOL CoordFrame CreateCoordFrame(const geom::Vec3List& atom_pos, 
                                                const geom::Vec3& cell_size=geom::Vec3(),
                                                const geom::Vec3& cell_angles=geom::Vec3());

  // these should really be either in the entity view interface or in mol/alg

  /*!
    This function returns a vector containing the atom indices of the atoms in an EntityView;
    it is used to accelerate the extraction of information from a trajectory
  */
  DLLEXPORT_OST_MOL void GetIndices(const EntityView& sele, std::vector<unsigned long>& indices);

  /*!
    This function returns a vector containing the atom masses of the atoms in an EntityView;
    it is used together with GetIndices to accelerate the extraction of RMSD from a trajectory
  */
  DLLEXPORT_OST_MOL void GetMasses(const EntityView& sele, std::vector<Real>& masses);

  //! conveniece for GetIndices and GetMasses in one call
  DLLEXPORT_OST_MOL void GetIndicesAndMasses(const EntityView& sele,
                                             std::vector<unsigned long>& indices,
                                             std::vector<Real>& masses);

  //! Writes the positions of all atoms in the EntityView into the provided vec3 list
  DLLEXPORT_OST_MOL void GetPositions(const EntityView& sele, std::vector<geom::Vec3>& ref_pos);
  
  //! Writes the indices of all atoms in the EntityView into the provided list
  DLLEXPORT_OST_MOL void GetCaIndices(const EntityView& segment, std::vector<unsigned long>& indices_ca);
  
  //! Writes the backbone indices of all residues in the EntityView into the provided list
  DLLEXPORT_OST_MOL void GetCaCONIndices(const EntityView& segment, 
                                         std::vector<unsigned long>& indices_ca, 
                                         std::vector<unsigned long>& indices_c,
                                         std::vector<unsigned long>& indices_o, 
                                         std::vector<unsigned long>& indices_n);


}}

#endif
