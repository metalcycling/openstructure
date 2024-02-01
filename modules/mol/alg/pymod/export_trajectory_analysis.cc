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
#include <boost/python.hpp>
using namespace boost::python;

#include <ost/mol/alg/trajectory_analysis.hh>

using namespace ost;
using namespace ost::mol::alg;

//"thin wrappers" for default parameters
/*
BOOST_PYTHON_FUNCTION_OVERLOADS(extractapos, ExtractAtomPosition, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(extractcmpos, ExtractCMPosition, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(extractdist, ExtractDistance, 3, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(extractang, ExtractAngle, 4, 5)
BOOST_PYTHON_FUNCTION_OVERLOADS(extractdih, ExtractDihedral, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(extractcmdist, ExtractCMDistance, 3, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(extractrmsd, ExtractRMSD, 3, 4)
*/
void export_TrajectoryAnalysis()
{
  def("AnalyzeAtomPos",&AnalyzeAtomPos, (arg("traj"), arg("atom"), arg("stride")=1));
  def("AnalyzeCenterOfMassPos",&AnalyzeCenterOfMassPos, (arg("traj"), arg("selection"), arg("stride")=1));
  def("AnalyzeDistanceBetwAtoms",&AnalyzeDistanceBetwAtoms, (arg("traj"), arg("atom"), arg("atom"), arg("stride")=1));
  def("AnalyzeAngle",&AnalyzeAngle, (arg("traj"), arg("atom"), arg("atom"), arg("atom"), arg("stride")=1));
  def("AnalyzeDihedralAngle",&AnalyzeDihedralAngle, (arg("traj"), arg("atom"), arg("atom"), arg("atom"), arg("atom"), arg("stride")=1));
  def("AnalyzeDistanceBetwCenterOfMass",&AnalyzeDistanceBetwCenterOfMass, (arg("traj"), arg("selection"), arg("selection"), arg("stride")=1));
  def("AnalyzeRMSD",&AnalyzeRMSD, (arg("traj"), arg("reference_view"), arg("selection"), arg("stride")=1));
  def("AnalyzeRMSF",&AnalyzeRMSF, (arg("traj"), arg("selection"), arg("first")=0, arg("last")=-1, arg("stride")=1));
  def("AnalyzeMinDistance", &AnalyzeMinDistance, (arg("traj"), arg("view1"), arg("view2"), arg("stride")=1));
  def("AnalyzeMinDistanceBetwCenterOfMassAndView", &AnalyzeMinDistanceBetwCenterOfMassAndView, (arg("traj"), arg("view_cm"), arg("view_atoms"), arg("stride")=1));
  def("AnalyzeAromaticRingInteraction", &AnalyzeAromaticRingInteraction, (arg("traj"), arg("view_ring1"), arg("view_ring2"), arg("stride")=1));
  def("AnalyzeAlphaHelixAxis", &AnalyzeAlphaHelixAxis, (arg("traj"), arg("protein_segment"), arg("directions"), arg("centers"), arg("stride")=1));
  def("AnalyzeBestFitLine", &AnalyzeBestFitLine, (arg("traj"), arg("protein_segment"), arg("directions"), arg("centers"), arg("stride")=1));
  def("AnalyzeBestFitPlane", &AnalyzeBestFitPlane, (arg("traj"), arg("protein_segment"), arg("normals"), arg("origins"), arg("stride")=1));
  def("AnalyzeHelicity", &AnalyzeHelicity, (arg("traj"), arg("protein_segment"), arg("stride")=1));
  def("CreateMeanStructure", &CreateMeanStructure, (arg("traj"), arg("selection"), arg("from")=0, arg("to")=-1, arg("stride")=1 ));
}
