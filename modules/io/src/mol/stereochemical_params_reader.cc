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
#include <boost/filesystem/fstream.hpp>
#include <ost/platform.hh>
#include <ost/io/stereochemical_params_reader.hh>


namespace ost { namespace io {

ost::mol::alg::StereoChemicalProps ReadStereoChemicalPropsFile(String filename, bool check){
  boost::filesystem::path loc(filename);
  boost::filesystem::ifstream infile(loc);
  if (!infile) {
    std::stringstream serr;
    serr << "Could not find parameter file '" << filename << "'";   
    throw ost::Error(serr.str());
  }
  std::vector<String> stereo_chemical_props;
  String line;
  while (std::getline(infile, line))
  {
    std::stringstream line_stream(line);
    stereo_chemical_props.push_back(line);
  }

  ost::mol::alg::StereoChemicalParams bond_table;
  ost::mol::alg::StereoChemicalParams angle_table;
  ost::mol::alg::ClashingDistances nonbonded_table;

  bond_table = ost::mol::alg::FillStereoChemicalParams("Bond", stereo_chemical_props, check);
  // Angles
  angle_table = ost::mol::alg::FillStereoChemicalParams("Angle", stereo_chemical_props, check);
  // Not bonded
  nonbonded_table = ost::mol::alg::FillClashingDistances(stereo_chemical_props, check);

  return ost::mol::alg::StereoChemicalProps(bond_table, angle_table, nonbonded_table);
}

ost::mol::alg::StereoChemicalProps ReadStereoChemicalPropsFile(bool check) {
  String filename = GetStereoChemicalPropsFile();
  return ReadStereoChemicalPropsFile(filename, check);
}

String GetStereoChemicalPropsFile() {
  String filename;
  filename = ost::GetSharedDataPath() + "/stereo_chemical_props.txt";
  return filename;
}

}} // ns

