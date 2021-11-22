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
#include <boost/python.hpp>
#include <ost/config.hh>

#include <ost/mol/alg/entity_to_density.hh>

using namespace boost::python;
using namespace ost;
using namespace ost::mol::alg;

//"thin wrappers" for default parameters
BOOST_PYTHON_FUNCTION_OVERLOADS(etd_rosetta, EntityToDensityRosetta, 4, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(etd_scattering, EntityToDensityScattering, 4, 6)

void export_entity_to_density()
{
    def("EntityToDensityRosetta",EntityToDensityRosetta,etd_rosetta());
    def("EntityToDensityScattering",EntityToDensityScattering,etd_scattering());

    enum_<DensityType>("DensityType")
    .value("HIGH_RESOLUTION", HIGH_RESOLUTION)
    .value("LOW_RESOLUTION", LOW_RESOLUTION)
    .export_values()
    ;
}

