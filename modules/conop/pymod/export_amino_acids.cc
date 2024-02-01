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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <ost/conop/amino_acids.hh>

using namespace boost::python;

using namespace ost::conop;

void export_AminoAcids()
{
  enum_<AminoAcid>("AminoAcid")
    .value("ALA", ALA)
    .value("ARG", ARG)
    .value("ASN", ASN)
    .value("ASP", ASP)
    .value("GLN", GLN) 
    .value("GLU", GLU)
    .value("LYS", LYS)
    .value("SER", SER)
    .value("CYS", CYS)
    .value("TYR", TYR)
    .value("THR", THR)
    .value("VAL", VAL)
    .value("ILE", ILE)
    .value("LEU", LEU)
    .value("GLY", GLY)
    .value("PRO", PRO)
    .value("MET", MET)
    .value("HIS", HIS)
    .value("PHE", PHE)
    .value("TRP", TRP)
    .value("XXX", XXX)
    .export_values()
  ;
  
  class_<AminoAcidSet>("AminoAcidSet", init<>())
    .def("Add", &AminoAcidSet::Add)
    .def("Remove", &AminoAcidSet::Remove)
    .def("Contains", &AminoAcidSet::Contains)
    .def("Empty", &AminoAcidSet::Empty)
    .def("CreatePolarSet", 
         &AminoAcidSet::CreatePolarSet).staticmethod("CreatePolarSet")
    .def("CreateApolarSet", 
         &AminoAcidSet::CreateApolarSet).staticmethod("CreateApolarSet")
    .def("CreateAromaticSet", 
         &AminoAcidSet::CreateAromaticSet).staticmethod("CreateAromaticSet")
    .def("CreateThreeStateSet", 
         &AminoAcidSet::CreateThreeStateSet).staticmethod("CreateThreeStateSet")
    .def("CreatePseudoSet", 
         &AminoAcidSet::CreatePseudoSet).staticmethod("CreatePseudoSet")
    .def("CreateCompleteSet", 
         &AminoAcidSet::CreateCompleteSet).staticmethod("CreateCompleteSet")
    .def("CreateSet", &AminoAcidSet::CreateSet).staticmethod("CreateSet")
    .def(self_ns::str(self))
  ;
  
  class_<AminoAcidAlphabet>("AminoAcidAlphabet", init<>())
    .def(vector_indexing_suite<AminoAcidAlphabet>())
  ;
  
  def("ResidueToAminoAcid",&ResidueToAminoAcid);
  def("ResidueNameToAminoAcid",&ResidueNameToAminoAcid);
  def("AminoAcidToResidueName",&AminoAcidToResidueName);
  def("OneLetterCodeToResidueName",&OneLetterCodeToResidueName);
  def("ResidueNameToOneLetterCode",&ResidueNameToOneLetterCode); 
  def("OneLetterCodeToAminoAcid",&OneLetterCodeToAminoAcid);
}

