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
  Author: Gerardo Tauriello
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/conop/amino_acids.hh>

using namespace ost;
using namespace ost::conop;

BOOST_AUTO_TEST_SUITE(conop);

BOOST_AUTO_TEST_CASE(test_amino_acid_to_residue_name) {
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(ALA), "ALA");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(ARG), "ARG");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(ASN), "ASN");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(ASP), "ASP");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(CYS), "CYS");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(GLU), "GLU");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(GLN), "GLN");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(GLY), "GLY");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(HIS), "HIS");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(ILE), "ILE");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(LEU), "LEU");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(LYS), "LYS");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(MET), "MET");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(PHE), "PHE");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(PRO), "PRO");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(SER), "SER");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(THR), "THR");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(TRP), "TRP");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(TYR), "TYR");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(VAL), "VAL");
  BOOST_CHECK_EQUAL(AminoAcidToResidueName(XXX), "UNK");
}

BOOST_AUTO_TEST_CASE(test_one_letter_code_to_residue_name) {
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('A'), "ALA");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('R'), "ARG");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('N'), "ASN");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('D'), "ASP");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('C'), "CYS");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('E'), "GLU");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('Q'), "GLN");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('G'), "GLY");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('H'), "HIS");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('I'), "ILE");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('L'), "LEU");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('K'), "LYS");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('M'), "MET");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('F'), "PHE");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('P'), "PRO");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('S'), "SER");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('T'), "THR");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('W'), "TRP");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('Y'), "TYR");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('V'), "VAL");
  // should also work if not uppercase!
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('a'), "ALA");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('r'), "ARG");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('n'), "ASN");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('d'), "ASP");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('c'), "CYS");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('e'), "GLU");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('q'), "GLN");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('g'), "GLY");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('h'), "HIS");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('i'), "ILE");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('l'), "LEU");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('k'), "LYS");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('m'), "MET");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('f'), "PHE");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('p'), "PRO");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('s'), "SER");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('t'), "THR");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('w'), "TRP");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('y'), "TYR");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('v'), "VAL");
  // unknown ones
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('B'), "UNK");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('O'), "UNK");
  BOOST_CHECK_EQUAL(OneLetterCodeToResidueName('Z'), "UNK");
}

BOOST_AUTO_TEST_CASE(test_one_letter_code_to_amino_acid) {
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('A'), ALA);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('R'), ARG);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('N'), ASN);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('D'), ASP);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('C'), CYS);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('E'), GLU);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('Q'), GLN);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('G'), GLY);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('H'), HIS);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('I'), ILE);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('L'), LEU);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('K'), LYS);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('M'), MET);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('F'), PHE);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('P'), PRO);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('S'), SER);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('T'), THR);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('W'), TRP);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('Y'), TYR);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('V'), VAL);
  // should also work if not uppercase!
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('a'), ALA);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('r'), ARG);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('n'), ASN);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('d'), ASP);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('c'), CYS);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('e'), GLU);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('q'), GLN);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('g'), GLY);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('h'), HIS);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('i'), ILE);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('l'), LEU);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('k'), LYS);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('m'), MET);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('f'), PHE);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('p'), PRO);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('s'), SER);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('t'), THR);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('w'), TRP);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('y'), TYR);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('v'), VAL);
  // unknown ones
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('B'), XXX);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('O'), XXX);
  BOOST_CHECK_EQUAL(OneLetterCodeToAminoAcid('Z'), XXX);
}

BOOST_AUTO_TEST_CASE(test_residue_name_to_one_letter_code) {
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("ALA"), 'A');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("ARG"), 'R');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("ASN"), 'N');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("ASP"), 'D');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("CYS"), 'C');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("GLU"), 'E');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("GLN"), 'Q');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("GLY"), 'G');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("HIS"), 'H');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("ILE"), 'I');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("LEU"), 'L');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("LYS"), 'K');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("MET"), 'M');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("PHE"), 'F');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("PRO"), 'P');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("SER"), 'S');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("THR"), 'T');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("TRP"), 'W');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("TYR"), 'Y');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("VAL"), 'V');
  // should also work if not uppercase!
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Ala"), 'A');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Arg"), 'R');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Asn"), 'N');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Asp"), 'D');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Cys"), 'C');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Glu"), 'E');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Gln"), 'Q');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Gly"), 'G');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("His"), 'H');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Ile"), 'I');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Leu"), 'L');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Lys"), 'K');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Met"), 'M');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Phe"), 'F');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Pro"), 'P');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Ser"), 'S');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Thr"), 'T');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Trp"), 'W');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Tyr"), 'Y');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("Val"), 'V');
  // non-standard amino acids
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("MSE"), 'X');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("SEC"), 'X');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("PYL"), 'X');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("PCA"), 'X');
  // unknown ones
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("MUH"), 'X');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("FOO"), 'X');
  BOOST_CHECK_EQUAL(ResidueNameToOneLetterCode("BAR"), 'X');
}

BOOST_AUTO_TEST_CASE(test_residue_name_to_amino_acid) {
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("ALA"), ALA);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("ARG"), ARG);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("ASN"), ASN);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("ASP"), ASP);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("CYS"), CYS);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("GLU"), GLU);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("GLN"), GLN);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("GLY"), GLY);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("HIS"), HIS);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("ILE"), ILE);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("LEU"), LEU);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("LYS"), LYS);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("MET"), MET);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("PHE"), PHE);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("PRO"), PRO);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("SER"), SER);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("THR"), THR);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("TRP"), TRP);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("TYR"), TYR);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("VAL"), VAL);
  // should also work if not uppercase!
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Ala"), ALA);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Arg"), ARG);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Asn"), ASN);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Asp"), ASP);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Cys"), CYS);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Glu"), GLU);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Gln"), GLN);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Gly"), GLY);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("His"), HIS);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Ile"), ILE);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Leu"), LEU);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Lys"), LYS);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Met"), MET);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Phe"), PHE);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Pro"), PRO);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Ser"), SER);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Thr"), THR);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Trp"), TRP);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Tyr"), TYR);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("Val"), VAL);
  // non-standard amino acids
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("MSE"), XXX);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("SEC"), XXX);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("PYL"), XXX);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("PCA"), XXX);
  // unknown ones
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("MUH"), XXX);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("FOO"), XXX);
  BOOST_CHECK_EQUAL(ResidueNameToAminoAcid("BAR"), XXX);
}

BOOST_AUTO_TEST_SUITE_END();
