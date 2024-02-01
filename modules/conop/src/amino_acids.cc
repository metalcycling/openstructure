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
#include "amino_acids.hh"

#include <boost/version.hpp>
#if BOOST_VERSION>=103800
#define BOOST_SPIRIT_USE_OLD_NAMESPACE
#include <boost/spirit/include/classic_symbols.hpp>
#else
#include <boost/spirit/symbols.hpp>
#endif


#include <ost/mol/chain_handle.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/log.hh>
/*
 Author: Marco Biasini
 */
namespace ost { namespace conop {

using namespace boost::spirit;

namespace {

struct AminoAcidKeys : public symbols<AminoAcid> {
  AminoAcidKeys() {
    add("ALA", ALA)
       ("ARG", ARG)
       ("ASN", ASN)
       ("ASP", ASP)
       ("GLN", GLN) 
       ("GLU", GLU)
       ("LYS", LYS)
       ("SER", SER)
       ("CYS", CYS)
       ("TYR", TYR)
       ("THR", THR)
       ("VAL", VAL)
       ("ILE", ILE)
       ("LEU", LEU)
       ("PRO", PRO)
       ("GLY", GLY)
       ("MET", MET)
       ("HIS", HIS)
       ("TRP", TRP)
       ("PHE", PHE);
  }
};
} // anon namespace
   
AminoAcid ResidueToAminoAcid(const mol::ResidueHandle& r) {
  const mol::ResidueKey& key=r.GetKey();
  // use boost::spirits symbol table implementation that is particularly well 
  // suited for this job as it only has to traverse the keys once.
  static AminoAcidKeys aa_keys = AminoAcidKeys();
  AminoAcid* aa=find(aa_keys, key.c_str());
  if (aa)
    return *aa;
    
  return XXX;
}

String AminoAcidToResidueName(AminoAcid aa)
{
  switch(aa) {
    case ALA:
      return "ALA";
    case ARG:
      return "ARG";
    case ASN:
      return "ASN";
    case ASP:
      return "ASP";
    case GLN:
      return "GLN";
    case GLU:
      return "GLU";
    case LYS:
      return "LYS";
    case SER:
      return "SER";
    case CYS:
      return "CYS";
    case TYR:  
      return "TYR";
    case TRP:
      return "TRP";
    case THR:
      return "THR";
    case VAL:
      return "VAL";
    case ILE:
      return "ILE";
    case MET:
      return "MET";
    case LEU:  
      return "LEU";
    case GLY:  
      return "GLY";
    case PRO:  
      return "PRO";
    case HIS:  
      return "HIS";
    case PHE:
      return "PHE";
    default:
      return "UNK";
  }
}

String OneLetterCodeToResidueName(char olc)
{
  char upper_olc=toupper(olc);
  switch(upper_olc) {
    case 'A':
      return "ALA";
    case 'R':
      return "ARG";
    case 'N':
      return "ASN";
    case 'D':
      return "ASP";
    case 'Q':
      return "GLN";
    case 'E':
      return "GLU";
    case 'K':
      return "LYS";
    case 'S':
      return "SER";
    case 'C':
      return "CYS";
    case 'Y':  
      return "TYR";
    case 'W':
      return "TRP";
    case 'T':
      return "THR";
    case 'V':
      return "VAL";
    case 'I':
      return "ILE";
    case 'M':
      return "MET";
    case 'L':  
      return "LEU";
    case 'G':  
      return "GLY";
    case 'P':  
      return "PRO";
    case 'H':  
      return "HIS";
    case 'F':
      return "PHE";
    default:
      return "UNK";
  }
}

AminoAcid ResidueNameToAminoAcid(String rn) 
{
  std::transform(rn.begin(),rn.end(),rn.begin(),toupper);
  static AminoAcidKeys aa_keys = AminoAcidKeys();
  AminoAcid* aa=find(aa_keys, rn.c_str());
  if (aa)
    return *aa;
    
  return XXX;
}

char ResidueNameToOneLetterCode(String rn)
{
  if (rn.empty()) return 'X';

  std::transform(rn.begin(),rn.end(),rn.begin(),toupper);

  switch(rn[0]){
    case 'A': {
      if (rn == "ALA") return 'A';
      if (rn == "ARG") return 'R';      
      if (rn == "ASN") return 'N';
      if (rn == "ASP") return 'D';
      break;
    }

    case 'C': {
      if (rn == "CYS") return 'C';
      break;;
    }
    
    case 'G': {
      if (rn == "GLN") return 'Q';
      if (rn == "GLU") return 'E';
      if (rn == "GLY") return 'G';
      break;
    }

    case 'H': {
      if (rn == "HIS") return 'H';
      break;
    }

    case 'I': {
      if (rn == "ILE") return 'I';
      break;
    }

    case 'L': {
      if (rn == "LEU") return 'L';
      if (rn == "LYS") return 'K';
      break;
    }

    case 'M': {
      if (rn == "MET") return 'M';
      break;
    }  

    case 'P': {
      if (rn == "PRO") return 'P';
      if (rn == "PHE") return 'F';
      break;
    } 

    case 'S': {
      if (rn == "SER") return 'S';
      break;
    }  

    case 'T': {
      if (rn == "TYR") return 'Y';
      if (rn == "TRP") return 'W';
      if (rn == "THR") return 'T';
      break;
    }  

    case 'V': {
      if (rn == "VAL") return 'V';
      break;
    }  

    default: break;
  }

  return 'X';
}

AminoAcid OneLetterCodeToAminoAcid(char olc)
{
  char upper_olc=toupper(olc);
  switch(upper_olc) {
    case 'A':
      return ALA;
    case 'R':   
      return ARG;
    case 'N':   
      return ASN;
    case 'D':   
      return ASP;
    case 'Q':   
      return GLN;
    case 'E':   
      return GLU;
    case 'K':   
      return LYS;
    case 'S':   
      return SER;
    case 'C':   
      return CYS;
    case 'Y':   
      return TYR;
    case 'W':   
      return TRP;
    case 'T':   
      return THR;
    case 'V':   
      return VAL;
    case 'I':   
      return ILE;
    case 'M':   
      return MET;
    case 'L':   
      return LEU;
    case 'G':   
      return GLY;
    case 'P':   
      return PRO;
    case 'H':   
      return HIS;
    case 'F':   
      return PHE;
    default:
      return XXX;
  }
}

AminoAcidSet AminoAcidSet::CreatePolarSet()
{
  AminoAcidSet polar;
  polar.Add(ARG);
  polar.Add(ASN);
  polar.Add(ASP);
  polar.Add(GLN);
  polar.Add(GLU);
  polar.Add(SER);
  polar.Add(HIS);
  polar.Add(LYS);
  polar.Add(THR);
  return polar;  
}

AminoAcidSet AminoAcidSet::CreateAromaticSet()
{
  AminoAcidSet aromatic;
  aromatic.Add(TRP);
  aromatic.Add(TYR);
  aromatic.Add(PHE);
  return aromatic;  
}

AminoAcidSet AminoAcidSet::CreateApolarSet()
{
  AminoAcidSet apolar;
  apolar.Add(ALA);
  apolar.Add(VAL);
  apolar.Add(LEU);
  apolar.Add(ILE);
  apolar.Add(GLY);
  apolar.Add(PRO);
  apolar.Add(CYS);
  apolar.Add(MET);
  return apolar;
}

AminoAcidSet AminoAcidSet::CreateSet(AminoAcid aa)
{
  AminoAcidSet one_aa_set;
  one_aa_set.Add(aa);
  return one_aa_set;
}

std::vector<AminoAcidSet> AminoAcidSet::CreateCompleteSet() {
  std::vector<AminoAcidSet> alphabet;

  alphabet.push_back(AminoAcidSet::CreateSet(ALA));
  alphabet.push_back(AminoAcidSet::CreateSet(ARG));
  alphabet.push_back(AminoAcidSet::CreateSet(ASN));
  alphabet.push_back(AminoAcidSet::CreateSet(ASP));
  alphabet.push_back(AminoAcidSet::CreateSet(GLN));
  alphabet.push_back(AminoAcidSet::CreateSet(GLU));
  alphabet.push_back(AminoAcidSet::CreateSet(LYS));
  alphabet.push_back(AminoAcidSet::CreateSet(SER));
  alphabet.push_back(AminoAcidSet::CreateSet(CYS));
  alphabet.push_back(AminoAcidSet::CreateSet(MET));
  alphabet.push_back(AminoAcidSet::CreateSet(TRP));
  alphabet.push_back(AminoAcidSet::CreateSet(TYR));
  alphabet.push_back(AminoAcidSet::CreateSet(THR));
  alphabet.push_back(AminoAcidSet::CreateSet(VAL));
  alphabet.push_back(AminoAcidSet::CreateSet(ILE));
  alphabet.push_back(AminoAcidSet::CreateSet(LEU));
  alphabet.push_back(AminoAcidSet::CreateSet(GLY));
  alphabet.push_back(AminoAcidSet::CreateSet(PRO));
  alphabet.push_back(AminoAcidSet::CreateSet(HIS));
  alphabet.push_back(AminoAcidSet::CreateSet(PHE));
  return alphabet;
}


std::vector<AminoAcidSet> AminoAcidSet::CreateThreeStateSet() {
  std::vector<AminoAcidSet> alphabet;

  alphabet.push_back(CreatePolarSet());
  alphabet.push_back(CreateApolarSet());
  alphabet.push_back(CreateAromaticSet());
  return alphabet;
}

std::vector<AminoAcidSet> AminoAcidSet::CreatePseudoSet() {
  std::vector<AminoAcidSet> alphabet;
  AminoAcidSet full_alphabet;

  full_alphabet.Add(ALA);
  full_alphabet.Add(VAL);
  full_alphabet.Add(LEU);
  full_alphabet.Add(ILE);
  full_alphabet.Add(GLY);
  full_alphabet.Add(PRO);
  full_alphabet.Add(CYS);
  full_alphabet.Add(MET);
  full_alphabet.Add(TRP);
  full_alphabet.Add(TYR);
  full_alphabet.Add(PHE);
  full_alphabet.Add(ARG);
  full_alphabet.Add(ASN);
  full_alphabet.Add(ASP);
  full_alphabet.Add(GLN);
  full_alphabet.Add(GLU);
  full_alphabet.Add(SER);
  full_alphabet.Add(HIS);
  full_alphabet.Add(LYS);
  full_alphabet.Add(THR);

  alphabet.push_back(full_alphabet);
  return alphabet;
}


AminoAcidSet::AminoAcidSet():
  bits_(0)
{
  
}

void AminoAcidSet::Add(AminoAcid amino_acid)
{
  bits_|=(1 << (amino_acid));
}

void AminoAcidSet::Remove(AminoAcid amino_acid)
{
  bits_&=~(1 << (amino_acid));
}

bool AminoAcidSet::Contains(AminoAcid amino_acid) const
{
  return (bits_ & (1 << (amino_acid)))!=0;
}

bool AminoAcidSet::Empty() const
{
  return bits_==0;
}

AminoAcidSet::Iterator AminoAcidSet::Begin() const
{
  // find first bit that is set
  int start=0;
  while (start<=XXX+1 && !(bits_ & (1 << start))) {
    ++start;
  }
  return AminoAcidSetIterator(bits_, start);
}

AminoAcidSet::Iterator AminoAcidSet::End() const
{
  return AminoAcidSetIterator(bits_, XXX+1);
}

bool AminoAcidSet::operator==(const AminoAcidSet& rhs) const
{
  return bits_==rhs.bits_;
}

bool AminoAcidSet::operator!=(const AminoAcidSet& rhs) const
{
  return !this->operator==(rhs);
}

std::ostream& operator<<(std::ostream& os, 
                         const AminoAcidSet& aa_set)
{
  os << "[";
  for (AminoAcidSet::Iterator i=aa_set.Begin(), e=aa_set.End(); i!=e; ++i) {
    os << " " << AminoAcidToResidueName(*i);    
  }
  os << " ]";
  return os;
}

}}
