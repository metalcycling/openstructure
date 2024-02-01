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
#include <ost/mol/mol.hh>

#include "sec_structure_segments.hh"

/* 
  Author: Marco Biasini
 */
 
namespace ost { namespace mol { namespace alg {


namespace {

class SecStrExtractor : public EntityVisitor {
public:
  SecStrExtractor(): index_(0) {}
  virtual bool VisitResidue(const ResidueHandle& r)
  {
    if (sec_.empty()) {
      sec_.push_back(SecStructureSegment(index_, index_, r.GetSecStructure()));
    } else if (sec_.back().ss_type!=r.GetSecStructure()) {
      sec_.push_back(SecStructureSegment(index_, index_, r.GetSecStructure()));
    } else {
      sec_.back().last+=1;
    }    
    index_+=1;
    return false;
  }
  
  const SecStructureSegments& GetList() const { return sec_; }
protected:
  SecStructureSegments  sec_;
  int                   index_;
};


// pretty bad case of code duplication.

class HelicalExtractor : public SecStrExtractor {
public:
  HelicalExtractor(): SecStrExtractor(), extend_(false)
  { }
  
  virtual bool VisitResidue(const ResidueHandle& r)
  {
    mol::SecStructure helical(mol::SecStructure::ALPHA_HELIX);
    if (r.GetSecStructure().IsHelical()) {
      if (sec_.empty()) {
        sec_.push_back(SecStructureSegment(index_, index_, helical));
        extend_=true;
      } else if (extend_) {
        sec_.back().last+=1;
      } else {
        sec_.push_back(SecStructureSegment(index_, index_, helical));
        extend_=true;        
      }
    } else {
      extend_=false;
    }    
    index_+=1;
    return false;
  }
private:
  bool          extend_;
};

class ExtendedExtractor : public SecStrExtractor {
public:
  ExtendedExtractor(): SecStrExtractor(), extend_(false)
  { }
  
  virtual bool VisitResidue(const ResidueHandle& r)
  {
    mol::SecStructure extended(mol::SecStructure::EXTENDED);
    if (r.GetSecStructure().IsExtended()) {
      if (sec_.empty()) {
        sec_.push_back(SecStructureSegment(index_, index_, extended));
        extend_=true;
      } else if (extend_) {
        sec_.back().last+=1;
      } else {
        sec_.push_back(SecStructureSegment(index_, index_, extended));
        extend_=true;        
      }
    } else {
      extend_=false;
    }        
    index_+=1;
    return false;
  }
private:
  bool          extend_;
};


template <typename E, typename C>
inline SecStructureSegments extract_ss(C c)
{
  E extractor;
  c.Apply(extractor);
  return extractor.GetList();
}

}

SecStructureSegments ExtractSecStructureSegments(const ChainHandle& chain)
{
  return extract_ss<SecStrExtractor, ChainHandle>(chain);
}

SecStructureSegments ExtractSecStructureSegments(const ChainView& chain)
{
  return extract_ss<SecStrExtractor, ChainView>(chain);
}
  
SecStructureSegments ExtractHelicalSegments(const ChainHandle& chain)
{
  return extract_ss<HelicalExtractor, ChainHandle>(chain);
}

SecStructureSegments ExtractHelicalSegments(const ChainView& chain)
{
  return extract_ss<HelicalExtractor, ChainView>(chain);
}


SecStructureSegments ExtractExtendedSegments(const ChainHandle& chain)
{
  return extract_ss<ExtendedExtractor, ChainHandle>(chain);
}

SecStructureSegments ExtractExtendedSegments(const ChainView& chain)
{
  return extract_ss<ExtendedExtractor, ChainView>(chain);
}



}}}
