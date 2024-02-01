//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#ifndef COMPOSITE_MASK_HH_
#define COMPOSITE_MASK_HH_

#include "mask_base.hh"

namespace ost { namespace img {

class OpBase
{
public:
  virtual ~OpBase() {}
  virtual OpBase* Clone()=0;
  virtual bool operator()(const MaskPtr& lhs,const  MaskPtr& rhs, const Vec2& v) const =0;
  virtual String GetName() const =0;
};

class AndOp:public OpBase
{
public:
  virtual OpBase* Clone(){return new AndOp();}
  virtual bool operator()(const MaskPtr& lhs,const  MaskPtr& rhs, const Vec2& v) const;
  virtual String GetName() const {return "AND";}
};

class OrOp:public OpBase
{
public:
  virtual OpBase* Clone(){return new OrOp();}
  virtual bool operator()(const MaskPtr& lhs,const  MaskPtr& rhs, const Vec2& v) const;
  virtual String GetName() const {return "OR";}
};

class XorOp:public OpBase
{
public:
  virtual OpBase* Clone(){return new XorOp();}
  virtual bool operator()(const MaskPtr& lhs,const  MaskPtr& rhs, const Vec2& v) const;
  virtual String GetName() const {return "XOR";}
};
typedef boost::shared_ptr<OpBase> OpPtr;


class DLLEXPORT_OST_IMG_BASE CompositeMask: public MaskBase
{
public:
  CompositeMask(const MaskPtr& lhs, const MaskPtr& rhs,const String& op);
  CompositeMask(const CompositeMask& m);

  String GetOperatorName() {
    return op_->GetName();
  }

  virtual ~CompositeMask();
  virtual MaskPtr Clone();
  virtual bool IsInside(const Vec2& v);
  
  virtual void Shift(const Vec2& v);
  virtual void Expand(Real d);
  virtual void Scale(Real d);

  // traverse mask hierarchy with a visitor
  virtual void Apply(MaskVisitor& v);


protected:
  static OpBase* create_op_(const String& name);
  MaskPtr lhs_;
  MaskPtr rhs_;
  OpBase* op_;
};

}}//ns


#endif /*COMPOSITE_MASK_HH_*/
