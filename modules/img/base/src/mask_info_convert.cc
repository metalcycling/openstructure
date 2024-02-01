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

#include <stack>

#include <ost/info/info_visitor.hh>
#include <ost/info/info.hh>
#include "mask_info_convert.hh"
#include "polygon_mask.hh"
#include "circle_mask.hh"
#include "extent_mask.hh"
#include "composite_mask.hh"
#include "inverted_mask.hh"
#include "spherical_mask.hh"
#include "mask_visitor.hh"

namespace ost{ namespace img {

namespace {



class PolygonInfoVisitor: public info::InfoConstVisitor
{
public:
  PolygonInfoVisitor():
    p_()
  {}

  Polygon2 GetPolygon() const {return p_;}

  virtual bool VisitGroup(const info::InfoGroup& group)
  {
    if(group.GetName()=="Polygon"){
      return true;
    } else if(group.GetName()=="Node"){
      p_.AddNode(geom::Vec2(group.GetItem("x").AsFloat(),
			    group.GetItem("y").AsFloat()));
      return false;
    }
    return true;
  }

  virtual void VisitGroupFinish(const info::InfoGroup& group)
  {}

private:
  Polygon2 p_;  
};

// this is used recursively
class MaskInfoVisitor: public info::InfoConstVisitor
{
public:
  MaskInfoVisitor():
    mask_list_()
  {}

  std::vector<MaskPtr> GetMaskList() const {return mask_list_;}

  virtual bool VisitGroup(const info::InfoGroup& group) 
  {
    String name=group.GetName();
    if(name=="Composite"){
      String opname=group.GetItem("Operator").GetValue();
      MaskInfoVisitor miv;
      group.Apply(miv,false); // subgroups
      std::vector<MaskPtr> mpl = miv.GetMaskList();
      mask_list_.push_back(MaskPtr(new CompositeMask(mpl[0],mpl[1],opname)));
      return false;
    }else if(name=="Invert"){
      MaskInfoVisitor miv;
      group.Apply(miv,false); // subgroups
      std::vector<MaskPtr> mpl = miv.GetMaskList();
      mask_list_.push_back(MaskPtr(new InvertedMask(mpl[0])));
      return false;
    }else if(name=="Circle"){
      Circle2 c;
      c.SetRadius(group.GetItem("Radius").AsFloat());
      c.SetCenter(geom::Vec2(group.GetItem("Center/x").AsFloat(),
			     group.GetItem("Center/y").AsFloat()));
      mask_list_.push_back(MaskPtr(new CircleMask(c)));
      return false;
    }else if(name=="Polygon"){
      PolygonInfoVisitor pvis;
      group.Apply(pvis);
      mask_list_.push_back(MaskPtr(new PolygonMask(pvis.GetPolygon())));
      return false;
    }else if(name=="Extent"){
      Extent e(Point(group.GetItem("Start/x").AsInt(),
		     group.GetItem("Start/y").AsInt()),
	       Point(group.GetItem("End/x").AsInt(),
		     group.GetItem("End/y").AsInt()));
      mask_list_.push_back(MaskPtr(new ExtentMask(e)));
      return false;
    }else if(name=="Sphere"){
      Sphere s;
      s.SetRadius(group.GetItem("Radius").AsFloat());
      s.SetOrigin(geom::Vec3(group.GetItem("Origin/x").AsFloat(),
			     group.GetItem("Origin/y").AsFloat(),
			     group.GetItem("Origin/z").AsFloat()));
      mask_list_.push_back(MaskPtr(new SphericalMask(s)));
      return false;
    }
    return true;    
  }

  virtual void VisitGroupFinish(const info::InfoGroup& group)
  {
  }

private:
  std::vector<MaskPtr> mask_list_;
};


// mask visitor for mask->info conversion
class MaskToInfoVisitor: public MaskVisitor
{
  typedef std::stack<info::InfoGroup> InfoGroupStack;
public:
  MaskToInfoVisitor(info::InfoGroup& g):
    igs_()
  {
    igs_.push(g);
  }
  virtual void VisitCircleMask(CircleMask& m) 
  {
    info::InfoGroup cg=igs_.top().CreateGroup("Circle");
    info::InfoGroup centerg=cg.CreateGroup("Center");
    centerg.CreateItem("x","").SetFloat(m.GetCenter()[0]);
    centerg.CreateItem("y","").SetFloat(m.GetCenter()[1]);
    cg.CreateItem("Radius","").SetFloat(m.GetRadius());
  }
  virtual void VisitCompositeMask(CompositeMask& m) 
  {
    info::InfoGroup compg=igs_.top().CreateGroup("Composite");
    compg.CreateItem("Operator",m.GetOperatorName());
    igs_.push(compg);
  }
  virtual void VisitCompositeMaskFinish(CompositeMask& m) 
  {
    igs_.pop();
  }
  virtual void VisitExtentMask(ExtentMask& m) 
  {
    info::InfoGroup extg=igs_.top().CreateGroup("Extent");
    info::InfoGroup startg=extg.CreateGroup("Start");
    startg.CreateItem("x","").SetInt(m.GetStart()[0]);
    startg.CreateItem("y","").SetInt(m.GetStart()[1]);
    info::InfoGroup endg=extg.CreateGroup("End");
    endg.CreateItem("x","").SetInt(m.GetEnd()[0]);
    endg.CreateItem("y","").SetInt(m.GetEnd()[1]);
  }
  virtual void VisitInvertedMask(InvertedMask& m) 
  {
    info::InfoGroup invertg=igs_.top().CreateGroup("Invert");
    igs_.push(invertg);
  }
  virtual void VisitInvertedMaskFinish(InvertedMask& m) 
  {
    igs_.pop();
  }
  virtual void VisitPolygonMask(PolygonMask& m) 
  {
    info::InfoGroup polyg=igs_.top().CreateGroup("Polygon");
    for(unsigned int nn=0;nn<m.GetNodeCount();++nn) {
      info::InfoGroup nodeg=polyg.CreateGroup("Node");
      geom::Vec2 node=m.GetNode(nn);
      nodeg.CreateItem("x","").SetFloat(node[0]);
      nodeg.CreateItem("y","").SetFloat(node[1]);
    }
  }
  virtual void VisitSphericalMask(SphericalMask& m) 
  {
    info::InfoGroup cg=igs_.top().CreateGroup("Sphere");
    info::InfoGroup centerg=cg.CreateGroup("Origin");
    centerg.CreateItem("x","").SetFloat(m.GetOrigin()[0]);
    centerg.CreateItem("y","").SetFloat(m.GetOrigin()[1]);
    centerg.CreateItem("z","").SetFloat(m.GetOrigin()[2]);
    cg.CreateItem("Radius","").SetFloat(m.GetRadius());
  }
private:
  InfoGroupStack igs_;
};

} // anon ns

MaskPtr DLLEXPORT InfoToMask(const info::InfoGroup& g)
{
  MaskInfoVisitor vis;
  g.Apply(vis);
  std::vector<MaskPtr> mlist = vis.GetMaskList();
  if(!mlist.empty()) {
    return mlist[0];
  }
  return MaskPtr();
}

void MaskToInfo(const MaskPtr& mptr ,info::InfoGroup& g)
{
  MaskToInfoVisitor vis(g);
  mptr->Apply(vis);
}

}} //ns
