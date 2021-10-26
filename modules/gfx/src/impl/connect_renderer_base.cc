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
  Authors: Ansgar Philippsen, Marco Biasini
 */
#include <ost/gfx/scene.hh>
#include <ost/gfx/impl/connect_renderer_base.hh>

namespace ost { namespace gfx { namespace impl {


namespace {

typedef enum {
  SEL_ATOM,
  SEL_BOND,
  SEL_TORSION
} SelectionMode;

class Selector: public mol::EntityVisitor {
public:
  Selector(const geom::Line3& l, float zlim, float dlim):
    mode(SEL_ATOM), // atom selection is default
    line_(l),
    dlim_(dlim),
    zmin_(zlim),
    atom_(),
    bond_()
  {}

  virtual bool VisitAtom(const mol::AtomHandle& atom) {
    if(mode==SEL_ATOM) {
      float d=geom::Distance(line_,atom.GetPos());
      if((fixed_radius && d<=dlim_) || 
         (!fixed_radius && d<=atom.GetRadius())) {
        geom::Vec3 v=Scene::Instance().Project(atom.GetPos());
        if(v[2]>0.0 && v[2]<zmin_) {
          zmin_=v[2];
          atom_=atom;
        }
      }
    }
    return true;
  }

  virtual bool VisitBond(const mol::BondHandle& bond) {
    if(mode==SEL_BOND) {
      float d =geom::Distance(line_,bond.GetPos());

      if(d<=dlim_) {
        geom::Vec3 v=Scene::Instance().Project(bond.GetPos());
        if (v[2]>0.0 && v[2]<zmin_) {
          zmin_=v[2];
          bond_=bond;
        }
      }
    }
    return true;
  }

  virtual bool VisitTorsion(const mol::TorsionHandle& torsion) {
    if(mode==SEL_TORSION) {
      float d =geom::Distance(line_,torsion.GetPos());
      if(d<=dlim_) {
        geom::Vec3 v=Scene::Instance().Project(torsion.GetPos());
        if (v[2]>0.0 && v[2]<zmin_) {
          zmin_=v[2];
          torsion_=torsion;
        }
      }
    }
    return true;
  }
  
  const mol::AtomHandle GetAtom() const {return atom_;}
  const mol::BondHandle GetBond() const {return bond_;}
  const mol::TorsionHandle GetTorsion() const {return torsion_;}

  void SetMode(SelectionMode m) {mode=m;}
public:
  bool fixed_radius;
private:
  SelectionMode mode; // atom, bond or torsion selection
  geom::Line3 line_;
  float dlim_,zmin_;
  mol::AtomHandle atom_;
  mol::BondHandle bond_;
  mol::TorsionHandle torsion_;
};



class ConnectBuilder: public mol::EntityViewVisitor {
public:
  ConnectBuilder(GfxView* v, bool af, bool bf):
    v_(v), af_(af), bf_(bf) {}

  virtual bool VisitAtom(const mol::AtomView& a) {
    if(af_) {
      v_->AddAtom(a);
    }
    return true;
  }

  virtual bool VisitBond(const mol::BondHandle& b) {
    if(bf_)
      v_->AddBond(b);
    return true;
  }

private:
  GfxView* v_;
  bool af_,bf_;
};

template <typename T1>
inline void apply_color_op(ConnectRendererBase* rend, GfxView* v, T1 get_col, const ColorOp& op)
{
  if ((op.GetMask() & MAIN_COLOR)==0) {
    return;
  }
  rend->UpdateViews();
  mol::EntityView view;
  if(op.IsSelectionOnly()){
    view = rend->GetEffectiveView().Select(op.GetSelection(),op.GetSelectionFlags());
  }
  else{
    view = op.GetView();
  }
  for(AtomEntryMap::iterator it=v->atom_map.begin();it!=v->atom_map.end();++it){
    if(view.ViewForHandle(it->second.atom)){
      it->second.color=get_col.ColorOfAtom(it->second.atom).second;
    }
  }
};

} // anon ns
  
ConnectRendererBase::ConnectRendererBase(): pick_radius_(0.0)
{
  
}

void ConnectRendererBase::PrepareRendering()
{
}

void ConnectRendererBase::UpdateViews()
{
  if (state_ & DIRTY_VIEW) {
    mol::EntityView view=this->GetEffectiveView();
    state_&=~DIRTY_VIEW;    
    view_.Clear();  
    if (view.IsValid()) {
      ConnectBuilder connect_builder(&view_, true, true);
      view.Apply(connect_builder);
    }    
  }
  if (sel_state_ & DIRTY_VIEW) {
    sel_state_&=~DIRTY_VIEW;
    sel_view_.Clear();    
    if (this->HasSelection()) {
      ConnectBuilder connect_builder(&sel_view_, true, true);
      sel_.Apply(connect_builder);
      for(AtomEntryMap::iterator it=sel_view_.atom_map.begin();
          it!=sel_view_.atom_map.end();++it) {
        it->second.color=this->GetSelectionColor();
      }
    }
  }
}

geom::AlignedCuboid ConnectRendererBase::GetBoundingBox() const
{
  assert((state_ & DIRTY_VIEW)==0);
  if (view_.atom_map.empty()) {
    throw Error("Can't calculate bounding box of empty renderer");
  }
  Real rmax=std::numeric_limits<Real>::max();
  Real rmin=-std::numeric_limits<Real>::max();  
  geom::Vec3 minc(rmax, rmax, rmax), maxc(rmin, rmin, rmin);
  for (AtomEntryMap::const_iterator it=view_.atom_map.begin();
       it!=view_.atom_map.end();++it) {
    minc=geom::Min(minc, it->second.atom.GetPos());
    maxc=geom::Max(maxc, it->second.atom.GetPos());
  }
  // todo base this on max radius...
  minc-=geom::Vec3(2.0,2.0,2.0);
  maxc+=geom::Vec3(2.0,2.0,2.0);
  return geom::AlignedCuboid(minc, maxc);
}

void ConnectRendererBase::Apply(const gfx::ByElementColorOp& op)
{
  apply_color_op(this,&view_,ByElementGetCol(),op);
  state_|=DIRTY_VA;
}

void ConnectRendererBase::Apply(const gfx::ByChainColorOp& op)
{
  apply_color_op(this,&view_,ByChainGetCol(op),op);
  state_|=DIRTY_VA;
}

void ConnectRendererBase::Apply(const gfx::UniformColorOp& op)
{
  apply_color_op(this,&view_,UniformGetCol(op.GetColor()),op);
  state_|=DIRTY_VA;
}

void ConnectRendererBase::Apply(const gfx::GradientLevelColorOp& op)
{
  apply_color_op(this,&view_,GradientLevelGetCol(op),op);
  state_|=DIRTY_VA;
}

void ConnectRendererBase::Apply(const gfx::EntityViewColorOp& op)
{
  apply_color_op(this,&view_,EntityViewGetCol(op),op);
  state_|=DIRTY_VA;
}

void ConnectRendererBase::Apply(const gfx::MapHandleColorOp& op)
{
  apply_color_op(this,&view_,MapHandleGetCol(op),op);
  state_|=DIRTY_VA;
}

void ConnectRendererBase::PickAtom(const geom::Line3& line, Real line_width,
                                   mol::AtomHandle& picked_atom)
{
  Real zlim=1.0;
  if (picked_atom.IsValid()) {
    geom::Vec3 pp=Scene::Instance().Project(picked_atom.GetPos());
    if (pp[2]>0.0 && pp[2]<zlim) {
       zlim=pp[2];
    }
  }
  float radius=pick_radius_==0.0 ? 0.1 : pick_radius_;
  Selector sel(line, zlim, radius);
  sel.fixed_radius=this->HasFixedPickRadius();
  sel.SetMode(SEL_ATOM);
  mol::EntityView v=this->GetEffectiveView();
  v.Apply(sel);
  if (sel.GetAtom().IsValid()) {
    picked_atom=sel.GetAtom();
  }
}


void ConnectRendererBase::PickBond(const geom::Line3& line, Real line_width,
                                   mol::BondHandle& picked_bond)
{
  Real zlim=1.0;
  if (picked_bond.IsValid()) {
    geom::Vec3 pp=Scene::Instance().Project(picked_bond.GetPos());
    if (pp[2]>0.0 && pp[2]<zlim) {
       zlim=pp[2];
    }
  }
  Selector sel(line, zlim, 0.4);
  sel.fixed_radius=true;
  sel.SetMode(SEL_BOND);
  mol::EntityView v=this->GetEffectiveView();
  v.Apply(sel);
  if (sel.GetBond().IsValid()) {
    picked_bond=sel.GetBond();
  }
}

}}}
