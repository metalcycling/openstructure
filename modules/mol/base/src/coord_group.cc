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
#include <ost/invalid_handle.hh>
#include <ost/integrity_error.hh>
#include <ost/log.hh>
#include <ost/geom/transform.hh>
#include "in_mem_coord_source.hh"
#include "view_op.hh"
#include "mol.hh"
#include "coord_group.hh"


namespace ost { namespace mol {

CoordGroupHandle CreateCoordGroup(const AtomHandleList& atomlist)
{
  return CoordGroupHandle(CoordSourcePtr(new InMemCoordSource(atomlist)));
}

CoordGroupHandle::CoordGroupHandle():
  source_()
{}

CoordGroupHandle::CoordGroupHandle(CoordSourcePtr source):
  source_(source)
{
  this->CheckValidity();
}

bool CoordGroupHandle::IsValid() const
{
  this->CheckValidity();  
  return source_.get()!=0;
}

EntityHandle CoordGroupHandle::GetEntity() const
{
  this->CheckValidity();  
  return source_->GetEntity();
}

uint CoordGroupHandle::GetAtomCount() const
{
  this->CheckValidity();  
  return source_->GetAtomCount();
}

uint CoordGroupHandle::GetFrameCount() const
{
  this->CheckValidity();  
  return source_->GetFrameCount();
}

float CoordGroupHandle::GetDelta() const
{
  this->CheckValidity();  
  return source_->GetFrameDelta();
}

void CoordGroupHandle::SetDelta(float d)
{
  this->CheckValidity();  
  source_->SetFrameDelta(d);
}

float CoordGroupHandle::GetStartTime() const
{
  this->CheckValidity();  
  return source_->GetStartTime();
}

void CoordGroupHandle::SetStartTime(float t)
{
  this->CheckValidity();  
  source_->SetStartTime(t);
}

void CoordGroupHandle::SetFrameCellSize(uint frame, const geom::Vec3& size)
{
  this->CheckValidity();  
  if (source_->IsMutable()) {
    source_->SetFrameCellSize(frame, size);
  } else {
    throw IntegrityError("Can't add set frame cell size in immutable CoordGroup");
  }
}

void CoordGroupHandle::SetFramePositions(uint frame, 
                                         const geom::Vec3List& clist)
{
  this->CheckValidity();  
  if (source_->IsMutable()) {
    source_->SetFramePositions(frame, clist);
  } else {
    throw IntegrityError("Can't add set frame positions in immutable CoordGroup");
  }
}

  
geom::Vec3List CoordGroupHandle::GetFramePositions(uint frame)
  {
    this->CheckValidity();
    return *(this->GetFrame(frame));
  }
    
  
void CoordGroupHandle::CopyFrame(uint frame)
{
  this->CheckValidity();
  source_->CopyFrame(frame);
}

CoordGroupHandle::operator bool() const
{
  return this->IsValid();
}

//void CoordGroupHandle::AddFrame(const std::vector<geom::Vec3>& clist)
  void CoordGroupHandle::AddFrame(const geom::Vec3List& clist)
  {
  this->CheckValidity();
  if (source_->IsMutable()) {
    source_->AddFrame(clist);    
  } else {
    throw IntegrityError("Can't add frame to immutable CoordGroup");
  }
}

  void CoordGroupHandle::AddFrame(const geom::Vec3List& clist, const geom::Vec3& cell_size, const geom::Vec3& cell_angles)
{
  this->CheckValidity();
  if (source_->IsMutable()) {
    source_->AddFrame(clist,cell_size,cell_angles);    
  } else {
    throw IntegrityError("Can't add frame to immutable CoordGroup");
  }
}

void CoordGroupHandle::AddFrames(const CoordGroupHandle& cg)
{
  this->CheckValidity();
  if (source_->IsMutable()) {
    if (cg.GetAtomCount()!=this->GetAtomCount()) {
      throw IntegrityError("Atom number don't match");
    }
    for (size_t i=0; i<cg.GetFrameCount(); ++i) {
      CoordFramePtr p = cg.GetFrame(i);
      source_->AddFrame(*p,p->GetCellSize(),p->GetCellAngles());
    }
  } else {
    throw IntegrityError("Can't add frame to immutable CoordGroup");
  }
}

void CoordGroupHandle::CheckValidity() const
{
  if (!source_) {
    throw InvalidHandle();
  }
}

void CoordGroupHandle::SetAtomPos(uint frame, AtomHandle atom, const geom::Vec3& pos)
{
  this->CheckValidity();
  if (source_->IsMutable()) {  
    source_->SetAtomPos(frame,atom,pos);
  } else {
    throw IntegrityError("Can't set atom position in immutable CoordGroup");
  }
}

geom::Vec3 CoordGroupHandle::GetAtomPos(uint frame, AtomHandle atom) const
{
  this->CheckValidity();
  return source_->GetAtomPos(frame,atom);
}

CoordFramePtr CoordGroupHandle::GetFrame(uint frame) const
{
  this->CheckValidity();
  return source_->GetFrame(frame);
}

CoordFrame CoordGroupHandle::GetFrame2(uint frame)
{
  this->CheckValidity();
  return *(source_->GetFrame(frame));
}

AtomHandleList CoordGroupHandle::GetAtomList() const
{
  this->CheckValidity();
  return source_->GetAtomList();
}

void CoordGroupHandle::Capture()
{
  this->CheckValidity();
  if (source_->IsMutable()) {
    source_->Capture();    
  } else {
    throw IntegrityError("Can't capture. CoordGroup is immutable");
  }  
}

void CoordGroupHandle::CaptureInto(int pos)
{
  this->CheckValidity();
  if (source_->IsMutable()) {
    source_->CaptureInto(pos);    
  } else {
    throw IntegrityError("Can't capture. CoordGroup is immutable");
  }  
}

void CoordGroupHandle::Capture(uint frame)
{
  this->CheckValidity();
  if (source_->IsMutable()) {
    source_->Capture(frame);    
  } else {
    throw IntegrityError("Can't capture. CoordGroup is immutable");
  }  
}

CoordGroupHandle CoordGroupHandle::Filter(const EntityView& selected, int first, int last, int stride) const
{
  this->CheckValidity();
  std::vector<unsigned long> indices;
  EntityHandle new_ent;  
  if (!selected.IsValid()) {
    indices.reserve(this->GetAtomCount());
    for (size_t i=0;i<this->GetAtomCount(); ++i) {
      indices.push_back(i);
    }
    new_ent=this->GetEntity().Copy();
  } else {
    AtomViewList atoms=selected.GetAtomList();
    indices.reserve(atoms.size());
    for (AtomViewList::const_iterator i=atoms.begin(), 
         e=atoms.end(); i!=e; ++i) {
      indices.push_back(i->GetIndex());
    }
    new_ent=CreateEntityFromView(selected, false);
    //we have to correct the atom indices in the new entity
    AtomHandleList new_atoms=new_ent.GetAtomList();
    for (int i=0,e=atoms.size();i!=e;++i){
      new_atoms[i].SetIndex(i);
    }
  }
  CoordGroupHandle filtered_cg=CreateCoordGroup(new_ent.GetAtomList());
  std::vector<geom::Vec3> vecs(indices.size());
  if (last==-1) last=this->GetFrameCount();
  if (stride<=0) stride=1;
  for (int i=first; i<last; i=i+stride) {
    LOG_INFO("Filtering frame " << i << "/" << last);
    CoordFramePtr frame=this->GetFrame(i);
    for (std::vector<unsigned long>::const_iterator 
         j=indices.begin(), e2=indices.end(); j!=e2; ++j) {
      vecs[j-indices.begin()]=(*frame)[*j];
    }
    filtered_cg.AddFrame(vecs,(*frame).GetCellSize(),(*frame).GetCellAngles());
  }
  return filtered_cg;
}

void CoordGroupHandle::ApplyTransform(const geom::Transform& tf)
{
  this->CheckValidity();
  if (source_->IsMutable()) {
    source_->ApplyTransform(tf);
  } else {
    throw IntegrityError("Cannot apply transform, CoordGroup is immutable");
  }  
} 
  
}} // ns
