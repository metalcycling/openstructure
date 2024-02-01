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
  Author: Stefan Scheuber
 */



#include <ost/mol/mol.hh>
#include <ost/mol/view_op.hh>

#include "sequence_row.hh"
#include "secstr_row.hh"

#include "sequence_view_object.hh"

#include <QtGui>
namespace ost { namespace gui {

const QString SequenceViewObject::properties_mode = "Highlight properties";
const QString SequenceViewObject::secondary_structure_mode = "Secondary structure";

AlignPropertiesPainter* SequenceViewObject::align_properties_painter = new AlignPropertiesPainter();
ConservationPainter* SequenceViewObject::conservation_painter = new ConservationPainter();
BackgroundPainter* SequenceViewObject::background_painter = new BackgroundPainter();
SeqSecStrPainter* SequenceViewObject::seq_secondary_structure_painter = new SeqSecStrPainter();
SeqSelectionPainter* SequenceViewObject::seq_selection_painter = new SeqSelectionPainter();
SeqTextPainter* SequenceViewObject::seq_text_painter = new SeqTextPainter();

SequenceViewObject::SequenceViewObject(seq::SequenceList& sequences, const QList<QString>& names, QObject *parent): BaseViewObject(parent), entity_(gfx::EntityP())
{
  this->Init();
  if(names.size() == sequences.GetCount()){
    for(int i=0; i<sequences.GetCount(); i++){
      seq::SequenceHandle seq = sequences[i];
      this->AddSequence(seq, names[i]);
    }
  }
  this->SetDisplayMode(properties_mode);
}

SequenceViewObject::SequenceViewObject(seq::SequenceHandle& sequence, const QString& name, QObject *parent): BaseViewObject(parent), entity_(gfx::EntityP())
{
  this->Init();
  this->AddSequence(sequence, name);
  this->SetDisplayMode(properties_mode);
}

SequenceViewObject::SequenceViewObject(mol::ChainView& chain, const QString& name, QObject *parent): BaseViewObject(parent), entity_(gfx::EntityP())
{
  this->Init();
  this->AddChain(chain, name);
  this->SetDisplayMode(properties_mode);
}

SequenceViewObject::SequenceViewObject(const gfx::EntityP& entity, QObject* parent): BaseViewObject(parent), entity_(entity)
{
  this->Init();
  mol::EntityView view =entity->GetView();
  for (mol::ChainViewList::const_iterator c=view.GetChainList().begin(),
       e1=view.GetChainList().end(); c!=e1; ++c) {
    mol::ChainView chain=*c;
    QString name = QString(entity->GetName().c_str());
    if (chain.GetName()!="" && chain.GetName()!=" ") {
      name= name + " ("+chain.GetName().c_str()+")";
    }
    this->AddChain(chain, name);
  }
  this->SetDisplayMode(secondary_structure_mode);
}

SequenceViewObject::SequenceViewObject(QObject* parent): BaseViewObject(parent), entity_(gfx::EntityP())
{
  this->Init();
  this->SetDisplayMode(properties_mode);
}

void SequenceViewObject::Init()
{
  this->AddDisplayMode(properties_mode);
  if(entity_){
    this->AddDisplayMode(secondary_structure_mode);
  }
}

void SequenceViewObject::AddSequence(seq::ConstSequenceHandle sequence, const QString& name)
{
  SequenceRow* new_row = new SequenceRow(name, sequence, this);
  new_row->InsertPainter(background_painter);
  new_row->InsertPainter(align_properties_painter);
  new_row->InsertPainter(seq_text_painter);
  new_row->InsertPainter(seq_selection_painter);
  rows_.append(new_row);
}

void SequenceViewObject::AddChain(mol::ChainView& chain, const QString& name)
{
  SecStrRow* new_row = new SecStrRow(name, chain, this);
  new_row->InsertPainter(background_painter);
  new_row->InsertPainter(seq_secondary_structure_painter);
  new_row->InsertPainter(seq_text_painter);
  new_row->InsertPainter(seq_selection_painter);
  rows_.append(new_row);
}

void SequenceViewObject::AttachGfxObject(gfx::EntityP& ent)
{
  entity_ = ent;
}

gfx::EntityP& SequenceViewObject::GetGfxObject()
{
  return entity_;
}

void SequenceViewObject::SetDisplayMode(const QString& mode)
{
  if(this->display_modes_.contains(mode) &&  mode != this->GetCurrentDisplayMode()){
    if(mode == properties_mode){
      for(int i=0 ; i<this->GetRowCount(); i++){
        BaseRow* row = this->GetRow(i);
        row->RemovePainter(seq_secondary_structure_painter);
        row->RemovePainter(conservation_painter);
        row->InsertPainter(align_properties_painter,1);
      }
    }
    else if(mode == secondary_structure_mode){
      for(int i=0 ; i<this->GetRowCount(); i++){
        BaseRow* row = this->GetRow(i);
        row->RemovePainter(align_properties_painter);
        row->RemovePainter(conservation_painter);
        row->InsertPainter(seq_secondary_structure_painter,1);
      }
    }
  }
  BaseViewObject::SetDisplayMode(mode);
}

QMap<int, QList<int> > SequenceViewObject::GetIndexesForView(const mol::EntityView& view)
{
  if(view.GetChainCount()==0){
    return QMap<int, QList<int> >();
  }
  else{
    QMap<int, QList<int> > selected_indexes;
    for(int i=0; i< rows_.size(); i++){
      if(SecStrRow* secstr_row = qobject_cast<SecStrRow*>(rows_[i])){
        mol::ChainView dst_chain=(secstr_row->GetChain());
        seq::ConstSequenceHandle seq = secstr_row->GetSequence();
        if (mol::ChainView src_chain=view.FindChain(dst_chain.GetName())) {
          // for each residue in the selection deduce index in sequence
          for (mol::ResidueViewList::const_iterator j=src_chain.GetResidueList().begin(),
             e2=src_chain.GetResidueList().end(); j!=e2; ++j) {
            mol::ResidueView dst_res=dst_chain.ViewForHandle(j->GetHandle());
            assert(dst_res.IsValid());
            int p=dst_res.GetIndex()+1;
            assert(p>=0 && p<=seq.GetLength());
            selected_indexes[i].append(p);
          }
        }
      }
    }
    return selected_indexes;
  }
}

QMap<int, QList<int> > SequenceViewObject::GetIndexesForSubject(const QString& subject, const QString& sequence_name)
{
  if(subject.size()==0){
    return QMap<int, QList<int> >();
  }
  QMap<int, QList<int> > selected_indexes;
  String subject_str = subject.toStdString();
  for(int i=0; i< rows_.size(); i++){
    if(SequenceRow* secstr_row = qobject_cast<SequenceRow*>(rows_[i])){
      if(sequence_name.size()==0 || sequence_name==secstr_row->GetName()){
        seq::ConstSequenceHandle seq = secstr_row->GetSequence();
        String seq_str=seq.GetString();
        size_t pos=0;
        size_t first=String::npos;
        while ((pos=seq_str.find(subject_str, pos))!=String::npos) {
          if (first==String::npos) {
            first=pos;
          }
          for(int j=0; j < subject.size(); j++){
            selected_indexes[i].append(pos+j+1);
          }
          pos+=subject.length();
        }
      }
    }
  }
  return selected_indexes;
}

}}
