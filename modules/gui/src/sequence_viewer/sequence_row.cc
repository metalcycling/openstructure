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

#include <ost/gfx/entity.hh>

#include "sequence_view_object.hh"
#include "sequence_row.hh"

#include <QtGui>
namespace ost { namespace gui {

SequenceRow::SequenceRow(const QString& name, seq::ConstSequenceHandle sequence, 
                         SequenceViewObject* parent) : 
  BaseRow(QFont("Courier",11),parent), name_(name), 
  name_font_(QFont("Courier",11)), sequence_(sequence)
{
  name_font_.setBold(true);
  name_font_.setItalic(true);
}

SequenceRow::SequenceRow(const QString& name, SequenceViewObject* parent) : 
  BaseRow(QFont("Courier",11),parent), name_(name), 
  name_font_(QFont("Courier",11))
{
  name_font_.setBold(true);
  name_font_.setItalic(true);
}

int SequenceRow::GetColumnCount() const
{
  return sequence_.GetLength()+1;
}

void SequenceRow::SetName(const QString& name)
{
  this->name_ = name;
}

const QString& SequenceRow::GetName() const
{
  return this->name_;
}

void SequenceRow::SetSequence(seq::ConstSequenceHandle sequence)
{
  this->sequence_ = sequence;
}

const seq::ConstSequenceHandle& SequenceRow::GetSequence() const
{
  return this->sequence_;
}

QVariant SequenceRow::GetData(int column, int role) const
{
  if(column<0 || column > sequence_.GetLength())return QVariant();

  if (role == Qt::ForegroundRole){
    return QColor(Qt::black);
  }

  if(column == 0) {
    if (role == Qt::DisplayRole){
      return QVariant(this->name_);
    }
    if (role == Qt::FontRole){
      return QVariant(this->name_font_);
    }
    if (role == Qt::TextAlignmentRole){
      return QVariant(Qt::AlignLeft|Qt::AlignVCenter);
    }
    if (role==Qt::ToolTipRole){
      return QVariant(this->name_);
    }
  }
  else if(column > 0) {
    if (role==Qt::DisplayRole) {
      return QVariant(QString(this->sequence_.GetOneLetterCode(column - 1)));
    }
    if (role==Qt::FontRole){
      return QVariant(this->GetFont());
    }
    if (role==Qt::SizeHintRole){
      return QVariant(this->GetCellSize());
    }
    if (role==Qt::ToolTipRole){
      return QVariant(QString(this->sequence_.GetOneLetterCode(column - 1)));
    }
  }
  return QVariant();
}

Qt::ItemFlags SequenceRow::Flags(int column) const
{
  if(column<0 || column >= this->GetColumnCount())return Qt::NoItemFlags;

  if(column==0){
    return Qt::ItemIsSelectable|Qt::ItemIsEnabled;
  }
  else if(column>0){
    return Qt::ItemIsSelectable|Qt::ItemIsEnabled;
  }
  return BaseRow::Flags(column);
}

void SequenceRow::DoubleClicked(int column)
{
  QSet<int> all;
  int seq_length = this->sequence_.GetLength();
  for(int i = 0; i < seq_length; i++){
    all.insert(i+1);
  }
  this->SetSelection(all,QSet<int>());
}

void SequenceRow::SetSelection(const QSet<int>& added, const QSet<int>& removed)
{
  SequenceViewObject* view_object = qobject_cast<SequenceViewObject*>(this->parent());
  if(view_object){
    if(gfx::EntityP entity = view_object->GetGfxObject()){
      mol::EntityView sel = entity->GetSelection();
      mol::EntityView view = this->sequence_.GetAttachedView().GetHandle().CreateEmptyView();

      QSetIterator<int> i(removed);
      while (i.hasNext()){
        int column = i.next()-1;
        if(column >= 0 && column < this->sequence_.GetLength()){
          if (mol::ResidueView rv = this->sequence_.GetResidue(this->sequence_.GetResidueIndex(column))) {
            view.AddResidue(rv, mol::ViewAddFlag::INCLUDE_ATOMS);
          }
        }
      }
      sel = mol::Difference(sel,view);
      view = this->sequence_.GetAttachedView().GetHandle().CreateEmptyView();
      i = QSetIterator<int>(added);
      while (i.hasNext()){
        int column = i.next()-1;
        if(column >= 0 && column < this->sequence_.GetLength()){
          if (mol::ResidueView rv=this->sequence_.GetResidue(this->sequence_.GetResidueIndex(column))) {
            view.AddResidue(rv, mol::ViewAddFlag::INCLUDE_ATOMS);
          }
        }
      }
      sel = mol::Union(sel,view);
      sel.AddAllInclusiveBonds();
      entity->SetSelection(sel);
    }
  }
}

}}
