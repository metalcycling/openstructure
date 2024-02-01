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

#include <ost/seq/aligned_region.hh>
#include <ost/seq/alg/conservation.hh>

#include "sequence_row.hh"
#include "secstr_row.hh"

#include "painter.hh"
#include "background_painter.hh"
#include "seq_secstr_painter.hh"
#include "seq_selection_painter.hh"
#include "seq_text_painter.hh"

#include "alignment_view_object.hh"

#include <QtGui>
namespace ost { namespace gui {

namespace {
QMap<QString,int> GetGroupMap(){
  QMap<QString,int> map;
  map["G"]=1;
  map["A"]=1;
  map["V"]=1;
  map["L"]=1;
  map["I"]=1;
  map["F"]=2;
  map["Y"]=2;
  map["W"]=2;
  map["C"]=3;
  map["M"]=3;
  map["S"]=4;
  map["T"]=4;
  map["K"]=5;
  map["R"]=5;
  map["H"]=5;
  map["D"]=6;
  map["E"]=6;
  map["N"]=6;
  map["Q"]=6;
  map["P"]=7;
  return map;
}

}

QMap<QString,int> AlignmentViewObject::group_map_ = GetGroupMap();


const QString AlignmentViewObject::conservation_mode_1 = "Highlight conservation 1";
const QString AlignmentViewObject::conservation_mode_2 = "Highlight conservation 2";

AlignmentViewObject::AlignmentViewObject(const seq::AlignmentHandle& alignment, QObject* parent): SequenceViewObject(parent), alignment_(alignment)
{
  for(int i=0; i<alignment.GetCount(); i++){
    seq::ConstSequenceHandle seq_handle = alignment.GetSequence(i);
    this->AddSequence(seq_handle, seq_handle.GetName().c_str());
  }

  std::vector<Real> values = seq::alg::Conservation(alignment,false);

  gradient_.SetColorAt(0,gfx::Color(0,0,1));
  gradient_.SetColorAt(0.5,gfx::Color(1,1,1));
  gradient_.SetColorAt(1,gfx::Color(1,0,0));
  //Calculate Conservation Mode 1
  for(unsigned int i=0; i<values.size(); i++){
    gfx::Color color = gradient_.GetColorAt(values[i]);
    conservation_1_[i] = QColor(color.Red()*255,color.Green()*255,color.Blue()*255,100);
  }

  //Calculate Conservation Mode 2
  for(int j=0; j<alignment.GetLength(); j++){
    int group = 0;
    QString element = "";
    for(int i=0; i<alignment.GetCount();i++){
      QString code = QString(alignment.GetOneLetterCode(i,j));
      if(element.size()==0){
        element = code;
      }
      else if (element != code){
        element = "  ";
      }
      if(group_map_.contains(code) && group>=0){
        if(group == 0) {
            group = group_map_[code];
        }
        else if(group_map_[code] != group){
          group = -1;
        }
      }
      else{
        group = -1;
      }

    }
    if(element.size()==1){
      conservation_2_[j] = QColor(175,175,175);
    }
    else if(group > 0){
      conservation_2_[j] = QColor(200,200,200);
    }
    else{
      conservation_2_[j] = Qt::transparent;
    }
  }
  this->AddDisplayMode(conservation_mode_1);
  this->AddDisplayMode(conservation_mode_2);
}


QVariant AlignmentViewObject::GetData(int row, int column, int role)
{
  if(column > 0){
    if(this->GetCurrentDisplayMode() == conservation_mode_1){
      if(role == Qt::UserRole+3 ){
        if(column -1 < conservation_1_.size()){
          return QVariant(conservation_1_[column-1]);
        }
          return QVariant(QColor(Qt::transparent));
      }

      if(role == Qt::ForegroundRole){
        if(column -1 < conservation_1_.size()){
          if(conservation_1_[column-1].red()>128){
            return QVariant(QColor(Qt::black));
          }
          else{
            return QVariant(QColor(Qt::white));
          }
        }
        return QVariant(QColor(Qt::transparent));
      }
    }
    else if(this->GetCurrentDisplayMode() == conservation_mode_2){
      if(role == Qt::UserRole+3 ){
        if(column -1 < conservation_2_.size()){
          return QVariant(conservation_2_[column-1]);
        }
        return QVariant(QColor(Qt::transparent));
      }
      if(role == Qt::ForegroundRole){
        if(column -1 < conservation_2_.size()){
            return QVariant(QColor(Qt::black));
        }
        return QVariant(QColor(Qt::transparent));
      }
    }
  }
  return BaseViewObject::GetData(row,column,role);
}

void AlignmentViewObject::SetDisplayMode(const QString& mode)
{
  if(this->display_modes_.contains(mode) && mode != this->GetCurrentDisplayMode()){
    if(mode == conservation_mode_1 || mode == conservation_mode_2){
      for(int i=0 ; i<this->GetRowCount(); i++){
        BaseRow* row = this->GetRow(i);
        row->RemovePainter(seq_secondary_structure_painter);
        row->RemovePainter(align_properties_painter);
        row->InsertPainter(conservation_painter,1);
      }
    }
  }
  SequenceViewObject::SetDisplayMode(mode);
}

const seq::AlignmentHandle& AlignmentViewObject::GetAlignment()
{
  return alignment_;
}

}}
