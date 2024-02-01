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


#include "title_row.hh"

#include "background_painter.hh"
#include "tick_painter.hh"

#include "sequence_model.hh"

#include <QMap>
#include <QMapIterator>

#include <QtGui>
namespace ost { namespace gui {

SequenceModel::SequenceModel(QObject *parent)
    : QAbstractTableModel(parent), max_columns(0)
{
  this->beginInsertRows(QModelIndex(),this->rowCount(),this->rowCount());
  BaseViewObject* title = new BaseViewObject(this);
  TitleRow* title_row = new TitleRow(title);
  Painter* p = new BackgroundPainter(title_row);
  title_row->InsertPainter(p);
  p = new TickPainter(title_row);
  title_row->InsertPainter(p);
  title->InsertRow(0,title_row);
  objects_.append(title);
  this->endInsertRows();
}

void SequenceModel::InsertSequence(QString& name, seq::SequenceHandle& seq)
{
  int cols = this->columnCount();
  int new_cols = seq.GetLength();
  this->beginInsertRows(QModelIndex(),this->rowCount(),this->rowCount());
  objects_.append(new SequenceViewObject(seq, name, this));
  if(new_cols > cols){
    this->max_columns = new_cols;
    this->beginInsertColumns(QModelIndex(), cols, new_cols);
    this->endInsertColumns();
  }
  this->endInsertRows();
}

void SequenceModel::InsertChain(QString& name, mol::ChainView& view) {
  int cols = this->columnCount();
  int new_cols = view.GetResidueCount();
  this->beginInsertRows(QModelIndex(),this->rowCount(),this->rowCount());
  objects_.append(new SequenceViewObject(view, name, this));
  if(new_cols > cols){
    this->max_columns = new_cols;
    this->beginInsertColumns(QModelIndex(), cols, new_cols);
    this->endInsertColumns();
  }
  this->endInsertRows();
}

void SequenceModel::InsertAlignment(const seq::AlignmentHandle& alignment){
  int cols = this->columnCount();
  int new_cols = alignment.GetLength();
  this->beginInsertRows(QModelIndex(),this->rowCount(),this->rowCount()+alignment.GetCount()-1);
  objects_.append(new AlignmentViewObject(alignment, this));
  if(new_cols > cols){
    this->max_columns = new_cols;
    this->beginInsertColumns(QModelIndex(), cols, new_cols);
    this->endInsertColumns();
  }
  this->endInsertRows();
}

void SequenceModel::InsertSequences(const QList<QString>& names, 
                                    seq::SequenceList& list){
  this->beginInsertRows(this->index(this->rowCount(),0),
                        this->rowCount(),this->rowCount()+list.GetCount());
  objects_.append(new SequenceViewObject(list, names, this));
  this->endInsertRows();
}

void SequenceModel::InsertGfxEntity(const gfx::EntityP& ent){
  mol::EntityView view=ent->GetView();
  int size = view.GetChainList().size();
  int cols = this->columnCount();
  this->beginInsertRows(QModelIndex(),this->rowCount(),this->rowCount()+size);
  SequenceViewObject* obj = new SequenceViewObject(ent, this);
  objects_.append(obj);
  int new_cols = obj->GetMaxColumnCount();
  if(new_cols > cols){
    this->max_columns = new_cols;
    this->beginInsertColumns(QModelIndex(), cols, new_cols);
    this->endInsertColumns();
  }
  this->endInsertRows();
}

void SequenceModel::RemoveGfxEntity(const gfx::EntityP& entity){
  if(SequenceViewObject* obj = this->GetItem(entity)){
    int index = this->GetGlobalRow(obj,0);
    this->beginRemoveRows(QModelIndex(),index,index+obj->GetRowCount()-1);
    int cols_before = this->columnCount();
    objects_.removeOne(obj);
    this->endRemoveRows();
    int cols = this->GetColumnCount();
    if(cols_before>cols){
      this->max_columns = cols;
      this->beginRemoveColumns(QModelIndex(), cols, cols_before);
      this->endRemoveColumns();
    }
  }
}

void SequenceModel::RemoveAlignment(const seq::AlignmentHandle& alignment){
  if(AlignmentViewObject* obj = this->GetItem(alignment)){
    int index = this->GetGlobalRow(obj,0);
    this->beginRemoveRows(QModelIndex(),index,index+obj->GetRowCount()-1);
    int cols_before = this->columnCount();
    objects_.removeOne(obj);
    this->endRemoveRows();
    int cols = this->GetColumnCount();
    if(cols_before>cols){
      this->max_columns = cols;
      this->beginRemoveColumns(QModelIndex(), cols, cols_before);
      this->endRemoveColumns();
    }
  }
}

int SequenceModel::GetColumnCount() const {
  int cols = 0;
  for(int i = 0; i<objects_.size();i++){
    int ob_cols = objects_[i]->GetMaxColumnCount();
    if(ob_cols>cols){
      cols = ob_cols;
    }
  }
  return cols;
}

SequenceViewObject* SequenceModel::GetItem(const gfx::EntityP& entity)
{
  if(entity != NULL){
    for (int i = 0 ; i< objects_.size(); i++){
      if(SequenceViewObject* seq_view_object = qobject_cast<SequenceViewObject*>(objects_[i])){
        if(entity == seq_view_object->GetGfxObject()){
          return seq_view_object;
        }
      }
    }
  }
  return NULL;
}

AlignmentViewObject* SequenceModel::GetItem(const seq::AlignmentHandle& alignment)
{
  for (int i = 0 ; i< objects_.size(); i++){
    if(AlignmentViewObject* alignment_object = qobject_cast<AlignmentViewObject*>(objects_[i])){
      if(alignment == alignment_object->GetAlignment()){
        return alignment_object;
      }
    }
  }
  return NULL;
}

const PainterList& SequenceModel::GetPainters(const QModelIndex& index) const
{
  QPair<int, BaseViewObject*> pair = this->GetRowWithItem(index);
  if(pair.second){
    return pair.second->GetRow(pair.first)->GetPainters();
  }
  return empty_painter_list_;
}

QPair<int, BaseViewObject*> SequenceModel::GetRowWithItem(int row) const {
  if(!objects_.isEmpty()){
    int rows = 0;
    int i = -1;
    int last_row = 0;
    while (rows <= row && i < objects_.size()){
      i++;
      last_row =objects_[i]->GetRowCount();
      rows += last_row;
    }
    int sub_index = row - (rows-last_row);
    return QPair<int, BaseViewObject*>(sub_index, objects_[i]);
  }
  return QPair<int, BaseViewObject*>(-1, NULL);
}

void SequenceModel::EmitRowChanged(int row)
{
  emit this->dataChanged(this->index(row, 0), 
                         this->index(row, this->columnCount()-1));
}

void SequenceModel::Clear()
{
  // remove everything but the title row. 
  this->beginRemoveRows(QModelIndex(), 1, this->rowCount());
  objects_.erase(objects_.begin()+1, objects_.end());
  this->endRemoveRows();
  this->beginRemoveColumns(QModelIndex(), 0, this->max_columns);
  this->max_columns=0;
  this->endRemoveColumns();
}

QPair<seq::AlignmentHandle, int> SequenceModel::GetAlignmentForRow(int row)
{
  QPair<int, BaseViewObject*> p=this->GetRowWithItem(row);
  if (p.second) {
    AlignmentViewObject* avo=dynamic_cast<AlignmentViewObject*>(p.second);
    if (avo) {
      return QPair<seq::AlignmentHandle, int>(avo->GetAlignment(), p.first);
    }    
  }
  return QPair<seq::AlignmentHandle, int>(seq::AlignmentHandle(), -1);
}


QPair<int, BaseViewObject*> SequenceModel::GetRowWithItem(const QModelIndex& index) const{
  return this->GetRowWithItem(index.row());
}

BaseViewObject* SequenceModel::GetItem(const QModelIndex& index) const
{
  return this->GetRowWithItem(index).second;
}

int SequenceModel::GetGlobalRow(BaseViewObject* obj, int row) const
{
  int glob_row = -1;
  int index = objects_.indexOf(obj);
  if(index >= 0){
    glob_row = 0;
    for(int i=0; i<index; i++){
      glob_row += objects_[i]->GetRowCount();
    }
    return glob_row + row;
  }
  return glob_row;
}

QModelIndexList SequenceModel::GetModelIndexes(gfx::EntityP& entity, 
                                               const mol::EntityView& view)
{
  QModelIndexList list;
  if(BaseViewObject* object = this->GetItem(entity)){
    if(SequenceViewObject* seq_view_object = qobject_cast<SequenceViewObject*>(object)){
      QMap<int, QList<int> > indexes = seq_view_object->GetIndexesForView(view);
      QMapIterator< int, QList<int> > i(indexes);
      while (i.hasNext()) {
        i.next();
        int row = this->GetGlobalRow(seq_view_object, i.key());
        const QList<int>& index_list = i.value();
        for(int i=0; i<index_list.size(); i++){
          list.append(this->index(row,index_list[i]));
        }
      }
    }
  }
  return list;
}

QModelIndexList SequenceModel::GetModelIndexes(const QString& subject, 
                                               const QString& sequence_name)
{
  QModelIndexList list;
  for (int i = 0; i<objects_.size(); i++){
    if(SequenceViewObject* seq_view_object = qobject_cast<SequenceViewObject*>(objects_[i])){
      QMap<int, QList<int> > indexes = seq_view_object->GetIndexesForSubject(subject,sequence_name);
      QMapIterator< int, QList<int> > iter(indexes);
      while (iter.hasNext()) {
        iter.next();
        int row = this->GetGlobalRow(seq_view_object, iter.key());
        const QList<int>& index_list = iter.value();
        for(int j=0; j<index_list.size(); j++){
          list.append(this->index(row,index_list[j]));
        }
      }
    }
  }
  return list;
}

void SequenceModel::SelectionChanged(const QItemSelection& sel, 
                                     const QItemSelection& desel)
{
  QMap<int,QPair<QSet<int>,QSet<int> > > sel_map;
  const QModelIndexList& sel_indexes = sel.indexes();
  for(int i =0; i< sel_indexes.size(); i++){
     sel_map[sel_indexes[i].row()].first.insert(sel_indexes[i].column());
  }

  const QModelIndexList& desel_indexes = desel.indexes();
  for(int i =0; i< desel_indexes.size(); i++){
     sel_map[desel_indexes[i].row()].second.insert(desel_indexes[i].column());
  }

  QMapIterator< int,QPair<QSet<int>,QSet<int> > > i(sel_map);
  while (i.hasNext()) {
    i.next();
    QPair<int, BaseViewObject*> item = this->GetRowWithItem(i.key());
    item.second->SetSelection(item.first,i.value().first, i.value().second);
  }
}

void SequenceModel::DoubleClicked(const QModelIndex& index)
{
  QPair<int, BaseViewObject*> item = this->GetRowWithItem(index);
  if(item.second){
    item.second->DoubleClicked(item.first,index.column());
  }
}

const QStringList& SequenceModel::GetDisplayModes()
{
  display_modes_.clear();
  QMap<QString,int> string_map;
  for (int i = 0; i<objects_.size(); i++){
    const QStringList& list = objects_[i]->GetDisplayModes();
    for(int j=0; j<list.size(); j++){
      if(!string_map.contains(list.at(j))){
        string_map.insert(list.at(j),1);
      }
      else {
        string_map[list.at(j)] = string_map[list.at(j)] + 1;
      }
    }
  }
  bool removed = false;
  QMapIterator<QString, int> i(string_map);
  while (i.hasNext()) {
    i.next();
    if(objects_.size()-1 <= i.value()){
      display_modes_.append(i.key());
    }
    else{
      removed = true;
    }
  }
  if(removed){
    display_modes_.insert(0," ");
  }
  return display_modes_;
}

const QStringList& SequenceModel::GetDisplayModes(const gfx::EntityP& entity)
{
  BaseViewObject* item = this->GetItem(entity);
  if(item){
    return item->GetDisplayModes();
  }
  else{
    return empty_string_list_;
  }
}

const QStringList& SequenceModel::GetDisplayModes(const seq::AlignmentHandle& alignment)
{
  BaseViewObject* item = this->GetItem(alignment);
  if(item){
    return item->GetDisplayModes();
  }
  else{
    return empty_string_list_;
  }
}

const QString& SequenceModel::GetCurrentDisplayMode()
{
  current_display_mode_.clear();
  for (int i = 0; i<objects_.size(); i++){
    const QString& mode = objects_[i]->GetCurrentDisplayMode();
    if(current_display_mode_.isEmpty()){
      current_display_mode_ = mode;
    }
    else if(current_display_mode_ != mode){
      current_display_mode_ = " ";
      break;
    }
  }
  return current_display_mode_;
}

const QString& SequenceModel::GetCurrentDisplayMode(const gfx::EntityP& entity)
{
  BaseViewObject* item = this->GetItem(entity);
  if(item){
    return item->GetCurrentDisplayMode();
  }
  else{
    return empty_string_;
  }
}

const QString& SequenceModel::GetCurrentDisplayMode(const seq::AlignmentHandle& alignment)
{
  BaseViewObject* item = this->GetItem(alignment);
  if(item){
    return item->GetCurrentDisplayMode();
  }
  else{
    return empty_string_;
  }
}

void SequenceModel::SetDisplayMode(const QString& mode)
{
  for (int i = 0; i<objects_.size(); i++){
    objects_[i]->SetDisplayMode(mode);
  }
}

void SequenceModel::SetDisplayMode(const gfx::EntityP& entity, const QString& mode)
{
  BaseViewObject* item = this->GetItem(entity);
  if(item){
    return item->SetDisplayMode(mode);
  }
}

void SequenceModel::SetDisplayMode(const seq::AlignmentHandle& alignment, const QString& mode)
{
  BaseViewObject* item = this->GetItem(alignment);
  if(item){
    return item->SetDisplayMode(mode);
  }
}

int SequenceModel::rowCount(const QModelIndex& parent) const
{
  int rows = 0;
  for (int i = 0; i<objects_.size(); i++){
    rows += objects_[i]->GetRowCount();
  }
  return rows;
}

int SequenceModel::columnCount(const QModelIndex& parent) const
{
  return max_columns;
}

QVariant SequenceModel::data(const QModelIndex& index, int role) const
{
  QPair<int, BaseViewObject*> item = this->GetRowWithItem(index);
  if(!item.second) return QVariant();
  QVariant data = item.second->GetData(item.first,index.column(),role);
  return data;
}

QVariant SequenceModel::headerData(int section, Qt::Orientation orientation,
                                    int role) const
{
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole) {
    return QVariant("");
  }
  return QVariant();
}

Qt::ItemFlags SequenceModel::flags(const QModelIndex& index) const
{
  QPair<int, BaseViewObject*> item = GetRowWithItem(index);
  if(item.second){
    return item.second->Flags(item.first, index.column());
  }
  return QAbstractItemModel::flags(index);
}

void SequenceModel::ZoomIn()
{
  for (int i = 0; i<objects_.size(); i++){
    objects_[i]->ZoomIn();
  }
}

void SequenceModel::ZoomOut()
{
  for (int i = 0; i<objects_.size(); i++){
    objects_[i]->ZoomOut();
  }
}

}}
