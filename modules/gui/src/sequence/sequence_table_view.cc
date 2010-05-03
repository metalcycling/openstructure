//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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

#include <QHeaderView>
#include <QScrollBar>
#include <QTableWidgetItem>

#include <iostream>
#include "sequence_table_view.hh"

namespace ost { namespace gui {


SequenceTableView::SequenceTableView(QAbstractItemModel * model)
{
  this->setModel(model);
  column_not_move_ = new QTableView(this);

  column_not_move_->setModel(this->model());
  column_not_move_->setFocusPolicy(Qt::NoFocus);
  column_not_move_->verticalHeader()->hide();
  column_not_move_->horizontalHeader()->setResizeMode(QHeaderView::Fixed);

  this->viewport()->stackUnder(column_not_move_);

  column_not_move_->setSelectionModel(this->selectionModel());
  for(int col=1; col<this->model()->columnCount(); col++){
   column_not_move_->setColumnHidden(col, true);
  }

  column_not_move_->setColumnWidth(0, this->columnWidth(0) );

  column_not_move_->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  column_not_move_->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
  column_not_move_->show();
  column_not_move_->setStyleSheet("QTableView { border: 0px;"
                                 "background-color: #8EDE21;"
                                 "selection-background-color: #999}"
                                 "QTableView::item{ border: none;"
                                 "padding: 0px; border-width: 0px; margin: 0px;}");
  column_not_move_->setShowGrid(false);

  column_not_move_->horizontalHeader()->setStyleSheet(
      "QHeaderView::section {"
         "padding-bottom: 0px;"
         "padding-top: 0px;"
         "padding-left: 0px;"
         "padding-right: 0px;"
         "margin: 0px;"
      "}"
    );

  this->setStyleSheet("QTableView {"
      "show-decoration-selected: 1;"
  "}"
  "QTableView::item {"
       "border: 1px solid #d9d9d9;"
      "border-left-color: transparent;"
      "border-right-color: transparent;"
  "}"
  "QTableView::item:hover {"
      "background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #e7effd, stop: 1 #cbdaf1);"
      "border: 1px solid #bfcde4;"
  "}"
  "QTableView::item:selected {"
      "border: 1px solid #567dbc;"
  "}"
  "QTableView::item:selected:active{"
      "background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #6ea1f1, stop: 1 #567dbc);"
  "}"
  "QTableView::item:selected:!active {"
      "background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, stop: 0 #6b9be8, stop: 1 #577fbf);"
  "}");
  this->setShowGrid(false);

  this->updateNotMoveColumn();

  this->setHorizontalScrollMode(ScrollPerPixel);
  this->setVerticalScrollMode(ScrollPerPixel);
  column_not_move_->setVerticalScrollMode(ScrollPerPixel);

  connect(this->horizontalHeader(),SIGNAL(sectionResized(int,int,int)), this, SLOT(ResizeWidth(int,int,int)));
  connect(this->verticalHeader(),SIGNAL(sectionResized(int,int,int)), this, SLOT(ResizeHeight(int,int,int)));

  connect(column_not_move_->verticalScrollBar(), SIGNAL(valueChanged(int)), this->verticalScrollBar(), SLOT(setValue(int)));
  connect(verticalScrollBar(), SIGNAL(valueChanged(int)), column_not_move_->verticalScrollBar(), SLOT(setValue(int)));
 }

void SequenceTableView::ResizeWidth(int index, int, int size)
{
  if(index == 0){
    column_not_move_->setColumnWidth(0,size);
    updateNotMoveColumn();
  }
}

void SequenceTableView::ResizeHeight(int index, int, int size)
{
  column_not_move_->setRowHeight(index, size);
}

void SequenceTableView::resizeEvent(QResizeEvent * event)
{
  QTableView::resizeEvent(event);
  updateNotMoveColumn();
}

QModelIndex SequenceTableView::moveCursor(CursorAction action, Qt::KeyboardModifiers modifiers)
{
  QModelIndex current = QTableView::moveCursor(action, modifiers);

  if(action == MoveLeft && current.column()>0
     && this->visualRect(current).topLeft().x() < column_not_move_->columnWidth(0) ){
   const int new_value = horizontalScrollBar()->value() + this->visualRect(current).topLeft().x() - column_not_move_->columnWidth(0);
   horizontalScrollBar()->setValue(new_value);
  }
  return current;
}

void SequenceTableView::scrollTo(const QModelIndex & index, ScrollHint hint){
  if(index.column()>0){
    QTableView::scrollTo(index, hint);
  }
}

void SequenceTableView::updateNotMoveColumn()
{
  int x = this->verticalHeader()->width()+this->frameWidth();
  int y = this->frameWidth();
  int w = this->columnWidth(0);
  int h = this->viewport()->height()+this->horizontalHeader()->height();
  column_not_move_->setGeometry(x,y,w,h);
}

void SequenceTableView::columnCountChanged(const QModelIndex& index, int old_count, int new_count){
  if(old_count > 1 && old_count <= new_count){
    for(int col=old_count-1; col<this->model()->columnCount(); col++){
      column_not_move_->setColumnHidden(col, true);
    }
  }
}

void SequenceTableView::resizeColumnsToContents(){
  QTableView::resizeColumnsToContents();
  column_not_move_->setColumnWidth(0,this->columnWidth(0));
  updateNotMoveColumn();
}

SequenceTableView::~SequenceTableView(){}

}}
