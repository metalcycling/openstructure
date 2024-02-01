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


#include "sequence_delegate.hh"
#include <QtGui>

namespace ost { namespace gui {

SequenceDelegate::SequenceDelegate(SequenceModel* seq_model, QObject *parent)
    : QItemDelegate(parent), seq_model_(seq_model)
{}

void SequenceDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option,
                         const QModelIndex &index) const
{
  const PainterList& painters = seq_model_->GetPainters(index);
  for(int i=0; i < painters.size(); i++){
    painters[i]->parent();
    if(painters[i]){
      painters[i]->Paint(painter, option, index);
    }
  }
}

}}
