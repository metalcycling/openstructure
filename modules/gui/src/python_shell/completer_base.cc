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
  Authors: Marco Biasini, Andreas Schenk
 */

#include "completer_base.hh"

#include <QStringList>
#include <QScrollBar>
#include <QAbstractItemView>
#include <QKeyEvent>
#include <QApplication>



namespace ost { namespace gui {



OstCompleterBase::OstCompleterBase(QObject * parent):
 QCompleter(parent)
{
  setCompletionMode(QCompleter::PopupCompletion);
  setCaseSensitivity(Qt::CaseSensitive);
  setModelSorting(QCompleter::CaseSensitivelySortedModel);
  setWrapAround(false);
  installEventFilter(popup());
}

OstCompleterBase::~OstCompleterBase()
{
}

void OstCompleterBase::complete(const QRect & rect,bool inline_completion)
{
  QString common_match=GetCommonMatch();
  if(common_match != completionPrefix() && inline_completion){
    emit activated(common_match);
    return;
  }
  QCompleter::complete(rect);
  popup()->setCurrentIndex(currentIndex());
  unsigned int popup_width=popup()->sizeHintForColumn(0)+ 
                           (popup()->verticalScrollBar()->isVisible() ? 
                            popup()->verticalScrollBar()->width():0);
  popup()->setFixedWidth(popup_width);
  popup()->move(widget()->mapToGlobal(rect.topLeft()+QPoint(0,-popup()->height())));
}

QString OstCompleterBase::GetCommonMatch()
{
  if(completionCount()==0){
    return completionPrefix();
  }
  QStringList matches;
  for(int i=0;i<completionCount();++i){
    setCurrentRow(i);
    matches.append(currentCompletion().remove(0,completionPrefix().length()));
  }
  QString common_match;
  bool break_flag=false;
  for(int i=0;i<matches[0].length() && break_flag==false;++i){
    QChar current_char=matches[0][i];
    for(int j=0;j<completionCount();++j){
      if(i>=matches[j].length() || matches[j][i]!=current_char){
        break_flag=true;
        break;
      }
    }
    if(!break_flag){
      common_match.append(current_char);
    }
  }
  setCurrentRow(0);
  return completionPrefix()+common_match;
}

bool OstCompleterBase::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::KeyPress && popup()->isVisible()) {
    
    QKeyEvent* keyevent = static_cast<QKeyEvent *>(event);
    if(!keyevent->text().isEmpty()){
      if(keyevent->text()[0].isLetterOrNumber()){
        emit recomplete(completionPrefix()+keyevent->text());
        return true;
      } else if(keyevent->key()==Qt::Key_Backspace){
        QString match=completionPrefix();
        match.chop(1);
        emit recomplete(match);
        return true;
      } else if (keyevent->key()!=Qt::Key_Return && 
              keyevent->key()!=Qt::Key_Enter && 
              keyevent->key()!=Qt::Key_Up &&
              keyevent->key()!=Qt::Key_Down) {
        popup()->close();
        QKeyEvent* new_event=new QKeyEvent(keyevent->type(),keyevent->key(),
                                           keyevent->modifiers(),
                                           keyevent->text(), 
                                           keyevent->isAutoRepeat(),
                                           keyevent->count());
        QApplication::postEvent(widget(), new_event);
        return false;
      }
    }
  } 
  return QCompleter::eventFilter(obj, event);
}

}}//ns

