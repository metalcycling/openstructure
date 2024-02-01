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

#include <ost/gui/python_shell/python_syntax_highlighter.hh>

#include "file_viewer.hh"

#include <QMainWindow>
#include <QSyntaxHighlighter>
#include <QTextEdit>
#include <QVBoxLayout>

namespace ost { namespace gui {

FileViewer::FileViewer(const QString& path, QWidget* parent):
 QWidget(parent)
{
  QLayout* layout = new QVBoxLayout(this);
  layout->setContentsMargins(0,0,0,0);
  text_viewer_ = new QTextEdit(this);
  text_viewer_->setReadOnly(true);

  QFile file(path);
  QTextDocument* doc = text_viewer_->document();
  highlighter_ = new PythonSyntaxHighlighter(doc);
  if (file.open(QFile::ReadOnly | QFile::Text))
     text_viewer_->setText(file.readAll());
  if(doc->blockCount()>1 && path.endsWith(".py")){
    for(int i= 0; i< doc->lineCount();i++){
      QTextBlock block = doc->findBlockByLineNumber(i);
      block.setUserState(BLOCKTYPE_CODE);
    }
    QTextBlock start = doc->findBlockByLineNumber(0);
    QTextBlock end = doc->findBlockByLineNumber(doc->lineCount()-1);
    doc->markContentsDirty(start.position(),end.position()+end.length()-1);
  }
  layout->addWidget(text_viewer_);
  this->setWindowTitle(path);
  this->resize(800, 800);

}

}}

