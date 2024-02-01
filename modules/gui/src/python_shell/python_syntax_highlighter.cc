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


#include "python_shell.hh"
#include "string_literal_positions.hh"
#include "python_syntax_highlighter.hh"
#include "python_context_parser.hh"

#include <iostream>

#include <QDebug>
#include <QTextCursor>

namespace ost { namespace gui {    

PythonSyntaxHighlighter::PythonSyntaxHighlighter(QTextDocument* parent):
  QSyntaxHighlighter(parent),
  formats_(),
  block_formats_()
{
  QTextCharFormat format;
  formats_[PythonToken::IDENTIFIER]=format;  
  formats_[PythonToken::END]=format;  
  format.setForeground(QBrush(QColor(150,0,0)));
  formats_[PythonToken::OPERATOR]=format;  
  format.setForeground(QBrush(QColor(100,100,100)));  
  formats_[PythonToken::COMMENT]=format;
  format.setForeground(QBrush(QColor(0,0,150)));
  formats_[PythonToken::GROUPING]=format;  
  format.setForeground(QBrush(QColor(35,107,142)));
  formats_[PythonToken::KEYWORD]=format;  
  format.setForeground(QBrush(QColor(0,34,100)));
  formats_[PythonToken::NUMBER]=format;    
  format.setForeground(QBrush(QColor(22,111,22)));
  format.setProperty(QTextFormat::UserProperty,QVariant(true));
  formats_[PythonToken::STRING_LITERAL]=format;  
  format.setForeground(QBrush(QColor(66,111,66)));
  formats_[PythonToken::STRING_DELIM]=format;
  QTextBlockFormat blockformat;
  blockformat.setBackground(QColor(255,245,235));
  block_formats_[BLOCKTYPE_ERROR]=blockformat;
  blockformat.setBackground(QColor(255,255,235));
  block_formats_[BLOCKTYPE_BLOCKEDIT]=blockformat;
  blockformat.setBackground(QColor(255,255,255));
  block_formats_[BLOCKTYPE_CODE]=blockformat;
  block_formats_[BLOCKTYPE_ACTIVE]=blockformat;
  blockformat.setForeground(QColor(0,32,64));
  blockformat.setBackground(QColor(245,245,255));
  block_formats_[BLOCKTYPE_OUTPUT]=blockformat;
  block_formats_[BLOCKTYPE_OUTPUT]=blockformat;

}

void PythonSyntaxHighlighter::highlightBlock(const QString& text_block) { 
  QTextCursor cursor(currentBlock());
  int base_state=-1;  
  if (currentBlockState() & BLOCKTYPE_OUTPUT) {
    base_state=BLOCKTYPE_OUTPUT;
  } else if (currentBlockState() & BLOCKTYPE_ERROR) {
    base_state=BLOCKTYPE_ERROR;
  } else if (currentBlockState() & BLOCKTYPE_ACTIVE) {
    base_state=BLOCKTYPE_ACTIVE;
  } else if (currentBlockState() & BLOCKTYPE_BLOCKEDIT) {
    base_state=BLOCKTYPE_BLOCKEDIT;
  } else if (currentBlockState() & BLOCKTYPE_CODE) {
    base_state=BLOCKTYPE_CODE;
  }
#if QT_VERSION < 0x040700
  cursor.setBlockFormat(block_formats_[static_cast<BlockType>(base_state)]);
#endif
  if (currentBlockState() & BLOCKTYPE_BLOCKEDIT || 
      currentBlockState() & BLOCKTYPE_ACTIVE || 
      currentBlockState() & BLOCKTYPE_CODE) {
    StringLiteralPositions* positions=new StringLiteralPositions();
    int blockposition=currentBlock().position();
    int string_state=0;
    int bs=previousBlockState()==-1 ? 0 : previousBlockState();
    if (bs & BLOCKTYPE_MULTILINE_SQ) {
      string_state=2;
    } else if (bs & BLOCKTYPE_MULTILINE_DQ) {
      string_state=1;
    }
    PythonTokenizer pt(text_block,string_state);
    PythonToken t;
    PythonToken last;
    int type_before_last=PythonToken::END;
    // Loop over all the tokens identified and color them appropriately. 
    while((t = pt.NextToken()).GetType() != PythonToken::END) {
      setFormat(t.GetRange().location,t.GetRange().length,formats_[t.GetType()]);
      if(t.GetType()==PythonToken::STRING_LITERAL){
        positions->Append(blockposition+t.GetRange().location,
                          blockposition+t.GetRange().location+
                          t.GetRange().length);
      }
      type_before_last=last.GetType();
      last=t;
    }
    if(last.GetType()==PythonToken::STRING_DELIM && 
       type_before_last!=PythonToken::STRING_LITERAL && 
       type_before_last!=PythonToken::STRING_DELIM){
      positions->Append(blockposition+last.GetRange().location+1,
                        blockposition+last.GetRange().location+1);
    }
    int new_state=base_state;
    if (pt.InString()) {
      if (pt.GetDelim()=="'''") {
        new_state|=BLOCKTYPE_MULTILINE_SQ;
      } else if (pt.GetDelim()=="\"\"\""){
        new_state|=BLOCKTYPE_MULTILINE_DQ;
      }
    }
    setCurrentBlockState(new_state);
    currentBlock().setUserData(positions);
  }
}

}} // ost::gui
