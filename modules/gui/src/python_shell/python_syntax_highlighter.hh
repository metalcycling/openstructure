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
#ifndef OST_GUI_PYTHON_SYNTAX_HIGHLIGHTER_HH
#define OST_GUI_PYTHON_SYNTAX_HIGHLIGHTER_HH
/*
  Authors: Marco Biasini, Andreas Schenk
 */


#include "python_shell_fw.hh"
#include "python_tokenizer.hh"

#include <QSyntaxHighlighter>
#include <QTextDocument>
namespace ost { namespace gui {

/** 
 * encodes the syntax rules of the python programming language
 */
class DLLEXPORT_OST_GUI PythonSyntaxHighlighter : public QSyntaxHighlighter {
public:
  PythonSyntaxHighlighter(QTextDocument* parent);
  
  /// \brief hightlight the given block of text. Called by QTextDocument upon 
  ///     need
  virtual void highlightBlock(const QString& text_block);
  

private:
  std::map<PythonToken::Type,QTextCharFormat>  formats_;
  std::map<BlockType,QTextBlockFormat> block_formats_;
  
};
 
}}

#endif //OST_GUI_PYTHON_GRAMMAR
