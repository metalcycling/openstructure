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


#ifndef OST_GUI_PYTHON_SHELL_WIDGET_HH
#define OST_GUI_PYTHON_SHELL_WIDGET_HH


#include "python_interpreter.hh"
#include "shell_history.hh"
#include  "python_shell_fw.hh"
#include "python_syntax_highlighter.hh"
#include <ost/gui/module_config.hh>

#include <QPlainTextEdit>
#include <QHash>
namespace ost { namespace gui {

// fw decl
class Gutter;
class PythonCompleter;
class PathCompleter;
class State;
class StateMachine;

typedef std::vector<GutterBlock> GutterBlockList;

class DLLEXPORT_OST_GUI PythonShellWidget : public QPlainTextEdit {

  typedef QHash<unsigned int,QTextBlock> OutputBlockList;
  Q_OBJECT

public:
  PythonShellWidget(QWidget* parent=NULL);
  GutterBlockList GetGutterBlocks(const QRect& rect);
  void SetTabWidth(int width);
  int GetTabWidth() const;
  QString GetCommand();
  QTextBlock GetEditStartBlock();

  int tab_width_;

signals:
  void Execute(const QString& command);
  void RequestCompletion(const QRect& rect,bool inline_completion);
  void SetCompletionPrefix(const QString& prefix);
  void RequestPathCompletion(const QRect& rect,bool inline_completion);
  void SetPathCompletionPrefix(const QString& prefix);

public slots:
  void InsertCompletion(const QString& completion);
  void InsertPathCompletion(const QString& completion);
  void AppendOutput(unsigned int id,const QString& output);
  void AppendError(unsigned int id,const QString& output);
  void OutputFinished(unsigned int id, bool error);
  void Complete(bool inline_completion=true);
  void Recomplete(const QString& completion);
  // slots for state machine
  void OnSingleLineStateEntered();
  void OnMultiLineActiveStateEntered();
  void OnMultiLineInactiveStateEntered();
  void OnHistoryUpStateEntered();
  void OnHistoryDownStateEntered();
  void OnExecuteStateEntered();
  void NewLineAtEnd();
  void NewLineAtCursor();
  void OnReadonlyEntered();
  void OnReadwriteEntered();
  void OnMixedToReadwrite();
protected slots:
  void handle_completion_();
  void handle_output_toggle_();
  void handle_wrap_to_function_();
  void handle_select_all_();
  void handle_select_all_rw_();
  void handle_delete_();
  void handle_backspace_();
  void handle_clear_all_();

protected:
  virtual void resizeEvent(QResizeEvent* event);
  virtual void showEvent(QShowEvent* event);
  virtual void insertFromMimeData( const QMimeData * source );
  void set_output_visible_(bool flag=true);
  void setup_readonly_state_machine_();
  void setup_state_machine_();
  void wrap_into_function_(const QString& command);
  void set_command_(const QString& command);
  void set_block_edit_mode_(BlockEditMode flag);
  BlockEditMode get_block_edit_mode_();
  void set_block_type_(const QTextBlock& start,const QTextBlock& end, BlockType type);

  PythonSyntaxHighlighter highlighter_;
  PythonCompleter* completer_;
  PathCompleter* path_completer_;
  BlockEditMode block_edit_mode_;
  ShellHistory history_;
  Gutter* gutter_;
  bool output_visible_;
  int completion_start_;
  int completion_end_;
  QTextBlock block_edit_start_;
  OutputBlockList output_blocks_;
  StateMachine* machine_;
  StateMachine* readonly_machine_;
  State* readwrite_state_;
  State* multiline_active_state_;
};

}}//ns

#endif
