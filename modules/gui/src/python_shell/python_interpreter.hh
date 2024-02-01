//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  python interpreter

  Authors: Andreas Schenk, Ansgar Philippsen, Marco Biasini
*/

#ifndef PYTHON_INTERPRETER_HH
#define PYTHON_INTERPRETER_HH


#include "python_interpreter_worker.hh"
#include "output_redirector.hh"
#include "main_thread_runner.hh"
#include <ost/gui/module_config.hh>

#include <QQueue>
#include <QMetaType>

namespace ost { namespace gui {
namespace bp = boost::python;

enum InterpreterStatus{
  STATUS_OK,
  STATUS_ERROR
};


typedef enum {
   CODE_BLOCK_COMPLETE=1,
   CODE_BLOCK_ERROR=2,
   CODE_BLOCK_INCOMPLETE=4
} CodeBlockStatus;

class DLLEXPORT_OST_GUI PythonInterpreter: public QObject
{
  Q_OBJECT
public:
  static PythonInterpreter& Instance();

  ~PythonInterpreter();
  bool IsRunning();
  void AppendModulePath(const QString& entry);
  
  void AppendCommandlineArgument(const QString& arg);
  
  bp::dict GetMainNamespace() const;
  bp::object GetMainModule() const;
  /// \brief Get code block status
  /// 
  /// Determines whether the command contains errors, is incomplete or ready
  /// to be run.
  CodeBlockStatus GetCodeBlockStatus(const QString& command);

  
  /// \brief stop execution of Python commands
  /// \sa Start
  void Stop();
  /// \brief start/continue execution of Python commands
  /// 
  /// \sa Stop
  void Start();
  /// \brief stop script execution
  void StopScript();
  
  

public slots:
  /// \brief execute python command
  unsigned int RunScript(const QString& script);
  unsigned int RunCommand(const QString& command);
signals:
  void Output(unsigned int id,const QString& output);
  void ErrorOutput(unsigned int id,const QString& output);
  void Finished(unsigned int id, bool error_state);
  void Exit();
  void WakeWorker();
protected:
  PythonInterpreter();
  bp::object main_module_;
  bp::dict main_namespace_;
  bool running_;
  bp::object compile_command_;
  PythonInterpreterWorker worker_;
};

}} // ns

#endif
