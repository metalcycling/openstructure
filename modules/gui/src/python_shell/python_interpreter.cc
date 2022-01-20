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
  python interpreter

  Authors: Andreas Schenk, Ansgar Philippsen, Marco Biasini
*/

#include "python_interpreter.hh"

#include <iostream>
#include <csignal>
#include <string>
#include <fstream>

#include <boost/filesystem/convenience.hpp>
#include <ost/log.hh>

#include <QFile>
#include <QDebug>
#include <QFile>
#include <QStringList>

namespace ost { namespace gui {

PythonInterpreter::PythonInterpreter():
  QObject(),
  main_module_(),
  main_namespace_(),
  running_(false),
  compile_command_(),
  worker_()
{
  connect(this,SIGNAL(WakeWorker()),&worker_,SLOT(Wake()),Qt::QueuedConnection);
  connect(&worker_,SIGNAL(Output(unsigned int,const QString&)),this,SIGNAL(Output(unsigned int,const QString&)));
  connect(&worker_,SIGNAL(ErrorOutput(unsigned int,const QString&)),this,SIGNAL(ErrorOutput(unsigned int,const QString&)));
  connect(&worker_,SIGNAL(Finished(unsigned int, bool)),this,SIGNAL(Finished(unsigned int, bool)));
  main_module_ = bp::import("__main__");
  main_namespace_ = bp::extract<bp::dict>(main_module_.attr("__dict__"));
  main_namespace_["os"]=bp::import("os");
  main_namespace_["builtins"]=bp::import("builtins");
  main_namespace_["keyword"]=bp::import("keyword");
  bp::object code=bp::import("code");
  compile_command_=code.attr("compile_command");
}

PythonInterpreter::~PythonInterpreter()
{
}

PythonInterpreter& PythonInterpreter::Instance()
{
  if(!Py_IsInitialized()){
    Py_InitializeEx(1);
  }
  static PythonInterpreter instance;
  return instance;
}



void PythonInterpreter::Stop()
{
  running_=false;
}

void PythonInterpreter::Start()
{
  running_=true;
    emit WakeWorker();
}

void PythonInterpreter::StopScript()
{
  PyErr_SetInterrupt();
}


unsigned int PythonInterpreter::RunScript(const QString& fname)
{
  QFile script(fname);
  if (!script.open(QIODevice::ReadOnly | QIODevice::Text)){
    LOG_ERROR("could not open " << fname.toStdString());
    return RunCommand("");
  }
  QString command=script.readAll();
  command.append('\n');
  command.append('\n');
  return RunCommand(command);
}

unsigned int PythonInterpreter::RunCommand(const QString& command)
{
  unsigned int id=worker_.AddCommand(command);
  if(running_){
    emit WakeWorker();
  }
  return id;
}

CodeBlockStatus PythonInterpreter::GetCodeBlockStatus(const QString& command)
{
  // move to pure c++ to avoid clash with running worker
  // move to worker class
  CodeBlockStatus status=CODE_BLOCK_COMPLETE;
/*#ifndef _MSC_VER
  SigIntHandler sh(sig_act_);
#endif*/
  QStringList lines=command.split("\n");  
  String cmd;  
  for (int i=0; i<lines.size(); ++i) {
    if (status==CODE_BLOCK_ERROR) {
      break;
    }
    bool complete=false;
    if (!(! lines[i].isEmpty() && (lines[i][0]=='\t' || lines[i][0]==' '))) {
      cmd.clear();
    }
    while (!complete && i<lines.size()) {
      cmd+=lines[i].toStdString();
      if (i<lines.size()-1) {
        cmd+="\n";
      }
      try {
        bp::object result=compile_command_(bp::str(cmd));
        try {
          String doc=bp::extract<String>(result.attr("__doc__"));
          complete=true;
          status=CODE_BLOCK_COMPLETE;
        } catch(bp::error_already_set&) {
          status=CODE_BLOCK_INCOMPLETE;
          ++i;          
          PyErr_Clear();
        }
      } catch(bp::error_already_set&) {
        status=CODE_BLOCK_ERROR;
        complete=true;
        PyErr_Clear();
      }
    }
  }
   return status;
}


void PythonInterpreter::AppendModulePath(const QString& entry)
{
  RunCommand(QString("sys.path.insert(0, os.path.expanduser(os.path.expandvars('")+ entry + QString("')))\n"));
}

void PythonInterpreter::AppendCommandlineArgument(const QString& arg)
{
  QString init_code("");
  static bool sys_argv_present=false;
  if(!sys_argv_present){
    sys_argv_present=true;
    init_code += "\nsys.argv=[]";
  }
  init_code += "\nsys.argv.append('" ;
  init_code += arg;
  init_code += "')\n";
  RunCommand(init_code);
}

bp::object PythonInterpreter::GetMainModule() const
{
  return main_module_;
}
bp::dict PythonInterpreter::GetMainNamespace() const
{
  bp::dict result=const_cast<bp::dict&>(main_namespace_).copy();
  bp::dict builtins_dict=bp::extract<bp::dict>(main_namespace_["__builtin__"].attr("__dict__"));
  result.update(builtins_dict);
  bp::list kwlist=bp::extract<bp::list>(main_namespace_["keyword"].attr("kwlist"));
  unsigned int kwlist_length=bp::len(kwlist);
  for(unsigned int i=0;i<kwlist_length;++i){
    result[kwlist[i]]=bp::str();
  }
  return result;
}


}} // ns

