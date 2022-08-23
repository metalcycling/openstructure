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


#include "python_interpreter_worker.hh"
#include "python_interpreter.hh"


namespace ost { namespace gui {

    namespace {
#ifndef _MSC_VER
  class SigIntHandler {
  public:
    SigIntHandler(const struct sigaction& act) {
      sigaction(SIGINT,&act,&old_);
    }
    ~SigIntHandler() {
      sigaction(SIGINT,&old_,0);
    }
  private:
    struct sigaction old_;
  };
#endif
}



PythonInterpreterWorker::PythonInterpreterWorker():
  exec_queue_(),
  command_id_(0),
  output_redirector_(new OutputRedirector),
  error_redirector_(new OutputRedirector),
  #ifndef _MSC_VER
  sig_act_(),
  #endif
  parse_expr_cmd_(),
  repr_(),
  main_namespace_(),
  current_id_(),
  awake_(false)
{
  parse_expr_cmd_=bp::import("__main__").attr("__builtins__").attr("compile");
  main_namespace_ = bp::extract<bp::dict>(bp::import("__main__").attr("__dict__"));
  repr_=bp::import("__main__").attr("__builtins__").attr("repr");
  #ifndef _MSC_VER
  sigaction(SIGINT,0,&sig_act_);
  #endif
  main_namespace_["OutputRedirector"] = bp::class_<OutputRedirector,boost::shared_ptr<OutputRedirector>, boost::noncopyable >("OutputRedirector", bp::init<>())
                                          .def("write", &OutputRedirector::Write);
  connect(output_redirector_.get(),SIGNAL(OnOutput(const QString&)),this,SLOT(handle_redirector_output(const QString&)));
  connect(error_redirector_.get(),SIGNAL(OnOutput(const QString&)),this,SLOT(handle_redirector_error(const QString&)));
  main_namespace_["sys"]=bp::import("sys");
  main_namespace_["sys"].attr("stderr") = error_redirector_;
  main_namespace_["sys"].attr("stdout") = output_redirector_;
}

PythonInterpreterWorker::~PythonInterpreterWorker()
{
  // we have to manually run the exit functions because we cannot use 
  // Py_Finalize due to some problems in boost python.
  run_command_(std::pair<unsigned int, QString>(0,
               "import atexit\n"
               "atexit._run_exitfuncs()\n"));
  //Py_Finalize();
}

void PythonInterpreterWorker::Wake()
{
  if (awake_) return;

  awake_=true;
  while (!exec_queue_.isEmpty()){
    std::pair<unsigned int, QString> pair=exec_queue_.dequeue();
    run_command_(pair);
  }
  awake_=false;
}

unsigned int PythonInterpreterWorker::AddCommand(const QString& command)
{
  exec_queue_.enqueue(std::pair<unsigned int, QString>(++command_id_,command)); //wrap around in command_id_ intended
  return command_id_;
}

void PythonInterpreterWorker::run_command_(std::pair<unsigned int,QString> pair)
{
  #ifndef _MSC_VER
  SigIntHandler sh(sig_act_);
  #endif
  try{
    // todo handle quit and exit without parantheses
    current_id_=pair.first;
    QString command=pair.second;
    if(is_simple_expression(command)){
      bp::object result=bp::eval(command.toStdString().c_str(),
                                 main_namespace_, main_namespace_);
      String rstring=bp::extract<String>(repr_(result));
      if(rstring!="None"){
        handle_redirector_output(QString::fromStdString(rstring)+"\n");
      }
    } else {
      bp::exec(command.toStdString().c_str(), main_namespace_, main_namespace_);
    }
    output_redirector_->Flush();
    error_redirector_->Flush();
    emit Finished(pair.first,true);
    return;
  }catch(bp::error_already_set&){
    if(PyErr_ExceptionMatches(PyExc_SystemExit)){
      PyErr_Clear();
      //emit Exit();
    } else{
      PyErr_Print();
    }
  }
  output_redirector_->Flush();
  error_redirector_->Flush();
  emit Finished(pair.first,false);
  return;
}

bool PythonInterpreterWorker::is_simple_expression(const QString& expr)
{
  try {
    parse_expr_cmd_(bp::str(expr.toStdString()), "<string>", "eval");
    return true;
  } catch(bp::error_already_set&) {
    PyErr_Clear();
    return false;
  }
}

void PythonInterpreterWorker::handle_redirector_output(const QString& output)
{
  emit Output(current_id_,output);
}

void PythonInterpreterWorker::handle_redirector_error(const QString& output)
{
  emit ErrorOutput(current_id_,output);
}

}} //ns

