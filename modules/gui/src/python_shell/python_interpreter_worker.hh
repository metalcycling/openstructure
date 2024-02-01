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


#ifndef PYTHON_INTERPRETER_WORKER_HH
#define PYTHON_INTERPRETER_WORKER_HH

// workaround for QTBUG-22829: https://bugreports.qt-project.org/browse/QTBUG-22829
#ifndef Q_MOC_RUN
  // https://stackoverflow.com/questions/23068700/embedding-python3-in-qt-5
  #pragma push_macro("slots")
  #undef slots
  #include <boost/python.hpp>
  #include <boost/shared_ptr.hpp>
  #pragma pop_macro("slots")
#endif
#include "output_redirector.hh"
#include <csignal>
#include <utility>
#include <QObject>
#include <QQueue>

namespace ost { namespace gui {
namespace bp = boost::python;

class PythonInterpreterWorker: public QObject
{
Q_OBJECT
public:
  PythonInterpreterWorker();
  ~PythonInterpreterWorker();
  unsigned int AddCommand(const QString& command);

signals:
  void Finished(unsigned int id, bool error_state);
  void Output(unsigned int id,const QString& output);
  void ErrorOutput(unsigned int id,const QString& output);

public slots:
  void Wake();

protected slots:
  void handle_redirector_output(const QString& output);
  void handle_redirector_error(const QString& output);

protected:
  bool is_simple_expression(const QString& expr);
  void run_command_(std::pair<unsigned int,QString> pair);
  QQueue<std::pair<unsigned int,QString> > exec_queue_;
  unsigned int command_id_;
  boost::shared_ptr<OutputRedirector>  output_redirector_;
  boost::shared_ptr<OutputRedirector>  error_redirector_;
#ifndef _MSC_VER
  struct sigaction sig_act_;
#endif
  bp::object parse_expr_cmd_;
  bp::object repr_;
  bp::dict main_namespace_;
  unsigned int current_id_;
  bool         awake_;
};

}} //ns

#endif // PYTHON_INTERPRETER_WORKER_HH
