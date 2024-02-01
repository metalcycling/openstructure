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
#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION<103400
#include <boost/python/detail/api_placeholder.hpp>
#endif

using namespace boost::python;

#include <ost/log.hh>
#include <ost/log_sink.hh>
#include <ost/dyn_cast.hh>

using namespace ost;

struct PyLogSink: public LogSink {

};

struct WrappedLogSink : public PyLogSink, public wrapper<PyLogSink> {
  WrappedLogSink(PyObject* self): self_(self)
  { }
  virtual void LogMessage(const String& message , int severity) 
  {
    call_method<void>(self_, "LogMessage", message, severity);
  }
  
  void LogMessageDefault(const String& message, int severity)
  {
  }
  PyObject* self_;
};



typedef boost::shared_ptr<WrappedLogSink> WrappedLogSinkPtr;

void push_verb(int n) 
{
  Logger::Instance().PushVerbosityLevel(n);
}

void pop_verb() 
{
  Logger::Instance().PopVerbosityLevel();
}

int get_verb()
{
  return Logger::Instance().GetVerbosityLevel();
}


void push_log_sink(LogSinkPtr sink)
{
  Logger::Instance().PushSink(sink);
}

void pop_log_sink()
{
  Logger::Instance().PopSink();
}

LogSinkPtr get_log_sink()
{
  return Logger::Instance().GetCurrentSink();
}

String args_to_string(tuple args, dict kwargs)
{
  std::stringstream ss;
  bool empty = true;
  for (size_t i = 0, l = len(args); i < l; ++i) {
    if (!empty) {
      ss << " ";
    }
    empty = false;
    String string_val;
    extract<String> tst(args[i]);
    if (tst.check()) {
      string_val = tst();
    } else {
      string_val = extract<String>(str(args[i])); // use boost python str
    }
    ss << string_val;
  }
  return ss.str();
}

object log_error(tuple args, dict kwargs) 
{
  LOG_ERROR(args_to_string(args, kwargs));
  return object();
}
object log_warning(tuple args, dict kwargs) 
{
  LOG_WARNING(args_to_string(args, kwargs));
  return object();
}
object log_script(tuple args, dict kwargs) 
{
  LOG_SCRIPT(args_to_string(args, kwargs));
  return object();  
}
object log_info(tuple args, dict kwargs) 
{
  LOG_INFO(args_to_string(args, kwargs));
  return object();  
}
object log_verbose(tuple args, dict kwargs) 
{
  LOG_VERBOSE(args_to_string(args, kwargs));
  return object();  
}
object log_debug(tuple args, dict kwargs) 
{
  LOG_DEBUG(args_to_string(args, kwargs));
  return object();  
}
object log_trace(tuple args, dict kwargs) 
{
  LOG_TRACE(args_to_string(args, kwargs));
  return object();  
}

void reset_sinks()
{
  Logger::Instance().ResetSinks();
}

void export_Logger()
{
  class_<LogSink, LogSinkPtr, boost::noncopyable>("_LogSink", no_init)
    .def("LogMessage", &LogSink::LogMessage)
  ;
  class_<PyLogSink, WrappedLogSinkPtr, bases<LogSink>,
         boost::noncopyable>("LogSink")
    .def("LogMessage", &WrappedLogSink::LogMessageDefault)
  ;
  class_<MultiLogSink, MultiLogSinkPtr, bases<LogSink>, 
         boost::noncopyable >("MultiLogSink", init<>())
    .def("AddSink",&MultiLogSink::AddSink)
    .def("RemoveSink",&MultiLogSink::RemoveSink)
    .def("LogMessage", &MultiLogSink::LogMessage)
  ;

  class_<FileLogSink, FileLogSinkPtr, bases<LogSink>, 
         boost::noncopyable >("FileLogSink", init<const String&>())
    .def("LogMessage", &FileLogSink::LogMessage)
  ;

  class_<StringLogSink, StringLogSinkPtr, bases<LogSink>,
         boost::noncopyable >("StringLogSink", init<>())
    .def("LogMessage", &StringLogSink::LogMessage)
    .def("GetLog", &StringLogSink::GetLog)
  ;

  enum_<Logger::LogLevel>("LogLevel")
    .value("Error", Logger::QUIET)
    .value("Warning", Logger::WARNING)
    .value("Script", Logger::SCRIPT)
    .value("Info", Logger::INFO)
    .value("Verbose", Logger::VERBOSE)
    .value("Debug", Logger::DEBUG)
    .value("Trace", Logger::TRACE)
  ;

  def("PushVerbosityLevel", push_verb);
  def("PopVerbosityLevel", pop_verb);
  def("GetVerbosityLevel", get_verb);
  def("PushLogSink", push_log_sink);
  def("GetCurrentLogSink", get_log_sink);
  def("PopLogSink", pop_log_sink);
  def("LogError", raw_function(log_error, 1));
  def("LogWarning",raw_function(log_warning, 1));
  def("LogInfo", raw_function(log_info, 1));
  def("LogScript", raw_function(log_script, 1));
  def("LogVerbose", raw_function(log_verbose, 1));
  def("LogDebug", raw_function(log_debug, 1));
  def("LogTrace", raw_function(log_trace, 1));
  
  // this relatively ugly construct is required to work around a problem with
  // the "ost" command-line interpreter. If we don't remove all the sinks from
  // the sink stack, we will get "Fatal Python error: PyEval_SaveThread: 
  // NULL tstate" upon exiting ost. I don't completely understand why, though.
  scope().attr("__dict__")["atexit"]=handle<>(PyImport_ImportModule("atexit"));

  def("_reset_sinks", &reset_sinks);
  object r=scope().attr("_reset_sinks");
  scope().attr("atexit").attr("register")(r);
}
