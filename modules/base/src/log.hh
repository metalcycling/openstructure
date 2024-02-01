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
#ifndef OST_LOG_HH
#define OST_LOG_HH

#include <ostream>
#include <stack>

#include <ost/log_sink.hh>
#include <ost/module_config.hh>

namespace ost {

typedef std::stack<LogSinkPtr> LogSinkStack;

  // singleton
class DLLEXPORT_OST_BASE Logger {
public:
  enum LogLevel {
    QUIET   =0,
    WARNING =1,
    SCRIPT  =2,
    INFO    =3,
    VERBOSE =4,
    DEBUG   =5,
    TRACE   =6
  };

  void PushVerbosityLevel(int level);
  void PopVerbosityLevel();
  void PushSink(LogSinkPtr& sink);
  void PushFile(const String& filename);
  //! DEPRECATED use PopSink() instead
  void PopFile();
  void PopSink();

  static Logger& Instance();
  LogSinkPtr GetCurrentSink() { return sink_stack_.top(); }
  int GetVerbosityLevel() const {return level_;}
  
  void ResetSinks() {
    while (sink_stack_.size()>1) {
      sink_stack_.pop();
    }
  }
protected:
  Logger();
  Logger(const Logger&);
  Logger& operator=(const Logger&);
  
private:
  int level_;
  std::stack<int> level_stack_;
  LogSinkStack sink_stack_;
};


#define OST_DO_LOGGING_(m, l) if (::ost::Logger::Instance().GetVerbosityLevel()>=l) {\
    std::stringstream tmp_s__;                                                 \
    tmp_s__ << m << std::endl;                                                 \
    ::ost::Logger::Instance().GetCurrentSink()->LogMessage(tmp_s__.str(), l);  \
  }

#define WARN_DEPRECATED(m) OST_DO_LOGGING_(m, ::ost::Logger::WARNING)
#define PUSH_VERBOSITY(n) ::ost::Logger::Instance().PushVerbosityLevel(n)
#define POP_VERBOSITY(n) ::ost::Logger::Instance().PopVerbosityLevel()

#define LOG_ERROR(m) OST_DO_LOGGING_(m, ::ost::Logger::QUIET)
#define LOG_WARNING(m) OST_DO_LOGGING_(m, ::ost::Logger::WARNING)
#define LOG_SCRIPT(m) OST_DO_LOGGING_(m, ::ost::Logger::SCRIPT)
#define LOG_INFO(m) OST_DO_LOGGING_(m, ::ost::Logger::INFO)
#define LOG_VERBOSE(m) OST_DO_LOGGING_(m, ::ost::Logger::VERBOSE)
#ifdef NDEBUG
#  define LOG_DEBUG(m)
#  define LOG_TRACE(m)
#else
#  define LOG_DEBUG(m) OST_DO_LOGGING_(m, ::ost::Logger::DEBUG)
#  define LOG_TRACE(m) OST_DO_LOGGING_(m, ::ost::Logger::TRACE)
#endif

}

#endif
