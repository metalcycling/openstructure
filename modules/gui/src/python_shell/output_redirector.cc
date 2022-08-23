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
  OutputRedirector

  Author: Andreas Schenk
*/

#include "output_redirector.hh"
#include <QTimer>


namespace ost { namespace gui {

OutputRedirector::OutputRedirector():
 QObject(),
 buffer_(),
 timer_()
{
}


void OutputRedirector::Write( String const& str )
{
  if(timer_.isValid()){
    if(timer_.elapsed()>1000){
      QString output = buffer_+QString::fromStdString(str);
      buffer_="";
      emit OnOutput(output);
      timer_.restart();
    }else{
      buffer_+=QString::fromStdString(str);
    }
  }else{
    buffer_+=QString::fromStdString(str);
    timer_.start();
    QTimer::singleShot(0, this, SLOT(Flush()));
  }
}

void OutputRedirector::Flush()
{
  if(! timer_.isValid()){
    return;
  }
  timer_=QElapsedTimer();
  QString output = buffer_;
  buffer_="";
  emit OnOutput(output);
}



}}//ns
