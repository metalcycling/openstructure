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

#include <stdarg.h>
#include <cassert>
#include <complex>
#include <ost/log.hh>

#include "tiff_util.hh"


namespace ost { namespace io { namespace detail {

typedef void (*TIFFWarningHandler)(const char*, const char*, va_list);

void tiff_warning_handler(const char *mod, const char* fmt, va_list ap)
{
  static char message[1024];
#ifdef _MSC_VER
  _vsnprintf(message,1024,fmt,ap);
#else
  vsnprintf(message,1024,fmt,ap);
#endif
  LOG_INFO(mod << ": " << message);
}

tsize_t CustomTIFFReadProcIStream(thandle_t thandle, tdata_t tdata, tsize_t tsize)
{
   std::istream* file= reinterpret_cast<std::istream*>(thandle);
   char* data= reinterpret_cast<char*>(tdata);
   file->read(data,tsize);
   return file->gcount();
}

tsize_t CustomTIFFReadProcOStream(thandle_t thandle, tdata_t tdata, tsize_t tsize)
{
  assert(false);
  return -1;
}

tsize_t CustomTIFFWriteProcIStream(thandle_t thandle, tdata_t tdata, tsize_t tsize)
{
  assert(false);
  return -1;
}

tsize_t CustomTIFFWriteProcOStream(thandle_t thandle, tdata_t tdata, tsize_t tsize)
{
   std::ostream* file= reinterpret_cast<std::ostream*>(thandle);
   char* data= reinterpret_cast<char*>(tdata);
   tsize_t before = file->tellp();
   file->write(data,tsize);
   tsize_t after = file->tellp();
   return after-before;
}

toff_t CustomTIFFSeekProcIStream(thandle_t thandle, toff_t toff, int dir)
{
  std::istream* stream= reinterpret_cast<std::istream*>(thandle);

  std::ios_base::seekdir sdir = std::ios_base::beg;

  if (dir == 1)
  {
    sdir = std::ios_base::cur;
  } else if (dir == 2 ) {
    sdir = std::ios_base::end;
  }

  stream->seekg(toff,
                sdir);
  return stream->rdstate();
}

toff_t CustomTIFFSeekProcOStream(thandle_t thandle, toff_t toff, int dir)
{
  std::ostream* stream= reinterpret_cast<std::ostream*>(thandle);

  std::ios_base::seekdir sdir = std::ios_base::beg;

  if (dir == 1)
  {
    sdir = std::ios_base::cur;
  } else if (dir == 2 ) {
    sdir = std::ios_base::end;
  }

  stream->seekp(toff,sdir);
  return stream->rdstate();
}

int CustomTIFFCloseProc(thandle_t thandle)
{
  return 0;
}

toff_t CustomTIFFSizeProcIStream(thandle_t thandle)
{
   std::istream* stream= reinterpret_cast<std::istream*>(thandle);
   toff_t curr_pos = stream->tellg();
   stream->seekg(0,std::ios::end);
   toff_t size = stream->tellg();
   stream->seekg(curr_pos,std::ios::beg);
   return size;
}

toff_t CustomTIFFSizeProcOStream(thandle_t thandle)
{
   std::ostream* stream= reinterpret_cast<std::ostream*>(thandle);
   toff_t curr_pos = stream->tellp();
   stream->seekp(0,std::ios::end);
   toff_t size = stream->tellp();
   stream->seekp(curr_pos,std::ios::beg);
   return size;
}

int CustomTIFFMapFileProc(thandle_t thandle, tdata_t* tdata, toff_t* toff)
{
  assert(false);
  return(0);
}

void CustomTIFFUnmapFileProc(thandle_t thandle, tdata_t tdata, toff_t toff)
{
  assert(false);
}

}}} //ns
