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
#ifndef OST_IMG_IO_TIFF_HH
#define OST_IMG_IO_TIFF_HH

#include <ost/stdint.hh>
#include <tiff.h>
#include <tiffio.h>
#include <boost/filesystem.hpp>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/io/img/image_format.hh>


namespace ost { namespace io { namespace detail {
/// \internal
typedef void (*TIFFWarningHandler)(const char*, const char*, va_list);

/// \internal
void tiff_warning_handler(const char *mod, const char* fmt, va_list ap);

/// \internal
struct tiff_warning_handler_wrapper {
  tiff_warning_handler_wrapper() {
    handler_ = TIFFSetWarningHandler(tiff_warning_handler);
  }

  ~tiff_warning_handler_wrapper() {
    TIFFSetWarningHandler(handler_);
  }

  TIFFWarningHandler handler_;
};

/// \internal
class complexint32:public std::complex<int32_t>{
public:

  operator std::complex<Real>()
  {
    return std::complex<Real>(real(),imag());
  }
};

/// \internal
class complexint16:public std::complex<int16_t>{
public:

  operator std::complex<Real>()
  {
    return std::complex<Real>(real(),imag());
  }
};

/// \internal
class complexint8:public std::complex<int8_t>{
public:

  operator std::complex<Real>()
  {
    return std::complex<Real>(real(),imag());
  }
};

/// \internal
tsize_t CustomTIFFReadProcIStream(thandle_t thandle, tdata_t  tdata, tsize_t tsize);
/// \internal
tsize_t CustomTIFFReadProcOStream(thandle_t thandle, tdata_t  tdata, tsize_t tsize);
/// \internal
tsize_t CustomTIFFWriteProcIStream(thandle_t thandle, tdata_t  tdata, tsize_t tsize);
/// \internal
tsize_t CustomTIFFWriteProcOStream(thandle_t thandle, tdata_t  tdata, tsize_t tsize);
/// \internal
toff_t CustomTIFFSeekProcIStream(thandle_t thandle, toff_t toff, int dir);
/// \internal
toff_t CustomTIFFSeekProcOStream(thandle_t thandle, toff_t toff, int dir);
/// \internal
int CustomTIFFCloseProc(thandle_t thandle);
/// \internal
toff_t CustomTIFFSizeProcIStream(thandle_t thandle);
/// \internal
toff_t CustomTIFFSizeProcOStream(thandle_t thandle);
/// \internal
int CustomTIFFMapFileProc(thandle_t thandle, tdata_t* tdata, toff_t* toff);
/// \internal
void CustomTIFFUnmapFileProc(thandle_t thandle, tdata_t tdata, toff_t toff);

/// \internal
template<typename IN_TYPE,typename OUT_TYPE, class IST>
void do_tiff_read(tdata_t buf,unsigned int rps, unsigned int width, IST* is,unsigned int& current_row, uint16_t spp)
{
  IN_TYPE* dp = static_cast<IN_TYPE*>(buf);
  for(uint r=0;r<rps;r++) {
    for(uint c=0;c<width;c++) {
      is->Value(img::Point(c,current_row))=static_cast<OUT_TYPE>(dp[(r*width+c)*spp]);
    }
    current_row++;
  }
}

/// \internal
template<typename IN_TYPE,typename OUT_TYPE, class IST>
void do_tiff_write(TIFF *tif, IST* is,unsigned int rowsperstrip,unsigned int width, unsigned int height, unsigned int strip,const  img::NormalizerPtr& nptr)
{
  uint datalength=rowsperstrip*width;
  if((strip+1)*rowsperstrip>height){
    datalength=(height-strip*rowsperstrip)*width;
  }
  OUT_TYPE* buf=new OUT_TYPE[datalength];
  img::Point start = is->GetExtent().GetStart();

  uint i=0;
  for(uint r=strip*rowsperstrip;r<(strip+1)*rowsperstrip && r<height;r++) {
    for(uint c=0;c<width;c++) {
      buf[i] = static_cast<OUT_TYPE>(nptr->Convert(is->Value(img::Point(c,r)+start)));
      ++i;
    }
  }
  TIFFWriteEncodedStrip(tif,strip,buf, sizeof(OUT_TYPE)*datalength);
  delete[] buf;
}

}}} // ns

#endif // TIFF_HH
