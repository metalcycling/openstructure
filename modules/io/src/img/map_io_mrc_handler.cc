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
  Authors: Ansgar Philippsen, Andreas Schenk, Valerio Mariani
*/
/*

MRC vs IPLT half-complex format

letters denote consecutive complex
values in the file, an appended *
stands for complex conjugate

  MRC                 IPLT


  6x6
a b c d . .       . . . . . .   -2
e f g h . .       . . . . . .   -1
i j k l . .       o*n*m n o p    0
m n o p . .       k*j*q r s t    1  q = i*, t = l*
q r s t . .       g*f*u v w x    2  u = e*, x = h*
u v w x . .       c*b*a*b c d    3

                 -2-1 0 1 2 3

  7x6
a b c d . . .     . . . . . . .  -2
e f g h . . .     . . . . . . .  -1
i j k l . . .     p*o*n*m n o p   0
m n o p . . .     l*k*j*q r s t   1
q r s t . . .     h*g*f*u v w x   2
u v w x . . .     d*c*b*a b c d   3

                 -3-2-1 0 1 2 3

*/

#include <cstdio>
#include <iostream>
#include <sstream>
#include <streambuf>

#include <boost/scoped_array.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/format.hpp>

#include <ost/base.hh>
#include <ost/log.hh>
#include <ost/base.hh>

#include <ost/img/image_state.hh>
#include <ost/img/alg/stat_accumulator.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/img/util.hh>
#include <ost/img/progress.hh>
#include <ost/img/alg/fractional_shift.hh>
#include <ost/img/alg/alg_shift.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/img/image_format.hh>
#include <ost/io/swap_util.hh>
#include <ost/io/converting_streams.hh>
#include <ost/io/io_utils.hh>


#include "map_io_mrc_handler.hh"

namespace ost { namespace io {

String MRC::FORMAT_STRING = "defined_mrc";


MRC::MRC(bool normalize_on_save, Subformat subformat,Endianess endianess_on_save, Format bit_depth):
  ImageFormatBase(FORMAT_STRING),
  subformat_(subformat),
  normalize_on_save_(normalize_on_save),
  endianess_on_save_(endianess_on_save),
  bit_depth_(bit_depth)
{
}
Format MRC::GetBitDepth() const
{
  return bit_depth_;
}

void MRC::SetBitDepth (Format bitdepth)
{
  bit_depth_ = bitdepth;
}


Endianess MRC::GetEndianessOnSave() const
{
  return endianess_on_save_;
}

void MRC::SetEndianessOnSave(Endianess end)
{
  endianess_on_save_ = end;
}

Subformat MRC::GetSubformat() const
{
  return subformat_;
}

void MRC::SetSubformat(Subformat subf)
{
  subformat_ = subf;
}

bool MRC::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void MRC::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}

CCP4::CCP4(bool normalize_on_save, Endianess endianess_on_save, Format bit_depth):
  MRC(normalize_on_save,MRC_NEW_FORMAT,endianess_on_save,bit_depth)
{}




namespace detail{

class ptristream : public std::istream
{
  class ptrinbuf : public std::streambuf
  {
  protected:
    char *      ptr;
    std::size_t len;

  public:
    ptrinbuf(char * _ptr, std::size_t _len) : ptr(_ptr), len(_len) {
      assert(ptr);
      if (*ptr && len == 0)
        len = std::strlen(ptr);

      setg(ptr,               // beginning of putback area
           ptr,               // read position
           ptr+len);          // end position
    }

  protected:
    virtual int_type underflow() {
      // is read position before end of buffer?
      if (gptr() < egptr())
        return traits_type::to_int_type(*gptr());
      else
        return EOF;
    }

    virtual pos_type seekoff(off_type off, ios_base::seekdir way,
                             ios_base::openmode mode =
                             ios_base::in | ios_base::out)
    {
      switch (way) {
      case std::ios::cur:
        setg(ptr, gptr()+off, ptr+len);
        break;
      case std::ios::beg:
        setg(ptr, ptr+off, ptr+len);
        break;
      case std::ios::end:
        setg(ptr, egptr()+off, ptr+len);
        break;

      default:
        assert(false);
        break;
      }
      return pos_type(gptr() - ptr);
    }
  };

protected:
  ptrinbuf buf;
public:
 ptristream(char * ptr, std::size_t len = 0)
    : std::istream(0), buf(ptr, len) {
    rdbuf(&buf);
  }
};


class header_base {
public:
  header_base():
    nc(0),nr(0),ns(0),
    mode(0),
    ncstart(0),nrstart(0),nsstart(0),
    nx(0),ny(0),nz(0),
    x(0.0),y(0.0),z(0.0),
    alpha(0.0),beta(0.0),gamma(0.0),
    mapc(1),mapr(2),maps(3),
    amin(0.0),amax(0.0),amean(0.0),
    ispg(0),
    nsymbt(0),
    nlabel(1),
    label()
  {
    for(unsigned int i=0;i<800;++i)
    {
      label[i]=' ';
    }
  }
  header_base(const img::ConstImageHandle& im, Format bit_depth):
    nc(),
    nr(static_cast<int>(im.GetExtent().GetSize().GetHeight())),
    ns(static_cast<int>(im.GetExtent().GetSize().GetDepth())),
    mode(),
    ncstart(im.GetExtent().GetStart()[0]),
    nrstart(im.GetExtent().GetStart()[1]),
    nsstart(im.GetExtent().GetStart()[2]),
    nx(static_cast<int>(im.GetExtent().GetSize().GetWidth())),
    ny(static_cast<int>(im.GetExtent().GetSize().GetHeight())),
    nz(static_cast<int>(im.GetExtent().GetSize().GetDepth())),
    x(),
    y(),
    z(),
    alpha(90.0),beta(90.0),gamma(90.0),
    mapc(1),mapr(2),maps(3),
    amin(0.0),amax(0.0),amean(0.0),
    ispg(0),
    nsymbt(0),
    nlabel(1),
    label()
  {
    if(im.GetType()==img::REAL){
      nc=static_cast<int>(im.GetExtent().GetSize().GetWidth());
      switch(bit_depth){
      case OST_BIT8_FORMAT:
        mode=0;
        break;
      case OST_BIT16_FORMAT:
        mode=1;
        break;
      case OST_FLOAT_FORMAT:
      case OST_DEFAULT_FORMAT:
        mode=2;
        break;
      default:
        throw(IOException("MRC/CCP4 export: Bit depth not supported."));
        break;
      }
      x=im.GetExtent().GetSize().GetWidth()*im.GetSpatialSampling()[0];
      y=im.GetExtent().GetSize().GetHeight()*im.GetSpatialSampling()[1];
      z=im.GetExtent().GetSize().GetDepth()*im.GetSpatialSampling()[2];
    }else{
      nc=static_cast<int>(im.GetExtent().GetSize().GetWidth()/2 +1);
      switch(bit_depth){
      case OST_BIT16_FORMAT:
        mode=3;
        break;
      case OST_FLOAT_FORMAT:
      case OST_DEFAULT_FORMAT:
        mode=4;
        break;
      default:
        throw(IOException("MRC/CCP4 export: Bit depth not supported."));
        break;
      }
      x=1.0;
      y=1.0;
      z=1.0;
    }
    ost::img::alg::StatAccumulator<1> acc;
    if(im.GetType()==img::REAL) {
      img::image_state::RealSpatialImageState *isr=dynamic_cast<img::image_state::RealSpatialImageState*>(im.ImageStatePtr().get());
      if(! isr){
        throw(IOException("MRC/CCP4 export: dynamic cast failed in header base."));
      }
      for(Real* ptr = isr->Data().GetData(); ptr<isr->Data().GetEnd(); ++ptr) {
        acc(*ptr);
      }
      amin=acc.GetMinimum();
      amax=acc.GetMaximum();
      amean=acc.GetMean();
    }
    for(unsigned int i=0;i<800;++i)
    {
      label[i]=' ';
    }
  }
  template<int CONVERSIONTYPE>
  void WriteCommonData(BinaryOStream<CONVERSIONTYPE>& out) const
  {
    out << nc << nr << ns;
    out << mode;
    out << ncstart << nrstart << nsstart;
    out << nx << ny << nz;
    out << x << y << z << alpha << beta << gamma;
    out << mapc << mapr << maps;
    out << amin << amax << amean;
    out << ispg;
    out << nsymbt;
  }
  template<int CONVERSIONTYPE>
  void WriteLabel(BinaryOStream<CONVERSIONTYPE>& out) const
  {
    out << nlabel;
    out.write(label,800);
  }
  template<int CONVERSIONTYPE>
  void ReadCommonData(BinaryIStream<CONVERSIONTYPE>& in)
  {
    in >> nc >> nr >> ns;
    in >> mode;
    in >> ncstart >> nrstart >> nsstart;
    in >> nx >> ny >> nz;
    in >> x >> y >> z >> alpha >> beta >> gamma;
    in >> mapc >> mapr >> maps;
    in >> amin >> amax >> amean;
    in >> ispg;
    in >> nsymbt;
  }
  template<int CONVERSIONTYPE>
  void ReadLabel(BinaryIStream<CONVERSIONTYPE>& in)
  {
    in >> nlabel;
    in.read(label,800);
  }
  void Print()
  {
    LOG_INFO("ncrs: " << nc << " " << nr << " " << ns);
    LOG_INFO("mode " << mode);
    LOG_INFO("nxyz: " << nx << " " << ny << " " << nz);
    LOG_INFO("cell: " << x << " " << y << " " << z);
    LOG_INFO(alpha << " " << beta << " " << gamma);
    LOG_INFO("order: " << mapc << mapr << maps);
    LOG_INFO("sg: " << ispg << " " << nsymbt);
    LOG_INFO("nlabel: " << nlabel);
  }

  int nc,nr,ns;
  int mode;
  int ncstart,nrstart,nsstart;
  int nx,ny,nz;
  float x,y,z;
  float alpha,beta,gamma;
  int mapc,mapr,maps;
  float amin,amax,amean;
  int ispg;
  int nsymbt;
  int nlabel;
  char label[800];
};

// mrc old style header
class mrc_header: public header_base
{
public:
  mrc_header():
    header_base(),
    xorigin(0.0),
    yorigin(0.0)
  {
  }
  mrc_header(const img::ConstImageHandle& im, Format bit_depth):
    header_base(im,bit_depth),
    xorigin(0.0), // todo determine origin
    yorigin(0.0)
  {
  }
  void Print()
  {
    LOG_INFO("mrc header");
    header_base::Print();
    LOG_INFO("ori: " << xorigin << " " << yorigin);
  }
  static int DetermineDataFormat( std::istream& f)
  {
    float x;
    unsigned int mapc;
    BinaryIStream<OST_LITTLE_ENDIAN> bs(f);
    bs.seekg(10*4,std::ios::beg); // seek to x
    bs >> x;
    bs.seekg(5*4,std::ios::cur); // seek to mapc
    bs >> mapc;
    bs.seekg(0,std::ios::beg); // seek to beginning
    if(mapc<1 || mapc>3) {
      // try byte swapping
      Convert<OST_BIG_ENDIAN,unsigned int>::FromIP(&mapc);
      if(mapc<1 || mapc>3) {
        throw(IOException("error reading mrc file (even tried byte-swapping"));
      }
      return OST_BIG_ENDIAN;
    }else if(x>0.0 || x<1.0) {
      Convert<OST_VAX_DATA,float>::FromIP(&x);
      if(x<1.0) {
        LOG_VERBOSE("suspicious floating Point value in mrc header");
        return OST_LITTLE_ENDIAN;
      }else{
        return OST_VAX_DATA;
      }
    }
    return OST_LITTLE_ENDIAN;
  }

  float xorigin,yorigin;
};

template<int CONVERSIONTYPE>
BinaryOStream<CONVERSIONTYPE>& operator<< (BinaryOStream<CONVERSIONTYPE>& out, const mrc_header& h )
{
  h.WriteCommonData(out);
  char dummy[116];
  for(int i=0;i<116;i++) dummy[i]=0;
  out.write(dummy,116);
  out << h.xorigin << h.yorigin;
  h.WriteLabel(out);
  return  out;
}

template<int CONVERSIONTYPE>
BinaryIStream<CONVERSIONTYPE>& operator>> (BinaryIStream<CONVERSIONTYPE>& in, mrc_header& h)
{
  char dummy[116];
  h.ReadCommonData(in);
  in.read(dummy,116);
  in >> h.xorigin >> h.yorigin;
  h.ReadLabel(in);

  if(h.nc<1) h.nc=1;
  if(h.nr<1) h.nr=1;
  if(h.ns<1) h.ns=1;

  // skip symm info, seek to start of actual data
  if(h.nsymbt>0){
	for (int counter=0;counter<h.nsymbt;++counter){
	  in.read(dummy,1);
	}
  }
  return in;
}



class ccp4_header: public header_base {
public:
  ccp4_header():
    header_base(),
    lskflag(),
    skwmat(),
    skwtrn(),
    ox(0.0),
    oy(0.0),
    oz(0.0),
    arms()
  {
  }
  ccp4_header(const img::ConstImageHandle& im, Format bit_depth):
    header_base(im,bit_depth),
    lskflag(),
    skwmat(),
    skwtrn(),
    ox(0.0),
    oy(0.0),
    oz(0.0),
    arms()
  {
    lskflag=0;
    for(unsigned int i=0;i<9;++i){
      skwmat[i]=0.0;
    }
    for(unsigned int i=0;i<3;++i){
      skwtrn[i]=0.0;
    }
    arms = 0.1;
  }

  static int DetermineDataFormat(std::istream& f)
  {
    char machst[4];
    // from the ccp4 documentation
    // The machine stamp is a 32-bit quantity containing a set of four `nibbles'
    // (half-bytes)---only half the space is used. Each nibble is a number
    // specifying the representation of (in C terms) double (d) , float (f),
    // int (i) and unsigned char (c) types. Thus each stamp is of the form
    // 0xdfic0000. The values for the floating point nibbles may be taken from the
    // list (following HDF):
    // 1        Big-endian ieee
    // 2        VAX
    // 3        Cray
    // 4        Little-endian ieee
    // 5        Convex native
    // 6        Fijitsu VP
    f.seekg(53*4,std::ios::beg); // seek to machine stamp
    f.read(machst,4);
    f.seekg(0,std::ios::beg); // seek to beginning
    char float_machst= machst[0] & 0x0f;
    if(float_machst == 1){
      LOG_DEBUG("CCP4Import: reading big endian data");
      return OST_BIG_ENDIAN;
    }else if(float_machst == 2){
      LOG_DEBUG("CCP4Import: reading vax data");
      return OST_VAX_DATA;
    }else if(float_machst == 4){
      LOG_DEBUG("CCP4Import: reading little endian data");
      return OST_LITTLE_ENDIAN;
    } else{
      LOG_VERBOSE("unknown machine stamp (" << (int)float_machst << "), using little endian");
      return OST_LITTLE_ENDIAN;
      //throw(IOException("CCP4Import: Cray, Convex native and Fijitsu VP formats are not supported."));
    }
  }

	
  void Print()
  {
    using boost::format;
    namespace bf = boost::filesystem;

    LOG_INFO("ccp4 header:");
    header_base::Print();
    LOG_INFO("arms: " << arms);
  }

  int lskflag;
  float skwmat[9];
  float skwtrn[3];
  float ox,oy,oz;
  float arms;
};




template <int CONVERSIONTYPE>
BinaryIStream<CONVERSIONTYPE>& operator>> (BinaryIStream<CONVERSIONTYPE>& in, ccp4_header& header)
{
  char dummy[48];
  header.ReadCommonData(in);
  in >> header.lskflag;
  in >> header.skwmat[0] >> header.skwmat[1] >> header.skwmat[2];
  in >> header.skwmat[3] >> header.skwmat[4] >> header.skwmat[5];
  in >> header.skwmat[6] >> header.skwmat[7] >> header.skwmat[8];
  in >> header.skwtrn[0] >> header.skwtrn[1] >> header.skwtrn[2];
  in.read(dummy,48);
  in >> header.ox >> header.oy >> header.oz;
  in.read(dummy,8);
  in >> header.arms;
  header.ReadLabel(in);

  if(header.nx<1) header.nx=1;
  if(header.ny<1) header.ny=1;
  if(header.nz<1) header.nz=1;

  // skip symm info, seek to start of actual data
  if(header.nsymbt>0) {
	for (int counter=0;counter<header.nsymbt;++counter){
      in.read(dummy,1);
	}
  }
  return in;
}
template <int CONVERSIONTYPE>
BinaryOStream<CONVERSIONTYPE>& operator<< (BinaryOStream<CONVERSIONTYPE>& out, ccp4_header& header)
{
  header.WriteCommonData(out);
  out << header.lskflag;
  out << header.skwmat[0] << header.skwmat[1] << header.skwmat[2];
  out << header.skwmat[3] << header.skwmat[4] << header.skwmat[5];
  out << header.skwmat[6] << header.skwmat[7] << header.skwmat[8];
  out << header.skwtrn[0] << header.skwtrn[1] << header.skwtrn[2];
  char dummy[48];
  for(unsigned int i=0;i<48;++i){
    dummy[i]=0;
  }
  out.write(dummy,48);
  out << header.ox << header.oy << header.oz;
  //space after MAP has to be written to file to produce valid new style format
  out.write("MAP ",4);
  char machst[4]={0x44,0x41,0x00,0x00};// little endian
  if(CONVERSIONTYPE==OST_VAX_DATA){
     machst[0]=0x22;
     machst[1]=0x41;
  }else if(CONVERSIONTYPE==OST_BIG_ENDIAN){
     machst[0]=0x11;
     machst[1]=0x11;
  }
  out.write(reinterpret_cast<char*>(&machst),4);
  out << header.arms;
  header.WriteLabel(out);
  return out;
}

using namespace ost::img;
using namespace ost;


template <typename B,int CONVERSIONTYPE>
void real_filler(img::image_state::RealSpatialImageState& isi,
                 BinaryIStream<CONVERSIONTYPE>& f,
                 const header_base& header)
{
  int mapc=header.mapc-1;
  int mapr=header.mapr-1;
  int maps=header.maps-1;

  img::Point pnt;
  char this_dummy; //create dummy variable to give to img::Progress as this
  img::Progress::Instance().Register(&this_dummy,header.ns*header.nr,100);
  for(int is=0;is<header.ns;++is) {
    pnt[maps]=header.nsstart+is;
    std::vector<B> buffer(header.nr*header.nc);
    f.read(&buffer[0],header.nr*header.nc);
    for(int ir=0;ir<header.nr;++ir) {
      pnt[mapr]=header.nrstart+ir;
      for(int ic=0;ic<header.nc;++ic) {
        pnt[mapc]=header.ncstart+ic;
        isi.Value(pnt) = img::Val2Val<B,Real>(buffer[header.nc*ir+ic]);
      }
      img::Progress::Instance().AdvanceProgress(&this_dummy);
    }
  }
  img::Progress::Instance().DeRegister(&this_dummy);
}

template <typename B,int CONVERSIONTYPE>
void complex_filler(img::image_state::ComplexHalfFrequencyImageState& isi,
                    BinaryIStream<CONVERSIONTYPE> & fhandle,
                    header_base& header)
{
  int mapc=header.mapc-1;
  int mapr=header.mapr-1;
  int maps=header.maps-1;
  B real, imag;
  if(header.nx & 0x1 || header.ny & 0x1) {
    throw ost::io::IOException("Only EVENxEVEN MRC half-complex format supported");
  }
  // always assume hcomplex (!)
  if((header.nc-1)*2 == header.nr) {
    char this_dummy; //create dummy variable to give to Progress as this
    Progress::Instance().Register(&this_dummy,header.ns*header.nr,100);
    Point p;
    for(int ss=0;ss<header.ns;++ss) {
      p[maps] = ss>header.ns/2 ? header.ns-ss : ss;
      int sr=header.nr;
      int cr=0;
      // first half
      for(;sr>header.nr/2;--sr,++cr) {
        p[mapr] = sr-header.nr/2;
        int sc=0;
        for(;sc<header.nc-1;sc++) {
          p[mapc]=-sc;
          // complex conjugate
          fhandle >> real >> imag;
          isi.Value(p)=Complex(Real(real),Real(imag));
          LOG_DEBUG(" " << p  << " " << isi.Value(p));
        }
        fhandle >> real >> imag;
        if(sr==header.nr) {
        // why set point (py,header.ny/2,pz)?
        //  isi.Value(Point(py,header.ny/2,pz))=scale*Complex(Real(real),Real(imag));
        //  LOG_DEBUG("+" << Point(py,header.ny/2,pz) << " <- " << Point(sx,cy,sz) << " " << " " << isi.Value(Point(py,header.ny/2,pz)));
          p[mapc]=p[header.mapr];
          p[mapr]=header.nr/2;
          isi.Value(p)=Complex(Real(real),-Real(imag));
          LOG_DEBUG("+" << p << " " << isi.Value(p));
        }
        Progress::Instance().AdvanceProgress(&this_dummy);
      }
      // second half, sr starts at header.nr/2
      for(;sr>0;--sr,++cr) {
        p[mapr]=header.nr/2-sr;
        int sc=0;
        for(;sc<header.nc-1;++sc) {
          p[mapc]=sc;
          fhandle >> real >> imag;
          isi.Value(p)=Complex(Real(real),-Real(imag));
          LOG_DEBUG(" " << p << " " << isi.Value(p));
        }
        fhandle >> real >> imag;
        p[mapc]=sc;
        isi.Value(p)=Complex(Real(real),Real(imag));
        LOG_DEBUG(" " << p << " " << isi.Value(p));
        Progress::Instance().AdvanceProgress(&this_dummy);
      }
      // set second half of r=0 and r=nr/2 line
      for(int sc=1;sc<header.nc-1;++sc) {
        p[mapc]=sc;
        p[mapr]=0;
        p[maps]=ss;
        isi.Value(Point(-p[0],p[1],p[2]))=std::conj(isi.Value(p));
        LOG_DEBUG(" " << Point(-p[0],p[1],p[2]) << " <- " << p << "**");
        p[mapr]=header.nr/2;
        isi.Value(p)=std::conj(isi.Value(Point(-p[0],p[1],p[2])));
        LOG_DEBUG(" " << p << " <- " << Point(-p[0],p[1],p[2]) << "**");
      }
    }
    Progress::Instance().DeRegister(&this_dummy);
  } else if(header.nc*2 == header.nr+1) {
    throw IOException("Half Complex Format NX = (NY+1)/2 not supported");
  } else {
    throw IOException("expected half-complex format, ie NX = NY/2+1 or NX = (NY+1)/2");
  }
}

template <typename B,int CONVERSIONTYPE>
void real_dumper(BinaryOStream<CONVERSIONTYPE>& f, header_base& header, const ConstImageHandle& image, const MRC& formatmrc)
{
  image_state::RealSpatialImageState *isr=dynamic_cast<image_state::RealSpatialImageState*>(image.ImageStatePtr().get());
  if(! isr){
    throw(IOException("MRC/CCP4 export: dynamic cast failed in real dumper."));
  }
  int mapc=header.mapc-1;
  int mapr=header.mapr-1;
  int maps=header.maps-1;

  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatmrc.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(image,formatmrc.GetMinimum(),formatmrc.GetMaximum());
  }

  img::Point pnt;
  char this_dummy; //create dummy variable to give to Progress as this
  Progress::Instance().Register(&this_dummy,header.ns*header.nr,100);
  for(int ss=0;ss<header.ns;++ss) {
    std::vector<B> buffer(header.nr*header.nc);
    pnt[maps]=header.nsstart+ss;
    for(int sr=0;sr<header.nr;++sr) {
      pnt[mapr]=header.nrstart+sr;
      for(int sc=0;sc<header.nc;++sc) {
        pnt[mapc]=header.ncstart+sc;
        buffer[header.nc*sr+sc]=img::Val2Val<Real,B>(norm.Convert(isr->Value(pnt)));
      }
      Progress::Instance().AdvanceProgress(&this_dummy);
    }
    f.write(&buffer[0],header.nr*header.nc);
  }
  Progress::Instance().DeRegister(&this_dummy);
}

template <typename B,int CONVERSIONTYPE>
void complex_dumper(BinaryOStream<CONVERSIONTYPE>& f,
                    header_base& header,
                    const ConstImageHandle& image,const MRC& formatmrc)
{
  int mapc=header.mapc-1;
  int mapr=header.mapr-1;
  int maps=header.maps-1;

  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatmrc.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(image,formatmrc.GetMinimum(),formatmrc.GetMaximum());
  }

  image_state::ComplexHalfFrequencyImageState *isc=dynamic_cast<image_state::ComplexHalfFrequencyImageState*>(image.ImageStatePtr().get());
  if(! isc){
    throw(IOException("MRC/CCP4 export: dynamic cast failed in complex dumper."));
  }
  char this_dummy; //create dummy variable to give to Progress as this
  Progress::Instance().Register(&this_dummy,header.ns*header.nr,100);
  img::Point pnt;
  for(int ss=0;ss<header.ns;++ss) {
    pnt[maps]=ss;
    // first half
    int sr=0;
    for(;sr<header.nr/2;++sr) {
      pnt[mapr]=header.nr/2-sr;
      for(int sc=0;sc<header.nc-1;++sc) {
        pnt[mapc]=-sc;
        Complex val = norm.Convert(isc->Value(pnt));
        f << static_cast<B>(val.real()) << static_cast<B>(val.imag());
        LOG_DEBUG(" " << pnt  << " " << val);
      }
      f << static_cast<B>(0.0) << static_cast<B>(0.0);
      Progress::Instance().AdvanceProgress(&this_dummy);
    }
    // second half
    for(;sr<header.nr;++sr) {
      pnt[mapr]=sr-header.nr/2;
      int sc=0;
      for(;sc<header.nc-1;++sc) {
        pnt[mapc]=sc;
        Complex  val =norm.Convert(isc->Value(pnt));
        f << static_cast<B>(val.real()) << static_cast<B>(-val.imag());
        LOG_DEBUG(" " << pnt << " " << val);
      }
      pnt[mapc]=sc;
      Complex  val = norm.Convert(conj(isc->Value(pnt)));
      f << static_cast<B>(val.real()) << static_cast<B>(-val.imag());
      LOG_DEBUG(" " << pnt << " " << val);
      Progress::Instance().AdvanceProgress(&this_dummy);
    }
  }
  Progress::Instance().DeRegister(&this_dummy);
}

template<class HEADER,int CONVERSIONTYPE>
void import_helper(img::MapHandle& image, std::istream& in,const MRC& formatmrc)
{
  BinaryIStream<CONVERSIONTYPE> f(in);
  HEADER header;
  f >> header;
  if(header.mode==5) header.mode=0;

  if(Logger::Instance().GetVerbosityLevel()>=4) {
    header.Print();
  }
  if(header.mode==3 || header.mode==4) {
    // always assume half-complex mode
    image.Reset(img::Size(header.nx,header.ny,header.nz),img::COMPLEX,img::HALF_FREQUENCY);
    if(img::image_state::ComplexHalfFrequencyImageState *cs=dynamic_cast<img::image_state::ComplexHalfFrequencyImageState*>(image.ImageStatePtr().get())) {
      if (header.mode==3) {
        detail::complex_filler<ushort,CONVERSIONTYPE>(*cs,f,header);
      } else if(header.mode==4) {
        detail::complex_filler<float,CONVERSIONTYPE>(*cs,f,header);
      }
    } else {
      throw IOException("internal error in MRC io: expected ComplexHalfFrequencyImageState");
    }
  } else if (header.mode==0 || header.mode==2 || header.mode==5 || header.mode==6 || header.mode==7) {
    img::Size msize;
    msize[header.mapc-1]=header.nc;
    msize[header.mapr-1]=header.nr;
    msize[header.maps-1]=header.ns;
    img::Point mstart;
    mstart[header.mapc-1]=header.ncstart;
    mstart[header.mapr-1]=header.nrstart;
    mstart[header.maps-1]=header.nsstart;
    image.Reset(img::Extent(mstart,msize),img::REAL,img::SPATIAL);
    if(header.x>0.0 && header.y >0.0 && header.z > 0){
      image.SetSpatialSampling(geom::Vec3(static_cast<Real>(header.x)/static_cast<Real>(header.nx),
                                          static_cast<Real>(header.y)/static_cast<Real>(header.ny),
                                          static_cast<Real>(header.z)/static_cast<Real>(header.nz)));
    }else{
      LOG_INFO("Suspicious cell dimensions found. Cannot set sampling.");
    }
    LOG_INFO("resulting image extent: " << image.GetExtent());
    if(img::image_state::RealSpatialImageState *rs=dynamic_cast<img::image_state::RealSpatialImageState*>(image.ImageStatePtr().get())) {
      if(header.mode==0) {
        detail::real_filler<int8_t,CONVERSIONTYPE>(*rs,f,header);
      } else if(header.mode==1 || header.mode==5) {
        detail::real_filler<int16_t,CONVERSIONTYPE>(*rs,f,header);
      } else if(header.mode==2) {
        detail::real_filler<float,CONVERSIONTYPE>(*rs,f,header);
      } else if(header.mode==6) {
        detail::real_filler<uint16_t,CONVERSIONTYPE>(*rs,f,header);
      } else if(header.mode==7) {
        detail::real_filler<int32_t,CONVERSIONTYPE>(*rs,f,header);
      }
    } else {
      throw IOException("internal error in MRC/CCP4 io: expected RealSpatialImageState");
    }
  } else {
    ::String message_string;
    std::ostringstream mesg(message_string);
    mesg << "MRC/CCP4 import: unknown data type: " << header.mode;
    throw IOException(mesg.str());
  }
}

template<class HEADER>
void import_endianess_switcher(img::MapHandle& image,
                             std::istream& f,
                             std::istream& header_str,
                             const MRC& formatmrc)
{
  switch(HEADER::DetermineDataFormat(header_str)){
  case OST_BIG_ENDIAN:
   import_helper<HEADER,OST_BIG_ENDIAN>(image,f,formatmrc);
   break;
  case OST_LITTLE_ENDIAN:
   import_helper<HEADER,OST_LITTLE_ENDIAN>(image,f,formatmrc);
   break;
  case OST_VAX_DATA:
   import_helper<HEADER,OST_VAX_DATA>(image,f,formatmrc);
   break;
  }
}

template<class HEADER,int CONVERSIONTYPE>
void export_helper(const img::MapHandle& image,
                                    std::ostream& out,
                                    const MRC& formatmrc)
{
  BinaryOStream<CONVERSIONTYPE> f(out);
  HEADER header(image,formatmrc.GetBitDepth());

  f << header;
  if(image.GetType()==img::REAL) {
    switch(formatmrc.GetBitDepth()){
    case OST_BIT8_FORMAT:
      detail::real_dumper<int8_t,CONVERSIONTYPE>(f,header,image,formatmrc);
      break;
    case OST_BIT16_FORMAT:
      detail::real_dumper<int16_t,CONVERSIONTYPE>(f,header,image,formatmrc);
      break;
    case OST_FLOAT_FORMAT:
    case OST_DEFAULT_FORMAT:
      detail::real_dumper<float,CONVERSIONTYPE>(f,header,image,formatmrc);
      break;
    default:
      throw(IOException("MRC/CCP4 export: Bit depth not supported."));
      break;
    }
  } else {
    if(image.GetDomain()==img::HALF_FREQUENCY){
      switch(formatmrc.GetBitDepth()){
        case OST_BIT16_FORMAT:
          detail::complex_dumper<int16_t,CONVERSIONTYPE>(f,header,image,formatmrc);
          break;
        case OST_FLOAT_FORMAT:
        case OST_DEFAULT_FORMAT:
          detail::complex_dumper<float,CONVERSIONTYPE>(f,header,image,formatmrc);
          break;
        default:
          throw(IOException("MRC/CCP4 export: Bit depth not supported."));
          break;
      }
    } else {
      throw(IOException("MRC/CCP4 export: full complex export not supported."));
    }
  }
  if(Logger::Instance().GetVerbosityLevel()>=4) {
    header.Print();
  }
}

template<class HEADER>
void export_endianess_switcher(const img::MapHandle& image,
                                std::ostream& f,
                                const MRC& formatmrc)
{
  switch(formatmrc.GetEndianessOnSave()){
  case OST_BIG_ENDIAN:
    export_helper<HEADER,OST_BIG_ENDIAN>(image,f,formatmrc);
    break;
  case OST_LITTLE_ENDIAN:
    export_helper<HEADER,OST_LITTLE_ENDIAN>(image,f,formatmrc);
    break;
  case OST_VAX_DATA:
    export_helper<HEADER,OST_VAX_DATA>(image,f,formatmrc);
    break;
  }
}

} //ns common

namespace bf = boost::filesystem;

void MapIOMrcHandler::Import(img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct )
{
  boost::filesystem::ifstream infile(loc, std::ios::binary);
  if(!infile)
  {
    throw IOException("could not open "+loc.string());
  }
  boost::iostreams::filtering_stream<boost::iostreams::input> in;
  if (detail::FilenameEndsWith(loc.string(),".map.gz")) {
    in.push(boost::iostreams::gzip_decompressor());
  }
  in.push(infile);
  in.read(reinterpret_cast<char*>(&header_),256);
  // seekg does not work on compressed streams
  // the only way to rewind is closing and reopening the stream
  // by popping out the file source and pushing it in again
  in.pop();
  infile.close();
  boost::filesystem::ifstream infile2(loc, std::ios::binary);
  if(!infile2)
  {
    throw IOException("could not open "+loc.string());
  }
  in.push(infile2);
  is_file_=true;
  filename_=loc.string();
  this->Import(sh,in,formatstruct);
  infile2.close();
}

void MapIOMrcHandler::Import(img::MapHandle& sh, std::istream& loc, const ImageFormatBase& formatstruct)
{
   MRC form;
   MRC& formatmrc = form;
   if (formatstruct.GetFormatString()==MRC::FORMAT_STRING) {
     formatmrc = formatstruct.As<MRC>();
   } else {
     assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
   }
   if (is_file_ == false){
     loc.read(reinterpret_cast<char*>(&header_),256);
     loc.seekg(0,std::ios::beg);
   }
   char* headerptr=(char*)&header_;
   boost::iostreams::filtering_streambuf<boost::iostreams::input> head_strbuf;
   std::istream head_str(&head_strbuf);
   head_strbuf.push(boost::iostreams::basic_array_source<char>(headerptr,sizeof(header_)));
   if (formatmrc.GetSubformat()==MRC_OLD_FORMAT) {
     LOG_DEBUG("mrc io: importing old style format");
     detail::import_endianess_switcher<detail::mrc_header>(sh,loc,head_str,formatmrc);
   } else if (formatmrc.GetSubformat()==MRC_NEW_FORMAT) {
     LOG_DEBUG("mrc io: importing new style format");
     detail::import_endianess_switcher<detail::ccp4_header>(sh,loc,head_str,formatmrc);
   } else if (is_file_ && (detail::FilenameEndsWith(filename_,".ccp4") || detail::FilenameEndsWith(filename_,".map") || detail::FilenameEndsWith(filename_,".map.gz"))) {
     LOG_DEBUG("mrc io: importing new style format");
     detail::import_endianess_switcher<detail::ccp4_header>(sh,loc,head_str,formatmrc);
   } else {
	 unsigned char header_content[256];
	 memcpy(&header_content[0],&header_,256*sizeof(char));
     if (MatchContent(header_content) == true) {
       LOG_DEBUG("mrc io: importing new style format");
       detail::import_endianess_switcher<detail::ccp4_header>(sh,loc,head_str,formatmrc);
     } else {
       LOG_DEBUG("mrc io: importing old style format");
       detail::import_endianess_switcher<detail::mrc_header>(sh,loc,head_str,formatmrc);
     }
  }
}

void MapIOMrcHandler::Export(const img::MapHandle& image,
                         const boost::filesystem::path& loc,const ImageFormatBase& formatstruct) const
{
  bf::ofstream outfile(loc, std::ios::binary);
  if(!outfile)
  {
    throw IOException("could not open "+loc.string());
  }
  is_file_=true;
  filename_=loc.string();
  Export(image,outfile,formatstruct);
  outfile.close();
}


void MapIOMrcHandler::Export(const img::MapHandle& sh, std::ostream& loc,const ImageFormatBase& formatstruct) const
{
  MRC form;
  MRC& formatmrc = form;
  if (formatstruct.GetFormatString()==MRC::FORMAT_STRING) {
    formatmrc = formatstruct.As<MRC>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }
  img::MapHandle sh2=sh;
  if(geom::Length2(sh.GetAbsoluteOrigin())>0.0) {
    /*
      CCP4/MRC format does not allow arbitrary absolute origin to be used, hence
      this correction routine to re-sample the map appropriately
    */
    LOG_VERBOSE("mrc io: applying absolute origin shift correction");
    geom::Vec3 frac_shift=geom::CompDivide(sh.GetAbsoluteOrigin(),sh.GetSpatialSampling());
    // apply fractional shift via Fourier space -> will wrap around!
    sh2=sh.Apply(img::alg::FractionalShift(frac_shift));
    // force origin to be zero now
    sh2.SetAbsoluteOrigin(geom::Vec3(0.0,0.0,0.0));
    // and correct for most of the wrap-around
    img::Point int_shift(frac_shift);
    sh2.ApplyIP(img::alg::Shift(-int_shift));
    sh2.SetSpatialOrigin(int_shift);
  }
  if (formatmrc.GetSubformat()==MRC_OLD_FORMAT) {
   LOG_DEBUG("mrc io: exporting old style format");
   detail::export_endianess_switcher<detail::mrc_header>(sh2,loc,formatmrc);
  } else {
   LOG_DEBUG("mrc io: exporting new style format");
   detail::export_endianess_switcher<detail::ccp4_header>(sh2,loc,formatmrc);
  }
}

bool MapIOMrcHandler::MatchContent(unsigned char* header)
{
  if(header[208]=='M' &&
     header[209]=='A' &&
     header[210]=='P'){
      return true;
  }
  return false;
}

bool MapIOMrcHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()=="defined_mrc") {
    return true;
  }
  return false;
}

bool MapIOMrcHandler::MatchSuffix(const String& loc)
{
  if(detail::FilenameEndsWith(loc,".mrc") || detail::FilenameEndsWith(loc,".map") ||
     detail::FilenameEndsWith(loc,".ccp4") || detail::FilenameEndsWith(loc,".map.gz")) {
     return true;
  }
  return false;
}


}} // ns
