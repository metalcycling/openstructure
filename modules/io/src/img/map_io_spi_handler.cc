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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <cstdio>
#include <fstream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/scoped_array.hpp>
#include <boost/format.hpp>
#include <ost/img/progress.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/log.hh>
#include <ost/io/converting_streams.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/swap_util.hh>
#include <ost/img/alg/stat.hh>

#include "map_io_spi_handler.hh"

namespace bf = boost::filesystem;

namespace ost { namespace io {

String Spider::FORMAT_STRING = "defined_spider";

Spider::Spider(bool normalize_on_save,Endianess endianness_on_save):
  ImageFormatBase(FORMAT_STRING),
  normalize_on_save_(normalize_on_save),
  endianness_on_save_(endianness_on_save)
{
}

Endianess Spider::GetEndianessOnSave() const
{
  return endianness_on_save_;
}

void Spider::SetEndianessOnSave(Endianess endianness_on_save)
{
  endianness_on_save_ = endianness_on_save;
}

bool Spider::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void Spider::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}
namespace detail {

using boost::format;

class spider_header {
public:
  spider_header():
    fNslice(),fNrow(),fNrec(),fNlabel(),fIform(),fImami(),fFmax(),fFmin(),fAv(),fSig(),
    fIhist(),fNcol(),fLabrec(),fIangle(),fPhi(),fTheta(),fPsi(),fXoff(),fYoff(),fZoff(),
    fScale(),fLabbyt(),fLenbyt(),ffIstack(),fFlag(),fMaxim(),fImgnum(),fLastindx(),
    fUn1(),fUn2(),Kangle(),fPhi1(),fTheta1(),fPsi1(),fPhi2(),fTheta2(),fPsi2(),
    fEmpty1(),fRes(),fEmpty2(),szIDat(),szITim(),szITit()
  {
    for(unsigned i=0;i<12;++i){
      fEmpty1[i]=0;
    }
    for(unsigned i=0;i<27;++i){
      fRes[i]=0;
    }
    for(unsigned i=0;i<135;++i){
      fEmpty2[i]=0;
    }
    for(unsigned i=0;i<12;++i){
      szIDat[i]=0;
    }
    for(unsigned i=0;i<8;++i){
      szITim[i]=0;
    }
    for(unsigned i=0;i<160;++i){
      szITit[i]=0;
    }
  }

  float fNslice;  // nr of slices (1 for an image, negative nr of slices indicates very old 2d image
  float fNrow;    // nr row per slice (Y)
  float fNrec;    // total nr of records (unused).
  float fNlabel;  // obsolete, unused
  float fIform;   // file type specifier:
                  // 1   (r)     2D image.
                  // 3   (r)     3D volume.
                  // -11     (fo)    2D Fourier, mixed radix odd.
                  // -12     (fe)    2D Fourier, mixed radix even.
                  // -21     (fo)    3D Fourier, mixed radix odd.
                  // -22     (fe)    3D Fourier, mixed radix even.
  float fImami;   // max/min flag. Is  Set at 0 at file creation
                  // and at 1 when max min have been computed and stored
  float fFmax;    // max
  float fFmin;    // min
  float fAv;      // avg
  float fSig;     // Std dev -1 if not computed previously
  float fIhist;   // flag indicating if histogram has been computed (not used in 3d files)
  float fNcol;    // nr pixel per line (Columns X)
  float fLabrec;  // nr label records in file header
  float fIangle;  // flag indicating that tilt angles have been filled
  float fPhi;     // euler: rot angle
  float fTheta;   // euler: tilt angle
  float fPsi;     // euler: psi tilt angle
  float fXoff;    // X transl
  float fYoff;    // y transl
  float fZoff;    // z transl
  float fScale;   // scale
  float fLabbyt;  // nr bytes in label
  float fLenbyt;  // record length in bytes
  float ffIstack; // stack flag 0=single 2D or 3D image. Other=stack
  float fFlag;    // obsolete, unused
  float fMaxim;   // used for stacks
  float fImgnum;  // used for stacks
  float fLastindx;// used for stacks
  float fUn1;     // unused
  float fUn2;     // unused
  float Kangle;   // flag inidcating additional rotations are present in
                   // the heaser 1=additional rotation 2=additional
                   // rotations that precedes the one stored in fPhi,
                   // fTheta, fPsi, fXoff, fYoff, Zoff;
  float fPhi1;    // euler: rot angle
  float fTheta1;  // euler: tilt angle
  float fPsi1;    // euler: psi tilt angle
  float fPhi2;    // euler: rot angle
  float fTheta2;  // euler: tilt angle
  float fPsi2;    // euler: psi tilt angle

  float fEmpty1[12]; //
  float fRes[27];    //
  float fEmpty2[135];//
  char szIDat[12];   // array sized eleven, date of creation (0 if absent)
  char szITim[8];    // array sized eight, time of creation (0 if absent)
  char szITit[160];  // array dimensioned 160, title
};

void prep_header(spider_header& header, const img::Size& size, const geom::Vec3& spacing, float min, float max, float mean, float std_dev)
{
  if (spacing[0] != spacing[1] ||
      spacing[1] != spacing[2] ||
      spacing[0] != spacing[2] )  {
    throw(ost::io:: IOException("This format does not support anisotropic pixel sampling"));
  }

  uint ncol = size[0];
  uint nrow = size[1];
  uint nslice = size[2];
  float spatial_sampling = spacing[0];
  time_t time_value;
  time(&time_value);
  String time_string=ctime(&time_value);
  String date_for_header=time_string.substr(8,2) + "-" + time_string.substr(4,3) + "-" +
                         time_string.substr(22,4) + " ";
  String time_for_header=time_string.substr(11,8);
  String title_for_header = "Generated by OST (www.openstructure.org)";

  header.fNslice = nslice;
  header.fNrow = nrow;
  header.fNrec = 1.0;
  if (nslice == 1) {
    header.fIform=1.0;
  } else {
    header.fIform=3.0;
  }
  header.fImami = 1.0;
  header.fFmax = max;
  header.fFmin = min;
  header.fAv = mean;
  header.fSig = std_dev;
  header.fNcol =  ncol;
  header.fLenbyt = ncol*4.0;  // record length in bytesS
  header.fLabrec = ceil(1024.0 / header.fLenbyt);  // nr label records in file header
  header.fLabbyt = header.fLabrec * header.fLenbyt;
  header.fIangle = 0.0;  // flag indicating that tilt angles have been filled
  header.fScale = spatial_sampling;   // scale
  header.ffIstack = 0.0; // stack flag 0=single 2D or 3D image. Other=stack
  header.Kangle = 0.0;
  for (int counter = 0; counter < 12; ++counter)  {
    header.szIDat[counter] = date_for_header[counter];
  }
  for (int counter = 0; counter < 8; ++counter)  {
    header.szITim[counter] = time_for_header[counter];
  }
  for (size_t counter = 0; counter < title_for_header.size(); ++counter)  {
    header.szITit[counter] = title_for_header[counter];
  }
}

void print_header(const spider_header& h)
{
  std::ostringstream s;
  s << "spider header: " << std::endl;
  s << format(" columns,rows,slices:  %1% %2% %3%") % h.fNslice  % h.fNcol % h.fNrow << std::endl;
  if (h.fImami == 1) {
    s << format(" man,min,avg: %1% %2% %3%") %  h.fFmax % h.fFmin % h.fAv << std::endl;
  }
  if (h.fSig != -1) {
    s << format(" std dev: %1%") % h.fSig << std::endl;
  }
  s << format(" scale: %1%") %  h.fScale << std::endl;
  s << format(" header size:  %1%") % h.fLabbyt;
  LOG_INFO(s.str());
}

int determine_byte_order(std::istream& f)
{
  float iform;
  BinaryIStream<OST_LITTLE_ENDIAN> bs(f);
  bs.seekg(4*sizeof(float),std::ios::beg); // seek to iform
  bs >> iform;
  bs.seekg(0,std::ios::beg); // seek to beginning
  if (iform<-22.0 ||iform>3.0 ||static_cast<int>(iform)==0) {
    // try byte swapping
    float  bs_iform=Convert<OST_BIG_ENDIAN,float>::From(iform);
    if (bs_iform<-22.0 ||bs_iform>3.0 ||static_cast<int>(bs_iform)==0) {
      throw(ost::io::IOException("error reading spider file (even tried byte-swapping"));
    }
    LOG_VERBOSE("Spider IO: reading big endian data");
    return OST_BIG_ENDIAN;
  }
  LOG_VERBOSE("Spider IO: reading little endian data");
  return OST_LITTLE_ENDIAN;
}

template<int CONVERSIONTYPE>
void header_filler(std::istream& in,  spider_header& header)
{
  BinaryIStream<CONVERSIONTYPE> f(in);
  f >> header.fNslice;  
  if(header.fNslice<1.0){
    header.fNslice=1.0; // fix for very old 2d images with Nslice=-1
  }
  f >> header.fNrow;  
  f >> header.fNrec;    
  f >> header.fNlabel; 
  f >> header.fIform;   
  f >> header.fImami;  
  f >> header.fFmax;  
  f >> header.fFmin;  
  f >> header.fAv;      
  f >> header.fSig;    
  f >> header.fIhist;  
  f >> header.fNcol;    
  f >> header.fLabrec;
  f >> header.fIangle; 
  f >> header.fPhi;   
  f >> header.fTheta; 
  f >> header.fPsi;    
  f >> header.fXoff;  
  f >> header.fYoff;  
  f >> header.fZoff;    
  f >> header.fScale;  
  f >> header.fLabbyt; 
  f >> header.fLenbyt;
  f >> header.ffIstack;
  f >> header.fFlag;   
  f >> header.fMaxim;  
  f >> header.fImgnum; 
  f >> header.fLastindx;
  f >> header.fUn1;   
  f >> header.fUn2;   
  f >> header.Kangle;  
  f >> header.fPhi1;    
  f >> header.fTheta1; 
  f >> header.fPsi1;    
  f >> header.fPhi2;   
  f >> header.fTheta2; 
  f >> header.fPsi2;   
  f.read(header.fEmpty1,12);
  f.read(header.fRes,27);   
  f.read(header.fEmpty2,135);
  f.read(header.szIDat,12);   
  f.read(header.szITim,8);  
  f.read(header.szITit,160);  


  if(header.fIform!=3.0  && header.fIform!=1.0) {
  throw(ost::io::IOException("spiderRead: can only read modes 1 (2D image) or 3 (3D volume)\n"));
  }

}

template<int CONVERSIONTYPE>
void header_dumper(std::ostream& out, const spider_header& header)
{
  BinaryOStream<CONVERSIONTYPE> f(out);
  f << header.fNslice;  
  f << header.fNrow;  
  f << header.fNrec;    
  f << header.fNlabel; 
  f << header.fIform;   
  f << header.fImami;  
  f << header.fFmax;  
  f << header.fFmin;  
  f << header.fAv;      
  f << header.fSig;    
  f << header.fIhist;  
  f << header.fNcol;    
  f << header.fLabrec;
  f << header.fIangle; 
  f << header.fPhi;   
  f << header.fTheta; 
  f << header.fPsi;    
  f << header.fXoff;  
  f << header.fYoff;  
  f << header.fZoff;    
  f << header.fScale;  
  f << header.fLabbyt; 
  f << header.fLenbyt;
  f << header.ffIstack;
  f << header.fFlag;   
  f << header.fMaxim;  
  f << header.fImgnum; 
  f << header.fLastindx;
  f << header.fUn1;   
  f << header.fUn2;   
  f << header.Kangle;  
  f << header.fPhi1;    
  f << header.fTheta1; 
  f << header.fPsi1;    
  f << header.fPhi2;   
  f << header.fTheta2; 
  f << header.fPsi2;   
  f.write(header.fEmpty1,12);
  f.write(header.fRes,27);   
  f.write(header.fEmpty2,135);
  f.write(header.szIDat,12);   
  f.write(header.szITim,8);  
  f.write(header.szITit,160);  
}


template <int CONVERSIONTYPE,typename B>
void real_filler(std::istream& in, const spider_header& header, img::ImageHandle& mh,double scale)
{
  BinaryIStream<CONVERSIONTYPE> f(in);

  // size to skip forward to
  int size =( static_cast<int>(header.fNcol) *  static_cast<int>(header.fLabrec) * 4); //Size of whole header
  f.seekg(size,std::ios::beg); // move position of original stream to beginning of data

  char this_dummy; //create dummy variable to give to Progress as this
  img::Progress::Instance().Register(&this_dummy,static_cast<int>(header.fNslice*header.fNrow),100);
  for(int w=0;w<static_cast<int>(header.fNslice);++w) {

   int rawp_index=0;

   for(int v=0;v<static_cast<int>(header.fNrow);++v) {
      for(int u=0;u<static_cast<int>(header.fNcol);++u) {
        float buffer;
        f >> buffer;
        mh.SetReal(img::Point(u,v,w),static_cast<double>(buffer));
        ++rawp_index;
      }
      img::Progress::Instance().AdvanceProgress(&this_dummy);
    }
  }
 img:: Progress::Instance().DeRegister(&this_dummy);
}

template <typename B >
void real_dumper(std::ostream& f,  const spider_header& header, const img::ImageHandle& mh,const img::alg::Normalizer& norm, bool swap_flag)
{
  int padding = header.fLabbyt-f.tellp();
  std::vector<char> buffer(padding,0);
  f.write(&buffer[0],padding);
  int slice_size=static_cast<int>(header.fNcol) * static_cast<int>(header.fNrow);
  boost::scoped_array<B> rawp(new B[slice_size]);

  char this_dummy; //create dummy variable to give to Progress as this
  img::Progress::Instance().Register(&this_dummy,static_cast<int>(header.fNslice*header.fNrow),100);
  for(int w=0;w< static_cast<int>(header.fNslice);++w) {
    int rawp_index=0;
    for(int v=0;v< static_cast<int>(header.fNrow);++v) {
      for(int u=0;u< static_cast<int>(header.fNcol);++u) {
        rawp[rawp_index]=static_cast<B>(norm.Convert(mh.GetReal(mh.GetSpatialOrigin() + img::Point(u,v,w))));;
        ++rawp_index;
      }
      img::Progress::Instance().AdvanceProgress(&this_dummy);
    }

    if(swap_flag)  {
      ost::io::swap_buf<B>(&rawp[0],slice_size);
    }

    f.write(reinterpret_cast<char*>(&rawp[0]),sizeof(B)*slice_size);
  }
  img::Progress::Instance().DeRegister(&this_dummy);
}

template<int CONVERSIONTYPE>
void ExportSpiderHelper(const img::MapHandle& mh2, std::ostream& outfile,const Spider& formatspider)
{

  bool swap_flag=false;

  img::alg::Stat stat;
  mh2.Apply(stat);

  if (CONVERSIONTYPE==OST_BIG_ENDIAN && OST_LOCAL_ENDIAN==OST_LITTLE_ENDIAN) {
    swap_flag = true;
  }

  if (CONVERSIONTYPE==OST_LITTLE_ENDIAN && OST_LOCAL_ENDIAN==OST_BIG_ENDIAN) {
    swap_flag = true;
  }

  detail::spider_header header;

  detail::prep_header(header,mh2.GetSize(),mh2.GetSpatialSampling(),stat.GetMinimum(),stat.GetMaximum(),stat.GetMean(),stat.GetStandardDeviation());

  header_dumper<CONVERSIONTYPE>(outfile, header);

  detail::print_header(header);

  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatspider.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(mh2,formatspider.GetMinimum(),formatspider.GetMaximum());
  }

  detail::real_dumper<float>(outfile,header,mh2,norm,swap_flag);
}

} // ns detail

void MapIOSpiHandler::Import(img::MapHandle& mh, const bf::path& loc,const ImageFormatBase& formatstruct)
{
  bf::ifstream infile(loc, std::ios::binary);
  if(!infile)
  {
    throw IOException("could not open "+loc.string());
  }
  Import(mh,infile,formatstruct);
  infile.close();
}

void MapIOSpiHandler::Import(img::MapHandle& mh, std::istream& infile,const ImageFormatBase& formatstruct)
{
  Spider form;
  Spider& formatspider = form;
  if (formatstruct.GetFormatString()==Spider::FORMAT_STRING) {
    formatspider = formatstruct.As<Spider>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  String line;
  detail::spider_header header;
  switch(detail::determine_byte_order(infile)){
  case OST_BIG_ENDIAN:
    detail::header_filler<OST_BIG_ENDIAN>(infile, header);
    break;
  case OST_LITTLE_ENDIAN:
    detail::header_filler<OST_LITTLE_ENDIAN>(infile, header);
    break;
  default:
    throw(IOException("unknown byte order"));
    break;
  }

  detail::print_header(header);

  img::Size msize;
  msize[0]=static_cast<unsigned int>(header.fNcol);
  msize[1]=static_cast<unsigned int>(header.fNrow);
  msize[2]=static_cast<unsigned int>(header.fNslice);

  mh.Reset(msize);

  LOG_INFO("density map extent: " << mh.GetExtent());


  switch(detail::determine_byte_order(infile)){
  case OST_BIG_ENDIAN:
    detail::real_filler<OST_BIG_ENDIAN,float>(infile,header,mh,1.0);
    break;
  case OST_LITTLE_ENDIAN:
    detail::real_filler<OST_LITTLE_ENDIAN,float>(infile,header,mh,1.0);
    break;
  default:
    throw(IOException("unknown byte order"));
    break;
  }

  if (header.fScale==0.0)  {
    mh.SetSpatialSampling(geom::Vec3(1.0,1.0,1.0));
    LOG_INFO("assuming 1 A spatial sampling");
  } else {
    Real fscale=static_cast<Real>(header.fScale);
    mh.SetSpatialSampling(geom::Vec3(fscale,fscale,fscale));
  }

  mh.SetAbsoluteOrigin(geom::Vec3(0,0,0));

  LOG_INFO("assuming origin at the first voxel");
  LOG_INFO("assuming absolute origin at (0,0,0)");

}

void MapIOSpiHandler::Export(const img::MapHandle& mh2,
                                  const bf::path& loc,const ImageFormatBase& formatstruct) const
{
  bf::ofstream outfile(loc, std::ios::binary);
  if(!outfile)
  {
    throw IOException("could not open "+loc.string());
  }
  this->Export(mh2,outfile,formatstruct);
  outfile.close();
}

void MapIOSpiHandler::Export(const img::MapHandle& mh2,
                                  std::ostream& outfile,const ImageFormatBase& formatstruct) const
{
  Spider form;
  Spider& formatspider = form;
  if (formatstruct.GetFormatString()==Spider::FORMAT_STRING) {
    formatspider = formatstruct.As<Spider>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  switch (formatspider.GetEndianessOnSave()) {
  case OST_BIG_ENDIAN:
    detail::ExportSpiderHelper<OST_BIG_ENDIAN>(mh2,outfile,formatspider);
    break;
  case OST_LITTLE_ENDIAN:
    detail::ExportSpiderHelper<OST_LITTLE_ENDIAN>(mh2,outfile,formatspider);
    break;
  case OST_VAX_DATA:
    throw IOException ("Spider IO: Vax data format is not supported for exporting");
    break;
  default:  
    detail::ExportSpiderHelper<OST_LOCAL_ENDIAN>(mh2,outfile,formatspider);
    break;
  }  
}

bool MapIOSpiHandler::MatchContent(unsigned char* header)
{
  return false;
}
bool MapIOSpiHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()=="defined_spider") {
    return true;
  }
  return false;
}
bool MapIOSpiHandler::MatchSuffix(const String& loc)
{
    if (detail::FilenameEndsWith(loc,".spi") ) {
      return true;
    }
    return false;
}

}} //namespaces
