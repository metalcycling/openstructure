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
#include <iostream>
#include <sstream>
#include <ost/message.hh>

#include <boost/shared_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <ost/log.hh>
#include <ost/img/util.hh>
#include <ost/img/image.hh>
#include <ost/img/image_state.hh>
#include <ost/img/progress.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/swap_util.hh>

#include "map_io_dm3_handler.hh"

namespace ost { namespace io {

String DM3::FORMAT_STRING = "defined_dm3";

DM3::DM3():
  ImageFormatBase(FORMAT_STRING)
{}

namespace bf = boost::filesystem;

namespace detail {

/*
  convert a tag entry from a void buffer that is of
  known type SOURCE to a given destination of type
  DEST
*/
template<typename SOURCE, typename DEST>
void value_convert(void* src_ptr, DEST* dest_ptr)
{
  SOURCE* r_src_ptr = reinterpret_cast<SOURCE*>(src_ptr);
  *dest_ptr = static_cast<DEST>(*r_src_ptr);
}

/*
  array version of tag_value_convert
*/
template<typename SOURCE, typename DEST>
void value_convert(void* src_ptr, DEST* dest_ptr, size_t count)
{
  SOURCE* r_src_ptr = reinterpret_cast<SOURCE*>(src_ptr);
  for(uint c=0;c<count;++c) {
    LOG_TRACE(" " << c << " / " << count << " .");
    DEST val=static_cast<DEST>(r_src_ptr[c]);
    LOG_TRACE(" ." << std::endl);
    dest_ptr[c] = val;
  }
}

template<typename DEST>
void tag_value_convert(int number_type, void* src_ptr, DEST* dest_ptr, size_t n=1)
{
  if(number_type==2) {
    value_convert<short,DEST>(src_ptr,dest_ptr,n);
  } else if(number_type==3) {
    value_convert<int,DEST>(src_ptr,dest_ptr,n);
  } else if(number_type==4) {
    value_convert<unsigned short,DEST>(src_ptr,dest_ptr,n);
  } else if(number_type==5) {
    value_convert<unsigned int,DEST>(src_ptr,dest_ptr,n);
  } else if(number_type==6) {
    value_convert<float,DEST>(src_ptr,dest_ptr,n);
  } else if(number_type==7) {
    value_convert<double,DEST>(src_ptr,dest_ptr,n);
  } else if(number_type==8) {
    value_convert<int,DEST>(src_ptr,dest_ptr,n);
  } else if(number_type==9) {
    value_convert<char,DEST>(src_ptr,dest_ptr,n);
  } else {
    throw IOException("I/O DM3: unsupported number type in tag_value_convert");
  }
}


template<typename SOURCE, typename DEST>
void img_data_convert2(void* src_ptr, img::RealSpatialImageState& isi)
{
  SOURCE* r_src_ptr = reinterpret_cast<SOURCE*>(src_ptr);
  int width=isi.GetSize()[0];
  int height=isi.GetSize()[1];
  for(int hi=0;hi<height;++hi) {
    for(int wi=0;wi<width;++wi) {
      isi.Value(img::Point(wi,hi))=static_cast<DEST>(r_src_ptr[hi*width+wi]);
    }
  }
}

template<typename DEST>
void img_data_convert(int number_type, void* src_ptr, img::RealSpatialImageState& isi)
{
  if(number_type==2) {
    img_data_convert2<short,DEST>(src_ptr,isi);
  } else if(number_type==3) {
    img_data_convert2<int,DEST>(src_ptr,isi);
  } else if(number_type==4) {
    img_data_convert2<unsigned short,DEST>(src_ptr,isi);
  } else if(number_type==5) {
    img_data_convert2<unsigned int,DEST>(src_ptr,isi);
  } else if(number_type==6) {
    img_data_convert2<float,DEST>(src_ptr,isi);
  } else if(number_type==7) {
    img_data_convert2<double,DEST>(src_ptr,isi);
  } else if(number_type==8) {
    img_data_convert2<int,DEST>(src_ptr,isi);
  } else if(number_type==9) {
    img_data_convert2<char,DEST>(src_ptr,isi);
  } else {
    throw IOException("I/O Dm3: unsupported number type in img_data_convert2");
  }
}

class DM3Collector
{
public:
  DM3Collector();

  // parse file
  void ParseFile(const bf::path& loc);
  void ParseStream(std::istream& is);
  img::ImageHandle RetrieveImage() const;

private:
  typedef std::vector<uint> TagTypeInfo;

  void parse_tag_entry(int depth,std::istream& fp);
  void parse_tag_group(int depth, const String& label,std::istream& fp);
  void parse_tag_type(int depth, const String& label,std::istream& fp);
  void parse_tag_type_single(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp);
  void parse_tag_type_String(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp);
  void parse_tag_type_array(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp);
  void parse_tag_type_simple_struct(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp);
  void parse_tag_type_array_struct(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp);
  void handle_tag_group(const String& name);
  void handle_single_tag(const String& name, int number_type, void* data);
  void handle_struct_tag(const String& name, std::vector<uint> number_type, std::vector<long> data);
  void handle_array_tag(const String& name, int number_type, boost::shared_array<uchar> data, uint size);
  void check_image();
  size_t file_read ( void * ptr, size_t size, size_t count, std::istream& stream );

  size_t type_sizeof(uint n);

  bool swap_tag_;
  bool swap_data_;

  std::vector<img::ImageHandle> image_handle_list_;
  boost::shared_array<uchar> image_data_;
  uint image_data_byte_size_;
  uint image_number_type_;
  int image_width_;
  int image_height_;
  bool expect_width_flag_;
  bool expect_height_flag_;
};

DM3Collector::DM3Collector():
  swap_tag_(true),
  swap_data_(false),
  image_handle_list_(),
  image_data_(),
  image_data_byte_size_(0),
  image_number_type_(0),
  image_width_(0),
  image_height_(0),
  expect_width_flag_(false),
  expect_height_flag_(false)
{}

size_t DM3Collector::file_read ( void * ptr, size_t size, size_t count, std::istream& stream )
{
    stream.read(reinterpret_cast<char*>(ptr),size*count);
    size_t s=stream.gcount();
    size_t pos=stream.tellg();
    img::Progress::Instance().AdvanceProgress(this,pos-0);
    return s;
}

void DM3Collector::ParseFile(const bf::path& loc)
{
  img::Progress::Instance().Register(this,bf::file_size(loc),1000);
  bf::ifstream ifile(loc, std::ios::binary);
  if(!ifile)
  {
    throw IOException("could not open "+loc.string());
  }
  LOG_DEBUG("Starting DM3 Parsing" << std::endl);
  this->ParseStream(ifile);
  ifile.close();
}

void DM3Collector::ParseStream(std::istream& fp)
{
  // top level entries
  uint version;
  file_read(&version,sizeof(uint),1,fp);
  uint bytecount;
  file_read(&bytecount,sizeof(uint),1,fp);
  uint endianess;
  file_read(&endianess,sizeof(uint),1,fp);

  LOG_INFO("version: " << version);
  LOG_INFO("bytecount: " << bytecount);
  LOG_INFO("endianess: " << endianess);

  parse_tag_group(0,"ROOT",fp);
  img::Progress::Instance().DeRegister(this);
  if(fp.eof()==0) {
    LOG_INFO("file not parsed completely");
  } else {
    LOG_INFO("file parsed completely");
  }

}

img::ImageHandle DM3Collector::RetrieveImage() const
{
  img::ImageHandle ih;
  int maxdim=0;
  for(std::vector<img::ImageHandle>::const_iterator it=image_handle_list_.begin();
      it!=image_handle_list_.end();++it) {
    if(it->GetExtent().GetVolume()>maxdim) {
      ih=*it;
      maxdim=ih.GetExtent().GetVolume();
    }
  }

  return ih;
}

void DM3Collector::handle_tag_group(const String& name)
{
  if(name=="Dimensions") {
    expect_width_flag_=true;
  }
}

void DM3Collector::handle_single_tag(const String& name, int number_type, void* data)
{
  if(expect_width_flag_) {
    tag_value_convert<int>(number_type,data,&image_width_);
    expect_width_flag_=false;
    expect_height_flag_=true;
  }
  if(expect_height_flag_) {
    tag_value_convert<int>(number_type,data,&image_height_);
    expect_width_flag_=false;
    expect_height_flag_=false;
    check_image();
  }
}

void DM3Collector::handle_struct_tag(const String& name, std::vector<uint> number_type, std::vector<long> data)
{
}

void DM3Collector::handle_array_tag(const String& name, int number_type, boost::shared_array<uchar> data, uint size)
{
  if(name=="Data") {
    image_data_=data;
    image_data_byte_size_=size;
    image_number_type_=number_type;
    check_image();
  }
}

void DM3Collector::check_image()
{
  if(image_data_ && image_width_>0 && image_height_>0) {
    LOG_VERBOSE("assembling image " << image_width_ << "x" << image_height_ << std::endl);

    // implicit row-major to colimn major conversion
    img::ImageHandle ih=CreateImage(img::Extent(img::Size(image_height_,image_width_)),img::REAL,img::SPATIAL);
    img::RealSpatialImageState* is = dynamic_cast<img::RealSpatialImageState*>(ih.ImageStatePtr().get());

    if(!is) {
      throw IOException("I/O Dm3: unexpected failure in dynamic_cast in check_image");
    }

    int data_size=image_width_*image_height_;

    if(data_size*type_sizeof(image_number_type_)==image_data_byte_size_) {

      img_data_convert<Real>(image_number_type_,image_data_.get(),*is);

      if(swap_data_)
      swap_buf<Real>(reinterpret_cast<Real*>(is->Data().GetData()),data_size);

      image_handle_list_.push_back(ih);
    } else {
      LOG_ERROR("mismatch in byte size of data (" 
                << data_size*type_sizeof(image_number_type_) << " != " 
                << image_data_byte_size_ << "), ignoring this image entry");
    }

    image_width_=0;
    image_height_=0;
    image_data_=boost::shared_array<uchar>();
  }
}

size_t DM3Collector::type_sizeof(uint n)
{
  if(n==2 || n==4) {
    return 2; // short or ushort
  } else if(n==3 || n==5) {
    return 4; // int or uint
  } else if(n==6) {
    return sizeof(float);
  } else if(n==7) {
    return sizeof(double);
  } else if(n==0) {
    throw IOException("I/O Dm3: unexpected number type 0");
  }
  // rest is one
  return 1;
}

void DM3Collector::parse_tag_group(int depth, const String& name,std::istream& fp)
{
  LOG_TRACE(String(depth,' ') << "parsing tag group [" << name << "]");

  uchar sorted;
  file_read(&sorted,sizeof(uchar),1,fp);

  uchar open;
  file_read(&open,sizeof(uchar),1,fp);

  LOG_TRACE(String(depth,' ') << "s,o: " << static_cast<int>(sorted) << "," 
            << static_cast<int>(open));

  uint tag_count;
  file_read(&tag_count,sizeof(uint),1,fp);
  if(swap_tag_) swap_uint(&tag_count,1);

  LOG_TRACE(String(depth,' ') << "tag_count: " << tag_count);

  handle_tag_group(name);

  for(uint tc=0; tc<tag_count;++tc) {
    parse_tag_entry(depth+1,fp);
  }
}

void DM3Collector::parse_tag_entry(int depth,std::istream& fp)
{
  LOG_TRACE(String(depth,' ') << "parsing tag entry");

  uchar id;
  file_read(&id,sizeof(uchar),1,fp);

  LOG_TRACE(String(depth,' ') << "id: " << static_cast<int>(id));

  if(!(id==20 || id==21)) {
    std::ostringstream msg("");
    msg << "I/O Dm3: unknown id [" << static_cast<int>(id) << "], expected 20 or 21";
    throw IOException(msg.str());
  }

  ushort nlabel;
  file_read(&nlabel,sizeof(ushort),1,fp);
  if(swap_tag_) swap_ushort(&nlabel,1);

  String label(nlabel,'\0');
  if(nlabel>0) {
    std::vector<char> buffer(nlabel);
    file_read(&buffer[0],sizeof(uchar),nlabel,fp);
    label=String(&buffer[0],nlabel);
  }

  LOG_TRACE(String(depth,' ') << "nlabel: " << nlabel << " label: '" 
            << label << "'" << std::endl);

  if(id==static_cast<uchar>(20)) {
    parse_tag_group(depth+1,label,fp);
  } else { // 21
    parse_tag_type(depth,label,fp);
  }

}

void DM3Collector::parse_tag_type(int depth, const String& name, std::istream& fp)
{
  LOG_TRACE(String(depth,' ') << "parsing tag type");

  char dummy[4];
  file_read(dummy,sizeof(char),4,fp);

  if(String(dummy,4)!=String("%%%%")) {
    throw IOException("I/O Dm3: invalid tag_type magic number");
  } else {
    LOG_TRACE(String(depth,' ') << "magic number %%%% found");
  }

  uint ntype;
  file_read(&ntype,sizeof(uint),1,fp);
  if(swap_tag_) swap_uint(&ntype,1);

  LOG_TRACE(String(depth,' ') << "ntype: " << ntype);

  std::vector<uint> tag_type_info(ntype);
  file_read(&tag_type_info[0],sizeof(uint),ntype,fp);
  if(swap_tag_) swap_uint(&tag_type_info[0],ntype);

  for(uint i=0;i<ntype;++i) {
    LOG_TRACE(String(depth,' ') << "type_info[" << i << "]: " << tag_type_info[i]);
  }

  if(ntype==1) { // simple type
    parse_tag_type_single(depth,name,tag_type_info,fp);
  } else if(ntype==2) { // String type
    parse_tag_type_String(depth,name,tag_type_info,fp);
  } else if(ntype==3) { // array of simple type
    parse_tag_type_array(depth,name,tag_type_info,fp);
  } else { // more complicated layout
    if(tag_type_info[0]==15) {
      parse_tag_type_simple_struct(depth,name,tag_type_info,fp);
    } else if(tag_type_info[0]==20) {
      parse_tag_type_array_struct(depth,name,tag_type_info,fp);
    }
  }
}


void DM3Collector::parse_tag_type_single(int depth, const String& name,TagTypeInfo& tag_type_info,std::istream& fp)
{
  LOG_TRACE(String(depth,' ') << "parsing simple tag [" << name << "]");

  // tag_type_info[0] contains number type
  uint number_type=tag_type_info[0];
  if(number_type>9) {
    std::ostringstream msg("");
    msg << "I/O Dm3: unexpected number type " << number_type << " for simple type";
    throw IOException(msg.str());
  } else {
    LOG_TRACE(String(depth,' ') << "reading simple value of number_type " << number_type);
  }

  uchar buffer[16]; // sufficient for all primitive types
  size_t bufsize = type_sizeof(number_type);
  file_read(buffer,bufsize,1,fp);
  // TODO: swap data
  double tmp;
  tag_value_convert<double>(number_type,buffer,&tmp);
  LOG_TRACE(String(depth,' ') << "value: " << tmp);

  handle_single_tag(name,number_type,buffer);
}

void DM3Collector::parse_tag_type_String(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp)
{
  LOG_TRACE(String(depth,' ') << "parsing String tag [" << name << "]");
  // tag_type_info[0] contains number type 18
  // tag_type_info[1] contains String length
  std::vector<char> buffer(tag_type_info[1]*2);
  file_read(&buffer[0],sizeof(char)*2,tag_type_info[1],fp); // unicode Strings
}

void DM3Collector::parse_tag_type_array(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp)
{
  LOG_TRACE(String(depth,' ') << "parsing array tag [" << name << "]");
  // tag_type_info[0] contains number type 20
  // tag_type_info[1] contains array number type
  // tag_type_info[2] contains array length
  uint number_type=tag_type_info[1];
  if(number_type>10) {
    std::ostringstream msg("");
    msg << "I/O Exception: unexpected number type " << number_type << " for array";
    throw IOException(msg.str());
  } else {
    LOG_TRACE(String(depth,' ') << "parsing array of number_type " << number_type);
  }

  size_t array_type_size = type_sizeof(number_type);

  LOG_TRACE(String(depth,' ') << "array byte size: " << array_type_size*tag_type_info[2]);

  //std::vector<uchar> data(array_type_size*tag_type_info[2]);
  boost::shared_array<uchar> array_data(new uchar[array_type_size*tag_type_info[2]]);
  file_read(array_data.get(),array_type_size,tag_type_info[2],fp);
  // TODO: swap data
  handle_array_tag(name,number_type,array_data,array_type_size*tag_type_info[2]);
}

void DM3Collector::parse_tag_type_simple_struct(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp)
{
  LOG_TRACE(String(depth,' ') << "parsing simple struct tag [" << name << "]");
  // tag_type_info[0] contains number type 15
  if(tag_type_info[0]!=15) {
    std::ostringstream msg("");
    msg << "I/O Dm3: expected number type [15], but found [" << tag_type_info[0] << "] in tag_type_info[0] for tag_type_simple_struct";
    throw IOException(msg.str());
  }

  // tag_type_info[2] contains number of entries in struct
  LOG_TRACE(String(depth,' ') << "reading in " << tag_type_info[2] << " values");
  std::vector<uint> number_types(tag_type_info[2]);
  std::vector<long> buffers(tag_type_info[2]);
  for(uint c=0;c<tag_type_info[2];++c) {
    number_types[c]=tag_type_info[c*2+4];
    if(number_types[c]>10) {
      std::ostringstream msg("");
      msg << "I/O Dm3: unexpected number type " << number_types[c] << " while parsing tag_type_struct";
      throw IOException(msg.str());
    } else {
      LOG_TRACE("reading value of number_type " << number_types[c]);
    }
    size_t bufsize = type_sizeof(number_types[c]);
    void* buffer = &buffers[c];
    file_read(buffer,bufsize,1,fp);
    double tmp;
    tag_value_convert<double>(number_types[c],buffer,&tmp);
    LOG_TRACE(String(depth,' ') << "value " << c << ": " << tmp);
    // TODO: swap data
  }
  handle_struct_tag(name,number_types,buffers);

}

void DM3Collector::parse_tag_type_array_struct(int depth, const String& name, TagTypeInfo& tag_type_info,std::istream& fp)
{
  LOG_DEBUG(String(depth,' ') << "parsing array struct tag [" << name << "]");
  // tag_type_info[0] contains number type 20
  if(tag_type_info[0]!=20) {
    std::ostringstream msg("");
    msg << "I/O Dm3: expected number type [20], but found [" << tag_type_info[0] << "] in tag_type_info[0] for tag_type_array_struct";
    throw IOException(msg.str());
  }

  // tag_type_info[size-1] contains length of array!
  uint array_size = tag_type_info[tag_type_info.size()-1];

  LOG_TRACE(String(depth,' ') << "looping over " << array_size 
            << " array elements");

  uchar buffer[16]; // sufficient for all primitive types
  for(uint ac=0;ac<array_size;++ac) {
    // tag_type_info[3] contains number of entries in struct
    LOG_TRACE(String(depth,' ') << "reading in " << tag_type_info[3] << " values for array entry ");
    for(uint c=0;c<tag_type_info[3];++c) {
      uint number_type=tag_type_info[c*2+5];
      if(number_type>10) {
        std::ostringstream msg;
        msg << "I/O Dm3: unexpected number type " << number_type 
            << " while parsing tag_type_struct";
        throw IOException(msg.str());
      } else {
        LOG_TRACE(String(depth, ' ') << "reading value of number_type " << number_type);
      }
      size_t bufsize = type_sizeof(number_type);
      file_read(buffer,bufsize,1,fp);
      // TODO: swap data
    }
  }

}

} // ns detail

void MapIODm3Handler::Import(img::MapHandle& mh, const bf::path& loc,const ImageFormatBase& formatstruct)
{
    detail::DM3Collector collector;
    collector.ParseFile(loc);
    mh=collector.RetrieveImage();
}

void MapIODm3Handler::Import(img::MapHandle& mh, std::istream& is,const ImageFormatBase& formatstruct)
{
  DM3 form;
  DM3& formatdat = form;
  if (formatstruct.GetFormatString()==DM3::FORMAT_STRING) {
    formatdat = formatstruct.As<DM3>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  detail::DM3Collector collector;
  collector.ParseStream(is);
  mh=collector.RetrieveImage();
}

void MapIODm3Handler::Export(const img::MapHandle& sh, const bf::path& loc, const ImageFormatBase& formatstruct) const
{
  throw IOException("This file format is not supported for export");
}

void MapIODm3Handler::Export(const img::MapHandle& sh, std::ostream& loc, const ImageFormatBase& formatstruct) const
{
  throw IOException("This file format is not supported for export");
}


bool MapIODm3Handler::MatchContent(unsigned char* header)
{
   if(
     header[0]==0x00 &&
     header[1]==0x00 &&
     header[2]==0x00 &&
     header[3]==0x03 &&
     (header[18]==0x14 || header[18]==0x15)
     ){
      return true;
  }
  return false;
}
bool MapIODm3Handler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()=="defined_dm3") {
    return true;
  }
  return false;
}
bool MapIODm3Handler::MatchSuffix(const String& loc)
{
	if(detail::FilenameEndsWith(loc,".dm3") ) {
      return true;
    }
    return false;
}

}} // ns
