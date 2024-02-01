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

#include <png.h>

#include <boost/filesystem/fstream.hpp>

#include <ost/img/normalizer_impl.hh>
#include <ost/img/raster_image/raster_image.hh>
#include <ost/img/image.hh>
#include <ost/img/progress.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/io/io_exception.hh>

#include "map_io_png_handler.hh"


namespace ost { namespace io {

String PNG::FORMAT_STRING="defined_png";

PNG::PNG(bool normalize_on_save):
    ImageFormatBase(FORMAT_STRING),
    normalize_on_save_(normalize_on_save)
{
}

bool PNG::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void PNG::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}

Real PNG::GetMaximum() const
{
  return 255.0;
}
Real PNG::GetMinimum() const
{
  return 0.0;
}

namespace detail {

void user_read_data(png_structp pngPtr, png_bytep data, png_size_t length) {
    png_voidp a = png_get_io_ptr(pngPtr);
    ((std::istream*)a)->read((char*)data, length);
}

void user_write_data(png_structp pngPtr, png_bytep data, png_size_t length) {
    png_voidp a = png_get_io_ptr(pngPtr);
    ((std::ostream*)a)->write((char*)data, length);
}

void user_flush_data(png_structp pngPtr) {
   png_voidp a = png_get_io_ptr(pngPtr);
   ((std::ostream*)a)->flush();
}


} // ns

void MapIOPngHandler::Import(img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct )
{
  boost::filesystem::ifstream infile(loc, std::ios::binary);
  if(!infile)
  {
    throw IOException("could not open "+loc.string());
  }
  this->Import(sh,infile,formatstruct);
  infile.close();
}


void MapIOPngHandler::Import(img::MapHandle& image, std::istream& f,const ImageFormatBase& formatstruct)
{
  png_structp png_ptr;
  png_infop info_ptr;
  png_uint_32 channels;
  int png_transforms;
  png_bytep *row_pointers;
  png_uint_32 row,column;
  png_byte *png_data;
  png_uint_32 width,height;
  std::ostringstream mesg;

  PNG form;
  PNG& formatpng = form;
  if (formatstruct.GetFormatString()==PNG::FORMAT_STRING) {
    formatpng = formatstruct.As<PNG>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }



  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
           NULL,NULL,NULL);
  if (png_ptr == NULL) {
    throw IOException("png: error creating read_struct");
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
    throw IOException("png: error creating info_struct");
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
    throw IOException("png: error setting jmpbuf");
  }

  png_set_read_fn(png_ptr,&f,detail::user_read_data);

  png_transforms=PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING;
  png_read_png(png_ptr, info_ptr, png_transforms, NULL);

  channels=png_get_channels(png_ptr,info_ptr);

  if(!(channels==1 || channels==2)) {
    throw IOException("png: error: expected greyscale image");
  }

  // from here on channels is always 1 or 2, but the alpha value is ignored

  width=(unsigned int)png_get_image_width(png_ptr, info_ptr);
  height=(unsigned int)png_get_image_height(png_ptr, info_ptr);

  // TODO 16 bit greyscale


  png_data= (png_byte *)malloc(width*height*sizeof(png_byte));

  // default origin is top-left
  image.Reset(img::Extent(img::Point(0,0),img::Size(width,height)), img::REAL, img::SPATIAL);
  img::image_state::RealSpatialImageState *isr=dynamic_cast<img::image_state::RealSpatialImageState*>(image.ImageStatePtr().get());;

  row_pointers = png_get_rows(png_ptr, info_ptr);

  img::Progress::Instance().Register(this,height,100);
  for(row=0;row<height;row++) {
    for(column=0;column<width;column++) {
      isr->Value(img::Index(column,row,0))=((Real)(row_pointers[row])[column*channels]);
    }
    img::Progress::Instance().AdvanceProgress(this);
  }
  img::Progress::Instance().DeRegister(this);
  free(png_data);
  png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);

}


void MapIOPngHandler::Export(const img::MapHandle& image,
                         const boost::filesystem::path& loc,const ImageFormatBase& formatstruct) const
{
  boost::filesystem::ofstream outfile(loc, std::ios::binary);
  if(!outfile)
  {
    throw IOException("could not open "+loc.string());
  }

  Export(image,outfile,formatstruct);
  outfile.close();
}


void MapIOPngHandler::Export(const img::MapHandle& image, std::ostream& f,const ImageFormatBase& formatstruct) const
{
  png_structp png_ptr;
  png_infop info_ptr;
  std::ostringstream mesg;

  PNG form;
  PNG& formatpng = form;
  if (formatstruct.GetFormatString()==PNG::FORMAT_STRING) {
    formatpng = formatstruct.As<PNG>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatpng.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(image,formatpng.GetMinimum(),formatpng.GetMaximum());
  }

  img::Size size=image.GetSize();

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
            NULL, NULL, NULL);
  if (png_ptr == NULL) {
    throw IOException("png: error creating write_struct");
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
    throw IOException("png: error creating info_struct");
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    throw IOException("png: error setting jmpbuf");
  }

  png_set_write_fn(png_ptr,&f,detail::user_write_data,detail::user_flush_data);

  int png_color_type = (image.GetType()==img::COMPLEX) ? PNG_COLOR_TYPE_RGB : PNG_COLOR_TYPE_GRAY;
  unsigned int bpp = 8;

  png_set_IHDR(png_ptr, info_ptr, size.GetWidth(), size.GetHeight(), bpp,
         png_color_type,
         PNG_INTERLACE_NONE,
         PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(png_ptr, info_ptr);

  png_bytep* row_pointers = new png_bytep[size.GetHeight()];

  if(image.GetType()==img::COMPLEX) {
    uint width=size.GetWidth();
    uint height=size.GetHeight();

    img::RasterImage ri(width,height);
    ri.Fill(image, 0, image.GetExtent().GetStart().ToVec3(), 0,norm.GetImplementation(), img::RasterImage::PHASECOLOR);

    png_byte* png_data = (png_byte *) ri.GetDataPtr();
    img::Progress::Instance().Register(this,height+1,100);
    for(uint i=0;i<height;i++){
      row_pointers[i]=&png_data[i*3*width];
      img::Progress::Instance().AdvanceProgress(this);
    }
    png_write_image(png_ptr, row_pointers);
    img::Progress::Instance().AdvanceProgress(this);
    img::Progress::Instance().DeRegister(this);

  } else {

    png_byte* png_data = new png_byte[size.GetSlab()];

    int byte_val;
    img::Point start=image.GetExtent().GetStart();
    uint height=size.GetHeight();
    uint width=size.GetWidth();
    for(uint j=0;j<height;j++) {
      for(uint i=0;i<width;i++) {
        img::Point loc=img::Point(i,j)+start;
        byte_val = static_cast<int>(norm.Convert(image.GetReal(loc)));
        byte_val = (byte_val<0) ? 0 : byte_val;
        byte_val = (byte_val>255) ? 255 : byte_val;
        png_data[j*width+i]=static_cast<png_byte>(byte_val);
      }
      row_pointers[j]=&png_data[j*width];
    }

    png_write_image(png_ptr, row_pointers);

    delete []png_data;
  }
  delete []row_pointers;

  png_write_end(png_ptr, info_ptr);

  png_destroy_write_struct(&png_ptr, &info_ptr);

}

bool MapIOPngHandler::MatchContent(unsigned char* header)
{
  if(header[0]==0x89 &&
     header[1]==0x50 &&
     header[2]==0x4E &&
     header[3]==0x47 &&
     header[4]==0x0D &&
     header[5]==0x0A &&
     header[6]==0x1A &&
     header[7]==0x0A){
      return true;
  }
  return false;
}
bool MapIOPngHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()=="defined_png") {
    return true;
  }
  return false;
}
bool MapIOPngHandler::MatchSuffix(const String& loc)
{
    if(detail::FilenameEndsWith(loc,".png") ) {
      return true;
    }
    return false;
}

}} // ns


