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
#ifndef DX_IO_tiff_IO_PLUGIN_OST_HH
#define DX_IO_tiff_IO_PLUGIN_OST_HH

/*
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <boost/logic/tribool.hpp>

#include "map_io_handler.hh"



struct tiff;
    
namespace ost { namespace io {

class DLLEXPORT_OST_IO TIF: public  ImageFormatBase
{

 public:

  TIF(boost::logic::tribool normalize_on_save = boost::logic::indeterminate, Format bit_depth = OST_DEFAULT_FORMAT, bool sign = false, bool phasecolor = false, int subimage = -1);

  Format GetBitDepth() const;
  void SetBitDepth ( Format bitdepth);

  bool GetSigned () const;
  void SetSigned (bool sign);

  bool GetPhasecolor () const;
  void SetPhasecolor (bool phasecolor);

  Format GetFormat () const;
  void SetFormat (Format bit_depth);

  boost::logic::tribool GetNormalizeOnSave() const;
  void SetNormalizeOnSave(boost::logic::tribool normalize_on_save);

  int GetSubimage() const;
  void SetSubimage(int subimage);


  Real GetMaximum() const;
  Real GetMinimum() const;

  static String FORMAT_STRING;


 protected:

  TIF(String format_string, boost::logic::tribool normalize_on_save = boost::logic::indeterminate, Format bit_depth = OST_DEFAULT_FORMAT, bool sign = false, bool phasecolor = false, int subimage = -1);
  void do_export(const img::MapHandle& image,tiff* tfile,const TIF& formattif) const;

 private:

  boost::logic::tribool normalize_on_save_;
  Format bit_depth_;
  bool signed_;
  bool phasecolor_;
  int subimage_;

};


class DLLEXPORT_OST_IO MapIOTiffHandler: public MapIOHandler
{
  public:
    /// \brief Map IO handler to read/write Tiff files
    ///
    /// This map IO handler reads and writes tiff-formatted files.
    virtual void Import(img::MapHandle& sh, const boost::filesystem::path& loc,
                        const ImageFormatBase& formatstruct);
    virtual void Import(img::MapHandle& sh, std::istream& loc,
                        const ImageFormatBase& formatstruct);
    virtual void Export(const img::MapHandle& sh, 
                        const boost::filesystem::path& loc,
                        const ImageFormatBase& formatstruct) const;
    virtual void Export(const img::MapHandle& sh, std::ostream& loc,
                        const ImageFormatBase& formatstruct) const;
    static bool MatchContent(unsigned char* header);
    static bool MatchType(const ImageFormatBase& type);
    static bool MatchSuffix(const String& loc);
    static bool ProvidesImport() { return true; }
    static bool ProvidesExport() { return true; }
    static String GetFormatName() { return String( "Tiff"); }
    static String GetFormatDescription() { return String("Tagged Image File Format"); }

  protected:

    tiff* open_subimage_file(const boost::filesystem::path& location,
                             const TIF& formattif);
    tiff* open_subimage_stream(std::istream& location,const TIF& formattif);
    void load_image_data(tiff* tfile, img::ImageHandle& image,
                         const TIF& formattif);
    virtual void do_export(const img::MapHandle& sh,tiff* tfile,
                           TIF& formatstruct) const;

};

typedef MapIOHandlerFactory<MapIOTiffHandler> MapIOTiffHandlerFactory;

}} // ns

#endif
