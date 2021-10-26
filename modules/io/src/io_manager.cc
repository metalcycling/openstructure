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
#include <ost/io/io_manager.hh>
#include <ost/io/mol/entity_io_pdb_handler.hh>
#include <ost/io/mol/entity_io_crd_handler.hh>
#include <ost/io/mol/entity_io_sdf_handler.hh>
#include <ost/io/mol/entity_io_mae_handler.hh>
#include <ost/io/mol/entity_io_mmcif_handler.hh>
#include <ost/io/seq/fasta_io_handler.hh>
#include <ost/io/seq/pir_io_handler.hh>
#include <ost/io/seq/promod_io_handler.hh>
#include <ost/io/seq/hhm_io_handler.hh>
#include <ost/io/seq/pssm_io_handler.hh>
#include <ost/io/mol/surface_io_msms_handler.hh>
#include <ost/io/seq/clustal_io_handler.hh>
#include  <ost/io/img/map_io_dx_handler.hh>
#include  <ost/io/img/map_io_spi_handler.hh>
#include  <ost/io/img/map_io_mrc_handler.hh>
#include  <ost/io/img/map_io_dm3_handler.hh>
#include  <ost/io/img/map_io_situs_handler.hh>
#include  <ost/io/img/map_io_tiff_handler.hh>
#include  <ost/io/img/map_io_png_handler.hh>
#include  <ost/io/img/map_io_dat_handler.hh>
#include  <ost/io/img/map_io_jpk_handler.hh>
#include  <ost/io/img/map_io_nanoscope_handler.hh>
#include  <ost/io/img/map_io_df3_handler.hh>
#include  <ost/io/img/map_io_ipl_handler.hh>

namespace ost { namespace io {

IOManager::IOManager()
{
  RegisterFactory(EntityIOHandlerFactoryBaseP(new EntityIOMMCIFHandlerFactory));
  RegisterFactory(EntityIOHandlerFactoryBaseP(new EntityIOPDBHandlerFactory));
  RegisterFactory(EntityIOHandlerFactoryBaseP(new EntityIOCRDHandlerFactory));
  RegisterFactory(EntityIOHandlerFactoryBaseP(new EntityIOSDFHandlerFactory));
  RegisterFactory(EntityIOHandlerFactoryBaseP(new EntityIOMAEHandlerFactory));
  RegisterFactory(SequenceIOHandlerFactoryBasePtr(new FastaIOHandlerFactory));
  RegisterFactory(SequenceIOHandlerFactoryBasePtr(new PirIOHandlerFactory));
  RegisterFactory(SequenceIOHandlerFactoryBasePtr(new ClustalIOHandlerFactory));
  RegisterFactory(SequenceIOHandlerFactoryBasePtr(new PromodIOHandlerFactory));
  RegisterFactory(ProfileIOHandlerFactoryBasePtr(new HhmIOHandlerFactory));
  RegisterFactory(ProfileIOHandlerFactoryBasePtr(new PssmIOHandlerFactory));
  RegisterFactory(SurfaceIOHandlerFactoryBasePtr(new SurfaceIOMSMSHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIODxHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIOSpiHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIOSitusHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIOMrcHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIODm3HandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIOTiffHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIOPngHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIOJpkHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIODatHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIONanoscopeHandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIODF3HandlerFactory));
  RegisterFactory(MapIOHandlerFactoryBasePtr(new MapIOIPLHandlerFactory));
}

IOManager& IOManager::Instance()
{
  static IOManager iom;
  return iom;
}

EntityIOHandlerP IOManager::FindEntityImportHandler(const String& filename,
                                                    const String& format)
{
  for(EntityIOHFList::const_iterator it=entity_iohf_list_.begin();
      it!=entity_iohf_list_.end();++it) {
    if( (*it)->ProvidesImport(filename,format)) {
      return (*it)->Create();
    }
  }

  throw IOUnknownFormatException("no suitable entity io handler found for "+filename);
}

EntityIOHandlerP IOManager::FindEntityExportHandler(const String& filename,
                                                    const String& format)
{
  for (EntityIOHFList::const_iterator it=entity_iohf_list_.begin();
       it!=entity_iohf_list_.end();++it) {
    if ((*it)->ProvidesExport(filename, format)) {
      return (*it)->Create();
    }
  }
  throw IOUnknownFormatException("no suitable handler found for "+filename);
}

SequenceIOHandlerPtr IOManager::FindAlignmentImportHandler(const String& filename,
                                                           const String& format) {
  for(AlignmentIOFList::const_iterator it=alignment_io_list_.begin();
      it!=alignment_io_list_.end();++it) {
    if( (*it)->ProvidesImport(filename,format)) {
      return (*it)->Create();
    }
  }
  throw IOUnknownFormatException("no suitable alignment io handler found for "+filename);
}

SequenceIOHandlerPtr IOManager::FindAlignmentExportHandler(const String& filename,
                                                           const String& format)
{
  for(AlignmentIOFList::const_iterator it=alignment_io_list_.begin();
      it!=alignment_io_list_.end();++it) {
    if( (*it)->ProvidesExport(filename,format)) {
      return (*it)->Create();
    }
  }
  throw IOUnknownFormatException("no suitable alignment io handler found for "+filename);
}

ProfileIOHandlerPtr IOManager::FindProfileImportHandler(const String& filename,
                                                        const String& format) {
  for(ProfileIOFList::const_iterator it=profile_io_list_.begin();
      it!=profile_io_list_.end();++it) {
    if( (*it)->ProvidesImport(filename,format)) {
      return (*it)->Create();
    }
  }
  throw IOUnknownFormatException("no suitable profile io handler found for "+filename);
}

ProfileIOHandlerPtr IOManager::FindProfileExportHandler(const String& filename,
                                                        const String& format)
{
  for(ProfileIOFList::const_iterator it=profile_io_list_.begin();
      it!=profile_io_list_.end();++it) {
    if( (*it)->ProvidesExport(filename,format)) {
      return (*it)->Create();
    }
  }
  throw IOUnknownFormatException("no suitable profile io handler found for "+filename);
}

SurfaceIOHandlerPtr IOManager::FindSurfaceImportHandler(const String& filename,
                                                        const String& format)
{
  for(SurfaceIOFList::const_iterator it=surface_io_list_.begin();
      it!=surface_io_list_.end();++it) {
    if( (*it)->ProvidesImport(filename,format)) {
      return (*it)->Create();
    }
  }

  throw IOUnknownFormatException("no suitable surface io handler found for "+filename);
}

MapIOHandlerPtr IOManager::FindMapImportHandlerFile(const boost::filesystem::path& loc,
                                                    const ImageFormatBase& formatstruct)
{
  if(formatstruct.GetFormatString()!="undefined" ){
    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesImport() && (*it)->MatchType(formatstruct)) {
        return (*it)->Create();
      }
    }
    throw IOUnknownFormatException("Unsupported type in FindMapImportHandle.");
  }else{
    String match_suf_string=loc.string();
    std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesImport() && (*it)->MatchSuffix(match_suf_string)) {
        return (*it)->Create();
      }
    }
    unsigned char header[256];
    boost::filesystem::ifstream infile(loc, std::ios::binary);
    if(!infile)
    {
      throw IOException("could not open "+loc.string());
    }
    infile.read(reinterpret_cast<char*>(&header),256);
    infile.close();

    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesImport() && (*it)->MatchContent(header)) {
        return (*it)->Create();
      }
    }
    throw IOUnknownFormatException("No file suffix given for " + loc.string()+", and could not detect automatically, please indicate file type.");
  }
  return MapIOHandlerPtr();  // removes warning
}

MapIOHandlerPtr IOManager::FindMapImportHandlerStream(std::istream& stream,
                                                        const ImageFormatBase& formatstruct)
{
  if(formatstruct.GetFormatString()!="undefined" ){
    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesImport() && (*it)->MatchType(formatstruct)) {
        return (*it)->Create();;
      }
    }
    throw IOUnknownFormatException("Unsuspported type in FindMapImportHandle.");
  }else{
    unsigned char header[256];
    stream.read(reinterpret_cast<char*>(&header),256);
    stream.seekg(0,std::ios::beg);

    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesImport() && (*it)->MatchContent(header)) {
        return (*it)->Create();;
      }
    }
    throw IOUnknownFormatException("No type given for stream, and could not detect automatically, please indicate file type.");
  }
  return MapIOHandlerPtr(); //removes warning
}

MapIOHandlerPtr IOManager::FindMapExportHandlerFile(const boost::filesystem::path& loc,
                                                    const ImageFormatBase& formatstruct)
{
  if(formatstruct.GetFormatString()!="undefined" ){
    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesExport() && (*it)->MatchType(formatstruct)) {
        return (*it)->Create();
      }
    }
    throw IOUnknownFormatException("Unsupported type in FindMapExportHandle.");
  }else{
    String filename=loc.string();
    String::size_type pos = filename.rfind(".");
    if (pos == String::npos){
      throw IOUnknownFormatException("No file suffix given for " + filename+", please indicate file type.");
    }
    String match_suf_string=loc.string();
    std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesExport() && (*it)->MatchSuffix(match_suf_string)) {
        return(*it)->Create();
      }
    }
    throw IOUnknownFormatException("could not find io-plugin for " + filename);
  }
  return MapIOHandlerPtr(); // removes the warning
}

MapIOHandlerPtr IOManager::FindMapExportHandlerStream(std::istream& stream, const ImageFormatBase& formatstruct)
{
  if(formatstruct.GetFormatString()!="undefined" ){
    for(MapIOFList::const_iterator it=map_io_list_.begin(); it!=map_io_list_.end();++it) {
      if((*it)->ProvidesExport() && (*it)->MatchType(formatstruct)) {
        return (*it)->Create();
      }
    }
    throw IOUnknownFormatException("Unsupported type in FindMapExportHandlerStream");
    } else {
    throw IOUnknownFormatException("No type given. Please indicate stream type");
  }
  return MapIOHandlerPtr(); // removes the warning
}


void IOManager::RegisterFactory(const MapIOHandlerFactoryBasePtr& f) 
{
  map_io_list_.push_back(f);
}

const MapIOFList& IOManager::GetAvailableMapHandler() const
{
  return map_io_list_;
}

void IOManager::RegisterFactory(const EntityIOHandlerFactoryBaseP& f) 
{
  entity_iohf_list_.push_back(f);
}

void IOManager::RegisterFactory(const SequenceIOHandlerFactoryBasePtr& f) 
{
  alignment_io_list_.push_back(f);
}

void IOManager::RegisterFactory(const ProfileIOHandlerFactoryBasePtr& f) 
{
  profile_io_list_.push_back(f);
}

void IOManager::RegisterFactory(const SurfaceIOHandlerFactoryBasePtr& f) 
{
  surface_io_list_.push_back(f);
}

const EntityIOHFList& IOManager::GetAvailableEntityHandler() const
{
  return entity_iohf_list_;
}

const AlignmentIOFList& IOManager::GetAvailableAlignmentHandler() const
{
  return alignment_io_list_;
}

const ProfileIOFList& IOManager::GetAvailableProfileHandler() const
{
  return profile_io_list_;
}

const SurfaceIOFList& IOManager::GetAvailableSurfaceHandler() const
{
  return surface_io_list_;
}

}} // ns
