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
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
using namespace boost::python;

#include <ost/export_helper/pair_to_tuple_conv.hh>
#include <ost/io/mol/io_profile.hh>
#include <ost/io/mol/mmcif_reader.hh>
#include <ost/io/mol/mmcif_info.hh>
#include <ost/io/mol/star_writer.hh>
#include <ost/io/mol/mmcif_writer.hh>
#include <ost/io/mmcif_str.hh>
using namespace ost;
using namespace ost::io;
using namespace ost::mol;

template<typename T>
boost::python::list VecToList(std::vector<T>& vec){
  boost::python::list l;
  for(typename std::vector<T>::iterator it=vec.begin();it!=vec.end();++it){
    l.append(*it);
  }
  return l;
}

boost::python::list WrapGetNames(MMCifInfo *p){
  std::vector<String> names = p->GetEntityBranchChainNames();
  return VecToList<String>(names);
}

boost::python::tuple WrapMMCifStringToEntity(const String& mmcif,
                                             const IOProfile& profile=IOProfile(),
                                             bool process=false) {
  std::tuple<mol::EntityHandle, MMCifInfo, ost::seq::SequenceList> res =
  MMCifStringToEntity(mmcif, profile, process);
  return boost::python::make_tuple(std::get<0>(res),
                                   std::get<1>(res),
                                   std::get<2>(res));
}

String WrapEntityToMMCifStringEnt(const ost::mol::EntityHandle& ent,
                                 const String& data_name,
                                 ost::conop::CompoundLibPtr compound_lib,
                                 bool mmcif_conform) {
  return EntityToMMCifString(ent, data_name, compound_lib,
                             mmcif_conform);
}

String WrapEntityToMMCifStringView(const ost::mol::EntityView& ent,
                                   const String& data_name,
                                   ost::conop::CompoundLibPtr compound_lib,
                                   bool mmcif_conform) {
  return EntityToMMCifString(ent, data_name, compound_lib,
                             mmcif_conform);
}

void WrapStarLoopAddData(StarWriterLoop& sl, const boost::python::list& l) {
  std::vector<StarWriterValue> v;
  for (int i = 0; i < boost::python::len(l); ++i){
    v.push_back(boost::python::extract<StarWriterValue>(l[i]));
  }
  sl.AddData(v);
}

void WrapStarWriterWrite(StarWriter& writer, const String& data_name,
                         const String& filename) {
  writer.Write(data_name, filename);
}

void WrapSetStructureHandle(MMCifWriter& writer,
                            const ost::mol::EntityHandle& ent,
                            ost::conop::CompoundLibPtr compound_lib,
                            bool mmcif_conform,
                            const std::vector<MMCifWriterEntity>& entity_info) {
  writer.SetStructure(ent, compound_lib, mmcif_conform, entity_info);
}

void WrapSetStructureView(MMCifWriter& writer,
                          const ost::mol::EntityView& ent,
                          ost::conop::CompoundLibPtr compound_lib,
                          bool mmcif_conform,
                          const std::vector<MMCifWriterEntity>& entity_info) {
  writer.SetStructure(ent, compound_lib, mmcif_conform, entity_info);
}

void export_mmcif_io()
{
  class_<MMCifReader, boost::noncopyable>("MMCifReader", init<const String&, EntityHandle&, const IOProfile&>())
    .def("Parse", &MMCifReader::Parse)
    .def("SetRestrictChains", &MMCifReader::SetRestrictChains)
    .def("SetReadCanonicalSeqRes", &MMCifReader::SetReadCanonicalSeqRes)
    .def("GetSeqRes", &MMCifReader::GetSeqRes)
    .def("GetInfo", make_function(&MMCifReader::GetInfo,
                                  return_value_policy<copy_const_reference>()))
    .add_property("restrict_chains",
                  make_function(&MMCifReader::GetRestrictChains,
                                return_value_policy<copy_const_reference>()),
                  &MMCifReader::SetRestrictChains)
    .add_property("seqres", &MMCifReader::GetSeqRes)
    .add_property("read_seqres", &MMCifReader::GetReadSeqRes, 
                  &MMCifReader::SetReadSeqRes)
    .add_property("info", make_function(&MMCifReader::GetInfo,
                                   return_value_policy<copy_const_reference>()))
    ;

  class_<StarWriterObject, boost::noncopyable>("StarWriterObject", no_init);

  class_<StarWriterValue>("StarWriterValue", no_init)
    .def("FromInt", &StarWriterValue::FromInt, (arg("int_val"))).staticmethod("FromInt")
    .def("FromFloat", &StarWriterValue::FromFloat, (arg("float_val"), arg("decimals"))).staticmethod("FromFloat")
    .def("FromString", &StarWriterValue::FromString, (arg("string_val"))).staticmethod("FromString")
    .def("GetValue", &StarWriterValue::GetValue, return_value_policy<copy_const_reference>())
  ;

  class_<StarWriterDataItem, bases<StarWriterObject> >("StarWriterDataItem", init<const String&, const String&, const StarWriterValue&>())
    .def("GetCategory", &StarWriterDataItem::GetCategory, return_value_policy<copy_const_reference>())
    .def("GetAttribute", &StarWriterDataItem::GetAttribute, return_value_policy<copy_const_reference>())
    .def("GetValue", &StarWriterDataItem::GetValue, return_value_policy<copy_const_reference>())
  ;

  class_<StarWriterLoopDesc, bases<StarWriterObject> >("StarWriterLoopDesc", init<const String&>())
    .def("GetCategory", &StarWriterLoopDesc::GetCategory, return_value_policy<copy_const_reference>())
    .def("GetSize", &StarWriterLoopDesc::GetSize)
    .def("Add", &StarWriterLoopDesc::Add, (arg("attribute")))
    .def("GetIndex", &StarWriterLoopDesc::GetIndex, (arg("attribute")))
  ;

  class_<StarWriterLoop, bases<StarWriterObject> >("StarWriterLoop", init<const StarWriterLoopDesc&>())
    .def("GetDesc", &StarWriterLoop::GetDesc, return_value_policy<reference_existing_object>())
    .def("GetN", &StarWriterLoop::GetN)
    .def("AddData", &WrapStarLoopAddData, (arg("data_list")))
  ;

  class_<StarWriter>("StarWriter", init<>())
    .def("Push", &StarWriter::Push, arg("star_writer_object"))
    .def("Write", &WrapStarWriterWrite, (arg("data_name"), arg("filename")))
  ;

  class_<MMCifWriterEntity>("MMCifWriterEntity", no_init)
    .def("FromPolymer", &MMCifWriterEntity::FromPolymer).staticmethod("FromPolymer")
    .def("AddHet", &MMCifWriterEntity::AddHet, (arg("rnum"), arg("mon_id")))
    .add_property("type", &MMCifWriterEntity::type)
    .add_property("poly_type", &MMCifWriterEntity::poly_type)
    .add_property("branch_type", &MMCifWriterEntity::branch_type)
    .add_property("is_poly", &MMCifWriterEntity::is_poly)
    .add_property("mon_ids", &MMCifWriterEntity::mon_ids)
    .add_property("seq_olcs", &MMCifWriterEntity::seq_olcs)
    .add_property("seq_can_olcs", &MMCifWriterEntity::seq_can_olcs)
    .add_property("asym_ids", &MMCifWriterEntity::asym_ids)
  ;

  class_<std::vector<MMCifWriterEntity> >("MMCifWriterEntityList", init<>())
    .def(vector_indexing_suite<std::vector<MMCifWriterEntity> >())
  ;

  class_<MMCifWriter, bases<StarWriter> >("MMCifWriter", init<>())
    .def("SetStructure", &WrapSetStructureHandle, (arg("ent"), arg("compound_lib"),
                                                   arg("mmcif_conform")=true,
                                                   arg("entity_info")=std::vector<MMCifWriterEntity>()))
    .def("SetStructure", &WrapSetStructureView, (arg("ent"), arg("compound_lib"),
                                                 arg("mmcif_conform")=true,
                                                 arg("entity_info")=std::vector<MMCifWriterEntity>()))
    .def("GetEntities", &MMCifWriter::GetEntities, return_value_policy<copy_const_reference>())
  ;

  enum_<MMCifInfoCitation::MMCifInfoCType>("MMCifInfoCType")
    .value("Journal", MMCifInfoCitation::JOURNAL)
    .value("Book", MMCifInfoCitation::BOOK)
    .value("Unknown", MMCifInfoCitation::UNKNOWN)
  ;
 
  class_<MMCifInfoCitation>("MMCifInfoCitation", init<>())
    .def("SetID", &MMCifInfoCitation::SetID)
    .def("GetID", &MMCifInfoCitation::GetID)
    .def("SetCAS", &MMCifInfoCitation::SetCAS)
    .def("GetCAS", &MMCifInfoCitation::GetCAS)
    .def("SetISBN", &MMCifInfoCitation::SetISBN)
    .def("GetISBN", &MMCifInfoCitation::GetISBN)
    .def("SetPublishedIn", &MMCifInfoCitation::SetPublishedIn)
    .def("GetPublishedIn", &MMCifInfoCitation::GetPublishedIn)
    .def("SetVolume", &MMCifInfoCitation::SetVolume)
    .def("GetVolume", &MMCifInfoCitation::GetVolume)
    .def("SetPageFirst", &MMCifInfoCitation::SetPageFirst)
    .def("GetPageFirst", &MMCifInfoCitation::GetPageFirst)
    .def("SetPageLast", &MMCifInfoCitation::SetPageLast)
    .def("GetPageLast", &MMCifInfoCitation::GetPageLast)
    .def("SetDOI", &MMCifInfoCitation::SetDOI)
    .def("GetDOI", &MMCifInfoCitation::GetDOI)
    .def("SetPubMed", &MMCifInfoCitation::SetPubMed)
    .def("GetPubMed", &MMCifInfoCitation::GetPubMed)
    .def("SetYear", &MMCifInfoCitation::SetYear)
    .def("GetYear", &MMCifInfoCitation::GetYear)
    .def("SetTitle", &MMCifInfoCitation::SetTitle)
    .def("GetTitle", &MMCifInfoCitation::GetTitle)
    .def("SetBookPublisher", &MMCifInfoCitation::SetBookPublisher)
    .def("GetBookPublisher", &MMCifInfoCitation::GetBookPublisher)
    .def("SetBookPublisherCity", &MMCifInfoCitation::SetBookPublisherCity)
    .def("GetBookPublisherCity", &MMCifInfoCitation::GetBookPublisherCity)
    .def("SetCitationType", &MMCifInfoCitation::SetCitationType)
    .def("SetCitationTypeJournal", &MMCifInfoCitation::SetCitationTypeJournal)
    .def("SetCitationTypeBook", &MMCifInfoCitation::SetCitationTypeBook)
    .def("SetCitationTypeUnknown", &MMCifInfoCitation::SetCitationTypeUnknown)
    .def("GetCitationType", &MMCifInfoCitation::GetCitationType)
    .def("IsCitationTypeJournal", &MMCifInfoCitation::IsCitationTypeJournal)
    .def("IsCitationTypeBook", &MMCifInfoCitation::IsCitationTypeBook)
    .def("IsCitationTypeUnknown", &MMCifInfoCitation::IsCitationTypeUnknown)
    .def("SetAuthorList", &MMCifInfoCitation::SetAuthorList)
    .def("GetAuthorList", make_function(&MMCifInfoCitation::GetAuthorList,
                                   return_value_policy<copy_const_reference>()))
    .add_property("id", &MMCifInfoCitation::GetID, &MMCifInfoCitation::SetID)
    .add_property("cas", &MMCifInfoCitation::GetCAS, &MMCifInfoCitation::SetCAS)
    .add_property("isbn", &MMCifInfoCitation::GetISBN,
                  &MMCifInfoCitation::SetISBN)
    .add_property("published_in", &MMCifInfoCitation::GetPublishedIn,
                  &MMCifInfoCitation::SetPublishedIn)
    .add_property("volume", &MMCifInfoCitation::GetVolume,
                  &MMCifInfoCitation::SetVolume)
    .add_property("page_first", &MMCifInfoCitation::GetPageFirst,
                  &MMCifInfoCitation::SetPageFirst)
    .add_property("page_last", &MMCifInfoCitation::GetPageLast,
                  &MMCifInfoCitation::SetPageLast)
    .add_property("doi", &MMCifInfoCitation::GetDOI, &MMCifInfoCitation::SetDOI)
    .add_property("pubmed", &MMCifInfoCitation::GetPubMed,
                  &MMCifInfoCitation::SetPubMed)
    .add_property("year", &MMCifInfoCitation::GetYear,
                  &MMCifInfoCitation::SetYear)
    .add_property("title", &MMCifInfoCitation::GetTitle,
                  &MMCifInfoCitation::SetTitle)
    .add_property("book_publisher", &MMCifInfoCitation::GetBookPublisher,
                  &MMCifInfoCitation::SetBookPublisher)
    .add_property("book_publisher_city",
                  &MMCifInfoCitation::GetBookPublisherCity,
                  &MMCifInfoCitation::SetBookPublisherCity)
    .add_property("citation_type", &MMCifInfoCitation::GetCitationType,
                  &MMCifInfoCitation::SetCitationType)
    .add_property("authors", make_function(&MMCifInfoCitation::GetAuthorList,
                                   return_value_policy<copy_const_reference>()),
                  &MMCifInfoCitation::SetAuthorList)
    .def("__eq__", &MMCifInfoCitation::operator==) 
    .def("__ne__", &MMCifInfoCitation::operator!=)
  ;

  class_<std::vector<MMCifInfoCitation> >("MMCifInfoCitationList", init<>())
    .def(vector_indexing_suite<std::vector<MMCifInfoCitation> >())
  ;


  class_<MMCifInfoTransOp, MMCifInfoTransOpPtr>("MMCifInfoTransOp", init<>())
    .def("SetID", &MMCifInfoTransOp::SetID)
    .def("GetID", &MMCifInfoTransOp::GetID)
    .def("SetType", &MMCifInfoTransOp::SetType)
    .def("GetType", &MMCifInfoTransOp::GetType)
    .def("SetVector", &MMCifInfoTransOp::SetVector)
    .def("GetVector", &MMCifInfoTransOp::GetVector)
    .def("SetMatrix", &MMCifInfoTransOp::SetMatrix)
    .def("GetMatrix", &MMCifInfoTransOp::GetMatrix)
    .add_property("id", &MMCifInfoTransOp::GetID,
                  &MMCifInfoTransOp::SetID)
    .add_property("type", &MMCifInfoTransOp::GetType,
                  &MMCifInfoTransOp::SetType)
    .add_property("translation", &MMCifInfoTransOp::GetVector,
                  &MMCifInfoTransOp::SetVector)
    .add_property("rotation", &MMCifInfoTransOp::GetMatrix,
                  &MMCifInfoTransOp::SetMatrix)
  ;

  class_<std::vector<MMCifInfoTransOp> >("MMCifInfoTransOpList", init<>())
    .def(vector_indexing_suite<std::vector<MMCifInfoTransOp> >())
  ;

  typedef std::vector<MMCifInfoTransOpPtr> MMCifInfoTransOpPtrList;
  class_<std::vector<MMCifInfoTransOpPtr> >("MMCifInfoTransOpPtrList", init<>())
    .def(vector_indexing_suite<std::vector<MMCifInfoTransOpPtr>, true >())
  ;

  class_<std::vector<MMCifInfoTransOpPtrList > >("MMCifInfoTransOpPtrListList",
                                                init<>())
    .def(vector_indexing_suite<std::vector<MMCifInfoTransOpPtrList >, true >())
  ;
  class_<MMCifInfoStructRef, MMCifInfoStructRefPtr>("MMCifInfoStructRef", no_init)
  	.add_property("id", make_function(&MMCifInfoStructRef::GetID, 
  				        return_value_policy<copy_const_reference>()))
  	.add_property("db_name", make_function(&MMCifInfoStructRef::GetDBName, 
  				        return_value_policy<copy_const_reference>()))
  	.add_property("db_id", make_function(&MMCifInfoStructRef::GetDBID, 
  				        return_value_policy<copy_const_reference>()))
  	.add_property("entity_id", make_function(&MMCifInfoStructRef::GetEntityID, 
  				        return_value_policy<copy_const_reference>()))
  	.add_property("db_access", make_function(&MMCifInfoStructRef::GetDBAccess, 
  				        return_value_policy<copy_const_reference>()))
  	.def("GetAlignedSeq", &MMCifInfoStructRef::GetAlignedSeq, arg("align_id"))
  	.def("GetAlignedSeqs", &MMCifInfoStructRef::GetAlignedSeqs)
  	.add_property("aligned_seqs", &MMCifInfoStructRef::GetAlignedSeqs)
 ; 
  class_<MMCifInfoStructRefSeq, MMCifInfoStructRefSeqPtr>("MMCifInfoStructRefSeq", no_init)
  	.add_property("align_id", make_function(&MMCifInfoStructRefSeq::GetID, 
  				        return_value_policy<copy_const_reference>()))
  	.add_property("chain_name", make_function(&MMCifInfoStructRefSeq::GetChainName, 
  				        return_value_policy<copy_const_reference>()))
  	.add_property("seq_begin", &MMCifInfoStructRefSeq::GetSeqBegin)
  	.add_property("seq_end", &MMCifInfoStructRefSeq::GetSeqEnd)
  	.add_property("db_begin", &MMCifInfoStructRefSeq::GetDBBegin)
  	.add_property("db_end", &MMCifInfoStructRefSeq::GetDBEnd)
  	.add_property("difs", make_function(&MMCifInfoStructRefSeq::GetDifs,
  				        return_value_policy<copy_const_reference>()))
  ;
  class_<MMCifInfoStructRefSeqDif, 
  	     MMCifInfoStructRefSeqDifPtr>("MMCifInfoStructRefSeqDif", no_init)
  	.add_property("details", make_function(&MMCifInfoStructRefSeqDif::GetDetails, 
  				        return_value_policy<copy_const_reference>()))
  	.add_property("seq_rnum", &MMCifInfoStructRefSeqDif::GetSeqRNum)
    .add_property("db_rnum", make_function(&MMCifInfoStructRefSeqDif::GetDBRNum,
                                           return_value_policy<copy_const_reference>()))
  ;

  typedef std::pair<int, int> IntPair;
  to_python_converter<IntPair, PairToTupleConverter<int, int> >();
  typedef std::vector<IntPair> VectorIntPair;
  class_<VectorIntPair>("VectorIntPair", init<>())
    .def(vector_indexing_suite<VectorIntPair, true>())
  ;

  class_<MMCifInfoBioUnit>("MMCifInfoBioUnit", init<>())
    .def("SetDetails", &MMCifInfoBioUnit::SetDetails)
    .def("GetDetails", &MMCifInfoBioUnit::GetDetails)
    .def("SetMethodDetails", &MMCifInfoBioUnit::SetMethodDetails)
    .def("GetMethodDetails", &MMCifInfoBioUnit::GetMethodDetails)
    .def("AddChain", &MMCifInfoBioUnit::AddChain)
    .def("SetChainList", &MMCifInfoBioUnit::SetChainList)
    .def("GetChainList", make_function(&MMCifInfoBioUnit::GetChainList,
                                   return_value_policy<copy_const_reference>()))
    .def("GetChainIntervalList",
         make_function(&MMCifInfoBioUnit::GetChainIntervalList,
                       return_value_policy<copy_const_reference>()))
    .def("AddOperations", &MMCifInfoBioUnit::AddOperations)
    .def("GetOperations", make_function(&MMCifInfoBioUnit::GetOperations,
                                   return_value_policy<copy_const_reference>()))
    .def("GetOperationsIntervalList",
         make_function(&MMCifInfoBioUnit::GetOperationsIntervalList,
                       return_value_policy<copy_const_reference>()))
    .def("SetID", &MMCifInfoBioUnit::SetID)
    .def("GetID", &MMCifInfoBioUnit::GetID)
    .add_property("details", &MMCifInfoBioUnit::GetDetails,
                  &MMCifInfoBioUnit::SetDetails)
    .add_property("method_details", &MMCifInfoBioUnit::GetMethodDetails,
                  &MMCifInfoBioUnit::SetMethodDetails)
    .add_property("chains", make_function(&MMCifInfoBioUnit::GetChainList,
                                   return_value_policy<copy_const_reference>()))
    .add_property("chainintervalls", make_function(
                                  &MMCifInfoBioUnit::GetChainIntervalList,
                                  return_value_policy<copy_const_reference>()))
    .add_property("operations", make_function(&MMCifInfoBioUnit::GetOperations,
                                   return_value_policy<copy_const_reference>()))
    .add_property("operationsintervalls", make_function(
                                  &MMCifInfoBioUnit::GetOperationsIntervalList,
                                  return_value_policy<copy_const_reference>()))
    .add_property("id", &MMCifInfoBioUnit::GetID, &MMCifInfoBioUnit::SetID)
  ;

  class_<MMCifInfoStructRefs>("MMCifInfoStructRefs", init<>())
  	.def(vector_indexing_suite<MMCifInfoStructRefs, true>())
  ;
  class_<MMCifInfoStructRefSeqs>("MMCifInfoStructRefSeqs", init<>())
  	.def(vector_indexing_suite<MMCifInfoStructRefSeqs, true>())
  ;
  class_<MMCifInfoStructRefSeqDifs>("MMCifInfoStructRefSeqDifs", init<>())
  	.def(vector_indexing_suite<MMCifInfoStructRefSeqDifs, true>())
  ;
  class_<std::vector<MMCifInfoBioUnit> >("MMCifInfoBioUnitList", init<>())
    .def(vector_indexing_suite<std::vector<MMCifInfoBioUnit> >())
  ;

  class_<MMCifInfoStructDetails>("MMCifInfoStructDetails", init<>())
    .def("SetEntryID", &MMCifInfoStructDetails::SetEntryID)
    .def("GetEntryID", &MMCifInfoStructDetails::GetEntryID)
    .def("SetTitle", &MMCifInfoStructDetails::SetTitle)
    .def("GetTitle", &MMCifInfoStructDetails::GetTitle)
    .def("SetCASPFlag", &MMCifInfoStructDetails::SetCASPFlag)
    .def("GetCASPFlag", &MMCifInfoStructDetails::GetCASPFlag)
    .def("SetDescriptor", &MMCifInfoStructDetails::SetDescriptor)
    .def("GetDescriptor", &MMCifInfoStructDetails::GetDescriptor)
    .def("SetMass", &MMCifInfoStructDetails::SetMass)
    .def("GetMass", &MMCifInfoStructDetails::GetMass)
    .def("SetMassMethod", &MMCifInfoStructDetails::SetMassMethod)
    .def("GetMassMethod", &MMCifInfoStructDetails::GetMassMethod)
    .def("SetModelDetails", &MMCifInfoStructDetails::SetModelDetails)
    .def("GetModelDetails", &MMCifInfoStructDetails::GetModelDetails)
    .def("SetModelTypeDetails", &MMCifInfoStructDetails::SetModelTypeDetails)
    .def("GetModelTypeDetails", &MMCifInfoStructDetails::GetModelTypeDetails)
    .add_property("entry_id", &MMCifInfoStructDetails::GetEntryID,
                  &MMCifInfoStructDetails::SetEntryID)
    .add_property("title", &MMCifInfoStructDetails::GetTitle,
                  &MMCifInfoStructDetails::SetTitle)
    .add_property("casp_flag", &MMCifInfoStructDetails::GetCASPFlag,
                  &MMCifInfoStructDetails::SetCASPFlag)
    .add_property("descriptor", &MMCifInfoStructDetails::GetDescriptor,
                  &MMCifInfoStructDetails::SetDescriptor)
    .add_property("mass", &MMCifInfoStructDetails::GetMass,
                  &MMCifInfoStructDetails::SetMass)
    .add_property("mass_method", &MMCifInfoStructDetails::GetMassMethod,
                  &MMCifInfoStructDetails::SetMassMethod)
    .add_property("model_details", &MMCifInfoStructDetails::GetModelDetails,
                  &MMCifInfoStructDetails::SetModelDetails)
    .add_property("model_type_details",
                  &MMCifInfoStructDetails::GetModelTypeDetails,
                  &MMCifInfoStructDetails::SetModelTypeDetails)
  ;

  class_<MMCifInfoObsolete>("MMCifInfoObsolete", init<>())
    .def("SetDate", &MMCifInfoObsolete::SetDate)
    .def("GetDate", &MMCifInfoObsolete::GetDate)
    .def("SetID", &MMCifInfoObsolete::SetID)
    .def("GetID", &MMCifInfoObsolete::GetID)
    .def("SetPDBID", &MMCifInfoObsolete::SetPDBID)
    .def("GetPDBID", &MMCifInfoObsolete::GetPDBID)
    .def("SetReplacedPDBID", &MMCifInfoObsolete::SetReplacedPDBID)
    .def("GetReplacedPDBID", &MMCifInfoObsolete::GetReplacedPDBID)
    .add_property("date", &MMCifInfoObsolete::GetDate,
                  &MMCifInfoObsolete::SetDate)
    .add_property("id", &MMCifInfoObsolete::GetID,
                  &MMCifInfoObsolete::SetID)
    .add_property("pdb_id", &MMCifInfoObsolete::GetPDBID,
                  &MMCifInfoObsolete::SetPDBID)
    .add_property("replace_pdb_id", &MMCifInfoObsolete::GetReplacedPDBID,
                  &MMCifInfoObsolete::SetReplacedPDBID)
  ;

  class_<MMCifInfoRevisions>("MMCifInfoRevisions", init<>())
    .def("SetDateOriginal", &MMCifInfoRevisions::SetDateOriginal)
    .def("GetDateOriginal", &MMCifInfoRevisions::GetDateOriginal)
    .def("AddRevision", &MMCifInfoRevisions::AddRevision,
         (arg("num"), arg("date"), arg("status"), arg("major")=-1,
          arg("minor")=-1))
    .def("GetSize", &MMCifInfoRevisions::GetSize)
    .def("GetDate", &MMCifInfoRevisions::GetDate)
    .def("GetNum", &MMCifInfoRevisions::GetNum)
    .def("GetStatus", &MMCifInfoRevisions::GetStatus)
    .def("GetMajor", &MMCifInfoRevisions::GetMajor)
    .def("GetMinor", &MMCifInfoRevisions::GetMinor)
    .def("GetLastDate", &MMCifInfoRevisions::GetLastDate)
    .def("GetLastMajor", &MMCifInfoRevisions::GetLastMajor)
    .def("GetLastMinor", &MMCifInfoRevisions::GetLastMinor)
    .def("GetFirstRelease", &MMCifInfoRevisions::GetFirstRelease)
    .add_property("date_original", &MMCifInfoRevisions::GetDateOriginal,
                  &MMCifInfoRevisions::SetDateOriginal)
    .add_property("first_release", &MMCifInfoRevisions::GetFirstRelease)
  ;

  class_<MMCifInfoEntityBranchLink>("MMCifInfoEntityBranchLink",
                                    init<mol::AtomHandle,
                                    mol::AtomHandle, unsigned char>())
    .def("GetAtom1", &MMCifInfoEntityBranchLink::GetAtom1)
    .def("GetAtom2", &MMCifInfoEntityBranchLink::GetAtom2)
    .def("GetBondOrder", &MMCifInfoEntityBranchLink::GetBondOrder)
    .def("ConnectBranchLink", &MMCifInfoEntityBranchLink::ConnectBranchLink)
    .def("SetAtom1", &MMCifInfoEntityBranchLink::SetAtom1)
    .def("SetAtom2", &MMCifInfoEntityBranchLink::SetAtom2)
    .def("SetBondOrder", &MMCifInfoEntityBranchLink::SetBondOrder)
    .def(self_ns::str(self))
    .add_property("atom1", &MMCifInfoEntityBranchLink::GetAtom1,
                  &MMCifInfoEntityBranchLink::SetAtom1)
    .add_property("atom2", &MMCifInfoEntityBranchLink::GetAtom2,
                  &MMCifInfoEntityBranchLink::SetAtom2)
    .add_property("bond_order", &MMCifInfoEntityBranchLink::GetBondOrder,
                  &MMCifInfoEntityBranchLink::SetBondOrder)
  ;

  class_<MMCifInfoEntityBranchLinkMap>("MMCifInfoEntityBranchLinkMap", init<>())
    .def(map_indexing_suite<MMCifInfoEntityBranchLinkMap>())
  ;

  class_<std::vector<MMCifInfoEntityBranchLink> >(
                                                 "MMCifInfoEntityBranchLinkList",
                                                 init<>())
    .def(vector_indexing_suite<std::vector<MMCifInfoEntityBranchLink> >())
    .def(self_ns::str(self))
  ;

  class_<MMCifEntityDesc>("MMCifEntityDesc", init<>())
   .add_property("type", &MMCifEntityDesc::type)
   .add_property("entity_type", &MMCifEntityDesc::entity_type)
   .add_property("entity_poly_type", &MMCifEntityDesc::entity_poly_type)
   .add_property("branched_type", &MMCifEntityDesc::branched_type)
   .add_property("details", &MMCifEntityDesc::details)
   .add_property("seqres", &MMCifEntityDesc::seqres)
   .add_property("mon_ids", &MMCifEntityDesc::mon_ids)
   .add_property("hetero_num", &MMCifEntityDesc::hetero_num)
   .add_property("hetero_ids", &MMCifEntityDesc::hetero_ids)
  ;

  class_<MMCifInfo>("MMCifInfo", init<>())
    .def("AddCitation", &MMCifInfo::AddCitation)
    .def("GetCitations", make_function(&MMCifInfo::GetCitations,
                                   return_value_policy<copy_const_reference>()))
    .def("AddBioUnit", &MMCifInfo::AddBioUnit)
    .def("GetBioUnits", make_function(&MMCifInfo::GetBioUnits,
                                   return_value_policy<copy_const_reference>()))
    .def("SetMethod", &MMCifInfo::SetMethod)
    .def("GetMethod", &MMCifInfo::GetMethod)
    .def("SetResolution", &MMCifInfo::SetResolution)
    .def("GetResolution", &MMCifInfo::GetResolution)
    .def("SetEMResolution", &MMCifInfo::SetEMResolution)
    .def("GetEMResolution", &MMCifInfo::GetEMResolution)
    .def("SetRFree", &MMCifInfo::SetRFree)
    .def("GetRFree", &MMCifInfo::GetRFree)
    .def("SetRWork", &MMCifInfo::SetRWork)
    .def("GetRWork", &MMCifInfo::GetRWork)
    .def("AddAuthorsToCitation", &MMCifInfo::AddAuthorsToCitation, (arg("id"), arg("list"), arg("fault_tolerant")=false))
    .def("AddOperation", &MMCifInfo::AddOperation)
    .def("GetOperations", make_function(&MMCifInfo::GetOperations,
                                   return_value_policy<copy_const_reference>()))
    .def("SetStructDetails", &MMCifInfo::SetStructDetails)
    .def("GetStructDetails", &MMCifInfo::GetStructDetails)
    .def("SetObsoleteInfo", &MMCifInfo::SetObsoleteInfo)
    .def("GetObsoleteInfo", &MMCifInfo::GetObsoleteInfo)
    .def("AddMMCifPDBChainTr", &MMCifInfo::AddMMCifPDBChainTr)
    .def("GetMMCifPDBChainTr", &MMCifInfo::GetMMCifPDBChainTr)
    .def("AddPDBMMCifChainTr", &MMCifInfo::AddPDBMMCifChainTr)
    .def("GetPDBMMCifChainTr", &MMCifInfo::GetPDBMMCifChainTr)
    .def("AddMMCifEntityIdTr", &MMCifInfo::AddMMCifEntityIdTr)
    .def("GetMMCifEntityIdTr", &MMCifInfo::GetMMCifEntityIdTr)
    .def("SetRevisionsDateOriginal", &MMCifInfo::SetRevisionsDateOriginal)
    .def("AddRevision", &MMCifInfo::AddRevision,
         (arg("num"), arg("date"), arg("status"), arg("major")=-1,
          arg("minor")=-1))
    .def("GetRevisions", &MMCifInfo::GetRevisions)
    .def("AddEntityBranchLink", &MMCifInfo::AddEntityBranchLink)
    .def("GetEntityBranchLinks", &MMCifInfo::GetEntityBranchLinks)
    .def("GetEntityBranchByChain", &MMCifInfo::GetEntityBranchByChain)
    .def("ConnectBranchLinks", &MMCifInfo::ConnectBranchLinks)
    .def("GetEntityBranchChainNames", &WrapGetNames)
    .def("GetEntityBranchChains", &MMCifInfo::GetEntityBranchChains)
    .def("SetEntityDesc", &MMCifInfo::SetEntityDesc)
    .def("GetEntityDesc", &MMCifInfo::GetEntityDesc, return_value_policy<copy_const_reference>())
    .def("GetEntityIds", &MMCifInfo::GetEntityIds)
    .def("GetEntityIdsOfType", &MMCifInfo::GetEntityIdsOfType)
    .add_property("citations", make_function(&MMCifInfo::GetCitations,
                                   return_value_policy<copy_const_reference>()))
    .add_property("biounits", make_function(&MMCifInfo::GetBioUnits,
                                   return_value_policy<copy_const_reference>()))
    .add_property("method", &MMCifInfo::GetMethod, &MMCifInfo::SetMethod)
    .add_property("resolution", &MMCifInfo::GetResolution,
                  &MMCifInfo::SetResolution)
    .add_property("em_resolution", &MMCifInfo::GetEMResolution,
                  &MMCifInfo::SetEMResolution)
    .add_property("r_free", &MMCifInfo::GetRFree, &MMCifInfo::SetRFree)
    .add_property("r_work", &MMCifInfo::GetRWork, &MMCifInfo::SetRWork)
    .add_property("operations", make_function(&MMCifInfo::GetOperations,
                                   return_value_policy<copy_const_reference>()))
    .add_property("struct_details", &MMCifInfo::GetStructDetails,
                  &MMCifInfo::SetStructDetails)
    .add_property("struct_refs", make_function(&MMCifInfo::GetStructRefs,
    			        return_value_policy<copy_const_reference>()))
    .add_property("obsolete", &MMCifInfo::GetObsoleteInfo,
                  &MMCifInfo::SetObsoleteInfo)
    .add_property("revisions", &MMCifInfo::GetRevisions)
 ;

  def("MMCifStrToEntity", &WrapMMCifStringToEntity, (arg("pdb_string"),
                                                     arg("profile")=IOProfile(),
                                                     arg("process")=false));

  def("EntityToMMCifString",  &WrapEntityToMMCifStringEnt, (arg("ent"),
                                                            arg("data_name"),
                                                            arg("compound_lib"),
                                                            arg("mmcif_conform")));

  def("EntityToMMCifString",  &WrapEntityToMMCifStringView, (arg("ent"),
                                                             arg("data_name"),
                                                             arg("compound_lib"),
                                                             arg("mmcif_conform")));
}
