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
  info wrapper

  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <ost/info/info.hh>
#include <ost/info/info_visitor.hh>
#include <ost/info/info_verbose_visitor.hh>
#include <ost/info/geom_info_conversion.hh>
using namespace boost::python;
using namespace ost;
using namespace ost::info;


InfoHandle (*CreateInfoPtr1)()=CreateInfo;
InfoHandle (*CreateInfoPtr2)(const String&)=CreateInfo;

void info_group_apply1a(InfoGroup* g, InfoVisitor& v)
{
  g->Apply(v);
}

void info_group_apply1b(InfoGroup* g, InfoVisitor& v, bool b)
{
  g->Apply(v,b);
}

void info_group_apply2a(InfoGroup* g, InfoConstVisitor& v)
{
  g->Apply(v);
}

void info_group_apply2b(InfoGroup* g, InfoConstVisitor& v, bool b)
{
  g->Apply(v,b);
}

void info_group_remove1(InfoGroup* g, const String& path )
{
  g->Remove(path);
}

void info_group_remove2(InfoGroup* g, InfoGroup& group )
{
  g->Remove(group);
}
void info_group_remove3(InfoGroup* g, const String& path, bool use_defaults )
{
  g->Remove(path,use_defaults);
}

void info_handle_apply1a(InfoHandle* h, InfoVisitor& v)
{
  h->Apply(v);
}

void info_handle_apply1b(InfoHandle* h, InfoVisitor& v, bool b)
{
  h->Apply(v,b);
}

void info_handle_apply2a(InfoHandle* h, InfoConstVisitor& v)
{
  h->Apply(v);
}

void info_handle_apply2b(InfoHandle* h, InfoConstVisitor& v, bool b)
{
  h->Apply(v,b);
}

void info_handle_remove1(InfoHandle* h, const String& path )
{
  h->Remove(path);
}

void info_handle_remove2(InfoHandle* h, InfoGroup& group )
{
  h->Remove(group);
}
void info_handle_remove3(InfoHandle* h, const String& path, bool use_defaults )
{
  h->Remove(path,use_defaults);
}


class InfoVisitorProxy : public InfoVisitor {
public:
  InfoVisitorProxy(PyObject *p)
      : self(p) {}

  InfoVisitorProxy(PyObject *p, const InfoVisitor& i)
      : InfoVisitor(i), self(p) {}
  
  virtual bool VisitGroup(InfoGroup& group)
  {
    return call_method<bool, InfoGroup>(self, "VisitGroup", group);
  }
  
  virtual void VisitItem(InfoItem& item)
  {
    call_method<void, InfoItem>(self, "VisitItem", item);
  }  
  
  bool VisitItemDefault(const InfoItem&) 
  {
    return true;
  }
  
  bool VisitGroupDefault(const InfoGroup&) 
  {
    return true;
  }  

  virtual void VisitGroupFinish(InfoGroup& group)
 {
    call_method<void, InfoGroup>(self, "VisitGroupFinish", group);
 }

 void VisitGroupFinishDefault(const InfoGroup&) 
 {
 }  
private:
  PyObject* self;
};

InfoItem (InfoGroup::*create_item_a)(const String&, const String&)=&InfoGroup::CreateItem;
InfoItem (InfoGroup::*create_item_b)(const String&, int)=&InfoGroup::CreateItem;
InfoItem (InfoGroup::*create_item_c)(const String&, bool)=&InfoGroup::CreateItem;
InfoItem (InfoGroup::*create_item_d)(const String&, Real)=&InfoGroup::CreateItem;


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getgroup_overloads, GetGroup, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(retrievegroup_overloads, RetrieveGroup, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(hasgroup_overloads, HasGroup, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getitem_overloads, GetItem, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(hasitem_overloads, HasItem, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(retrieveitem_overloads, RetrieveItem, 1, 2)


BOOST_PYTHON_MODULE(_ost_info)
{

  enum_<Type>("ItemType")
    .value("STRING",IT_STRING)
    .value("INT",IT_INT)
    .value("FLOAT",IT_FLOAT)
    .value("BOOL",IT_BOOL)
    .value("VECTOR",IT_VECTOR)
    ;

  void (InfoItem::*set_type1)(int type)=&InfoItem::SetType;
  
  class_<InfoItem>("InfoItem",no_init)
    .def("GetParent",&InfoItem::GetParent)
    .def("GetName",&InfoItem::GetName)
    .add_property("name",&InfoItem::GetName)
    .def("GetValue",&InfoItem::GetValue)
    .def("SetValue",&InfoItem::SetValue)
    .add_property("value",&InfoItem::GetValue,&InfoItem::SetValue)
    .def("GetAttribute",&InfoItem::GetAttribute)
    .def("SetAttribute",&InfoItem::SetAttribute)
    .def("HasAttribute",&InfoItem::HasAttribute)
    .def("GetAttributeList",&InfoItem::GetAttributeList)
    .def("AsFloat",&InfoItem::AsFloat)
    .def("AsInt",&InfoItem::AsInt)
    .def("AsBool",&InfoItem::AsBool)
    .def("AsVector",&InfoItem::AsVector)
    .def("SetType",set_type1)
    .def("GetType",&InfoItem::GetType)
    .def("SetFloat",&InfoItem::SetFloat)
    .def("SetInt",&InfoItem::SetInt)
    .def("SetBool",&InfoItem::SetBool)
    .def("SetVector",&InfoItem::SetVector)
    .add_property("attribute",&InfoItem::GetAttribute,&InfoItem::SetAttribute)
    .def("GetComment",&InfoItem::GetComment)
    ;
  class_<InfoGroupList>("InfoGroupList", no_init)
    .def(vector_indexing_suite<InfoGroupList>())
  ;
  class_<InfoItemList>("InfoItemList", no_init)
    .def(vector_indexing_suite<InfoItemList>())
  ;  
  class_<InfoGroup>("InfoGroup",no_init)
    .def("GetParent",&InfoGroup::GetParent)
    .def("GetName",&InfoGroup::GetName)
    .def("SetName",&InfoGroup::SetName)
    .add_property("name",&InfoGroup::GetName,&InfoGroup::SetName)
    .def("GetGroup",&InfoGroup::GetGroup,getgroup_overloads())
    .def("GetGroups", &InfoGroup::GetGroups)
    .def("GetItems", &InfoGroup::GetItems)
    .def("CreateGroup",&InfoGroup::CreateGroup)
    .def("RetrieveGroup",&InfoGroup::RetrieveGroup,retrievegroup_overloads())
    .def("HasGroup",&InfoGroup::HasGroup,hasgroup_overloads())
    .def("GetItem",&InfoGroup::GetItem,getitem_overloads())
    .def("CreateItem",create_item_a)
    .def("CreateItem",create_item_b)
    .def("CreateItem",create_item_c)
    .def("CreateItem",create_item_d)            
    .def("HasItem",&InfoGroup::HasItem,hasitem_overloads())
    .def("RetrieveItem",&InfoGroup::RetrieveItem,retrieveitem_overloads())
    .def("Remove",info_group_remove1)
    .def("Remove",info_group_remove2)
    .def("Remove",info_group_remove3)
    .def("GetAttribute",&InfoGroup::GetAttribute)
    .def("SetAttribute",&InfoGroup::SetAttribute)
    .def("HasAttribute",&InfoGroup::HasAttribute)
    .def("GetAttributeList",&InfoGroup::GetAttributeList)
    .def("SetTextData", &InfoGroup::SetTextData)
    .def("Apply",info_group_apply1a)
    .def("Apply",info_group_apply1b)
    .def("Apply",info_group_apply2a)
    .def("Apply",info_group_apply2b)
    .def("GetTextData",&InfoGroup::GetTextData)
    .def("GetComment",&InfoGroup::GetComment)
    .def("GetPath",&InfoGroup::GetPath)
    ;

  class_<InfoHandle>("InfoHandle",no_init)
    .def("Import",&InfoHandle::Import)
    .def("Export",&InfoHandle::Export)
    .def("Root",&InfoHandle::Root)
    .def("AddDefault",&InfoHandle::AddDefault)
    .def("Copy",&InfoHandle::Copy)
    .def("HasDefaultGroup",&InfoHandle::HasDefaultGroup)
    .def("GetDefaultGroup",&InfoHandle::GetDefaultGroup)
    .def("HasDefaultItem",&InfoHandle::HasDefaultItem)
    .def("GetDefaultItem",&InfoHandle::GetDefaultItem)
    .def("GetParent",&InfoHandle::GetParent)
    .def("GetName",&InfoHandle::GetName)
    .def("SetName",&InfoHandle::SetName)
    .add_property("name",&InfoHandle::GetName,&InfoHandle::SetName)
    .def("GetGroup",&InfoHandle::GetGroup,getgroup_overloads())
    .def("CreateGroup",&InfoHandle::CreateGroup)
    .def("RetrieveGroup",&InfoHandle::RetrieveGroup,retrievegroup_overloads())
    .def("HasGroup",&InfoHandle::HasGroup,hasgroup_overloads())
    .def("GetItem",&InfoHandle::GetItem,getitem_overloads())
    .def("CreateItem",&InfoHandle::CreateItem)
    .def("HasItem",&InfoHandle::HasItem,hasitem_overloads())
    .def("RetrieveItem",&InfoHandle::RetrieveItem,retrieveitem_overloads())
    .def("Remove",info_handle_remove1)
    .def("Remove",info_handle_remove2)
    .def("Remove",info_handle_remove3)
    .def("GetAttribute",&InfoHandle::GetAttribute)
    .def("SetAttribute",&InfoHandle::SetAttribute)
    .def("HasAttribute",&InfoHandle::HasAttribute)
    .def("GetAttributeList",&InfoHandle::GetAttributeList)
    .def("Apply",info_handle_apply1a)
    .def("Apply",info_handle_apply1b)
    .def("Apply",info_handle_apply2a)
    .def("Apply",info_handle_apply2b)
    .def("GetTextData",&InfoHandle::GetTextData)
    ;

  def("CreateInfo",CreateInfoPtr1);
  def("CreateInfo",CreateInfoPtr2);
  def("LoadInfo",&LoadInfo);

  class_<InfoPath>("InfoPath",init<const String&>())
    .def("GetList",&InfoPath::GetList)
    .def("IsRelative",&InfoPath::IsRelative)
    .def(self_ns::str(self))
    ;

  def("GetFloatInfoItem",GetFloatInfoItem);
  def("GetIntInfoItem",GetIntInfoItem);
  def("GetBoolInfoItem",GetBoolInfoItem);
  def("GetVectorInfoItem",GetVectorInfoItem);
  def("GetStringInfoItem",GetStringInfoItem);

  def("SetFloatInfoItem",SetFloatInfoItem);
  def("SetIntInfoItem",SetIntInfoItem);
  def("SetBoolInfoItem",SetBoolInfoItem);
  def("SetVectorInfoItem",SetVectorInfoItem);
  def("SetStringInfoItem",SetStringInfoItem);

  def("Mat2ToInfo", &Mat2ToInfo);
  def("Mat3ToInfo", &Mat3ToInfo);  
  def("Mat4ToInfo", &Mat4ToInfo);  
  
  def("Mat2FromInfo", &Mat2FromInfo);
  def("Mat3FromInfo", &Mat3FromInfo);  
  def("Mat4FromInfo", &Mat4FromInfo);  
  implicitly_convertible<String,InfoPath>();

  class_<InfoVisitor, InfoVisitorProxy>("InfoVisitor",init<>())
    .def("VisitGroup", &InfoVisitorProxy::VisitGroupDefault)
    .def("VisitItem", &InfoVisitorProxy::VisitItemDefault)
    .def("VisitGroupFinish", &InfoVisitorProxy::VisitGroupFinishDefault)
  ;

  class_<VerboseInfoVisitor, bases<InfoVisitor> >("VerboseInfoVisitor",init<>())
  ;

  def("TransformToInfo", &TransformToInfo);
  def("TransformFromInfo", &TransformFromInfo);

}
