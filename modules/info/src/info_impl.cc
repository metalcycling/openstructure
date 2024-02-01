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
  low level info implementation
  
  Author: Ansgar Philippsen
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include <boost/filesystem/operations.hpp>

#include <QtCore/QFile>

#include "info_error.hh"
#include "info_impl.hh"
#include "item_type.hh"
#include "item_type_cast.hh"

namespace ost { namespace info { namespace detail {

namespace bf = boost::filesystem;

namespace {

// typecasting stuff

template <typename IN>
String to_String(const IN& in);

template<>
String to_String(const geom::Vec3& in)
{
  std::ostringstream str;
  str << in;
  return str.str();
}

template <typename IN>
String to_String(const IN& in) {
  std::ostringstream str;
  str<<in;
  return str.str();
}

String QS(const QString& s)
{
  return s.toStdString();
}

QString QS(const String& s)
{
  return QString(s.c_str());
}

} // anon ns

InfoImpl::InfoImpl():
  document_("EMDataInfo"),
  def_list_()
{
  document_.appendChild(document_.createElement("EMDataInfo"));
}

InfoImpl::InfoImpl(bool dummy,const String& text):
  document_("EMDataInfo"),
  def_list_()
{
  if (!document_.setContent(QString::fromStdString(text))) {
    throw(InfoError("QDomDocument.setContent failed"));
  }
}

InfoImpl::InfoImpl(const String& fname):
  document_("EMDataInfo"),
  def_list_()
{
  Import(fname);
}

RootPtr InfoImpl::Copy() const
{
  RootPtr nrvo(new InfoImpl());
  nrvo->document_=document_.cloneNode().toDocument();
  nrvo->def_list_=def_list_;
  return nrvo;
}

void InfoImpl::Import(const String& fname)
{
  document_=QDomDocument("EMDataInfo");
  QFile file(QS(fname));
  if(!file.exists()){
    file.close();
    throw(InfoError("File: "+fname+" doesn't exist."));
  }
  QString errormessage;
  int errorline,errorcolumn;
  if (!document_.setContent(&file,&errormessage,&errorline,&errorcolumn)) {
    file.close();
    if(file.error()){
      switch(file.error()){
      case QFile::ReadError:
        throw(InfoError("Error while reading file: "+fname));
        break;
      case QFile::OpenError:
        throw(InfoError("File: "+fname+" could not be opened."));
        break;
      case QFile::PermissionsError:
        throw(InfoError("Permission error accessing file: "+fname));
        break;
      default:
        throw(InfoError("Unspecified error while accessing file: "+fname));
        break;
      }
    }else{
      std::ostringstream stream;
      stream << "Error while accessing file " << fname << ": " << errormessage.toStdString() << " in line " << errorline << " column " << errorcolumn;
      throw(InfoError(stream.str()));
    }
  }
  file.close();
}

void InfoImpl::Export(const String& fname)
{
  QFile out(QS(fname));
  out.open(QIODevice::WriteOnly);
  out.write(document_.toByteArray());
}

EleImpl InfoImpl::Root()
{
  return EleImpl(document_.documentElement());
}


// ele impl
EleImpl::EleImpl(const QDomElement& e):
  ele_(e),
  type_(IT_STRING),
  String_repr_(""),
  int_repr_(0),
  float_repr_(0.0),
  bool_repr_(false),
  vec_repr_(0.0,0.0,0.0)
{
}

String EleImpl::GetAttribute(const String& name) const
{
  if(!HasAttribute(name)) {
    throw InfoError("Requested attribute '" + name + "' not found");
  }

  return QS(ele_.attribute(QS(name)));
}

void EleImpl::SetAttribute(const String& name, const String& value)
{
  ele_.setAttribute(QS(name),QS(value));
}

bool EleImpl::HasAttribute(const String& name) const
{
  return ele_.hasAttribute(QS(name));
}

void EleImpl::RemoveAttribute(const String& name)
{
  if(ele_.hasAttribute(QS(name))) {
    ele_.removeAttribute(QS(name));
  }
}

String EleImpl::GetName() const
{
  return QS(ele_.tagName());
}

void EleImpl::SetName(const String& n)
{
  ele_.setTagName(QS(n));
}

String EleImpl::GetPath() const
{
  QDomNode p = ele_.parentNode();
  if(p.isNull()|| GetName()=="EMDataInfo") {
    return "/";
  } else {
    return EleImpl(p.toElement()).GetPath()+"/"+GetName();
  }
}

QDomElement& EleImpl::GetElement()
{
  return ele_;
}

const QDomElement& EleImpl::GetElement() const
{
  return ele_;
}

ElePtr EleImpl::GetSub(const String& name) const
{
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isElement()) {
      QDomElement en = n.toElement();
      if(en.tagName()==QString(QS(name))) {
        return ElePtr(new EleImpl(en));
      }
    }
  }
  throw InfoError("Requested group '" + name + "' not found");
}

EleList EleImpl::GetSubs(const String& name) const
{
  EleList eles;
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isElement()) {
      QDomElement en = n.toElement();
      if(en.tagName()==QString(QS(name))) {
        eles.push_back(ElePtr(new EleImpl(en)));
      }
    }
  }
  return eles;
}

ElePtr EleImpl::CreateSub(const String& name)
{
  QDomDocument doc = ele_.ownerDocument();
  QDomElement new_ele = doc.createElement(QS(name));
  ele_.appendChild(new_ele);
  return ElePtr(new EleImpl(new_ele));
}

bool EleImpl::HasSub(const String& name) const
{
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isElement()) {
      QDomElement en = n.toElement();
      if(en.tagName()==QString(QS(name))) {
        return true;
      }
    }
  }
  return false;
}

void EleImpl::RemoveSub(const String& name)
{
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isElement()) {
      QDomElement en = n.toElement();
      if(en.tagName()==QString(QS(name))) {
	ele_.removeChild(en);
      }
    }
  }
}

void EleImpl::RemoveSub(const ElePtr& ele)
{
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isElement()) {
      QDomElement en = n.toElement();
      if(en==ele->GetElement()) {
    	  ele_.removeChild(en);
      }
    }
  }
}

std::vector<ElePtr> EleImpl::GetSubList() const
{
  std::vector<ElePtr> epl;
  
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isElement()) {
      QDomElement ne = n.toElement();
      epl.push_back(ElePtr(new EleImpl(ne)));
    }
  }

  return epl;
}

ElePtr EleImpl::GetParent() const
{
  QDomNode p = ele_.parentNode();
  if(p.isNull()) {
    return ElePtr(new EleImpl(QDomElement()));
  }
  return ElePtr(new EleImpl(p.toElement()));
}

String EleImpl::GetTextData() const
{
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isText()) {
      return QS(n.toText().data());
    }
  }
  return "";
}

void EleImpl::SetTextData(const String& td)
{
  bool succ=false;
  QDomNodeList list = ele_.childNodes();
  for(int i=0;i<list.length();++i) {
    QDomNode n = list.item(i);
    if(n.isText()) {
      n.toText().setData(QS(td));
      succ=true;
    }
  }
  if(!succ) {
    QDomText tn=ele_.ownerDocument().createTextNode(QS(td));
    ele_.appendChild(tn);
  }
}


String EleImpl::GetComment() const
{
  QDomNode previous=ele_.previousSibling();
  if(previous.isComment()){
    return previous.toComment().data().toStdString();
  }
  return String();
}

Type EleImpl::GetType() const 
{
  return type_;
}

namespace {
String TypeToString(Type t) {
  switch(t) {
  case IT_INT: return "int";
  case IT_FLOAT: return "float";
  case IT_BOOL: return "bool";
  case IT_VECTOR: return "vector";
  case IT_STRING: break; // to make the compiler happy
  }
  return "String";
}
}

void EleImpl::SetType(Type t) 
{
  SetAttribute("type",TypeToString(t));
  type_=t;
}


void EleImpl::SetNewType(int t)
{
  SetInfoItemNewType(*this,t);
}

std::vector<String> EleImpl::GetAttributeList() const
{
  std::vector<String> result;
  QDomNamedNodeMap nmap = ele_.attributes();
  for(int i=0;i<nmap.length();++i) {
    QDomNode n = nmap.item(i);
    if(n.isAttr()) {
      QDomAttr attr = n.toAttr();
      if(attr.name()!="type" && attr.name()!="value") {
        // exlude internal attributes
        result.push_back(attr.name().toStdString());
      }
    }
  }
  return result;
}

bool EleImpl::operator==(const EleImpl& ref) const
{
  return ele_==ref.ele_;
}

void EleImpl::SetStringRepr(const String& repr, bool settype)
{
  String_repr_=repr;
  SetAttribute("value",repr);
  if(settype) SetType(IT_STRING);
}

String EleImpl::GetStringRepr() const 
{
  return String_repr_;
}

void EleImpl::SetIntRepr(int repr) 
{
  int_repr_=repr;
  SetStringRepr(to_String<int>(repr),false);
  SetType(IT_INT);
}

int EleImpl::GetIntRepr() const 
{
  return int_repr_;
}

void EleImpl::SetFloatRepr(float repr)
{
  float_repr_=repr;
  SetStringRepr(to_String<Real>(repr),false);
  SetType(IT_FLOAT);
}

float EleImpl::GetFloatRepr() const 
{
  return float_repr_;
}

void EleImpl::SetBoolRepr(bool repr) 
{
  bool_repr_=repr;
  SetStringRepr(to_String<bool>(repr),false);
  SetType(IT_BOOL);
}

bool EleImpl::GetBoolRepr() const 
{
  return bool_repr_;
}

void EleImpl::SetVecRepr(const geom::Vec3& repr) 
{
  vec_repr_=repr;
  SetStringRepr(to_String<geom::Vec3>(repr),false);
  SetType(IT_VECTOR);
}

geom::Vec3 EleImpl::GetVecRepr() const 
{
  return vec_repr_;
}


}}} // ns
