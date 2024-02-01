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
  Collection of template specializations for 
  info item type interconversion

  Author: Johan Hebert, Ansgar Philippsen
*/

#ifndef OST_INFO_TYPE_CAST_HH
#define OST_INFO_TYPE_CAST_HH

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace ost { namespace info { namespace detail {

namespace {

struct ItemCastError
{};

template <typename IN, typename OUT>
OUT do_cast(const IN& in) {
  try {
    return boost::lexical_cast<OUT>(in);
  } catch (boost::bad_lexical_cast &) {
    throw ItemCastError();
  }
  return OUT(); // should never get here, but makes some compilers happy
}

// this is the to-be-specialized template function
// -> implementation detail
template <int INTYPE, int OUTTYPE>
void set_new_type(EleImpl& item);

template <> void set_new_type<IT_STRING,IT_STRING>(EleImpl& item)
{
  //item.SetStringRepr(item.GetStringRepr(),true);
}

template <> void set_new_type<IT_STRING,IT_INT>(EleImpl& item)
{
  try {
    item.SetIntRepr(do_cast<String, int>(item.GetStringRepr()));
  } catch (ItemCastError&) {
    item.SetIntRepr(0);
  }
}

template <> void set_new_type<IT_STRING,IT_FLOAT>(EleImpl& item)
{
  try {
    item.SetFloatRepr(do_cast<String, float>(item.GetStringRepr()));
  } catch (ItemCastError&) {
    item.SetFloatRepr(0.0);
  }
}

template <> void set_new_type<IT_STRING,IT_BOOL>(EleImpl& item)
{
  const int Strings = 3;
  String trueStrings[] = {"true", "1", "yes"};
  String falseStrings[] = {"false", "0", "no"};

  String value = item.GetStringRepr();

  for(int i=0;i<Strings;++i) {
    if(value == trueStrings[i]) {
      item.SetBoolRepr(true);
      i = Strings;
    } else if(value == falseStrings[i]) {
      item.SetBoolRepr(false);
      i = Strings;
    } else {
    }
  }
}

template <> void set_new_type<IT_STRING,IT_VECTOR>(EleImpl& item)
{
  static String num("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?");
  boost::regex expression("\\[ *("+num+") *, *("+num+") *, *("+num+") *\\]");
  boost::cmatch what;

  if(boost::regex_match(item.GetStringRepr().c_str(), what, expression)) { 
    item.SetVecRepr(geom::Vec3(do_cast<String, Real>(what[1]),
			       do_cast<String, Real>(what[3]),
			       do_cast<String, Real>(what[5])));
  } else {
    item.SetVecRepr(geom::Vec3(0.0, 0.0, 0.0));
  }
}

template <> void set_new_type<IT_INT,IT_STRING>(EleImpl& item)
{
  try {
    item.SetStringRepr(do_cast<int, String>(item.GetIntRepr()),true);
  } catch (ItemCastError&) {
    item.SetStringRepr("",true);
  }
}

template <> void set_new_type<IT_INT,IT_INT>(EleImpl& item)
{
  //item.SetIntRepr(item.GetIntRepr());
}

template <> void set_new_type<IT_INT,IT_FLOAT>(EleImpl& item)
{
  try {
    item.SetFloatRepr(do_cast<int, float>(item.GetIntRepr()));
  } catch (ItemCastError&) {
    item.SetFloatRepr(0.0);
  }
}

template <> void set_new_type<IT_INT,IT_BOOL>(EleImpl& item)
{
  if(item.GetIntRepr() == 0) {
    item.SetBoolRepr(false);
  } else {
    item.SetBoolRepr(true);
  }
}

template <> void set_new_type<IT_INT,IT_VECTOR>(EleImpl& item)
{
  try {
    item.SetVecRepr(geom::Vec3(do_cast<int, Real>(item.GetIntRepr()), 0.0, 0.0));
  } catch (ItemCastError&) {
    item.SetVecRepr(geom::Vec3(0.0, 0.0, 0.0));
  }
}

template <> void set_new_type<IT_FLOAT,IT_STRING>(EleImpl& item)
{
  try {
    item.SetStringRepr(do_cast<float, String>(item.GetFloatRepr()),true);
  } catch (ItemCastError&) {
    item.SetStringRepr("",true);
  }
}

template <> void set_new_type<IT_FLOAT,IT_INT>(EleImpl& item)
{
  item.SetIntRepr(static_cast<int>(round(item.GetFloatRepr())));
}

template <> void set_new_type<IT_FLOAT,IT_FLOAT>(EleImpl& item)
{
  //item.SetFloatRepr(item.GetFloatRepr());
}

template <> void set_new_type<IT_FLOAT,IT_BOOL>(EleImpl& item)
{
  if(item.GetFloatRepr() == 0.0) {
    item.SetBoolRepr(false);
  } else {
    item.SetBoolRepr(true);
  }
}

template <> void set_new_type<IT_FLOAT,IT_VECTOR>(EleImpl& item)
{
    item.SetVecRepr(geom::Vec3(item.GetFloatRepr(), 0.0, 0.0));
}

template <> void set_new_type<IT_BOOL,IT_STRING>(EleImpl& item)
{
  try {
    item.SetStringRepr(do_cast<bool, String>(item.GetBoolRepr()),true);
  } catch (ItemCastError&) {
    item.SetStringRepr("",true);
  }
}

template <> void set_new_type<IT_BOOL,IT_INT>(EleImpl& item)
{
  try {
    item.SetIntRepr(do_cast<bool, int>(item.GetBoolRepr()));
  } catch (ItemCastError&) {
    item.SetIntRepr(0);
  }
}

template <> void set_new_type<IT_BOOL,IT_FLOAT>(EleImpl& item)
{
  try {
    item.SetFloatRepr(do_cast<bool, float>(item.GetBoolRepr()));
  } catch (ItemCastError&) {
    item.SetFloatRepr(0.0);
  }
}

template <> void set_new_type<IT_BOOL,IT_BOOL>(EleImpl& item)
{
  //item.SetBoolRepr(item.GetBoolRepr());
}

template <> void set_new_type<IT_BOOL,IT_VECTOR>(EleImpl& item)
{
  try {
    item.SetVecRepr(geom::Vec3(do_cast<bool, Real>(item.GetBoolRepr()), 0.0, 0.0));
  } catch (ItemCastError&) {
    item.SetVecRepr(geom::Vec3(0.0, 0.0, 0.0));
  }
}

template <> void set_new_type<IT_VECTOR,IT_STRING>(EleImpl& item)
{
  std::ostringstream str;
  str << item.GetVecRepr();
  item.SetStringRepr(str.str(),true);
}

template <> void set_new_type<IT_VECTOR,IT_INT>(EleImpl& item)
{
  try {
    item.SetIntRepr(do_cast<Real, int>(item.GetVecRepr()[0]));
  } catch (ItemCastError&) {
    item.SetIntRepr(0);
  }
}

template <> void set_new_type<IT_VECTOR,IT_FLOAT>(EleImpl& item)
{
    item.SetFloatRepr(item.GetVecRepr()[0]);
}

template <> void set_new_type<IT_VECTOR,IT_BOOL>(EleImpl& item)
{
  if(item.GetVecRepr()[0] == 0.0) {
    item.SetBoolRepr(false);
  } else {
    item.SetBoolRepr(true);
  }
}

template <> void set_new_type<IT_VECTOR,IT_VECTOR>(EleImpl& item)
{
  //item.SetVecRepr(item.GetVecRepr());
}

} // anon ns

#define ITEM_TYPE_CAST_EVAL(INTYPE) \
    if(new_type==IT_STRING) { \
      set_new_type< INTYPE ,IT_STRING>(item); \
    } else if(new_type==IT_INT) { \
      set_new_type< INTYPE ,IT_INT>(item); \
    } else if(new_type==IT_FLOAT) { \
      set_new_type< INTYPE ,IT_FLOAT>(item); \
    } else if(new_type==IT_BOOL) { \
      set_new_type< INTYPE ,IT_BOOL>(item); \
    } else if(new_type==IT_VECTOR) { \
      set_new_type< INTYPE ,IT_VECTOR>(item); \
    }

//! this is the public interface
void SetInfoItemNewType(EleImpl& item, int new_type)
{
  Type in_type=item.GetType();
  if(in_type==IT_STRING) {
    ITEM_TYPE_CAST_EVAL(IT_STRING)
  } else if(in_type==IT_INT) {
    ITEM_TYPE_CAST_EVAL(IT_INT)
  } else if(in_type==IT_FLOAT) {
    ITEM_TYPE_CAST_EVAL(IT_FLOAT)
  } else if(in_type==IT_BOOL) {
    ITEM_TYPE_CAST_EVAL(IT_BOOL)
  } else if(in_type==IT_VECTOR) {
    ITEM_TYPE_CAST_EVAL(IT_VECTOR)
  }
}

void SetInfoItemNewType(EleImpl& item, Type new_type)
{
  Type in_type=item.GetType();
  if(in_type==IT_STRING) {
    ITEM_TYPE_CAST_EVAL(IT_STRING)
  } else if(in_type==IT_INT) {
    ITEM_TYPE_CAST_EVAL(IT_INT)
  } else if(in_type==IT_FLOAT) {
    ITEM_TYPE_CAST_EVAL(IT_FLOAT)
  } else if(in_type==IT_BOOL) {
    ITEM_TYPE_CAST_EVAL(IT_BOOL)
  } else if(in_type==IT_VECTOR) {
    ITEM_TYPE_CAST_EVAL(IT_VECTOR)
  }
}

}}} // ns

#endif
