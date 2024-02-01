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
#ifndef OST_GUI_TOOL_OPTION_HH
#define OST_GUI_TOOL_OPTION_HH

/*
   Author: Marco Biasini 
 */

#include <limits>
#include <vector>
#include <iostream>


#include <ost/gui/module_config.hh>
#include <QObject>


namespace ost { namespace gui {

class DLLEXPORT_OST_GUI ToolOption : public QObject {
  Q_OBJECT
public:
  typedef enum {
    INT, FLOAT, ENUM, BUTTON
  } Type;  
protected:
  ToolOption(const String& key, const String& verbose_name, Type type);
public:
  const String& GetKey() const;
  const String& GetVerboseName() const;  
  Type GetType() const;
private:
  String key_;
  String verbose_name_;
  Type   type_;
};




template <typename T, ToolOption::Type C>
class DLLEXPORT ToolOptionNum : public ToolOption {
public:
  ToolOptionNum(const String& key, const String& verbose_name, T default_v, 
                T min_value=std::numeric_limits<T>::min(),
                T max_value=std::numeric_limits<T>::max()):
    ToolOption(key, verbose_name, C), value_(default_v),
    default_(default_v), min_value_(min_value), max_value_(max_value)
  {
    
  }
                
  T GetDefault() const { return default_; }
  T GetUpperLimit() const { return max_value_; }
  T GetLowerLimit() const { return min_value_; }
  T GetValue() const { return value_; }
  bool SetValue(T value)
  {
    if (value>=min_value_ && value<=max_value_) {
      value_=value;
      return true;
    }
    return false;
  }
private:  
  T value_;
  T default_;
  T min_value_;
  T max_value_;
};

typedef ToolOptionNum<int, ToolOption::INT> ToolOptionInt;
typedef ToolOptionNum<float, ToolOption::FLOAT> ToolOptionFloat;


#if !defined(_MSC_VER)
  extern template class ToolOptionNum<int,ToolOption::INT>;
  extern template class ToolOptionNum<float,ToolOption::FLOAT>;
#endif



class DLLEXPORT_OST_GUI ToolOptionEnum : public ToolOption {
public:
  struct Tuple {
    String text;
    int    tag;
  };
  typedef std::vector<Tuple>::const_iterator ConstIterator;
public:
  ToolOptionEnum(const String& key, const String& verbose_name);
  
  void Add(const String& text, int tag);
  
  void SetIndex(int index) { index_=index; }
  
  int GetIndex() const { return index_; }
  int GetValue() const { return tuples_[index_].tag; }
  
  ConstIterator Begin() const { return tuples_.begin(); }
  ConstIterator End() const { return tuples_.end(); }  
  size_t Size() { return tuples_.size(); }

private:  
  std::vector<Tuple>  tuples_;
  int                 index_;
};

class DLLEXPORT_OST_GUI ToolOptionButton : public ToolOption {
public:
  ToolOptionButton(const String& key,
                   const String& verbose_name,
                   QObject* receiver,
                   const char *slot_method);
  const char* GetSlotMethod() const { return slot_method_; }
  QObject* GetReceiver() const { return receiver_; }
private:
  const char *slot_method_;
  QObject* receiver_;
};


typedef std::vector<ToolOption*> ToolOptionList;

}}

#endif
