//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2023 by the OpenStructure authors
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
#ifndef OST_IO_STAR_WRITER_HH
#define OST_IO_STAR_WRITER_HH

#include <map>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/shared_ptr.hpp>
#include <ost/string_ref.hh>
#include <ost/io/io_exception.hh>

namespace{
  // float to string with specified number of decimals
  void fts(Real f, int decimals, String& s) {
    char data[20];
    size_t len;
    switch(decimals){
      case 0:
        len = std::snprintf(data, sizeof(data), "%.0f", f);
        break;
      case 1:
        len = std::snprintf(data, sizeof(data), "%.1f", f);
        break;
      case 2:
        len = std::snprintf(data, sizeof(data), "%.2f", f);
        break;
      case 3:
        len = std::snprintf(data, sizeof(data), "%.3f", f);
        break;
      case 4:
        len = std::snprintf(data, sizeof(data), "%.4f", f);
        break;
      case 5:
        len = std::snprintf(data, sizeof(data), "%.5f", f);
        break;
      case 6:
        len = std::snprintf(data, sizeof(data), "%.6f", f);
        break;
      default:
        throw ost::io::IOException("Max decimals in float conversion: 6");
    }
  
    if(len < 0 || len > 20) {
      throw ost::io::IOException("float conversion failed");
    }
    s.assign(data, len);
  }
}

namespace ost { namespace io {

class StarWriterObject;
class StarWriterValue;
class StarWriterDataItem;
class StarWriterLoopDesc;
class StarWriterLoop;
typedef boost::shared_ptr<StarWriterObject> StarWriterObjectPtr;
typedef boost::shared_ptr<StarWriterValue> StarWriterValuePtr;
typedef boost::shared_ptr<StarWriterDataItem> StarWriterDataItemPtr;
typedef boost::shared_ptr<StarWriterLoopDesc> StarWriterLoopDescPtr;
typedef boost::shared_ptr<StarWriterLoop> StarWriterLoopPtr;


class DLLEXPORT_OST_IO StarWriterObject {
public:
  virtual ~StarWriterObject() { }
  virtual void ToStream(std::ostream& s) = 0;
};


class DLLEXPORT_OST_IO StarWriterValue{
public:
  static StarWriterValue FromInt(int int_value) {
    StarWriterValue value;
    value.value_ = std::to_string(int_value);
    return value;
  } 
  static StarWriterValue FromFloat(Real float_value, int decimals) {
    StarWriterValue value;
    fts(float_value, decimals, value.value_);
    return value;  
 }
  static StarWriterValue FromString(const String& string_value) {
    StarWriterValue value;
    // cases we still need to deal with:
    // - special characters in strings (put in quotation marks)
    // - long strings (semicolon based syntax)
    // see https://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html
    bool has_space = false;
    for(char c: string_value) {
      if(isspace(c)) {
        has_space = true;
        break;
      }
    }
    if(string_value == "") {
      value.value_ = ".";
    } else if(has_space) {
      value.value_ = "'" + string_value + "'";
    }
    else {
      value.value_ = string_value;
    }
    return value;
  }
  const String& GetValue() const { return value_; }
private:
// force construction through static members
StarWriterValue() { }
String value_;
};


class DLLEXPORT_OST_IO StarWriterDataItem : public StarWriterObject {
public:
  StarWriterDataItem(const String& category, const String& attribute, 
                     const StarWriterValue& value): category_(category),
                                                    attribute_(attribute),
                                                    value_(value) { }
  virtual void ToStream(std::ostream& s) {
    s << category_ << '.' << attribute_ << ' ' << value_.GetValue() << std::endl;
  }
  const String& GetCategory() const { return category_; }
  const String& GetAttribute() const { return attribute_; }
  const StarWriterValue& GetValue() const { return value_; }
private:
  String category_;
  String attribute_;
  StarWriterValue value_;
};


class DLLEXPORT_OST_IO StarWriterLoopDesc : public StarWriterObject {
public:
  StarWriterLoopDesc(const String& category): category_(category) { }
  
  int GetIndex(const String& attribute) const {
    std::map<String, int>::const_iterator i=index_map_.find(attribute);
    return i==index_map_.end() ? -1 : i->second;
  }

  void Add(const String& attribute) {
    index_map_.insert(std::make_pair(attribute, index_map_.size()));
  }

  size_t GetSize() const  {
    return index_map_.size();
  }

  virtual void ToStream(std::ostream& s) {
    std::vector<std::pair<int, String> > tmp;
    for(auto it = index_map_.begin(); it != index_map_.end(); ++it) {
      tmp.push_back(std::make_pair(it->second, it->first));
    }
    std::sort(tmp.begin(), tmp.end());
    for(auto it = tmp.begin(); it != tmp.end(); ++it) {
      s << category_ << "." << it->second << std::endl;
    }
  }

  const String& GetCategory() const { return category_; }
private:
  String                category_;
  std::map<String, int> index_map_;
};


class DLLEXPORT_OST_IO StarWriterLoop: public StarWriterObject {
public:

  StarWriterLoop(const StarWriterLoopDesc& desc): desc_(desc) { }

  const StarWriterLoopDesc& GetDesc() { return desc_; }

  void AddData(const std::vector<StarWriterValue>& data) {
    if(data.size() != desc_.GetSize()) {
      throw ost::io::IOException("Invalid data size when adding to StarLoop");
    }
    data_.insert(data_.end(), data.begin(), data.end());
  }

  const std::vector<StarWriterValue>& GetData() { return data_; }

  int GetN() {
    return data_.size() / desc_.GetSize();
  }

  virtual void ToStream(std::ostream& s) {
    if(data_.empty()) {
      return; // skip loop, including header
    }
    s << "loop_" << std::endl;
    desc_.ToStream(s);
    int desc_size = desc_.GetSize();
    for(size_t i = 0; i < data_.size(); ++i) {
      s << data_[i].GetValue();
      if((i+1) % desc_size == 0) {
        s << std::endl;
      } else {
        s << ' ';
      }
    }
  }

private:
  StarWriterLoopDesc desc_;
  std::vector<StarWriterValue> data_;
};


class DLLEXPORT_OST_IO StarWriter {
public:
  StarWriter(std::ostream& stream);
  StarWriter(const String& filename);
  virtual ~StarWriter() { }

  void Push(StarWriterObjectPtr obj) { categories_to_write_.push_back(obj); }
  void Write(const String& data_name);
private:
  String filename_;
  std::ofstream fstream_;
  boost::iostreams::filtering_stream<boost::iostreams::output> stream_;
  std::vector<StarWriterObjectPtr> categories_to_write_;
};

}} // ns

#endif
