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

  StarWriterValue() { }

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
   
    if(string_value == "") {
      value.value_ = "?";
    } else {
      // string requires quotes if any of the following is True
      // information from https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax
      // * space in string
      // * any string that starts with any of the following strings
      //   * _
      //   * #
      //   * $
      //   * '
      //   * "
      //   * [
      //   * ]
      //   * ;
      //   * data_ (case insensitive)
      //   * save_ (case insensitive)
      // * any string that is equal to any of the following reserved words 
      //   * loop_ (case insensitive)
      //   * stop_ (case insensitive)
      //   * global_ (case insensitive)
      bool needs_quotes = false;

      // space in string
      for(char c: string_value) {
        if(isspace(c)) {
          needs_quotes = true;
          break;
        }
      }

      // any string that starts with any of the special single characters
      if(!needs_quotes) {
        switch(string_value[0]) {
          case '_': {
            needs_quotes = true;
            break;
          }
          case '#': {
            needs_quotes = true;
            break;
          }
          case '$': {
            needs_quotes = true;
            break;
          }
          case '\'': {
            needs_quotes = true;
            break;
          }
          case '\"': {
            needs_quotes = true;
            break;
          }
          case '[': {
            needs_quotes = true;
            break;
          }
          case ']': {
            needs_quotes = true;
            break;
          }
          case ';': {
            needs_quotes = true;
            break;
          }
        }
      }

      // any string that starts with any of the special multi character thingies
      if(!needs_quotes && string_value.size() >= 5 && string_value[4] == '_') {
        // need to do case insensitive checking
        if((string_value[0] == 'd' || string_value[0] == 'D') &&
           (string_value[1] == 'a' || string_value[1] == 'A') &&
           (string_value[2] == 't' || string_value[2] == 'T') &&
           (string_value[3] == 'a' || string_value[3] == 'A')) {
            needs_quotes = true;
        }
        if((string_value[0] == 's' || string_value[0] == 'S') &&
           (string_value[1] == 'a' || string_value[1] == 'A') &&
           (string_value[2] == 'v' || string_value[2] == 'V') &&
           (string_value[3] == 'e' || string_value[3] == 'E')) {
            needs_quotes = true;
        }
      }

      // any string that is exactly one of the reserved words
      if(!needs_quotes && string_value.size() == 5 && string_value[4] == '_') {
        // need to do case insensitive checking
        if((string_value[0] == 'l' || string_value[0] == 'L') &&
           (string_value[1] == 'o' || string_value[1] == 'O') &&
           (string_value[2] == 'o' || string_value[2] == 'O') &&
           (string_value[3] == 'p' || string_value[3] == 'P')) {
            needs_quotes = true;
        }
        if((string_value[0] == 's' || string_value[0] == 'S') &&
           (string_value[1] == 't' || string_value[1] == 'T') &&
           (string_value[2] == 'o' || string_value[2] == 'O') &&
           (string_value[3] == 'p' || string_value[3] == 'P')) {
            needs_quotes = true;
        }
      }

      if(!needs_quotes && string_value.size() == 7 && string_value[6] == '_') {
        // need to do case insensitive checking
        if((string_value[0] == 'g' || string_value[0] == 'G') &&
           (string_value[1] == 'l' || string_value[1] == 'L') &&
           (string_value[2] == 'o' || string_value[2] == 'O') &&
           (string_value[3] == 'b' || string_value[3] == 'B') &&
           (string_value[4] == 'a' || string_value[4] == 'A') &&
           (string_value[5] == 'l' || string_value[5] == 'L')) {
            needs_quotes = true;
        }
      }

      if(needs_quotes) {
        value.value_ = "\"" + string_value + "\"";
      } else {
        value.value_ = string_value;
      }
    }
    return value;
  }
  const String& GetValue() const { return value_; }
private:
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
  StarWriter() { }
  virtual ~StarWriter() { }

  void Push(StarWriterObjectPtr obj) { categories_to_write_.push_back(obj); }

  void Write(const String& data_name, const String& filename);
  void Write(const String& data_name, std::ostream& stream);

private:
  std::vector<StarWriterObjectPtr> categories_to_write_;
};

}} // ns

#endif
