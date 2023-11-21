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
#ifndef OST_IO_STAR_BASE_HH
#define OST_IO_STAR_BASE_HH

#include <map>
#include <algorithm>
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


typedef enum {
  STAR_DIAG_WARNING,
  STAR_DIAG_ERROR
} StarDiagType;


class DLLEXPORT_OST_IO StarObject {
public:
  virtual ~StarObject() { }
  virtual void ToStream(std::ostream& s) = 0;
};

class DLLEXPORT_OST_IO StarDataItem : public StarObject{
public:
  StarDataItem(const StringRef& category, const StringRef& name, 
                const StringRef& value): 
    category_(category), name_(name), value_(value)
  { }

  virtual void ToStream(std::ostream& s) {
    s << category_ << '.' << name_ << ' ' << value_ << std::endl;
  }

  const StringRef& GetCategory() const { return category_; }
  const StringRef& GetName() const { return name_; }
  const StringRef& GetValue() const { return value_; }
private:
  StringRef category_;
  StringRef name_;
  StringRef value_;
};


class DLLEXPORT_OST_IO StarDataItemDO : public StarObject{
  // basically a copy of StarDataItem but with data ownership
  // Problem with StringRef: It's purely pointer based
  // => data needs to reside elsewhere
public:
  StarDataItemDO(const String& category, const String& name, 
                 const String& value): 
    category_(category), name_(name)
  {
    // cases we still need to deal with:
    // - special characters in strings (put in quotation marks)
    // - long strings (semicolon based syntax)
    // see https://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html
    if(value == "") {
      value_ = ".";
    } else {
      value_ = value;
    }
  }

  StarDataItemDO(const String& category, const String& name, 
                 Real value, int decimals): 
    category_(category), name_(name)
  {
    fts(value, decimals, value_);
  }

  StarDataItemDO(const String& category, const String& name, 
                 int value): 
    category_(category), name_(name)
  {
    value_ = std::to_string(value);
  }

  virtual void ToStream(std::ostream& s) {
    s << category_ << '.' << name_ << ' ' << value_ << std::endl;
  }

  const String& GetCategory() const { return category_; }
  const String& GetName() const { return name_; }
  const String& GetValue() const { return value_; }
private:
  String category_;
  String name_;
  String value_;
};

class DLLEXPORT_OST_IO StarDataCategoryDO : public StarObject {
public:
  StarDataCategoryDO(const String& category): category_(category) {}

  void Add(const StarDataItemDO& data_item) {
    if(data_item.GetCategory() != category_) {
      throw ost::io::IOException("category mismatch");
    }
    data_items_.push_back(data_item);
  }

  virtual void ToStream(std::ostream& s) {
    for(auto it = data_items_.begin(); it != data_items_.end(); ++it) {
      it->ToStream(s);
    }
  }

private:
  String category_;
  std::vector<StarDataItemDO> data_items_;
};

class DLLEXPORT_OST_IO StarLoopDesc : public StarObject {
public:
  StarLoopDesc():
    category_("")
  { }
  
  int GetIndex(const String& name) const
  {
    std::map<String, int>::const_iterator i=index_map_.find(name);
    return i==index_map_.end() ? -1 : i->second;
  }
  
  void SetCategory(const StringRef& category)
  {
    category_=category.str();
  }

  void SetCategory(const String& category)
  {
    category_=category;
  }  

  void Add(const StringRef& name)
  {
    index_map_.insert(std::make_pair(name.str(), index_map_.size()));
  }
  void Add(const String& name)
  {
    index_map_.insert(std::make_pair(name, index_map_.size()));
  }
  size_t GetSize() const 
  {
    return index_map_.size();
  }
  void Clear()
  {
    category_.clear();
    index_map_.clear();
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

class DLLEXPORT_OST_IO StarLoopDataItemDO{
public:

  StarLoopDataItemDO(const String& value) {
    // cases we still need to deal with:
    // - special characters in strings (put in quotation marks)
    // - long strings (semicolon based syntax)
    // see https://mmcif.wwpdb.org/docs/tutorials/mechanics/pdbx-mmcif-syntax.html
    if(value == "") {
      value_ = ".";
    } else {
      value_ = value;
    }
  }

  StarLoopDataItemDO(Real value, int decimals) {
    fts(value, decimals, value_);
  }

  StarLoopDataItemDO(int value) {
    value_ = std::to_string(value);
  }

  const String& GetValue() const { return value_; }

  virtual void ToStream(std::ostream& s) {
    s << value_;
  }

private:
  String value_;
};

class DLLEXPORT_OST_IO StarLoop: public StarObject {
public:

  StarLoop() { }

  StarLoop(const StarLoopDesc& desc): desc_(desc) { }

  void SetDesc(const StarLoopDesc& desc) {
    if(!data_.empty()) {
      throw ost::io::IOException("Can only set new StarLoop desc in "
                                 "in empty loop");
    }
    desc_ = desc;
  }

  void AddData(const std::vector<StarLoopDataItemDO>& data) {
    if(data.size() != desc_.GetSize()) {
      throw ost::io::IOException("Invalid data size when adding to StarLoop");
    }
    data_.insert(data_.end(), data.begin(), data.end());
  }

  int GetN() {
    return data_.size() / desc_.GetSize();
  }

  virtual void ToStream(std::ostream& s) {
    s << "loop_" << std::endl;
    desc_.ToStream(s);
    int desc_size = desc_.GetSize();
    for(size_t i = 0; i < data_.size(); ++i) {
      data_[i].ToStream(s);
      if((i+1) % desc_size == 0) {
        s << std::endl;
      } else {
        s << ' ';
      }
    }
  }

private:
  StarLoopDesc desc_;
  std::vector<StarLoopDataItemDO> data_;
};
 
}} // ns

#endif
