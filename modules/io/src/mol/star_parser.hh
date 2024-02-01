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
#ifndef OST_IO_STAR_PARSER_HH
#define OST_IO_STAR_PARSER_HH


/*
  Author: Marco Biasini
 */
#include <boost/iostreams/filtering_stream.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <ost/string_ref.hh>
#include <ost/io/module_config.hh>

namespace ost { namespace io {


typedef enum {
  STAR_DIAG_WARNING,
  STAR_DIAG_ERROR
} StarDiagType;


class DLLEXPORT_OST_IO StarDataItem {
public:
  StarDataItem(const StringRef& category, const StringRef& name, 
                const StringRef& value): 
    category_(category), name_(name), value_(value)
  { }
  const StringRef& GetCategory() const { return category_; }
  const StringRef& GetName() const { return name_; }
  const StringRef& GetValue() const { return value_; }
private:
  StringRef category_;
  StringRef name_;
  StringRef value_;
};

class DLLEXPORT_OST_IO StarLoopDesc {
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
  


  void Add(const StringRef& name)
  {
    index_map_.insert(std::make_pair(name.str(), index_map_.size()));
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

  const String& GetCategory() const { return category_; }
private:
  String                category_;
  std::map<String, int> index_map_;
};

/// \brief parser for the STAR file format
/// 
/// \section star_format STAR format description
/// 
/// In essence, STAR files are structured text files. A few simple control 
/// structures control the interpretation of the content.
/// 
/// \li Commented lines start with a '#'.
/// \li a data section starts with data_<id> and ends, whenever another data
///     block starts or at the end of the file.
/// \li data-items consist of a name-value pair. Names start with 
///     an underscore and are made up of a category and an identifier. The 
///     category and the identifier are separated by a '.'. The name is followed 
///     by a simple or quoted value.
/// \li lists of data items are specified using the loop_ keyword, followed by 
///     a list of column names and values.
class DLLEXPORT_OST_IO StarParser {
public:
  /// \brief create a StarParser
  ///
  /// \param stream input stream
  /// \param items_as_row if true, data-items are first gathered (identifier as
  ///                    header, values as row) and then parsed like a loop
  ///                    (OnBeginLoop(), OnDataRow(), OnEndLoop())
  explicit StarParser(std::istream& stream, bool items_as_row=false);
  explicit StarParser(const String& filename, bool items_as_row=false);
  virtual ~StarParser() { }
// callback interface
public:
  /// \brief called when a loop is encountered
  /// 
  /// \return When true, for each row in the loop, OnDataRow() is returned. 
  ///     If false, the loop is skipped.
  virtual bool OnBeginLoop(const StarLoopDesc& header) { return true; }
  /// \brief invoked when leaving a loop
  /// 
  /// OnEndLoop() is only invoked, when OnBeginLoop() returned true.
  virtual void OnEndLoop() { }
  /// \brief invoked when a data row in a loop is encountered.
  /// \param header describes the row format
  /// \param columns contains the data columns
  virtual void OnDataRow(const StarLoopDesc& header, 
                         const std::vector<StringRef>& columns) 
  {
  }
  /// \brief invoked when a data item is encountered
  virtual void OnDataItem(const StarDataItem& item) { }
  
  /// \brief invoked when a new data control structure is encountered
  /// \return When true, for each data-item belonging to the data block, 
  ///     OnDataItem() is invoked. For every loop, OnBeginLoop() is called.
  virtual bool OnBeginData(const StringRef& data_name) { return true; }
  
  /// \brief called when leaving a datasection. Will only be invoked when 
  ///     OnBeginData() returned true.
  virtual void OnEndData() { }

  /// \brief try to convert a value to Real, on failure raise an exception.
  ///
  /// \param data value to be converted
  /// \param name to be included in the message
  ///
  /// \return converted value
  Real TryGetReal(const StringRef& data, const String& name) const;

  /// \brief try to convert a value to float, on failure raise an exception.
  ///
  /// \param data value to be converted
  /// \param name to be included in the message
  ///
  /// \return converted value
  float TryGetFloat(const StringRef& data, const String& name) const;

  /// \brief try to convert a value to float, on failure raise an exception.
  ///
  /// \param data value to be converted
  /// \param name to be included in the message
  /// \param may_fail decides if an exception is raised (false) or not (true)
  ///
  /// \return converted value
  std::pair<bool, float> TryGetFloat(const StringRef& data,
                                     const String& name,
                                     bool may_fail) const;

  /// \brief try to convert a value to integer, on failure raise an exception.
  ///
  /// \param data value to be converted
  /// \param name to be included in the message
  ///
  /// \return converted value
  int TryGetInt(const StringRef& data, const String& name) const;

  /// \brief try to convert a value to integer, exception can be turned off.
  ///
  /// \param data value to be converted
  /// \param name to be included in the message
  /// \param may_fail decides if an exception is raised (false) or not (true)
  ///
  /// \return pair with value and indicator if conversion worked
  std::pair<bool, int> TryGetInt(const StringRef& data,
                                 const String& name,
                                 bool may_fail) const;

  /// \brief try to convert a value to bool, on failure raise an exception.
  ///
  /// \param data value to be converted
  /// \param name to be included in the message
  ///
  /// \return converted value
  bool TryGetBool(const StringRef& data, const String& name) const;

  /// \brief format diagnostic and returns it as a string.
  String FormatDiagnostic(StarDiagType type, const String& message,
                          int line=-1) const;

  void SetFilename(const String& filename)
  {
    filename_ = filename;
  }

  /// \brief retrieve the line, the parser is currently working on
  int GetCurrentLinenum() const
  {
    return line_num_;
  }
public:
  void Parse();
  
public:
  static bool SplitLine(const StringRef& line, 
                        std::vector<StringRef>& parts, bool clear=true);
private:
  void ParseLoop();
  /// \brief Calls the loop parsing functions on the last data item fetched to
  ///        be read as loop
  void ParseLastDataItemRow();
  /// \brief Calls functions for parsing a data item, either as singleton
  ///        (default) or loop.
  void ParseDataItemOrRow(StarDataItem& item);
  /// \brief If enabled, calls ParseLastDataItemRow()
  void ParseEndDataItemRow();
private:
  /// \brief read next line, replacing the current line
  bool NextLine(StringRef& str)
  {
    if (std::getline(stream_, current_line_)) {
      str=StringRef(current_line_.data(), current_line_.length());
      ++line_num_;
      has_current_line_=true;
      return true;
    }
    return false;
  }
  /// \brief return current line, or if no current line is set, read next line
  bool GetLine(StringRef& ref)
  {
    if (has_current_line_) {
      ref=StringRef(current_line_.data(), current_line_.size());
      return true;        
    }
    return this->NextLine(ref);
  }
  
  void ConsumeLine() 
  {
    assert(has_current_line_);
    has_current_line_=false;
  }

  // the reason for having this function are problems when generating 
  // StringRefs that contain pointers to strings in a vector.
  // When the vector grows, the memory of the strings might get 
  // reallocated and the StringRefs point into nothing.
  // This function Prepares the stringrefs given a constant
  // vector and calls OnDataRow to minimize the change in interface.
  void CallOnDataRow(const StarLoopDesc& header, 
                     const std::vector<String>& columns);

  void ParseDataItemIdent(const StringRef ident,
                          StringRef& cat, StringRef& name);
  void ParseGlobal();
  void ParseData();
  void ParseDataItem();
  void DiagnoseUnknown();
  bool ParseMultilineValue(String& value, bool skip=false);
  std::ifstream fstream_;
  boost::iostreams::filtering_stream<boost::iostreams::input> stream_;
  String        filename_;
  int           line_num_;
  bool          has_current_line_;
  String        current_line_;
  bool          items_as_row_;
  StarLoopDesc  items_row_header_;
  bool          file_open_;
  std::vector<String> items_row_values_;
};
 
}}

#endif
