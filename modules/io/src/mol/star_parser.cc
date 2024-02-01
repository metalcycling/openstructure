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

/*
  Author: Marco Biasini
 */
#include <boost/iostreams/filter/gzip.hpp>

#include <cassert>
#include <sstream>
#include <ost/log.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/mol/star_parser.hh>

namespace ost { namespace io {

StarParser::StarParser(std::istream& stream, bool items_as_row):
  filename_("<stream>"), line_num_(0),
  has_current_line_(false), current_line_(),
  items_row_header_(), file_open_(true),
  items_row_values_()
{
  items_as_row_ = items_as_row;
  
  if (!stream) {
    file_open_ = false;
  }

  stream_.push(stream);
}

StarParser::StarParser(const String& filename, bool items_as_row):
  fstream_(filename.c_str()), filename_(filename),
  line_num_(0), has_current_line_(false), current_line_(),
  items_row_header_(), file_open_(true),
  items_row_values_()
{
  items_as_row_=items_as_row;
  if (filename.length() >= 3 &&
      filename.substr(filename.length() - 3) == ".gz") {
    stream_.push(boost::iostreams::gzip_decompressor());
  }

  stream_.push(fstream_);

  if (!fstream_) {
    file_open_ = false;
  }
}

String StarParser::FormatDiagnostic(StarDiagType type, const String& message,
                                    int line) const
{
  std::stringstream ss;
  ss << filename_ << ":";
  if (line!=-1) {
    ss << line << ": ";
  } else {
    ss << " ";
  }
  switch (type) {
    case STAR_DIAG_ERROR:
      ss << "error: ";
      break;
    case STAR_DIAG_WARNING:
      ss << "warning: ";
      break;
  }
  ss << message;
  return ss.str();
}

Real StarParser::TryGetReal(const StringRef& data, const String& name) const
{
  std::pair<bool, Real> value = data.to_float();
  if (!value.first) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                         "Expecting real number for " +
                                               name + ", found '" + data.str() +
                                             "' instead.", line_num_));
  }
  return value.second;
}

float StarParser::TryGetFloat(const StringRef& data, const String& name) const
{
  std::pair<bool, float> value = data.to_float();
  if (!value.first) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                         "Expecting floating point value for " +
                                               name + ", found '" + data.str() +
                                             "' instead.", line_num_));
  }
  return value.second;
}

std::pair<bool, float> StarParser::TryGetFloat(const StringRef& data,
                                               const String& name,
                                               bool may_fail) const
{
  std::pair<bool, float> value = data.to_float();
  if (!value.first) {
    if (!may_fail) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                         "Expecting floating point value for " +
                                               name + ", found '" + data.str() +
                                             "' instead.", line_num_));
    }
    else {
      LOG_WARNING(this->FormatDiagnostic(STAR_DIAG_WARNING,
                                         "Expecting floating point value for " +
                                         name + ", found '" + data.str() +
                                         "' instead.", line_num_));
    }
  }
  return value;
}

int StarParser::TryGetInt(const StringRef& data, const String& name) const
{
  std::pair<bool, int> value = data.to_int();
  if (!value.first) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                         "Expecting integer value for " +
                                               name + ", found '" + data.str() +
                                             "' instead.", line_num_));
  }
  return value.second;
}

std::pair<bool, int> StarParser::TryGetInt(const StringRef& data,
                                           const String& name,
                                           bool may_fail) const
{
  std::pair<bool, int> value = data.to_int();
  if (!value.first) {
    if (!may_fail) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                               "Expecting integer value for " +
                                               name + ", found '" + data.str() +
                                               "' instead.", line_num_));
    } else {
      LOG_WARNING(this->FormatDiagnostic(STAR_DIAG_WARNING,
                                         "Expecting integer value for " +
                                         name + ", found '" + data.str() +
                                         "' instead.", line_num_));
    }
  }
  return value;
}

bool StarParser::TryGetBool(const StringRef& data, const String& name) const
{
  if (data.length() == 1) {
    if (data[0] == 'Y' || data[0] == 'y') {
      return true;
    } else if (data[0] == 'N' || data[0] == 'n') {
      return false;
    }
  }

  throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                           "Expecting Boolean (Y/N) value for "+
                                           name + ", found '" + data.str() +
                                           "' instead.", line_num_));
}

bool StarParser::SplitLine(const StringRef& line, 
                            std::vector<StringRef>& parts, bool clear)
{
  if (clear) {
    parts.clear();    
  }
  const char* s=line.begin();
  while (s!=line.end()) {
    while (isspace(*s)) {
      ++s;
      if (s==line.end()) {
        return true;
      }
    }
    if (*s=='\'' || *s=='"') {
      char delim=*s;
      const char* start=++s;
      while (s!=line.end()) {
        ++s;
        if (delim==*(s-1) && (s==line.end() || isspace(*s))) {
          break;
        }
      }
      parts.push_back(StringRef(start, s-start-1));
    } else {
      const char* start=s;
      while (s!=line.end() && !isspace(*s)) {
        ++s;
      }
      parts.push_back(StringRef(start, s-start));
    }
  }
  return true;
}

bool StarParser::ParseMultilineValue(String& value, bool skip)
{
  std::stringstream valuebuf;
  StringRef line;
  if (!this->GetLine(line)) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                               "Unexpected end of input",
                                               line_num_));
  }
  valuebuf << line.substr(1);
  bool found_semicolon = false;
  while (this->NextLine(line)) {
    StringRef tline=line.rtrim();
    if (!tline.empty() && tline[0]==';') {
      found_semicolon = true;
      break;
    }
    if (!skip) {
      valuebuf << tline << "\n";      
    }
  }
  if (!found_semicolon) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                             "Unterminated multiline value",
                                             line_num_));
  }
  if (!skip) {
    value=valuebuf.str();
  }
  this->ConsumeLine();
  return true;
}

void StarParser::ParseLoop()
{
  StringRef line;
  int prefix_len=0;
  int leave=false;
  StarLoopDesc header;
  this->ConsumeLine();
  while (this->GetLine(line)) {
    StringRef tline=line.rtrim();
    if (tline.empty()) {
      this->ConsumeLine();
      continue;
    }
    switch(tline[0]) {
      case '_':
        if (prefix_len==0) {
          prefix_len=tline.find('.')-tline.begin();
          header.SetCategory(tline.substr(1, prefix_len-1));
        } else {
          if (tline[prefix_len] != '.' ||
              StringRef(header.GetCategory().data(),
                        header.GetCategory().size())!=tline.substr(1,
                                                                prefix_len-1)) {
            throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                                   "Change of category in loop",
                                                     line_num_));
          }
        }
        header.Add(tline.substr(prefix_len+1));
        this->ConsumeLine();
        break;
      case '#':
        this->ConsumeLine();
        break;        
      default:
        leave=true;
        break;
    }
    if (leave) {
      break;
    }
  }
  bool process_rows=this->OnBeginLoop(header);
  std::vector<String> tmp_values;
  // optimized for the common case where all values are present on the same 
  // line.
  while (this->GetLine(line)) {
    StringRef tline=line.rtrim();
    if (tline.empty()) {
      this->ConsumeLine();
      continue;
    }
    switch (tline[0]) {
      case '#':
        this->ConsumeLine();
        break;
      case ';':
        if (process_rows) {
          tmp_values.push_back(String());
          this->ParseMultilineValue(tmp_values.back());
          if (tmp_values.size()==header.GetSize()) {
            this->CallOnDataRow(header, tmp_values);
            tmp_values.clear();
          }          
        } else {
          String s;
          this->ParseMultilineValue(s, true);
        }
        break;
      case '_':
        return;
      case 'd':
          if (tline.length()>=5 && StringRef("data_", 5)==tline.substr(0, 5)) {
            return;
          }
      case 'l':
        if (StringRef("loop_", 5)==tline) {
          return;
        }
      default:
        if (process_rows) {
          std::vector<StringRef> split_elements;
          StarParser::SplitLine(tline, split_elements);
          for(uint i = 0; i < split_elements.size(); ++i) {
            tmp_values.push_back(split_elements[i].str());
          }
          if (tmp_values.size()==header.GetSize()) {
            this->CallOnDataRow(header, tmp_values);           
            tmp_values.clear();
          }
        }
        this->ConsumeLine();
        break;
    }    
  }
  if (process_rows) {
    this->OnEndLoop();    
  }
}

void StarParser::ParseLastDataItemRow()
{
  if (items_row_header_.GetCategory().size() > 0) {
    if (this->OnBeginLoop(items_row_header_)) {
      this->CallOnDataRow(items_row_header_, items_row_values_);
      this->OnEndLoop();
    }
    items_row_values_.clear();
    items_row_header_.Clear();
  }
}

void StarParser::ParseDataItemOrRow(StarDataItem& item)
{
  if (items_as_row_) {
    // header
    if (StringRef(items_row_header_.GetCategory().c_str(),
                 items_row_header_.GetCategory().size())!=item.GetCategory()) {
      this->ParseLastDataItemRow();
      // set category for new section
      items_row_header_.SetCategory(item.GetCategory());
    }
    
    // row
    items_row_header_.Add(item.GetName());
    items_row_values_.push_back(item.GetValue().str());

  } else {
    this->OnDataItem(item);
  }
}

void StarParser::ParseEndDataItemRow()
{
  if (items_as_row_) {
    this->ParseLastDataItemRow();
  }
}

void StarParser::ParseDataItem()
{
  StringRef line;
  this->GetLine(line);
  // optimize for common case when name/value are present on the same line. 
  // We don't have to allocate any additional strings in that case.
  std::vector<StringRef> nv;
  StarParser::SplitLine(line, nv);
  if (nv.size()==1) {
    // remember identifier. 
    String identifier=line.trim().str();
    String value;

    while (this->NextLine(line)) {
      StringRef tline=line.rtrim();
      if (tline.empty()) {
        this->ConsumeLine();
        continue;
      }
      if (tline[0]==';') {
        this->ParseMultilineValue(value);
      } else {
        value=line.str();
        StarParser::SplitLine(StringRef(value.data(), value.length()), 
                              nv, false);
        if (nv.size()!=2) {
          throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                "More than 1 value for data item "+ identifier,
                                                   line_num_));
        }
        this->ConsumeLine();        
      }
      break;
    }
    if (value.empty()) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                               "Unexpected end of input",
                                               line_num_));
    }
    StringRef id_ref(identifier.data(), identifier.size());
    StringRef cat;
    StringRef name;
    StringRef value_ref=StringRef(value.data(),
                                  value.length()).trim();
    this->ParseDataItemIdent(id_ref, cat, name);
    StarDataItem data_item(cat, name, value_ref);
    this->ParseDataItemOrRow(data_item);
  } else {
    if (nv.size()!=2) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                 "More than 1 value for data item "+ line.str(),
                                               line_num_));
    }
    StringRef cat;
    StringRef name;
    this->ParseDataItemIdent(nv[0], cat, name);
    StarDataItem data_item(cat, name, nv[1]);
    this->ParseDataItemOrRow(data_item);
    this->ConsumeLine();
  }
}

void StarParser::ParseDataItemIdent(const StringRef ident,
                                    StringRef& cat, StringRef& name)
{
    StringRef::const_iterator i=ident.find('.');
    if (i == ident.end()) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                           "Invalid data-item identifier '" + ident.str() + "'",
                                               line_num_));
    }
    cat=ident.substr(1, i-ident.begin()-1);
    name=ident.substr(i-ident.begin()+1);
}

void StarParser::ParseData()
{
  StringRef line;
  this->GetLine(line);
  StringRef data_id=line.rtrim().substr(5);
  bool skip=!this->OnBeginData(data_id);
  this->ConsumeLine();
  while (this->GetLine(line)) {
    StringRef tline=line.rtrim();
    if (tline.empty()) {
      this->ConsumeLine();
      continue;
    }
    switch (tline[0]) {
      case '_':
        if (skip) {
          this->ConsumeLine();
        } else {
          this->ParseDataItem();          
        }
        break;
      case 'd':
        if (tline.length()>=5 && StringRef("data_", 5)==tline.substr(0, 5)) {
          this->ParseEndDataItemRow();
          this->OnEndData();
          return;
        }
      case ';':
        if (skip) {
          String s;
          this->ParseMultilineValue(s, true);
        }
        break;
      case 'l':
        if (tline==StringRef("loop_", 5)) {
          this->ParseEndDataItemRow();
          this->ParseLoop();
          break;
        }
        this->DiagnoseUnknown();
        break;
      case '#':
        // consume comment
        this->ConsumeLine();
        break;
      default:
        return;
    }
  }
  this->ParseEndDataItemRow();
  this->OnEndData();
}

void StarParser::DiagnoseUnknown()
{
  std::stringstream ss;
  StringRef line;
  this->GetLine(line);
  ss << "unknown control structure '"<< line.rtrim() << "'";
  throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR, ss.str(),
                                           line_num_));
}

void StarParser::ParseGlobal()
{
  assert(0 && "unhandled global section");
}

void StarParser::Parse()
{
  if (!file_open_) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                             "Failed to open file '" +
                                             filename_ + "'!"));
  }
  StringRef line;
  std::stringstream ss;
  while (this->GetLine(line)) {
    StringRef tline=line.rtrim();
    if (tline.empty()) {
      this->ConsumeLine();
      continue;
    }
    switch (tline[0]) {
      case 'd':
        if (tline.length()>=5 && StringRef("data_", 5)==tline.substr(0, 5)) {
          this->ParseData();
        } else {
          throw IOException("Missing 'data_' control structure");
        }
        break;
      case 'g':
        if (tline.length()>=7 && StringRef("global_", 7)==tline.substr(0, 7)) {
          this->ParseGlobal();
        } else {
          throw IOException("Missing 'global_' control structure");
        }
        break;
      case '#':
        // consume comment
        this->ConsumeLine();
        break;
      default:
        throw IOException("Missing 'data_' control structure");
        break;
    }
  }
}

void StarParser::CallOnDataRow(const StarLoopDesc& header,
                               const std::vector<String>& columns) {

  std::vector<StringRef> string_refs;
  string_refs.reserve(columns.size());
  for(uint i = 0; i < columns.size(); ++i) {
    string_refs.push_back(StringRef(columns[i].data(), 
                                    columns[i].length()).trim());
  }
  this->OnDataRow(header, string_refs);
}

}}
