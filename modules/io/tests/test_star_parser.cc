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
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <math.h>
#include <ost/mol/mol.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/mol/star_parser.hh>



using namespace ost;
using namespace ost::io;


class DataItemTestParser: public StarParser {
public:
  DataItemTestParser(std::istream& stream): StarParser(stream)
  { }

  DataItemTestParser(const String& filename): StarParser(filename)
  { }
  virtual void OnDataItem(const StarDataItem& item)
  {
    BOOST_CHECK_EQUAL(item.GetCategory().str(), "data-item");
    if (item.GetName()==StringRef("s1", 2)) {
      s1=item.GetValue().str();
    }      
    if (item.GetName()==StringRef("s2", 2)) {
      s2=item.GetValue().str();
    }
    if (item.GetName()==StringRef("s3", 2)) {
      s3=item.GetValue().str();
    }
    if (item.GetName()==StringRef("s4", 2)) {
      s4=item.GetValue().str();
    }    
  }
  virtual bool OnBeginLoop(const StarLoopDesc& header)
  {
    return true;
  }
  virtual void OnDataRow(const StarLoopDesc& header, 
                         const std::vector<StringRef>& columns) 
  {
    BOOST_CHECK_EQUAL(header.GetCategory(), "loop");
    BOOST_CHECK_EQUAL(columns.size(), size_t(4));
  }
  String s1;
  String s2;
  String s3;
  String s4;  
};

class LoopTestParser: public StarParser {
public:
  LoopTestParser(std::istream& stream): StarParser(stream), count_(0)
  { }
  virtual bool OnBeginLoop(const StarLoopDesc& header)
  {
    BOOST_CHECK_EQUAL(header.GetCategory(), "loop");
    return true;
  }
  virtual void OnDataRow(const StarLoopDesc& header, 
                         const std::vector<StringRef>& columns) 
  {
    BOOST_CHECK_EQUAL(header.GetCategory(), "loop");
    BOOST_CHECK_EQUAL(columns.size(), size_t(3));
    lines[count_][0]=columns[0].str();
    lines[count_][1]=columns[1].str();
    lines[count_][2]=columns[2].str();
    count_++;            
  }
  String lines[5][3];
private:
  int count_;
};

class MultiDataParser : public StarParser {
public:
  MultiDataParser(std::istream& stream): StarParser(stream), visit_one(false),
     visit_two(false)
  { }
  
  virtual bool OnBeginData(const StringRef& data)
  {
    if (data==StringRef("one", 3)) {
      visit_one=true;
    }
    if (data==StringRef("two", 3)) {
      visit_two=true;
    }
    return true;
  }
  
  bool visit_one;
  bool visit_two;
};

class HardLoopTestParser : public StarParser {
public:
  HardLoopTestParser(std::istream& stream): StarParser(stream), cur_char_('A')
  { }
  
  virtual bool OnBeginLoop(const StarLoopDesc& header)
  {
    return true;
  }  

  virtual void OnDataRow(const StarLoopDesc& header, 
                         const std::vector<StringRef>& columns) 
  {
    BOOST_CHECK_EQUAL(columns[0][0], cur_char_);
    BOOST_CHECK_EQUAL(columns[0].size(), static_cast<size_t>(1));
    ++cur_char_;
    BOOST_CHECK_EQUAL(columns[1][0], cur_char_);
    BOOST_CHECK_EQUAL(columns[1].size(), static_cast<size_t>(1));
    ++cur_char_;
    BOOST_CHECK_EQUAL(columns[2][0], cur_char_);
    BOOST_CHECK_EQUAL(columns[2].size(), static_cast<size_t>(1));    
    ++cur_char_;
    BOOST_CHECK_EQUAL(columns[3][0], cur_char_);
    BOOST_CHECK_EQUAL(columns[3].size(), static_cast<size_t>(1));    
    ++cur_char_;
    BOOST_CHECK_EQUAL(columns[4][0], cur_char_);
    BOOST_CHECK_EQUAL(columns[4].size(), static_cast<size_t>(1));    
    ++cur_char_;
    BOOST_CHECK_EQUAL(columns[5][0], cur_char_);
    BOOST_CHECK_EQUAL(columns[5].size(), static_cast<size_t>(1));
    ++cur_char_;
  }
  
  virtual void OnEndLoop()
  {
    BOOST_CHECK_EQUAL(cur_char_, 'S');
  }
  char cur_char_;
};

class ItemsAsRowTestParser : public StarParser {
public:
  ItemsAsRowTestParser(std::istream& stream): StarParser(stream, true),
                                              category("")
  { }

  virtual bool OnBeginLoop(const StarLoopDesc& header)
  {
    if ((header.GetCategory() == "data-item1") ||
        (header.GetCategory() == "data-item2") ||
        (header.GetCategory() == "data-item3") ||
        (header.GetCategory() == "data-item4")) {
      return true;
    }
    return false;
  }  

  virtual void OnDataRow(const StarLoopDesc& header, 
                         const std::vector<StringRef>& columns) 
  {
    BOOST_CHECK_EQUAL(columns.size(), size_t(4));
    category = header.GetCategory();
    if (header.GetCategory() == "data-item1") {
      s1 = columns[header.GetIndex("s1")].str();
      s2 = columns[header.GetIndex("s2")].str();
      s3 = columns[header.GetIndex("s3")].str();
      s4 = columns[header.GetIndex("s4")].str();
      return;
    }
    if (header.GetCategory() == "data-item2") {
      i1 = columns[header.GetIndex("i1")].to_int().second;
      i2 = columns[header.GetIndex("i2")].to_int().second;
      i3 = columns[header.GetIndex("i3")].to_int().second;
      i4 = columns[header.GetIndex("i4")].to_int().second;
      return;
    }
    if (header.GetCategory() == "data-item3") {
      f1 = columns[header.GetIndex("f1")].to_float().second;
      f2 = columns[header.GetIndex("f2")].to_float().second;
      f3 = columns[header.GetIndex("f3")].to_float().second;
      f4 = columns[header.GetIndex("f4")].to_float().second;
      return;
    }
    if (header.GetCategory() == "data-item4") {
      s5 = columns[header.GetIndex("s1")].str();
      s6 = columns[header.GetIndex("s2")].str();
      s7 = columns[header.GetIndex("s3")].str();
      s8 = columns[header.GetIndex("s4")].str();
      return;
    }
  }

  void OnEndLoop()
  {
    if (category == "data-item3") {
      f1 = ceil(f1*2);
      f2 = ceil(f2*2);
      f3 = ceil(f3*2);
      f4 = ceil(f4*2);
    }
  }

  String category;

  String s1;
  String s2;
  String s3;
  String s4;
  String s5;
  String s6;
  String s7;
  String s8;

  int i1;
  int i2;
  int i3;
  int i4;

  float f1;
  float f2;
  float f3;
  float f4;
};

BOOST_AUTO_TEST_SUITE( io );

BOOST_AUTO_TEST_CASE(star_split1)
{
  const char* line1="     ";
  StringRef r1(line1, strlen(line1));
  std::vector<StringRef> parts;
  BOOST_REQUIRE(StarParser::SplitLine(r1, parts));
  BOOST_REQUIRE(parts.size()==0);
}

BOOST_AUTO_TEST_CASE(star_split2)
{
  const char* line2="1 2";
  std::vector<StringRef> parts;  
  StringRef r2(line2, strlen(line2));
  BOOST_REQUIRE(StarParser::SplitLine(r2, parts));
  BOOST_REQUIRE(parts.size()==2);
  BOOST_CHECK_EQUAL(parts[0], StringRef("1", 1));
  BOOST_CHECK_EQUAL(parts[1], StringRef("2", 1));
}

BOOST_AUTO_TEST_CASE(star_split3)
{
  const char* line3="'1 2'  ";
  std::vector<StringRef> parts;  
  StringRef r3(line3, strlen(line3));
  BOOST_REQUIRE(StarParser::SplitLine(r3, parts));
  BOOST_REQUIRE(parts.size()==1);
  BOOST_CHECK_EQUAL(parts[0], StringRef("1 2", 3));
}

BOOST_AUTO_TEST_CASE(star_data_item)
{
  BOOST_TEST_MESSAGE("  Running star_data_item tests...");
  std::ifstream s("testfiles/data-item.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_NO_THROW(star_p.Parse());
  BOOST_CHECK_EQUAL(star_p.s1, "a");
  BOOST_CHECK_EQUAL(star_p.s2, "a b c");
  BOOST_CHECK_EQUAL(star_p.s3, "a\nb\nc");
  BOOST_CHECK_EQUAL(star_p.s4, "a'b");
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(format_diag_stream)
{
  BOOST_TEST_MESSAGE("  Running format_diag_stream tests...");
  std::ifstream s("testfiles/data-item.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_EQUAL(star_p.FormatDiagnostic(STAR_DIAG_WARNING, "bad", -1),
                    "<stream>: warning: bad");
  BOOST_CHECK_EQUAL(star_p.FormatDiagnostic(STAR_DIAG_ERROR, "really bad", -1),
                    "<stream>: error: really bad");
  BOOST_CHECK_EQUAL(star_p.FormatDiagnostic(STAR_DIAG_ERROR, "bad", 55),
                    "<stream>:55: error: bad");
  star_p.SetFilename("testname");
  BOOST_CHECK_EQUAL(star_p.FormatDiagnostic(STAR_DIAG_ERROR, "bad", 55),
                    "testname:55: error: bad");
}

BOOST_AUTO_TEST_CASE(format_diag_filename)
{
  BOOST_TEST_MESSAGE("  Running format_diag_filename tests...");
  DataItemTestParser star_p("testfiles/data-item.cif");
  BOOST_CHECK_EQUAL(star_p.FormatDiagnostic(STAR_DIAG_WARNING, "bad", -1),
                    "testfiles/data-item.cif: warning: bad");
  BOOST_CHECK_EQUAL(star_p.FormatDiagnostic(STAR_DIAG_ERROR, "really bad", -1),
                    "testfiles/data-item.cif: error: really bad");
  BOOST_CHECK_EQUAL(star_p.FormatDiagnostic(STAR_DIAG_ERROR, "bad", 55),
                    "testfiles/data-item.cif:55: error: bad");
}

BOOST_AUTO_TEST_CASE(star_multi)
{
  std::ifstream s("testfiles/multi-data.cif");
  MultiDataParser star_p(s);
  star_p.Parse();
  BOOST_CHECK(star_p.visit_one);
  BOOST_CHECK(star_p.visit_two);
}

BOOST_AUTO_TEST_CASE(star_multiline_loop)
{
  std::ifstream s("testfiles/multiline-loop.cif");
  HardLoopTestParser star_p(s);
  star_p.Parse();
}
  
BOOST_AUTO_TEST_CASE(star_loop)
{
  BOOST_TEST_MESSAGE("  Running star_loop tests...");
  std::ifstream s("testfiles/loop.cif");
  LoopTestParser star_p(s);
  star_p.Parse();
  BOOST_CHECK_EQUAL(star_p.lines[0][0], "d");
  BOOST_CHECK_EQUAL(star_p.lines[0][1], "2");
  BOOST_CHECK_EQUAL(star_p.lines[0][2], "3");
  
  BOOST_CHECK_EQUAL(star_p.lines[1][0], "4 x");
  BOOST_CHECK_EQUAL(star_p.lines[1][1], "5");
  BOOST_CHECK_EQUAL(star_p.lines[1][2], "6");
  
  BOOST_CHECK_EQUAL(star_p.lines[2][0], "7");
  BOOST_CHECK_EQUAL(star_p.lines[2][1], "8");
  BOOST_CHECK_EQUAL(star_p.lines[2][2], "9");

  BOOST_CHECK_EQUAL(star_p.lines[3][0], "10");
  BOOST_CHECK_EQUAL(star_p.lines[3][1], "11");
  BOOST_CHECK_EQUAL(star_p.lines[3][2], "12");
    
  BOOST_CHECK_EQUAL(star_p.lines[4][0], "13");
  BOOST_CHECK_EQUAL(star_p.lines[4][1], "14");
  BOOST_CHECK_EQUAL(star_p.lines[4][2], "15");
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(star_items_as_row)
{
  BOOST_TEST_MESSAGE("  Running star_items_as_row tests...");
  std::ifstream s("testfiles/items-as-row.cif");
  ItemsAsRowTestParser star_p(s);
  BOOST_TEST_MESSAGE("          checking parsing abilities...");
  star_p.Parse();
  BOOST_CHECK_EQUAL(star_p.s1, "a");
  BOOST_CHECK_EQUAL(star_p.s2, "a b c");
  BOOST_CHECK_EQUAL(star_p.s3, "a\nb\nc");
  BOOST_CHECK_EQUAL(star_p.s4, "a'b");
  BOOST_CHECK_EQUAL(star_p.i1, 1);
  BOOST_CHECK_EQUAL(star_p.i2, 2);
  BOOST_CHECK_EQUAL(star_p.i3, 3);
  BOOST_CHECK_EQUAL(star_p.i4, 4);
  BOOST_CHECK_EQUAL(ceil(star_p.f1), ceil(2.22));
  BOOST_CHECK_EQUAL(ceil(star_p.f2), ceil(2.44));
  BOOST_CHECK_EQUAL(ceil(star_p.f3), ceil(2.66));
  BOOST_CHECK_EQUAL(ceil(star_p.f4), ceil(2.88));

  BOOST_CHECK_EQUAL(star_p.s5, "1.11");
  BOOST_CHECK_EQUAL(star_p.s6, "1.22");
  BOOST_CHECK_EQUAL(star_p.s7, "1.33");
  BOOST_CHECK_EQUAL(star_p.s8, "1.44");
  BOOST_TEST_MESSAGE("          done.");
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(star_broken_data)
{
  std::ifstream s("testfiles/broken_data.cif");
  LoopTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_broken_global)
{
  std::ifstream s("testfiles/broken_global.cif");
  LoopTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_missing_data)
{
  std::ifstream s("testfiles/missing_data.cif");
  LoopTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_unterminated_dataitem)
{
  std::ifstream s("testfiles/unterminated_dataitem.cif");
  LoopTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_incomplete_dataitem)
{
  std::ifstream s("testfiles/incomplete_dataitem.cif");
  LoopTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_singleline_multiple_values)
{
  std::ifstream s("testfiles/singleline_multivalue_dataitem.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_multiline_multiple_values)
{
  std::ifstream s("testfiles/multiline_multivalue_dataitem.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_multiline_invalid_ident)
{
  std::ifstream s("testfiles/multiline_invalid_ident.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_singleline_invalid_ident)
{
  std::ifstream s("testfiles/singleline_invalid_ident.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_loop_category_change)
{
  std::ifstream s("testfiles/loop_category_change.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_loop_category_change_inplace)
{
  std::ifstream s("testfiles/loop_category_change_inplace.cif");
  DataItemTestParser star_p(s);
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
}

BOOST_AUTO_TEST_CASE(star_try_float_conversions)
{
  BOOST_TEST_MESSAGE("  Running star_try_float_conversions tests...");
  std::ifstream s("testfiles/loop_category_change_inplace.cif");
  DataItemTestParser star_p(s);
  StringRef data = StringRef("1.5", 3);
  BOOST_CHECK_CLOSE(star_p.TryGetFloat(data,"positive float test"),1.5f,0.001f);
  data = StringRef("foo", 3);
  BOOST_CHECK_THROW(star_p.TryGetFloat(data, "negative float test"),
                    IOException);
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(star_try_int_conversions)
{
  BOOST_TEST_MESSAGE("  Running star_try_int_conversions tests...");
  std::ifstream s("testfiles/loop_category_change_inplace.cif");
  DataItemTestParser star_p(s);
  StringRef data = StringRef("101", 3);
  BOOST_CHECK_EQUAL(star_p.TryGetInt(data, "positive int test"), 101);
  data = StringRef("foo", 3);
  BOOST_CHECK_THROW(star_p.TryGetInt(data, "negative int test"),
                    IOException);
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(star_try_bool_conversions)
{
  BOOST_TEST_MESSAGE("  Running star_try_bool_conversions tests...");
  std::ifstream s("testfiles/loop_category_change_inplace.cif");
  DataItemTestParser star_p(s);
  StringRef data = StringRef("Y", 1);
  BOOST_CHECK(star_p.TryGetBool(data, "positive bool test ("+data.str()+")"));
  data = StringRef("y", 1);
  BOOST_CHECK(star_p.TryGetBool(data, "positive bool test ("+data.str()+")"));
  data = StringRef("N", 1);
  BOOST_CHECK(!star_p.TryGetBool(data, "positive bool test ("+data.str()+")"));
  data = StringRef("n", 1);
  BOOST_CHECK(!star_p.TryGetBool(data, "positive bool test ("+data.str()+")"));
  data = StringRef("J", 1);
  BOOST_CHECK_THROW(star_p.TryGetInt(data,
                                     "negative bool test ("+data.str()+")"),
                    IOException);
  data = StringRef("Foo", 3);
  BOOST_CHECK_THROW(star_p.TryGetInt(data,
                                     "negative bool test ("+data.str()+")"),
                    IOException);
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_CASE(star_wrong_filename)
{
  BOOST_TEST_MESSAGE("  Running star_wrong_filename tests...");
  DataItemTestParser star_p("testfiles/doesnotexist.foo");
  BOOST_CHECK_THROW(star_p.Parse(), IOException);
  BOOST_TEST_MESSAGE("  done.");
}

BOOST_AUTO_TEST_SUITE_END();

