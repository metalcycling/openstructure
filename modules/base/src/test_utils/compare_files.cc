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

#include <iostream>
#include <fstream>
#include "compare_files.hh"

namespace ost {
bool compare_files(const String& test, const String& gold_standard)
{
  std::ifstream test_stream(test.c_str());
  if (!test_stream) {
    std::cerr << "output file '" << test << "' doesn't exist" << std::endl;
    return false;
  }
  std::ifstream gold_stream(gold_standard.c_str());
  if (!gold_stream) {
    std::cerr << "gold standard file '" << gold_standard 
              << "' doesn't exist" << std::endl;
    return false;
  }
  String test_line, gold_line;
  while (true) {
    bool test_read = static_cast<bool>(std::getline(test_stream, test_line));
    bool gold_read = static_cast<bool>(std::getline(gold_stream, gold_line));
    if (!test_read && !gold_read) {
      // nothing to read anymore in any of the files
      return true;
    }
    if (gold_read && !test_read) {
      std::cerr << gold_standard << " contains additional line(s):"
                << std::endl << gold_line << std::endl;
      return false;
    }
    if (test_read && !gold_read) {
      std::cerr << test << " contains additional line(s):"
                << std::endl << test_line << std::endl;
      return false;
    }
    if (gold_line!=test_line) {
      std::cerr << "line mismatch:" << std::endl << "test: " << test_line
                << std::endl << "gold: " << gold_line << std::endl;
      return false;
    }
  }
  return true;
}

}
