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
  Author: Ansgar Philippsen
*/

#include <ost/gfx/color.hh>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using boost::unit_test_framework::test_suite;

using namespace ost;
using namespace ost::gfx;

static const float tolerance=1e-4;

namespace {
  bool compare_colors(const Color& c1, const Color& c2, float tol=1e-6) {
    return std::abs<float>(c1.GetRed()-c2.GetRed()<tol) &&
      std::abs<float>(c1.GetGreen()-c2.GetGreen()<tol) &&
      std::abs<float>(c1.GetBlue()-c2.GetBlue()<tol) &&
      std::abs<float>(c1.GetAlpha()-c2.GetAlpha()<tol) &&
      std::abs<float>(c1.GetHue()-c2.GetHue()<tol) &&
      std::abs<float>(c1.GetSat()-c2.GetSat()<tol) &&
      std::abs<float>(c1.GetVal()-c2.GetVal()<tol);
  }
}

#define FIXED_BOOST_WARN_CLOSE(F1,F2,F3) BOOST_WARN_CLOSE(static_cast<float>(F1), static_cast<float>(F2), static_cast<float>(F3))

BOOST_AUTO_TEST_SUITE(gfx)

BOOST_AUTO_TEST_CASE(default_color)
{
  Color c;
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetAlpha(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
}

BOOST_AUTO_TEST_CASE(set_rgb)
{
  Color c;
  c = RGB(1.0f,0.0f,0.0f); // red
  FIXED_BOOST_WARN_CLOSE(c.GetAlpha(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),0.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = RGB(1.0f,1.0f,0.0f); // yellow
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),1.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = RGB(0.0f,1.0f,0.0f); // green
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),2.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = RGB(0.0f,1.0f,1.0f); // cyan
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),3.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = RGB(0.0f,0.0f,1.0f); // blue
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),4.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = RGB(1.0f,0.0f,1.0f); // purple
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),5.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
}

BOOST_AUTO_TEST_CASE(set_hsv)
{
  Color c;
  c = HSV(0.0f/6.0,1.0f,1.0f); // red
  FIXED_BOOST_WARN_CLOSE(c.GetAlpha(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),0.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(0.0f/6.0,0.5,1.0f); // light red
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),0.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(1.0f/6.0,1.0f,1.0f); // yellow
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),1.0f/6.0,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(2.0/6.0,1.0f,1.0f); // green
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),2.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(2.0/6.0,0.5,1.0f); // light green
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),2.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(3.0/6.0,1.0f,1.0f); // cyan
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),3.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(4.0/6.0,1.0f,1.0f); // blue
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),4.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(4.0/6.0,0.5,1.0f); // light blue
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),4.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),0.5f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
  c = HSV(5.0/6.0,1.0f,1.0f); // purple
  FIXED_BOOST_WARN_CLOSE(c.GetRed(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetGreen(),0.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetBlue(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetHue(),5.0f/6.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetSat(),1.0f,tolerance);
  FIXED_BOOST_WARN_CLOSE(c.GetVal(),1.0f,tolerance);
}

BOOST_AUTO_TEST_CASE(set_char)
{
  BOOST_WARN(compare_colors(RGBb(255,0,0),RGB(1.0f,0.0f,0.0f)));
  BOOST_WARN(compare_colors(RGBAb(0,127,0,127),RGBA(0.0f,0.5f,0.0f,0.5f)));
}

BOOST_AUTO_TEST_CASE(set_int)
{
  BOOST_WARN(compare_colors(RGBi(65535,65535,0),RGB(1.0f,1.0f,0.0f)));
  BOOST_WARN(compare_colors(RGBAi(0,32767,0,32767),RGBA(0.0f,0.5f,0.0f,0.5f)));
}

BOOST_AUTO_TEST_SUITE_END()
