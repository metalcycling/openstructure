/////////////////////////////////////////////////////////
// IPLT - Image Processing Library & Toolbox
// Copyright (c) 2003-2007 University of Basel
//
// The usage terms and conditions of this software 
// are governed by the GNU General Public License, 
// as described in the accompanying text file LICENSE.
// Please review this information carefully and contact
// info@iplt.org should questions arise.
/////////////////////////////////////////////////////////

/*
  Author: Ansgar Philippsen
*/

#ifndef GEOM_HELPER_HH
#define GEOM_HELPER_HH

#include <cstdlib>
#include <cmath>

#include <boost/random.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>
using boost::unit_test::test_suite;

#include <ost/geom/geom.hh>

using namespace geom;


inline bool dequal(Real a, Real b)
{
  return std::fabs(a-b)<EPSILON;
}
inline
bool match(const Vec2& v, Real x, Real y)
{
  return std::fabs(v[0]-x)<EPSILON && std::fabs(v[1]-y)<EPSILON;
}

inline
bool match(const Vec3& v, Real x, Real y, Real z)
{
  return std::fabs(v[0]-x)<EPSILON && std::fabs(v[1]-y)<EPSILON && std::fabs(v[2]-z)<EPSILON;
}

inline
bool match(const Vec4& v, Real x, Real y, Real z, Real w)
{
  return std::fabs(v[0]-x)<EPSILON && std::fabs(v[1]-y)<EPSILON && std::fabs(v[2]-z)<EPSILON && std::fabs(v[3]-w)<EPSILON;
}

inline  
bool match(const Mat2& m, Real i00, Real i01,  Real i10, Real i11)
{
  return std::fabs(m(0,0)-i00)<EPSILON &&
    std::fabs(m(0,1)-i01)<EPSILON &&
    std::fabs(m(1,0)-i10)<EPSILON &&
    std::fabs(m(1,1)-i11)<EPSILON;
}

inline
bool match(const Mat3& m, Real i00, Real i01, Real i02, Real i10, Real i11, Real i12, Real i20, Real i21, Real i22)
{
  return std::fabs(m(0,0)-i00)<EPSILON &&
    std::fabs(m(0,1)-i01)<EPSILON &&
    std::fabs(m(0,2)-i02)<EPSILON &&
    std::fabs(m(1,0)-i10)<EPSILON &&
    std::fabs(m(1,1)-i11)<EPSILON &&
    std::fabs(m(1,2)-i12)<EPSILON &&
    std::fabs(m(2,0)-i20)<EPSILON &&
    std::fabs(m(2,1)-i21)<EPSILON &&
    std::fabs(m(2,2)-i22)<EPSILON;

}
  
inline
bool match(const Mat4& m, Real i00, Real i01, Real i02, Real i03, Real i10, Real i11, Real i12, Real i13, Real i20, Real i21, Real i22, Real i23, Real i30, Real i31, Real i32, Real i33)
{
  static Real EPSILON=1e-10;
  return std::fabs(m(0,0)-i00)<EPSILON &&
    std::fabs(m(0,1)-i01)<EPSILON &&
    std::fabs(m(0,2)-i02)<EPSILON &&
    std::fabs(m(0,3)-i03)<EPSILON &&
    std::fabs(m(1,0)-i10)<EPSILON &&
    std::fabs(m(1,1)-i11)<EPSILON &&
    std::fabs(m(1,2)-i12)<EPSILON &&
    std::fabs(m(1,3)-i13)<EPSILON &&
    std::fabs(m(2,0)-i20)<EPSILON &&
    std::fabs(m(2,1)-i21)<EPSILON &&
    std::fabs(m(2,2)-i22)<EPSILON &&
    std::fabs(m(2,3)-i23)<EPSILON &&
    std::fabs(m(3,0)-i30)<EPSILON &&
    std::fabs(m(3,1)-i31)<EPSILON &&
    std::fabs(m(3,2)-i32)<EPSILON &&
    std::fabs(m(3,3)-i33)<EPSILON;
}

namespace {
  boost::mt19937 RandomGenerator;
  boost::uniform_01<boost::mt19937> UniformRandom(RandomGenerator);
}

inline  
Vec3 rnd_vec3()
{
  return Vec3(UniformRandom(),UniformRandom(),UniformRandom());
}

inline  
Mat3 rnd_mat3()
{
  return Mat3(UniformRandom(),UniformRandom(),UniformRandom(),
	      UniformRandom(),UniformRandom(),UniformRandom(),
	      UniformRandom(),UniformRandom(),UniformRandom());
}

inline  
Mat2 rnd_mat2()
{
  return Mat2(UniformRandom(),UniformRandom(),
	          UniformRandom(),UniformRandom());
}



#endif
