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
#ifndef OST_GEOM_VEC_MAT_PREDICATES_HH
#define OST_GEOM_VEC_MAT_PREDICATES_HH
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION<105900
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <ost/geom/geom.hh>


/*
  This file contains several predicates to check for equality of vectors and 
  matrices. Usage:
  
  BOOST_CHECK(vec2_is_close(geom::Vec2(1,2), geom::Vec2(1, 2)));
 */

template<class V>
boost::test_tools::predicate_result vec_is_close(const V& v1, const V& v2, 
                                                 Real tolerance, 
                                                 unsigned int dim)
{
  std::string labels[]={"x","y","z","w"};
  bool flag=true;
  boost::test_tools::predicate_result res( false );
#if BOOST_VERSION<105900
  boost::test_tools::close_at_tolerance<Real> close_test(boost::test_tools::percent_tolerance(tolerance));
#else
  boost::math::fpc::close_at_tolerance<Real> close_test(boost::math::fpc::percent_tolerance(tolerance));
#endif
  for(unsigned int i=0;i<dim;++i){
    if(v1[i]==0.0){
      if(!boost::test_tools::check_is_small(v2[i],tolerance)){
       flag=false;
       res.message() << "Value for " << labels[i] << " differs: ( " << v1[i] << " != " << v2[i] << " ). ";
      }
    }else if (v2[i]==0.0){
      if(!boost::test_tools::check_is_small(v1[i],tolerance)){
       flag=false;
       res.message() << "Value for "<< labels[i] << " differs: ( " << v1[i] << " != " << v2[i] << " ). ";
      }
    }else{
      if(!close_test(v1[i],v2[i])){
       flag=false;
       res.message() << "Value for " << labels[i] << " differs: ( " << v1[i] << " != " << v2[i] << " ). ";
      }
    }
  }
  if(flag){
    return true;
  }else{
    return res;
  }
}  
  
boost::test_tools::predicate_result vec2_is_close(const geom::Vec2& v1, 
                                                  const geom::Vec2& v2, 
                                                  Real tolerance=geom::EPSILON)
{
  return vec_is_close<geom::Vec2>(v1,v2,tolerance,2);
}

boost::test_tools::predicate_result vec3_is_close(const geom::Vec3& v1, 
                                                  const geom::Vec3& v2, 
                                                  Real tolerance=geom::EPSILON)
{
  return vec_is_close<geom::Vec3>(v1,v2,tolerance,3);
}

boost::test_tools::predicate_result vec4_is_close(const geom::Vec4& v1, 
                                                  const geom::Vec4& v2, 
                                                  Real tolerance=geom::EPSILON)
{
  return vec_is_close<geom::Vec4>(v1,v2,tolerance,4);
}

template<class M>
boost::test_tools::predicate_result mat_is_close(const M& m1, const M& m2, 
                                                 Real tolerance, 
                                                 unsigned int dim)
{
  bool flag=true;
  boost::test_tools::predicate_result res( false );
#if BOOST_VERSION<105900
  boost::test_tools::close_at_tolerance<Real> close_test(boost::test_tools::percent_tolerance(tolerance));
#else
  boost::math::fpc::close_at_tolerance<Real> close_test(boost::math::fpc::percent_tolerance(tolerance));
#endif
  for(unsigned int i=0;i<dim;++i){
    for(unsigned int j=0;j<dim;++j){
      if(m1(i,j)==0.0){
        if(!boost::test_tools::check_is_small(m2(i,j),tolerance)){
         flag=false;
         res.message() << "Value for (" << i << "," << j << ") differs: (" << m1(i,j) << "!=" << m2(i,j) << "). ";
        }
      }else if (m2(i,j)==0.0){
        if(!boost::test_tools::check_is_small(m1(i,j),tolerance)){
         flag=false;
         res.message() << "Value for (" << i << "," << j << ") differs: (" << m1(i,j) << "!=" << m2(i,j) << "). ";
        }
      }else{
        if(!close_test(m1(i,j),m2(i,j))){
         flag=false;
         res.message() << "Value for (" << i << "," << j << ") differs: (" << m1(i,j) << "!=" << m2(i,j) << "). ";
        }
      }
    }
  }
  if(flag){
    return true;
  }else{
    return res;
  }
}  

boost::test_tools::predicate_result mat2_is_close(const geom::Mat2& m1, 
                                                  const geom::Mat2& m2, 
                                                  Real tolerance=geom::EPSILON)
{
  return mat_is_close<geom::Mat2>(m1,m2,tolerance,2);
}

boost::test_tools::predicate_result mat3_is_close(const geom::Mat3& m1, 
                                                  const geom::Mat3& m2, 
                                                  Real tolerance=geom::EPSILON)
{
  return mat_is_close<geom::Mat3>(m1,m2,tolerance,3);
}

boost::test_tools::predicate_result mat4_is_close(const geom::Mat4& m1, 
                                                  const geom::Mat4& m2, 
                                                  Real tolerance=geom::EPSILON)
{
  return mat_is_close<geom::Mat4>(m1,m2,tolerance,4);
}

#endif
