//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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

#ifndef IMG_IMAGE_STATE_DISPATCH_HH
#define IMG_IMAGE_STATE_DISPATCH_HH

#include <ost/message.hh>
#include <string>

#include <ost/message.hh>

#include "image_state_impl.hh"
#include "image_state_def.hh"

namespace ost { namespace img { namespace image_state { 

class DLLEXPORT DispatchException: public Error {
public:
  DispatchException(const String& m):
    Error(m)
  {}
};

template <typename T, class D>
class ImageStateImpl;

namespace dispatch {

////////////////////////////////////////////
// unary RTTI dispatch 
// used instead of visitor pattern

//! In-place unary dispatch
/*!
  The template parameter must be a class with the
  following specifications:

  template <typename T, class D>
  struct FNC {
    void operator()(ImageStateImpl<T,D>* isi);
  };

  This allows partial specialization to be used
*/
template <template <typename T, class D> class FNC>
struct unary_dispatch_ip {

#define ISI_UNARY_DISPATCH_IP_CALL(TYPE,fnc,isb) \
  else if(TYPE * is = dynamic_cast<TYPE *>(isb)) fnc(is)

  void operator()(ImageStateBase* isb) {
    if(false); 
    ISI_UNARY_DISPATCH_IP_CALL(RealSpatialImageState,call,isb);
    ISI_UNARY_DISPATCH_IP_CALL(ComplexSpatialImageState,call,isb);
    ISI_UNARY_DISPATCH_IP_CALL(RealFrequencyImageState,call,isb);
    ISI_UNARY_DISPATCH_IP_CALL(ComplexFrequencyImageState,call,isb);
    ISI_UNARY_DISPATCH_IP_CALL(ComplexHalfFrequencyImageState,call,isb);
    ISI_UNARY_DISPATCH_IP_CALL(WordSpatialImageState,call,isb);
    else throw DispatchException("unknown image state in unary_dispatch_ip");
  }

#undef ISI_UNARY_DISPATCH_IP_CALL
  
  template <typename T, class D>
  void call(ImageStateImpl<T,D>* isi) {
    static FNC<T,D> fnc;
    fnc(isi);
  }

};


//! Out-of-place unary dispatch
/*!
  The template parameter must be a class with the
  following specifications:

  template <typename T, class D>
  struct FNC {
    ImageStateBasePtr operator()(const ImageStateImpl<T,D>* isi);
  };

  This allows partial specialization to be used
*/
template <template <typename T, class D> class FNC>
struct unary_dispatch_op {

#define ISI_UNARY_DISPATCH_OP_CALL(TYPE,fnc,isb) \
  else if(const TYPE * is = dynamic_cast<const TYPE *>(isb)) return fnc(is)

  ImageStateBasePtr operator()(const ImageStateBase* isb) {
    if(false); 
    ISI_UNARY_DISPATCH_OP_CALL(RealSpatialImageState,call,isb);
    ISI_UNARY_DISPATCH_OP_CALL(ComplexSpatialImageState,call,isb);
    ISI_UNARY_DISPATCH_OP_CALL(RealFrequencyImageState,call,isb);
    ISI_UNARY_DISPATCH_OP_CALL(ComplexFrequencyImageState,call,isb);
    ISI_UNARY_DISPATCH_OP_CALL(ComplexHalfFrequencyImageState,call,isb);
    ISI_UNARY_DISPATCH_OP_CALL(WordSpatialImageState,call,isb);
    else throw DispatchException("unknown image state in unary_dispatch_op");
  }

#undef ISI_UNARY_DISPATCH_OP_CALL

  template <typename T, class D>
  ImageStateBasePtr call(const ImageStateImpl<T,D>* isi) {
    static FNC<T,D> fnc;
    return fnc(isi);
  }
};

////////////////////////////////////////////
// binary RTTI dispatch

//! in-place binary dispatch
/*!
  the template parameter is expected to be a class with
  the following specifications:

  template<typename T1, class D1, typename T2, class D2>
  struct FNC {
    void operator()(ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs);
  };

  Note that the first parameter, lhs, is non-const; it is expected to be modified
  by the function object.

  Note also that partial specialization is possible.
*/
template <template <typename T1, class D1, typename T2, class D2> class FNC>
struct binary_dispatch_ip {

#define ISI_BINARY_DISPATCH_IP_CALL_LEFT(TYPE,fnc,lhs,rhs) \
  else if(TYPE * is = dynamic_cast<TYPE *>(lhs)) fnc(is,rhs)

  void operator()(ImageStateBase* lhs, const ImageStateBase* rhs) {
    if(false); 
    ISI_BINARY_DISPATCH_IP_CALL_LEFT(RealSpatialImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_LEFT(ComplexSpatialImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_LEFT(RealFrequencyImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_LEFT(ComplexFrequencyImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_LEFT(ComplexHalfFrequencyImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_LEFT(WordSpatialImageState,left,lhs,rhs);
    else {
      throw DispatchException("unknown lhs image state in binary_dispatch_ip");
    }
  }

#undef ISI_BINARY_DISPATCH_IP_CALL_LEFT

#define ISI_BINARY_DISPATCH_IP_CALL_RIGHT(TYPE,fnc,lhs,rhs) \
  else if(const TYPE * is = dynamic_cast<const TYPE *>(rhs)) fnc(lhs,is)

  template <typename T1, class D1>
  void left(ImageStateImpl<T1,D1>* lhs, const ImageStateBase* rhs) {
    if(false); 
    ISI_BINARY_DISPATCH_IP_CALL_RIGHT(RealSpatialImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_RIGHT(ComplexSpatialImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_RIGHT(RealFrequencyImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_RIGHT(ComplexFrequencyImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_RIGHT(ComplexHalfFrequencyImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_IP_CALL_RIGHT(WordSpatialImageState,right,lhs,rhs);
    else {
      throw DispatchException("unknown rhs image state in binary_dispatch_ip");
    }
  }
  
#undef ISI_BINARY_DISPATCH_IP_CALL_RIGHT

  template <typename T1, class D1, typename T2, class D2>
  void right(ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs) {
    static FNC<T1,D1,T2,D2> fnc;
    fnc(lhs,rhs);
  }
};



//! out-of-place binary dispatch
/*!
  the template parameter is expected to be a class with
  the following specifications:

  template<typename T1, class D1, typename T2, class D2>
  struct FNC {
    ImageStateBasePtr operator()(const ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs);
  };

  The resulting image state is returned.

  Note that partial specialization is possible.
*/
template <template <typename T1, class D1, typename T2, class D2> class FNC>
struct binary_dispatch_op {

#define ISI_BINARY_DISPATCH_OP_CALL_LEFT(TYPE,fnc,lhs,rhs) \
  else if(const TYPE * is = dynamic_cast<const TYPE *>(lhs)) return fnc(is,rhs)

  ImageStateBasePtr operator()(const ImageStateBase* lhs, const ImageStateBase* rhs) {
    if(false); 
    ISI_BINARY_DISPATCH_OP_CALL_LEFT(RealSpatialImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_LEFT(ComplexSpatialImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_LEFT(RealFrequencyImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_LEFT(ComplexFrequencyImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_LEFT(ComplexHalfFrequencyImageState,left,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_LEFT(WordSpatialImageState,left,lhs,rhs);
    else throw DispatchException("unknown lhs image state in binary_dispatch_op");
  }

#undef ISI_BINARY_DISPATCH_OP_CALL_LEFT

#define ISI_BINARY_DISPATCH_OP_CALL_RIGHT(TYPE,fnc,lhs,rhs) \
  else if(const TYPE * is = dynamic_cast<const TYPE *>(rhs)) return fnc(lhs,is)

  template <typename T1, class D1>
  ImageStateBasePtr left(const ImageStateImpl<T1,D1>* lhs, const ImageStateBase* rhs) {
    if(false); 
    ISI_BINARY_DISPATCH_OP_CALL_RIGHT(RealSpatialImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_RIGHT(ComplexSpatialImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_RIGHT(RealFrequencyImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_RIGHT(ComplexFrequencyImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_RIGHT(ComplexHalfFrequencyImageState,right,lhs,rhs);
    ISI_BINARY_DISPATCH_OP_CALL_RIGHT(WordSpatialImageState,right,lhs,rhs);
    else throw DispatchException("unknown rhs image state in binary_dispatch_op");
  }

#undef ISI_BINARY_DISPATCH_OP_CALL_RIGHT

  template <typename T1, class D1, typename T2, class D2>
  ImageStateBasePtr right(const ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs) {
    static FNC<T1,D1,T2,D2> fnc;
    return fnc(lhs,rhs);
  }


};



/*
 application of dispatch code

 template<typename T1, class D1, typename T2, class D2>
 struct fnc_add_ip {
   void operator()(ImageStateImpl<T1,D1>* lhs, const ImageStateImpl<T2,D2>* rhs);
 };

 // specialization possible, such as
 template<typename T1, typename T2, class D>
 struct fnc_add_ip<T1,D,T2,D> {
   void operator()(ImageStateImpl<T1,D>* lhs, const ImageStateImpl<T2,D>* rhs)
   {
     // same domain!
   }
 };

 typedef binary_dispatch_ip<add_ip> dispatch_add_ip;

 void f()
 {
   ImageStateBase *is1,*is2;
   dispatch_add_ip(is1,is2);
 }

*/

}}}} // ns

#endif
