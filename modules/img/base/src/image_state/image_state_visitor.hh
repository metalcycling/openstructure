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

#ifndef IMG_IMAGE_STATE_VISITOR_HH
#define IMG_IMAGE_STATE_VISITOR_HH

#include <ost/img/data_algorithm.hh>
#include "image_state_impl.hh"
#include "image_state_def.hh"

namespace ost { namespace img { namespace image_state {

//! ImageStateVisitor general exception
class DLLEXPORT InvalidImageStateVisitor: public Error {
public:
  InvalidImageStateVisitor(const String& m): Error(m) {}
};

/*
  one-time definition of the constructor adapters, allowing
  zero to 10 ctor parameters to be automagically used. There
  is probably a recursive way to do this more elegantly...
*/

#define IMAGE_STATE_VISITOR_CTOR_ADAPTERS(CLASS) \
  /* no param */                                 \
  CLASS (): FNC() {}                             \
  /* 1 param */                                  \
  template <class P0>                            \
  CLASS (const P0& p0):                          \
    FNC(p0) {}                                   \
  /* 2 params */                                 \
  template <class P0,                            \
            class P1>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1):                          \
    FNC(p0,p1) {}                                \
  /* 3 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2):                          \
    FNC(p0,p1,p2) {}                             \
  /* 4 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3):                          \
    FNC(p0,p1,p2,p3) {}                          \
  /* 5 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4):                          \
    FNC(p0,p1,p2,p3,p4) {}                       \
  /* 6 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4,                            \
            class P5>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4,                           \
         const P5& p5):                          \
    FNC(p0,p1,p2,p3,p4,p5) {}                    \
  /* 7 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4,                            \
            class P5,                            \
            class P6>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4,                           \
         const P5& p5,                           \
         const P6& p6):                          \
    FNC(p0,p1,p2,p3,p4,p5,p6) {}                 \
  /* 8 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4,                            \
            class P5,                            \
            class P6,                            \
            class P7>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4,                           \
         const P5& p5,                           \
         const P6& p6,                           \
         const P7& p7):                          \
    FNC(p0,p1,p2,p3,p4,p5,p6,p7) {}              \
  /* 9 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4,                            \
            class P5,                            \
            class P6,                            \
            class P7,                            \
            class P8>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4,                           \
         const P5& p5,                           \
         const P6& p6,                           \
         const P7& p7,                           \
         const P8& p8):                          \
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8) {}           \
  /* 10 params */                                \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4,                            \
            class P5,                            \
            class P6,                            \
            class P7,                            \
            class P8,                            \
            class P9>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4,                           \
         const P5& p5,                           \
         const P6& p6,                           \
         const P7& p7,                           \
         const P8& p8,                           \
         const P9& p9):                          \
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9) {}        \
  /* 11 params */                                \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4,                            \
            class P5,                            \
            class P6,                            \
            class P7,                            \
            class P8,                            \
            class P9,                            \
            class PA>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4,                           \
         const P5& p5,                           \
         const P6& p6,                           \
         const P7& p7,                           \
         const P8& p8,                           \
         const P9& p9,                           \
         const PA& pa):                          \
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa) {}     \
  /* 12 params */                                \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3,                            \
            class P4,                            \
            class P5,                            \
            class P6,                            \
            class P7,                            \
            class P8,                            \
            class P9,                            \
            class PA,                            \
            class PB>                            \
  CLASS (const P0& p0,                           \
	 const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3,                           \
         const P4& p4,                           \
         const P5& p5,                           \
         const P6& p6,                           \
         const P7& p7,                           \
         const P8& p8,                           \
         const P9& p9,                           \
         const PA& pa,                           \
         const PB& pb):                          \
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb) {}


//! non-modifying image state visitor base class
class DLLEXPORT_OST_IMG_BASE ImageStateNonModVisitorBase {
public:
  virtual ~ImageStateNonModVisitorBase() {}
  virtual void Visit_(const RealSpatialImageState& is) = 0;
  virtual void Visit_(const ComplexSpatialImageState& is) = 0;
  virtual void Visit_(const RealFrequencyImageState& is) = 0;
  virtual void Visit_(const ComplexFrequencyImageState& is) = 0;
  virtual void Visit_(const ComplexHalfFrequencyImageState& is) = 0;
  virtual void Visit_(const WordSpatialImageState& is) = 0;
};


//! non-modifying image state visitor
/*!
  Expects a class with a public VisitState() method, having the following signature:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    void VisitState(const img::image_state::ImageStateImpl<T,D>& isi);
  };

  \endcode

  Full specialization is supported, but not partial one, since the
  member function is templated, not the overall class.

  To use the visitor itself:

  \code

  // definition
  typedef img::image_state::ConstVisitor<MyVisitorFnc> MyVisitor;
  // creating an instance
  MyVisitor my_visitor;

  // applying it to some image
  some_image_handle.StateApply(my_visitor);

  \endcode

  MyVisitorFnc will be a base class of MyVisitor, so any additional
  interface components remain fully accessible! Also, due to some
  template and compiler cleverness, the constructor paramers for 
  MyVisitorFnc are automagically available in MyVisitor. Using the
  example above:

  \code

  class MyVisitorFnc {
  public:
    MyVisitorFnc(); // default ctor
    MyVisitorFnc(int a, int b); //initialize with two ints
    template<typename T, class D>
    void VisitState(const img::image_state::ImageStateImpl<T,D>& isi);
  };

  typedef img::image_state::ConstVisitor<MyVisitorFnc> MyVisitor;

  void f() {
    MyVisitor my_visitor1; // using default ctor
    MyVisitor my_visitor2(17,29); // using other ctor

  \endcode

  See also common algorithms, e.g. img::alg::Stat
*/
template <class FNC>
class DLLEXPORT ImageStateNonModVisitor: public FNC, public ImageStateNonModVisitorBase {
public:
  
  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateNonModVisitor)

  // base interface
  virtual void Visit_(const RealSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(const ComplexSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(const RealFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(const ComplexFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(const ComplexHalfFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(const WordSpatialImageState& is) {T_Visit_(is);}

private:
  template <typename T, class D>
  void T_Visit_(const ImageStateImpl<T,D>& isi) {
    this->VisitState(isi);
  }
};


//! in-place modifying image state visitor base class
class DLLEXPORT_OST_IMG_BASE ImageStateModIPVisitorBase {
public:
  ImageStateModIPVisitorBase(const String& name="") {}
  virtual ~ImageStateModIPVisitorBase() {}
  virtual void Visit_(RealSpatialImageState& is) = 0;
  virtual void Visit_(ComplexSpatialImageState& is) = 0;
  virtual void Visit_(RealFrequencyImageState& is) = 0;
  virtual void Visit_(ComplexFrequencyImageState& is) = 0;
  virtual void Visit_(ComplexHalfFrequencyImageState& is) = 0;
  virtual void Visit_(WordSpatialImageState& is) = 0;
};


//! in-place modifying image state visitor
/*!
  Expects a class with a public Visit() method, having the following signature:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    void VisitState(img::image_state::ImageStateImpl<T,D>& isi);
  };

  \endcode

  The image state passed as a parameter will be modified.

  See img::image_state::ImageStateNonModVisitor for other important remarks

*/
template <class FNC>
class DLLEXPORT ImageStateModIPVisitor: public FNC, public ImageStateModIPVisitorBase {
public:
  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateModIPVisitor)

  virtual void Visit_(RealSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(ComplexSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(RealFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(ComplexFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(ComplexHalfFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(WordSpatialImageState& is) {T_Visit_(is);}

private:
  template <typename T, class D>
  void T_Visit_(ImageStateImpl<T,D>& isi) {
    this->VisitState(isi);
  }
};


//! in-place modifying image state const visitor base class
class DLLEXPORT_OST_IMG_BASE ImageStateConstModIPVisitorBase {
public:
  ImageStateConstModIPVisitorBase(const String& name="") {}
  virtual ~ImageStateConstModIPVisitorBase() {}
  virtual void Visit_(RealSpatialImageState& is) const = 0;
  virtual void Visit_(ComplexSpatialImageState& is) const = 0;
  virtual void Visit_(RealFrequencyImageState& is) const = 0;
  virtual void Visit_(ComplexFrequencyImageState& is) const = 0;
  virtual void Visit_(ComplexHalfFrequencyImageState& is) const = 0;
  virtual void Visit_(WordSpatialImageState& is) const = 0;
};


//! in-place modifying image state visitor
/*!
  Expects a class with a public Visit() method, having the following signature:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    void VisitState(img::image_state::ImageStateImpl<T,D>& isi) const;
  };

  \endcode

  The image state passed as a parameter will be modified. Due to the const'nes
  of the VisitState method, this routine may not set member variables, but
  the visitor may be used as a temporary object.

  See img::image_state::ImageStateNonModVisitor for other important remarks

*/
template <class FNC>
class DLLEXPORT_OST_IMG_BASE ImageStateConstModIPVisitor: public FNC, public ImageStateConstModIPVisitorBase {
public:
  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateConstModIPVisitor)

  virtual void Visit_(RealSpatialImageState& is) const {T_Visit_(is);}
  virtual void Visit_(ComplexSpatialImageState& is) const {T_Visit_(is);}
  virtual void Visit_(RealFrequencyImageState& is) const {T_Visit_(is);}
  virtual void Visit_(ComplexFrequencyImageState& is) const {T_Visit_(is);}
  virtual void Visit_(ComplexHalfFrequencyImageState& is) const {T_Visit_(is);}
  virtual void Visit_(WordSpatialImageState& is) const {T_Visit_(is);}

private:
  template <typename T, class D>
  void T_Visit_(ImageStateImpl<T,D>& isi) const {
    this->VisitState(isi);
  }
};


//! out-of-place modifying image state visitor base class
class DLLEXPORT_OST_IMG_BASE ImageStateModOPVisitorBase {
public:
  ImageStateModOPVisitorBase(const String& name="") {}
  virtual ~ImageStateModOPVisitorBase() {}
  virtual ImageStateBasePtr Visit_(const RealSpatialImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(const ComplexSpatialImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(const RealFrequencyImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(const ComplexFrequencyImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(const ComplexHalfFrequencyImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(const WordSpatialImageState& is) = 0;
};

//! out-of-place modifying image state visitor
/*!
  Expects a class with a public VisitState() method, having the following signature:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    ImageStateBasePtr Visit(const img::image_state::ImageStateImpl<T,D>& isi);
  };

  \endcode

  The modified image state will be returned as a new object

  See img::image_state::ImageStateNonModVisitor for other important remarks

*/
template <class FNC>
class DLLEXPORT ImageStateModOPVisitor: public FNC, public ImageStateModOPVisitorBase {
public:

  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateModOPVisitor)

  virtual ImageStateBasePtr Visit_(const RealSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const RealFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexHalfFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const WordSpatialImageState& is) {return T_Visit_(is);}

private:
  template <typename T, class D>
  ImageStateBasePtr T_Visit_(const ImageStateImpl<T,D>& isi) {
    return this->VisitState(isi);
  }
};


//! out-of-place modifying image state const visitor base class
class DLLEXPORT_OST_IMG_BASE ImageStateConstModOPVisitorBase {
public:
  ImageStateConstModOPVisitorBase(const String& name="") {}
  virtual ~ImageStateConstModOPVisitorBase() {}
  virtual ImageStateBasePtr Visit_(const RealSpatialImageState& is) const = 0;
  virtual ImageStateBasePtr Visit_(const ComplexSpatialImageState& is) const = 0;
  virtual ImageStateBasePtr Visit_(const RealFrequencyImageState& is) const = 0;
  virtual ImageStateBasePtr Visit_(const ComplexFrequencyImageState& is) const = 0;
  virtual ImageStateBasePtr Visit_(const ComplexHalfFrequencyImageState& is) const = 0;
  virtual ImageStateBasePtr Visit_(const WordSpatialImageState& is) const = 0;
};

//! out-of-place modifying image state visitor
/*!
  Expects a class with a public VisitState() method, having the following signature:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    ImageStateBasePtr Visit(const img::image_state::ImageStateImpl<T,D>& isi) const;
  };

  \endcode

  The modified image state will be returned as a new object. Due to the const'nes of
  the VisitState method, no member variables may be modified, but the visitor can be
  used as an temporary object.

  See img::image_state::ImageStateNonModVisitor for other important remarks

*/
template <class FNC>
class DLLEXPORT_OST_IMG_BASE ImageStateConstModOPVisitor: public FNC, public ImageStateConstModOPVisitorBase {
public:

  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateConstModOPVisitor)

  virtual ImageStateBasePtr Visit_(const RealSpatialImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexSpatialImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const RealFrequencyImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexFrequencyImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexHalfFrequencyImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const WordSpatialImageState& is) const {return T_Visit_(is);}

private:
  template <typename T, class D>
  ImageStateBasePtr T_Visit_(const ImageStateImpl<T,D>& isi) const {
    return this->VisitState(isi);
  }
};

//! morphing image state visitor base class
class DLLEXPORT_OST_IMG_BASE ImageStateMorphVisitorBase {
public:
  ImageStateMorphVisitorBase(const String& name="") {}
  virtual ~ImageStateMorphVisitorBase() {}
  virtual ImageStateBasePtr Visit_(RealSpatialImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(ComplexSpatialImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(RealFrequencyImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(ComplexFrequencyImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(ComplexHalfFrequencyImageState& is) = 0;
  virtual ImageStateBasePtr Visit_(WordSpatialImageState& is) = 0;
};

//! morphing image state visitor
/*!
  This visitor will use the memory area located by the passed in
  image state, and re-use it in a new one, hence 'morphing' the
  state from one to another. The calling image handle will delete
  the input image state after the VisitState() call returns.

  A class with a public VisitState method of the following signature
  is expected as the template parameter:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    ImageStateBasePtr Visit(img::image_state::ImageStateImpl<T,D>& isi);
  };

  \endcode
*/
template <class FNC>
class DLLEXPORT_OST_IMG_BASE ImageStateMorphVisitor: public FNC, public ImageStateMorphVisitorBase {
public:

  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateMorphVisitor)

  virtual ImageStateBasePtr Visit_(RealSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(ComplexSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(RealFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(ComplexFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(ComplexHalfFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(WordSpatialImageState& is) {return T_Visit_(is);}

private:
  template <typename T, class D>
  ImageStateBasePtr T_Visit_(ImageStateImpl<T,D>& isi) {
    return this->VisitState(isi);
  }
};

#undef IMAGE_STATE_VISITOR_CTOR_ADAPTERS

}}} // ns

#endif
