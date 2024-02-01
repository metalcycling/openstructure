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

#ifndef IMG_IMAGE_STATE_ALGORITHM_HH
#define IMG_IMAGE_STATE_ALGORITHM_HH

#include <ost/img/image.hh>
#include <ost/img/data_algorithm.hh>

#include "image_state_visitor.hh"

namespace ost { namespace img { namespace image_state {

/*
  one-time definition of the constructor adapters, allowing
  zero to 12 ctor parameters to be automagically used. There
  is probably a recursive way to do this more elegantly...

  this version includes a call to a base class to allow a name
  to be set
*/

#define IMAGE_STATE_VISITOR_CTOR_ADAPTERS(CLASS, BASE ) \
  /* no param */                                 \
  CLASS (): FNC(), BASE (FNC::GetAlgorithmName()) {}                             \
  /* 1 param */                                  \
  template <class P0>                            \
  CLASS (const P0& p0):                          \
    FNC(p0), BASE (FNC::GetAlgorithmName()) {}                                   \
  /* 2 params */                                 \
  template <class P0,                            \
            class P1>                            \
  CLASS (const P0& p0,                           \
         const P1& p1):                          \
    FNC(p0,p1), BASE (FNC::GetAlgorithmName()) {}                                \
  /* 3 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2>                            \
  CLASS (const P0& p0,                           \
         const P1& p1,                           \
         const P2& p2):                          \
    FNC(p0,p1,p2), BASE (FNC::GetAlgorithmName()) {}                             \
  /* 4 params */                                 \
  template <class P0,                            \
            class P1,                            \
            class P2,                            \
            class P3>                            \
  CLASS (const P0& p0,                           \
         const P1& p1,                           \
         const P2& p2,                           \
         const P3& p3):                          \
    FNC(p0,p1,p2,p3), BASE (FNC::GetAlgorithmName()) {}                          \
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
    FNC(p0,p1,p2,p3,p4), BASE (FNC::GetAlgorithmName()) {}                       \
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
    FNC(p0,p1,p2,p3,p4,p5), BASE (FNC::GetAlgorithmName()) {}                    \
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
    FNC(p0,p1,p2,p3,p4,p5,p6), BASE (FNC::GetAlgorithmName()) {}                 \
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
    FNC(p0,p1,p2,p3,p4,p5,p6,p7), BASE (FNC::GetAlgorithmName()) {}              \
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
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8), BASE (FNC::GetAlgorithmName()) {}           \
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
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9), BASE (FNC::GetAlgorithmName()) {} \
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
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa), BASE (FNC::GetAlgorithmName()) {} \
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
    FNC(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb), BASE (FNC::GetAlgorithmName()) {}


//! non-modifying image state visitor plus non-mod algorithm
/*!
  This convenience hybrid class combines the functionality of
  the NonModAlgorithm and the ImageStateNonModVisitor by implementing
  both interfaces. The template parameter is expected to be a class
  offering VisitState (see img::image_state::ImageStateNonModVisitor ),
  plus a static GetAlgorithmName() method:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    void VisitState(const img::image_state::ImageStateImpl<T,D>& isi);

    static String GetAlgorithmName();
  };

  \endcode

  The second interface method Visit(const ConstImageHandle& h) is handled
  internally by re-routing it automatically to the ImageStateVisitor. 

  For comments about the constructor and remaining interface of FNC, see
  img::image_state::ImageStateNonModVisitor
*/
template <class FNC>
class TEMPLATE_EXPORT ImageStateNonModAlgorithm: public FNC,  public NonModAlgorithm, public ImageStateNonModVisitorBase {
public:
  
  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateNonModAlgorithm, NonModAlgorithm )

  // base interface
  virtual void Visit_(const RealSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(const ComplexSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(const RealFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(const ComplexFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(const ComplexHalfFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(const WordSpatialImageState& is) {T_Visit_(is);}

  // const algorithm interface
  virtual void Visit(const ConstImageHandle& h) {
    h.StateApply(*this);
  }
  
private:
  /*
    must be defined at this point to ensure proper instantiation of
    ImageStateNonModAlgorithm<X>
  */
  template <typename T, class D>
  void T_Visit_(const ImageStateImpl<T,D>& isi) {
    this->VisitState(isi);
  }

};


//! in-place modifying image state visitor plus ip algorithm
/*!
  offers a combined image state visitor as well as image handle algorithm. The template
  parameter is expected to be a class offering VisitState (see 
  img::image_state::ImageStateModIPVisitor ), plus a static GetAlgorithmName() method:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    void VisitState(img::image_state::ImageStateImpl<T,D>& isi);

    static String GetAlgorithmName();
  };

  \endcode

*/
template <class FNC>
class TEMPLATE_EXPORT ImageStateModIPAlgorithm: public FNC, public ModIPAlgorithm, public ImageStateModIPVisitorBase {
public:
  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateModIPAlgorithm, ModIPAlgorithm)

  virtual void Visit_(RealSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(ComplexSpatialImageState& is) {T_Visit_(is);}
  virtual void Visit_(RealFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(ComplexFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(ComplexHalfFrequencyImageState& is) {T_Visit_(is);}
  virtual void Visit_(WordSpatialImageState& is) {T_Visit_(is);}

  // mod ip algorithm interface
  virtual void Visit(ImageHandle& h) {
    h.StateApplyIP(*this);
  }

private:
  template <typename T, class D>
  void T_Visit_(ImageStateImpl<T,D>& isi) {
    this->VisitState(isi);
  }
};

//! in-place modifying image state const visitor plus ip algorithm
/*!
  offers a combined image state const visitor as well as image handle algorithm. The template
  parameter is expected to be a class offering VisitState (see 
  img::image_state::ImageStateConstModIPVisitor ), plus a static GetAlgorithmName() method:

  \code

  class MyVisitorFnc {
  public:
    template<typename T, class D>
    void VisitState(img::image_state::ImageStateImpl<T,D>& isi) const;

    static String GetAlgorithmName();
  };

  \endcode

*/
template <class FNC>
class TEMPLATE_EXPORT ImageStateConstModIPAlgorithm: public FNC, 
             public ConstModIPAlgorithm,
             public ImageStateConstModIPVisitorBase {
public:
  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateConstModIPAlgorithm, ConstModIPAlgorithm)

  virtual void Visit_(RealSpatialImageState& is) const {T_Visit_(is);}
  virtual void Visit_(ComplexSpatialImageState& is) const {T_Visit_(is);}
  virtual void Visit_(RealFrequencyImageState& is) const {T_Visit_(is);}
  virtual void Visit_(ComplexFrequencyImageState& is) const {T_Visit_(is);}
  virtual void Visit_(ComplexHalfFrequencyImageState& is) const {T_Visit_(is);}
  virtual void Visit_(WordSpatialImageState& is) const {T_Visit_(is);}

  // mod ip const algorithm interface
  virtual void Visit(ImageHandle& h) const {
    h.StateApplyIP(*this);
  }

private:
  template <typename T, class D>
  void T_Visit_(ImageStateImpl<T,D>& isi) const {
    this->VisitState(isi);
  }
};

//! out-of-place modifying image state visitor plus op algorithm
template <class FNC>
class TEMPLATE_EXPORT ImageStateModOPAlgorithm: public FNC, public ModOPAlgorithm, public ImageStateModOPVisitorBase {
public:

  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateModOPAlgorithm, ModOPAlgorithm)

  virtual ImageStateBasePtr Visit_(const RealSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const RealFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexHalfFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const WordSpatialImageState& is) {return T_Visit_(is);}

  // mod op algorithm interface
  virtual ImageHandle Visit(const ConstImageHandle& h) {
    return h.StateApply(*this);
  }

private:
  template <typename T, class D>
  ImageStateBasePtr T_Visit_(const ImageStateImpl<T,D>& isi) {
    return this->VisitState(isi);
  }
};

//! out-of-place modifying image state const visitor plus op algorithm
template <class FNC>
class TEMPLATE_EXPORT ImageStateConstModOPAlgorithm: public FNC, 
                                                     public ConstModOPAlgorithm,
                                                     public ImageStateConstModOPVisitorBase {
public:

  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateConstModOPAlgorithm, ConstModOPAlgorithm)

  virtual ImageStateBasePtr Visit_(const RealSpatialImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexSpatialImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const RealFrequencyImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexFrequencyImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const ComplexHalfFrequencyImageState& is) const {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(const WordSpatialImageState& is) const {return T_Visit_(is);}

  // mod op algorithm interface
  virtual ImageHandle Visit(const ConstImageHandle& h) const {
    return h.StateApply(*this);
  }

private:
  template <typename T, class D>
  ImageStateBasePtr T_Visit_(const ImageStateImpl<T,D>& isi) const {
    return this->VisitState(isi);
  }
};


//! morph visitor plus ip algorithm
template <class FNC>
class TEMPLATE_EXPORT ImageStateMorphAlgorithm: public FNC, public ModIPAlgorithm, public ImageStateMorphVisitorBase {
public:

  IMAGE_STATE_VISITOR_CTOR_ADAPTERS(ImageStateMorphAlgorithm, ModIPAlgorithm)

  virtual ImageStateBasePtr Visit_(RealSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(ComplexSpatialImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(RealFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(ComplexFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(ComplexHalfFrequencyImageState& is) {return T_Visit_(is);}
  virtual ImageStateBasePtr Visit_(WordSpatialImageState& is) {return T_Visit_(is);}

  // mod ip algorithm interface
  virtual void Visit(ImageHandle& h) {
    h.StateApply(*this);
  }

private:
  template <typename T, class D>
  ImageStateBasePtr T_Visit_(ImageStateImpl<T,D>& isi) {
    return this->VisitState(isi);
  }
};

#undef IMAGE_STATE_VISITOR_CTOR_ADAPTERS

}}} // ns

#endif
