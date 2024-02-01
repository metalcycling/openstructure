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
  basic image handle, manages shared instances of an image

  offers the Image interface

  Author: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_IMAGE_HANDLE_H
#define IMG_IMAGE_HANDLE_H

#include <boost/shared_ptr.hpp>

#include <ost/img/module_config.hh>
#include <ost/message.hh>
#include <ost/base.hh>

#include "data.hh"
#include "observable.hh"
#include "extent_iterator.hh"

#include "image_state/image_state_visitor_fw.hh"
#include "image_state/image_state_base_fw.hh"

#include <ost/img/module_config.hh>

namespace test_image {
  void test_Observer();
}

namespace ost { namespace img {

using image_state::ImageStateBase;
using image_state::ImageStateBasePtr;
using image_state::ImageStateNonModVisitorBase;
using image_state::ImageStateModIPVisitorBase;
using image_state::ImageStateConstModIPVisitorBase;
using image_state::ImageStateModOPVisitorBase;
using image_state::ImageStateConstModOPVisitorBase;
using image_state::ImageStateMorphVisitorBase;

// fw declarations
class Point;
class Extent;
class Size;
class NonModAlgorithm;
class ModIPAlgorithm;
class ConstModIPAlgorithm;
class ModOPAlgorithm;
class ConstModOPAlgorithm;
class DataObserver;

//! Exception thrown upon invalid image handle usage
class DLLEXPORT InvalidImageHandle: public Error {
public:
  InvalidImageHandle(): Error("Attempt to use invalid ImageHandle") {}
};

// fw
class ConstImageHandle;

//! Manage shared instances of images
/*!
  An ImageHandle is used to access an underlying implementation
  of an image. ImageHandle uses reference counting, ie copy and assignement
  return a new handle which points to the same image. This makes
  it cheap to pass ImageHandle around as value. The resources allocated
  by the underlying image implementation will be freed as soon as the
  reference count drops to zero.

  ImageHandle derives from Data and thus implements its interface. It extends 
  it to include image specific functionality. ImageHandle does shared_ptr 
  management of an ImageImpl instance. The observer functionality is implemented
  in parallel using shared_ptr as well, so that the set of observers is 
  preserved among the shared ImageHandles.

  Creating a new image is done with the CreateImage() function:

  \code

  ImageHandle ih = CreateImage(Extent(Point(0,0),Size(10,10)));

  \endcode

  To import and export an image from or to a file, the IOManager
  must be used.

  All image manipulation routines are provided as algorithms, which 
  are applied with the Apply() and ApplyIP() methods.
*/
class DLLEXPORT_OST_IMG_BASE ImageHandle: public Data {
  friend ImageHandle DoCreateImage(const Extent& e, DataType type, DataDomain dom);

  // for internal observer checking
  friend void test_image::test_Observer();

  // for access to Sampling
  friend class ConstImageHandle;

  typedef Observable DataObservable;
  typedef boost::shared_ptr<DataObservable> ObsPtr;

  typedef boost::shared_ptr<ImageStateBasePtr> StatePtrPtr;

public:
  //! Creates an empty, ie invalid handle
  /*!
    The default ctor is provided for e.g. using ImageHandle in STL containers. Any
    attempt to use the interface of an empty ImageHandle will throw an 
    InvalidImageHandle exception.
  */
  ImageHandle();

  //! Initialize with an existing handle
  /*!
    The newly created handle will share the image with the existing one
  */
  ImageHandle(const ImageHandle& h);

  //! Assign another image based on handle
  /*!
    The reference count the the previously shared images is
    decreased and the one of the assignee is increased
  */
  ImageHandle& operator=(const ImageHandle& h);

  //@{
  //! return a handle that manages a new instance of the image
  /*!
    Similar to the Clone concept, the ImageHandle that is
    returned points to a duplicated image instance
  */
  ImageHandle Copy(bool cc=true) const;

  //! return a handle to a new image containing a partial copy
  /*!
    The returning domain will match the type (ie REAL or COMPLEX), 
    but the domain will always be spatial
  */
  ImageHandle Extract(const Extent& e) const;

  //! pastes given data to this image, using data extent
  /*!
    Effectively replaces a part of the image with the values
    retrieved from the given data, either another image handle or 
    a function. The region that is replaced in this image is given
    by the extent of the data.
  */
  void Paste(const Data& d);

  //! sets content of image to that of the given one
  /*!
    In contrast to the assignement operator, which
    affects only the handle mechanism, the Set method
    assigns the content of the passwd handle to this
    one. Other image handle pointing to the same image
    content will therefore also be affected.
  */
  void Set(const ImageHandle& h);

  //! swap guts with another handle
  void Swap(ImageHandle& h);

  //! Reset image state to given size, type and domain
  /*!
    Essentially the same as creating a new ImageHandle with CreateImage()
    and assigning it to an existing handle. After a call to reset,
    all values are zero.

    This routine exhibits only weak exception safety, as the existing state
    is first cleared to conserve memory usage.
  */
  void Reset(const Extent &e, DataType type=REAL,DataDomain dom=SPATIAL);

  //! Returns true if handle is valid
  bool IsValid() const;

  //! Returns size in bytes of allocated memory
  long MemSize() const;
  //@}


  /*! @name Properties interface
   */
  //@{
  //! Return data type: \c REAL or \c COMPLEX
  virtual DataType GetType() const;

  //! Return domain: \c SPATIAL, \c FREQUENCY or \c HALF_FREQUENCY
  virtual DataDomain GetDomain() const;

  //! retrieve Extent
  /*!
    The extent that is returned here is the logical extent.
  */
  virtual Extent GetExtent() const;

  //! return spatial origin of the image
  /*!
    The spatial origin is independent on the current domain. It affects the 
    FFT operations by introducing a origin-dependent phase shift for each
    value in the frequency domain.
  */
  virtual Point GetSpatialOrigin() const;

  //! set spatial origin of image
  /*
    For images in the frequency domain, a modification of the spatial origin
    usually has drastic consequences, since all values are affected by a phase shift.
  */
  void SetSpatialOrigin(const Point& p);

  //! set origin in center
  /*
    See Extent for difference in even and odd dimensions. Only
    has an effect for spatial images, since frequency data always
    has origin in center.

    May have effect on values if image is in frequency domain.
  */
  void CenterSpatialOrigin();
  //@}

  //! Get absolute coordinate for the spatial origin of the image
  /*
    This method returns the aboslute coordinate of the real space spatial
    origin even when the map is in frequency space
  */
  Vec3 GetAbsoluteOrigin() const;

  //! Sets the absolute coordinate for the spatial origin of the image
  /*
    This method sets the aboslutes coordinate of the real space spatial
    origin even when the map is in frequency space
  */
  void SetAbsoluteOrigin(const Vec3& c) ;

  //! Index to coordinate conversion
  /*
    Returns the physical dimensions of a given Point based on
    the sampling of the current domain and the absolute coordinates
    of the origin of the image
  */
  Vec3 IndexToCoord(const Point &p) const;

  //! Convert coordinates back to (fractional) indices
  /*
    Returns the location within the Data that corresponds to
    the given physical coordinate, bases on the sampling of
    the current domain and the absolute coordinates
    of the origin of the image
  */
  Vec3 CoordToIndex(const Vec3& c) const;

    //! Index to coordinate conversion
  /*
    Returns the physical dimensions of a Vec containing fractional Point
    coordinates. The returned value is based on
    the sampling of the current domain and the absolute coordinates
    of the origin of the image
  */
  Vec3 FractionalIndexToCoord(const Vec3 &p) const;

  /*! @name Value access
    All value access with these getters/setters is boundary-checked;
    as a consequence, reading outside of the defined region will 
    return 0, and writing outside is silently ignored.
  */
  //@{
  //! Retrieve value at given Point as a \c Real
  virtual Real GetReal(const Point& p) const;

  //! Set explicit real value
  /// \sa \ref spoke_pattern.py "Spoke Pattern Example"
  void SetReal(const Point &p, Real r);

  //! Retrieve value at given Point as a \c Complex
  virtual Complex GetComplex(const Point& p) const;

  //! Set explicit complex value
  void SetComplex(const Point &p, const Complex& c);
  //@}

  //! Get interpolated real value from 3D,2D or 1D vec
  virtual Real GetIntpolReal(const Vec3 &v) const;
  virtual Real GetIntpolReal(const Vec2 &v) const;
  virtual Real GetIntpolReal(const Real &d) const;

  //! Get interpolated complex value from 3D,2D or 1D vec
  virtual Complex GetIntpolComplex(const Vec3 &v) const;
  virtual Complex GetIntpolComplex(const Vec2 &v) const;
  virtual Complex GetIntpolComplex(const Real &d) const;



  /*! @name DataObserver interface
    Allows (de)registration of a DataObserver as well as
    broadcasting a notification to all attached observers.
  */
  //! Attach observer
  virtual void Attach(DataObserver *o) const;

  //! Detach observer
  virtual void Detach(DataObserver *o) const;

  //! Notify each observer that state has changed
  virtual void Notify() const;
  virtual void Notify(const Extent& e) const;
  virtual void Notify(const Point& p) const;


  /*! @name DataAlgorithm interface
    A collection of methods that form the main interaction 
    between the underlying image and algorithms derived from
    one of the abstract algorithm classes.
  */
  //@{

  //! Apply a const algorithm to this image
  /*!
    Calls the implemented NonModAlgorithm::Visit(const ImageHandle&) method
  */
  void Apply(NonModAlgorithm& a) const;

  //! Apply a const algorithm to this image
  /*!
    This does exactly the same as ImageHandle::Apply(NonModAlgorithm &)
    and is provided for syntax consistency.
   */
  void ApplyIP(NonModAlgorithm& a) const;

  //! Apply in-place Algorithm in-place
  void ApplyIP(ModIPAlgorithm& a);

  //! Apply in-place Algorithm, return new image with modifications
  /*!
    The modifications will be done on a copy of this image, which
    is then returned.
  */
  ImageHandle Apply(ModIPAlgorithm& a) const;

  //! Apply in-place const Algorithm in-place
  void ApplyIP(const ConstModIPAlgorithm& a);

  //! Apply in-place const Algorithm, return new image with modifications
  /*!
    The modifications will be done on a copy of this image, which
    is then returned.
  */
  ImageHandle Apply(const ConstModIPAlgorithm& a) const;

  //! Apply out-of-place Algorithm in-place
  /*!
    The new image returned by the out-of-place algorithm will
    replace the current one. This replacement will of
    course reflect in all ImageHandle s that point to the same
    Image.
  */
  void ApplyIP(ModOPAlgorithm& a);

  //! Apply out-of-place Algorithm
  /*!
    The image modified by the out-of-place algorithm is returned.
  */
  ImageHandle Apply(ModOPAlgorithm& a) const;

  //! Apply out-of-place const Algorithm in-place
  /*!
    The new image returned by the out-of-place algorithm will
    replace the current one. This replacement will of
    course reflect in all ImageHandle s that point to the same
    Image.
  */
  void ApplyIP(const ConstModOPAlgorithm& a);

  //! Apply out-of-place const Algorithm
  /*!
    The image modified by the out-of-place algorithm is returned.
  */
  ImageHandle Apply(const ConstModOPAlgorithm& a) const;
  //@}


  /*! @name StateVisitor interface
    Allows access to underlying image state via ImageStateVisitor
    pattern.
   */
  //@{
  //! Apply non-modifying state visitor
  void StateApply(ImageStateNonModVisitorBase& v) const;

  //! Applies in-place visitor in-place
  void StateApplyIP(ImageStateModIPVisitorBase& v);
  //! Applies in-place visitor out-of-place
  ImageHandle StateApply(ImageStateModIPVisitorBase& v) const;
  //! Applies in-place const visitor in-place
  void StateApplyIP(const ImageStateConstModIPVisitorBase& v);
  //! Applies in-place const visitor out-of-place
  ImageHandle StateApply(const ImageStateConstModIPVisitorBase& v) const;
  //! Applies out-of-place visitor in-place
  void StateApplyIP(ImageStateModOPVisitorBase& v);
  //! Applies out-of-place visitor out-of-place
  ImageHandle StateApply(ImageStateModOPVisitorBase& v) const;
  //! Applies out-of-place const visitor in-place
  void StateApplyIP(const ImageStateConstModOPVisitorBase& v);
  //! Applies out-of-place const visitor out-of-place
  ImageHandle StateApply(const ImageStateConstModOPVisitorBase& v) const;
  //! Applies morph visitor in-place
  void StateApplyIP(ImageStateMorphVisitorBase& v);
  //! Applies morph visitor out-of-place
  ImageHandle StateApply(ImageStateMorphVisitorBase& v) const;
  //! Get ExtentIterator
  ExtentIterator GetIterator() const ;

  //@}

  /*! @name Comparison
    Returns \c true if this handle and reference handle
    point to the same image, \c false otherwise
  */
  //@{
  bool operator==(const ImageHandle& ih) const;
  bool operator!=(const ImageHandle& ih) const;
  //@}

  /*! @name Member operators
   */
  //@{
  ImageHandle& operator+=(Real v);
  ImageHandle& operator+=(const Complex& v);
  ImageHandle& operator-=(Real v);
  ImageHandle& operator-=(const Complex& v);
  ImageHandle& operator*=(Real v);
  ImageHandle& operator*=(const Complex& v);
  ImageHandle& operator/=(Real v);
  ImageHandle& operator/=(const Complex& v);

  void Reciproce();

  ImageHandle& operator+=(const ConstImageHandle& h);
  ImageHandle& operator-=(const ConstImageHandle& h);
  ImageHandle& operator*=(const ConstImageHandle& h);
  ImageHandle& operator/=(const ConstImageHandle& h);

  //@}

  /*! @name low-level access
   */
  //@{
  ImageStateBasePtr& ImageStatePtr();
  const ImageStateBasePtr& ImageStatePtr() const;
  //@}
  
protected:
  //! Initialization based on a StatePtrPtr
  ImageHandle(const StatePtrPtr& spp);

  //! see Data::Sampling()
  virtual PixelSampling& Sampling();
  //! see Data::Sampling()
  virtual const PixelSampling& Sampling() const;

  //! used by CreateImage factory function
  static ImageHandle Create(const Extent& e, DataType type, DataDomain dom);

private:
  // this is all there is to it: two shared_ptr<> instances
  StatePtrPtr impl_;
  ObsPtr obs_;

  bool equal(const ImageHandle& ih) const;
};


//! read-only ImageHandle implementation
/*!
  ConstImageHandle offers only the read-only (const) interface of
  ImageHandle. Any sort of out-of-place manipulation will return
  the result as a non-const ImageHandle again!

  ConstImageHandle provides two different 'upcasting' characteristics.
  First, it is derived from ConstData, and may thus be used as such.
  Second, it provides an implicit conversion to <tt>const Data&</tt>.
*/
/*
  In contrast to ConstData and Data, which use is-a inheritance,
  ConstImageHandle implements its functionality using an instance
  of ImageHandle as a private member variable, and an implicit
  copy ctor and assignement op. This avoids Diamond-of-Death 
  virtual inheritance.

  On the downside, the complete ConstData interface needs to be
  retyped and forwarded to the ImageHandle instance. It is a reference
  to this instance that is returned for the 'upcasting' to const Data&
*/
class DLLEXPORT_OST_IMG_BASE ConstImageHandle: public ConstData {
  friend class ImageHandle;
public:
  //! pseudo 'upcasting' to <tt>const Data&</tt>
  operator const Data& () const {return handle_;}

  //! Empty, ie invalid handle
  ConstImageHandle();

  //! Initialize with an existing handle
  /*!
    The newly created handle will share the image with the existing one
  */
  ConstImageHandle(const ConstImageHandle& h);

  //! Initialize with an existing non-const handle
  /*!
    Provides implicit ImageHandle -> ConstImageHandle conversion. The newly 
    created const handle will still share the image with the non-const one.
  */
  ConstImageHandle(const ImageHandle& h);

  //! Assign another image based on handle
  /*!
    The reference count the the previously shared images is
    decreased and the one of the assignee is increased
  */
  ConstImageHandle& operator=(const ConstImageHandle& h);

  //! Assign another image based on a non-const handle
  /*!
    Implicit conversion of non-const to const image handle.
  */
  ConstImageHandle& operator=(const ImageHandle& h);

  //@{
  //! see ImageHandle::Copy(bool)
  ImageHandle Copy(bool cc=true) const;

  //! see ImageHandle::Copy(const Extent&)
  ImageHandle Extract(const Extent& e) const;

  //! Returns true if handle is valid
  bool IsValid() const;

  //! Returns size in bytes of allocated memory
  long MemSize() const;
  //@}

  /*! @name Properties interface
   */
  //@{
  //! see ImageHandle::GetType()
  DataType GetType() const;

  //! see ImageHandle::GetDomain()
  DataDomain GetDomain() const;

  //! see ImageHandle::GetExtent()
  Extent GetExtent() const;

  //! see ImageHandle::GetSpatialOrigin()
  Point GetSpatialOrigin() const;
  //@}

  //! see ImageHandle::GetAbsoluteOrigin()
  Vec3 GetAbsoluteOrigin() const;

  //! see ImageHandle::IndexToCoord()
  Vec3 IndexToCoord(const Point &p) const;

  //! see ImageHandle::CoordToIndex()
  Vec3 CoordToIndex(const Vec3& c) const;

  //! See ImageHandle::FractionalIndexToCoord()
  Vec3 FractionalIndexToCoord(const Vec3 &p) const;

  /*! @name Value access
  */
  //@{
  //! see ImageHandle::GetReal()
  Real GetReal(const Point& p) const;

  //! see ImageHandle::GetComplex()
  Complex GetComplex(const Point& p) const;
  //@}

  //! Get interpolated real value from 3D,2D or 1D vec
  virtual Real GetIntpolReal(const Vec3 &v) const;
  virtual Real GetIntpolReal(const Vec2 &v) const;
  virtual Real GetIntpolReal(const Real &d) const;

  //! Get interpolated complex value from 3D,2D or 1D vec
  virtual Complex GetIntpolComplex(const Vec3 &v) const;
  virtual Complex GetIntpolComplex(const Vec2 &v) const;
  virtual Complex GetIntpolComplex(const Real &d) const;


  /*! @name DataObserver interface
  */
  //! see ImageHandle::Attach()
  void Attach(DataObserver *o) const;

  //! see ImageHandle::Detach()
  void Detach(DataObserver *o) const;

  //! see ImageHandle::Notify()
  void Notify() const;
  void Notify(const Extent& e) const;
  void Notify(const Point& p) const;


  /*! @name DataAlgorithm interface
    Only the out-of-place functionality is available
  */
  //@{

  //! see ImageHandle::Apply(NonModAlgorithm&)
  void Apply(NonModAlgorithm& a) const;
  void ApplyIP(NonModAlgorithm& a) const;

  //! see ImageHandle::Apply(ModIPAlgorithm&)
  ImageHandle Apply(ModIPAlgorithm& a) const;

  //! see ImageHandle::Apply(const ConstModIPAlgorithm&)
  ImageHandle Apply(const ConstModIPAlgorithm& a) const;

  //! see ImageHandle::Apply(ModOPAlgorithm)
  ImageHandle Apply(ModOPAlgorithm& a) const;

  //! see ImageHandle::Apply(const ConstModOPAlgorithm)
  ImageHandle Apply(const ConstModOPAlgorithm& a) const;
  //@}

  /*! @name StateVisitor interface
    Only the out-of-place functionality is available
   */
  //@{
  //! see ImageHandle::StateApply(ImageStateNonModVisitorBase&)
  void StateApply(ImageStateNonModVisitorBase& v) const;

  //! see ImageHandle::StateApply(ImageStateModIPVisitorBase&)
  ImageHandle StateApply(ImageStateModIPVisitorBase& v) const;

  //! see ImageHandle::StateApply(ImageStateConstModIPVisitorBase&)
  ImageHandle StateApply(const ImageStateConstModIPVisitorBase& v) const;

  //! see ImageHandle::StateApply(ImageStateModOPVisitorBase&)
  ImageHandle StateApply(ImageStateModOPVisitorBase& v) const;

  //! see ImageHandle::StateApply(ImageStateConstModOPVisitorBase&)
  ImageHandle StateApply(const ImageStateConstModOPVisitorBase& v) const;

  //! see ImageHandle::StateApply(ImageStateMorphVisitorBase&)
  ImageHandle StateApply(ImageStateMorphVisitorBase& v) const;
  //! Get ExtentIterator
  ExtentIterator GetIterator() const ;
  //@}

  /*! @name Comparison
  */
  //@{
  bool operator==(const ConstImageHandle& ih) const;
  bool operator==(const ImageHandle& ih) const;
  bool operator!=(const ConstImageHandle& ih) const;
  bool operator!=(const ImageHandle& ih) const;
  //@}

  const ImageStateBasePtr& ImageStatePtr() const;

protected:

  //! see ImageHandle::Sampling()
  PixelSampling& Sampling();
  //! see ImageHandle::Sampling()
  const PixelSampling& Sampling() const;

private:
  // actual implementation
  ImageHandle handle_;

  bool equal(const ConstImageHandle& h) const;
  bool equal(const ImageHandle& h) const;
};


/*! @name Global operators
  These operators are implemented in terms of the member
  operators of ImageHandle, and due to the implicit
  convertability of ImageHandle to ConstImageHandle, they
  may be applied to both handle variants.
*/

// the other comparison ops are handled as member operators
DLLEXPORT_OST_IMG_BASE bool operator==(const ImageHandle& lhs, const ConstImageHandle& rhs);
DLLEXPORT_OST_IMG_BASE bool operator!=(const ImageHandle& lhs, const ConstImageHandle& rhs);

DLLEXPORT_OST_IMG_BASE ImageHandle operator+(const ConstImageHandle& h, Real v);
DLLEXPORT_OST_IMG_BASE ImageHandle operator+(const ConstImageHandle& h, const Complex& v);
DLLEXPORT_OST_IMG_BASE ImageHandle operator+(Real v, const ConstImageHandle& h);
DLLEXPORT_OST_IMG_BASE ImageHandle operator+(const Complex& v, const ConstImageHandle& h);
DLLEXPORT_OST_IMG_BASE ImageHandle operator-(const ConstImageHandle& h, Real v);
DLLEXPORT_OST_IMG_BASE ImageHandle operator-(const ConstImageHandle& h, const Complex& v);
DLLEXPORT_OST_IMG_BASE ImageHandle operator-(Real v, const ConstImageHandle& h);
DLLEXPORT_OST_IMG_BASE ImageHandle operator-(const Complex& v, const ConstImageHandle& h);
DLLEXPORT_OST_IMG_BASE ImageHandle operator*(const ConstImageHandle& h, Real v);
DLLEXPORT_OST_IMG_BASE ImageHandle operator*(const ConstImageHandle& h, const Complex& v);
DLLEXPORT_OST_IMG_BASE ImageHandle operator*(Real v, const ConstImageHandle& h);
DLLEXPORT_OST_IMG_BASE ImageHandle operator*(const Complex& v, const ConstImageHandle& h);
DLLEXPORT_OST_IMG_BASE ImageHandle operator/(const ConstImageHandle& h, Real v);
DLLEXPORT_OST_IMG_BASE ImageHandle operator/(const ConstImageHandle& h, const Complex& v);

DLLEXPORT_OST_IMG_BASE ImageHandle operator+(const ConstImageHandle& h1, const ConstImageHandle& h2);
DLLEXPORT_OST_IMG_BASE ImageHandle operator-(const ConstImageHandle& h1, const ConstImageHandle& h2);
DLLEXPORT_OST_IMG_BASE ImageHandle operator*(const ConstImageHandle& h1, const ConstImageHandle& h2);
DLLEXPORT_OST_IMG_BASE ImageHandle operator/(const ConstImageHandle& h1, const ConstImageHandle& h2);


}} // namespace img

#endif
