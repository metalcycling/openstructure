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
  Author:Andreas Schenk
*/

#ifndef IMG_IMAGE_LIST_HH
#define IMG_IMAGE_LIST_HH


#include <ost/img/image.hh>

namespace ost { namespace img {



//! Manages a collection of images
/*!
  The image list, as its name suggests, contains a collection
  of images. The interface offers a collective load/save functionality
  as well as a collective Apply, which corresponds to a subset of the
  Apply methods in the img::ImageHandle interface.

*/
class DLLEXPORT_OST_IMG_BASE ImageList: public  std::vector<ImageHandle> {
public:
  ImageList();
  ImageList(const ImageList& il);
  ImageList(ImageList::iterator start,ImageList::iterator end);
  ImageHandle GetGallery(unsigned int columns=0,bool border=false);
  ImageHandle GetImageStack();


  /*! @name Algorithm apply
  */
  //@{

  //! Apply a const algorithm to each image in the list
  /*!
    see also img::ImageHandle::Apply(NonModAlgorithm&)
  */
  void Apply(NonModAlgorithm& a) const;

  //! Apply a const algorithm to each image in the list
  /*!
    see also img::ImageHandle::Apply(NonModAlgorithm&)
   */
  void ApplyIP(NonModAlgorithm& a) const;

  //! Apply in-place Algorithm in-place to each image in the list
  /*!
    see also img::ImageHandle::ApplyIP(ModIPAlgorithm&)
  */
  void ApplyIP(ModIPAlgorithm& a);

  //! Apply in-place Algorithm to each image in the list, return new list
  /*!
    The new list that is returned contains a copy of each image in
    the original list, on which the given algorithm has been applied
    out-of-place

    See also img::ImageHandle::Apply(ModIPAlgorithm&)
  */
  ImageList Apply(ModIPAlgorithm& a) const;

  //! Apply in-place const Algorithm in-place
  void ApplyIP(const ConstModIPAlgorithm& a);

  ImageList Apply(const ConstModIPAlgorithm& a) const;

  void ApplyIP(ModOPAlgorithm& a);

  ImageList Apply(ModOPAlgorithm& a) const;

  void ApplyIP(const ConstModOPAlgorithm& a);

  ImageList Apply(const ConstModOPAlgorithm& a) const;
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
  ImageList StateApply(ImageStateModIPVisitorBase& v) const;
  //! Applies in-place const visitor in-place
  void StateApplyIP(const ImageStateConstModIPVisitorBase& v);
  //! Applies in-place const visitor out-of-place
  ImageList StateApply(const ImageStateConstModIPVisitorBase& v) const;
  //! Applies out-of-place visitor in-place
  void StateApplyIP(ImageStateModOPVisitorBase& v);
  //! Applies out-of-place visitor out-of-place
  ImageList StateApply(ImageStateModOPVisitorBase& v) const;
  //! Applies out-of-place const visitor in-place
  void StateApplyIP(const ImageStateConstModOPVisitorBase& v);
  //! Applies out-of-place const visitor out-of-place
  ImageList StateApply(const ImageStateConstModOPVisitorBase& v) const;
  //! Applies morph visitor in-place
  void StateApplyIP(ImageStateMorphVisitorBase& v);
  //! Applies morph visitor out-of-place

  ImageList StateApply(ImageStateMorphVisitorBase& v) const;
  ImageList& operator+=(Real v);
  ImageList& operator+=(const Complex& v);
  ImageList& operator-=(Real v);
  ImageList& operator-=(const Complex& v);
  ImageList& operator*=(Real v);
  ImageList& operator*=(const Complex& v);
  ImageList& operator/=(Real v);
  ImageList& operator/=(const Complex& v);
  ImageList& operator+=(const ConstImageHandle& h);
  ImageList& operator-=(const ConstImageHandle& h);
  ImageList& operator*=(const ConstImageHandle& h);
  ImageList& operator/=(const ConstImageHandle& h);
private:
  Size get_max_size_();
  DataType get_data_type_();

};

}}

#endif
