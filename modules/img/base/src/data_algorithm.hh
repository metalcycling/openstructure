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
  Data algorithm abstract base classes

  Author: Ansgar Philippsen
*/

#ifndef IMG_DATA_ALG_H
#define IMG_DATA_ALG_H

#include <ost/img/module_config.hh>
#include "image_fw.hh"

/*
  TODO: 

  add info object as parameter in base class, this guarantees that every
  algorithm has access to an info object for logging
*/

namespace ost { namespace img {

//! Algorithm base class
/*!
  Provide common functionality for all algorithms. Only some sort of naming
  implemented.
*/
class DLLEXPORT_OST_IMG_BASE AlgorithmBase {
 public:
  //! Return name of algorithm
  const String& GetName() const;

  virtual ~AlgorithmBase();

 protected:
  //! Constructor must be initialized with algorithm name
  AlgorithmBase(const String& n);

  AlgorithmBase(const AlgorithmBase& a);
  AlgorithmBase& operator=(const AlgorithmBase& b);

 private:
  String name_;
};

//! Non-modifying Algorithm, abstract base class
class DLLEXPORT_OST_IMG_BASE NonModAlgorithm: public AlgorithmBase
{
public:
  //! visitor implementation for images
  virtual void Visit(const ConstImageHandle& i) = 0;

protected:
  //! ctor, must be initialized with name
  NonModAlgorithm(const String& name);
  //! copy ctor
  NonModAlgorithm(const NonModAlgorithm& a);
  //! assignement op
  NonModAlgorithm& operator=(const NonModAlgorithm& a);
};


//! In-place modification algorithm
/*!
  The Visit method is not const! As a consequence, algorithms
  of this type cannot be used as temporary objects. If the
  algorithm does not need to store values beyond the duraction
  of the Visit call, it is advisable to use ConstModIPAlgorithm
*/
class DLLEXPORT_OST_IMG_BASE ModIPAlgorithm: public AlgorithmBase
{
public:
  virtual void Visit(ImageHandle& ih) = 0;
protected:
  ModIPAlgorithm(const String& name);
  ModIPAlgorithm(const ModIPAlgorithm& a);
  ModIPAlgorithm& operator=(const ModIPAlgorithm& a);
};


//! In-place modification const algorithm
/*!
  The const refers to the visit method that is const; as
  a consequence, the algorithm may be used as a temporary
  object, but is not able to modify its members from
  within this method.
*/
class DLLEXPORT_OST_IMG_BASE ConstModIPAlgorithm: public AlgorithmBase
{
public:
  virtual void Visit(ImageHandle& ih) const = 0;
protected:
  ConstModIPAlgorithm(const String& name);
  ConstModIPAlgorithm(const ConstModIPAlgorithm& a);
  ConstModIPAlgorithm& operator=(const ConstModIPAlgorithm& a);
};


//! out-of-place modification algorithm
/*!
  The Visit method is not const! As a consequence, algorithms
  of this type cannot be used as temporary objects. If the
  algorithm does not need to store values beyond the duraction
  of the Visit call, it is advisable to use ConstModOPAlgorithm
*/
class DLLEXPORT_OST_IMG_BASE ModOPAlgorithm: public AlgorithmBase
{
public:
  virtual ImageHandle Visit(const ConstImageHandle& ih) = 0;

protected:
  ModOPAlgorithm(const String& name);
  ModOPAlgorithm(const ModOPAlgorithm& a);
  ModOPAlgorithm& operator=(const ModOPAlgorithm& a);
};

//! out-of-place modification const algorithm
/*!
  The const refers to the visit method that is const; as
  a consequence, the algorithm may be used as a temporary
  object, but is not able to modify its members from
  within this method.
*/
class DLLEXPORT_OST_IMG_BASE ConstModOPAlgorithm: public AlgorithmBase
{
public:
  virtual ImageHandle Visit(const ConstImageHandle& ih) const=  0;

protected:
  ConstModOPAlgorithm(const String& name);
  ConstModOPAlgorithm(const ConstModOPAlgorithm& a);
  ConstModOPAlgorithm& operator=(const ConstModOPAlgorithm& a);
};

}} // namespace

#endif
