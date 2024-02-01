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
#ifndef OST_ENTITY_BASE_HT
#define OST_ENTITY_BASE_HT

#include <ost/mol/module_config.hh>
#include <ost/mol/impl/entity_impl_fw.hh>

#include <ost/generic_property.hh>

namespace ost { namespace mol {

/// \brief definition of EntityBase
class DLLEXPORT_OST_MOL EntityBase: 
  public GenericPropContainer<EntityBase> {
public:
  friend class ConstGenericPropContainer<EntityBase>;
  EntityBase(const impl::EntityImplPtr& impl);
  EntityBase();
  
  ///\brief Get entity name.
  ///
  /// \sa      EntityBase::SetName
  const String& GetName() const;

  /// \brief Set entity name
  ///
  /// \param ent_name is the new entity name.
  ///
  /// \sa      EntityBase::GetName
  void SetName(const String& ent_name);

public:
  
  /// \name Handle validity
  //@{
  /// \brief check validity of handle
  /// 
  /// check, whether the entity handle points to a valid entity.
  /// \note It is an error to use any method other than #IsValid, #Impl and 
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const { return this->IsValid(); }
  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const { return Impl().get()!=0; }
  impl::EntityImplPtr& Impl();
  const impl::EntityImplPtr& Impl() const;
protected:
  
  GenericPropContainerImpl* GpImpl();
  
  const GenericPropContainerImpl* GpImpl() const;
  
  void CheckValidity() const;  
private:
  impl::EntityImplPtr   impl_;
};

DLLEXPORT_OST_MOL std::ostream& operator<<(std::ostream& os, 
                                            const EntityBase& ent);
}} // ns

#endif // OST_ENTITY_BASE_HH

