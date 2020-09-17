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
#ifndef OST_CHAIN_BASE_HH
#define OST_CHAIN_BASE_HH

#include <ost/mol/module_config.hh>
#include <ost/mol/impl/chain_impl_fw.hh>

#include <ost/mol/chain_type.hh>
#include <ost/generic_property.hh>

namespace ost { namespace mol {

/// \brief base class for ChainHandle and ChainView
///
/// This class implements the common functionality of ChainHandle and ChainView.
///
/// Like \ref ResidueHandle "residues" and \ref AtomHandle "atoms", chains are
/// generic property containers and a such able to store user-defined properties.
///
/// Chains are linear chains of \ref ResidueHandle "residues". Peptide chains are
/// ordered from N- to C-terminus. 
class DLLEXPORT_OST_MOL ChainBase: 
  public GenericPropContainer<ChainBase> {

public: // constructors
  ChainBase();
  ChainBase(const impl::ChainImplPtr& impl);
public:
  
  /// \name Handle validity
  //@{
  /// \brief check validity of handle
  /// 
  /// check, whether the chain handle points to a valid handle.
  /// \note It is an error to use any method other than #IsValid, #Impl() and 
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const { return this->IsValid(); }
  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const { return Impl().get()!=0; }
  //@}
  friend class ConstGenericPropContainer<ChainBase>;
  String GetName() const;

  /// \brief Get the type of a chain.
  ///
  /// \return chain type of ChainType
  ChainType GetType() const;

  /// \brief Get information about a chain.
  ///
  /// \return description
  String GetDescription() const;

  const impl::ChainImplPtr& Impl() const {
    return impl_;
  }

  /// \brief whether the chain is a polymer
  ///
  /// True if one of IsPolysaccharide(), IsPolynucleotide(), IsPolypeptide() is 
  /// true or the chain is of type CHAINTYPE_POLYMER.
  bool IsPolymer() const;
  
  /// \brief whether the chain is a polysaccharide
  bool IsPolysaccharide() const;

  /// \brief whether the chain is an oligsaccharide (branched mmCIF entity)
  bool IsOligosaccharide() const;

  /// \brief whether the chain is a polypeptide
  bool IsPolypeptide() const;
  
  /// \brief whether the chain is a polynucleotide
  bool IsPolynucleotide() const;
  
  impl::ChainImplPtr& Impl() {
    return impl_;
  }
protected:
  GenericPropContainerImpl* GpImpl();
  
  const GenericPropContainerImpl* GpImpl() const;
    
  void CheckValidity() const;
private:
  impl::ChainImplPtr   impl_;
};

DLLEXPORT_OST_MOL std::ostream& operator<<(std::ostream& os, 
                                            const ChainBase& chain);
}} // ns

#endif // OST_CHAIN_BASE_HH

