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

/*
  Author: Stefan Scheuber
 */



#include <ost/mol/mol.hh>
#include <ost/mol/view_op.hh>

#include <ost/gfx/entity.hh>

#include "secstr_row.hh"

#include <QtGui>
namespace ost { namespace gui {

SecStrRow::SecStrRow(const QString& name, mol::ChainView& chain, SequenceViewObject* parent) : SequenceRow(name,parent)
{ this->SetChain(chain); }

void SecStrRow::SetSequence(seq::ConstSequenceHandle sequence)
{
  //Do nothing
}

void SecStrRow::SetChain(mol::ChainView& chain)
{
  String seq_str;
  seq_str.reserve(chain.GetResidueCount());
  for (mol::ResidueViewList::const_iterator r=chain.GetResidueList().begin(),
       e2=chain.GetResidueList().end(); r!=e2; ++r) {
    mol::ResidueView res=*r;
    seq_str.append(1, res.GetOneLetterCode());
  }
  if (!seq_str.empty()) {
    seq::SequenceHandle sequence=seq::CreateSequence(this->GetName().toStdString(), seq_str);
    mol::EntityView v_one_chain=chain.GetEntity().GetHandle().CreateEmptyView();
    v_one_chain.AddChain(chain, mol::ViewAddFlag::INCLUDE_ALL);
    sequence.AttachView(v_one_chain);

    mol::alg::SecStructureSegments sec = mol::alg::ExtractSecStructureSegments(chain);
    secstr_ = QVarLengthArray<mol::SecStructure>(chain.GetResidueCount());
    for (mol::alg::SecStructureSegments::iterator i=sec.begin(),
         e=sec.end(); i!=e; ++i) {
      mol::alg::SecStructureSegment s=*i;
      for(int i = s.first; i <= s.last ;i++){
        secstr_[i] = s.ss_type;
      }
    }
    this->chain_ = chain;
    SequenceRow::SetSequence(sequence);
  }
}

const mol::ChainView& SecStrRow::GetChain() const
{
  return this->chain_;
}

QVariant SecStrRow::GetData(int column, int role) const
{
  if(column > 0 && column <= this->GetSequence().GetLength()){
    if (role==Qt::ToolTipRole){
      return QVariant(QString(this->chain_.GetResidueByIndex(column-1).GetQualifiedName().c_str()));
    }
    if (role==Qt::UserRole){
      QVariant variant;
      variant.setValue(secstr_);
      return variant;
    }
    if (role==Qt::UserRole+1){
      return QVariant(this->GetFontSize());
    }
  }
  return SequenceRow::GetData(column, role);
}

void SecStrRow::DoubleClicked(int column)
{
  if(column>0){
    column-=1;
    const QVarLengthArray<mol::SecStructure>& sec = this->secstr_;
    if(! sec.isEmpty()>0 && column < sec.size()){
      const mol::SecStructure& src_str = sec[column];
      QVarLengthArray<bool> src_type(3);
      src_type[0] = src_str.IsHelical();
      src_type[1] = src_str.IsExtended();
      src_type[2] = src_str.IsCoil();
      int i = column;
      QSet<int> cols_to_add;
      mol::SecStructure col_str = sec[i];
      while(i >= 0 && (col_str = sec[i])){
        if(src_type[0] == col_str.IsHelical()
            && src_type[1] == col_str.IsExtended()
            && src_type[2] == col_str.IsCoil()){
        cols_to_add.insert(i+1);
        --i;
        }
        else{break;}
      }
      i = column + 1;
      if(i < sec.size()){
        while(i < sec.size() && (col_str = sec[i])){
          if(src_type[0] == col_str.IsHelical()
              && src_type[1] == col_str.IsExtended()
              && src_type[2] == col_str.IsCoil()){
          cols_to_add.insert(i+1);
          ++i;
          }
          else{
            break;
          }
        }
      }
      this->SetSelection(cols_to_add, QSet<int>());
    }
  }
  else{
    SequenceRow::DoubleClicked(column);
  }
}

}}
