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
#ifndef OST_SEQUENCE_VIEWER_SEQUENCE_VIEW_OBJECT
#define OST_SEQUENCE_VIEWER_SEQUENCE_VIEW_OBJECT

/*
  Author: Stefan Scheuber
 */

#ifndef Q_MOC_RUN

#include <ost/mol/entity_handle.hh>

#include <ost/gfx/entity.hh>

#include <ost/seq/sequence_list.hh>

#include "align_properties_painter.hh"
#include "conservation_painter.hh"
#include "background_painter.hh"
#include "seq_secstr_painter.hh"
#include "seq_selection_painter.hh"
#include "seq_text_painter.hh"

#include "base_view_object.hh"
#endif

namespace ost { namespace gui {

class SequenceViewObject : public BaseViewObject
{
  Q_OBJECT

public:
  SequenceViewObject(seq::SequenceList& sequences, const QList<QString>& names, QObject* parent = 0);
  SequenceViewObject(seq::SequenceHandle& sequence, const QString& name, QObject* parent = 0);
  SequenceViewObject(mol::ChainView& chain, const QString& name, QObject* parent = 0);
  SequenceViewObject(const gfx::EntityP& entity, QObject* parent = 0);
  SequenceViewObject(QObject* parent = 0);

  void Init();

  void AddSequence(seq::ConstSequenceHandle sequence, const QString& name=QString());
  void AddChain(mol::ChainView& chain, const QString& name=QString());

  void AttachGfxObject(gfx::EntityP& ent);
  gfx::EntityP& GetGfxObject();

  virtual void SetDisplayMode(const QString& mode);

  QMap<int, QList<int> > GetIndexesForView(const mol::EntityView& view);
  QMap<int, QList<int> > GetIndexesForSubject(const QString& subject, const QString& sequence_name=QString());

protected:
  static const QString properties_mode;
  static const QString secondary_structure_mode;

  static AlignPropertiesPainter* align_properties_painter;
  static ConservationPainter* conservation_painter;
  static BackgroundPainter* background_painter;
  static SeqSecStrPainter* seq_secondary_structure_painter;
  static SeqSelectionPainter* seq_selection_painter;
  static SeqTextPainter* seq_text_painter;

private:
  gfx::EntityP entity_;


};


}}

#endif
