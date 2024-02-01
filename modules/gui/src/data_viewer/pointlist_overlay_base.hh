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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_POINTLIST_OVERLAY_BASE_HH
#define IMG_POINTLIST_OVERLAY_BASE_HH



#include <ost/img/point.hh>

#include <ost/img/point.hh>
#include <ost/gui/module_config.hh>
#include <ost/gui/data_viewer/overlay_base.hh>

#include "strategies_fw.hh"

#include <QColor>
#include <QAction>
#include <QDialog>
#include <QLineEdit>
#include <QColorDialog>
#include <QPushButton>
#include <QVBoxLayout>


namespace ost { namespace img { namespace gui {

class DLLEXPORT_OST_GUI PointlistOverlayBaseSettings: public QDialog
{
  Q_OBJECT;
 public:
  PointlistOverlayBaseSettings(const QColor& ac, const QColor& pc, int ssiz, int sstr, QWidget* p);

 public slots:
  void OnActiveColor();
  void OnPassiveColor();
  void OnSymbolSize(int v);
  void OnSymbolStrength(int v);

 public:
  QColor active_color;
  QColor passive_color;
  int symbol_size;
  int symbol_strength;

 protected:

  QVBoxLayout* main_layout_;
  QPushButton* active_color_b_;
  QPushButton* passive_color_b_;
};

class DLLEXPORT_OST_GUI PointlistOverlayBase: public Overlay
{
  Q_OBJECT;
public:
  PointlistOverlayBase(const String& name);

  // overlay interface
  virtual void OnDraw(QPainter& pnt,  DataViewerPanel* dvp, bool is_active) {};
  virtual void OnMenuEvent(QAction* e);
  virtual QMenu* GetMenu();

  // own virtual interface
  virtual void DrawPointList(QPainter& pnt, DataViewerPanel* dvp, const QColor& col ,const std::vector<QPoint>& pointlist);
  virtual void DrawVariableSizePointList(QPainter& pnt, DataViewerPanel* dvp, const QColor& col ,const std::vector<std::pair<QPoint,double> >& pointlist);

  //
  unsigned int GetSymbolSize() const;
  void SetSymbolSize(unsigned int symbolsize);
  unsigned int GetSymbolShape() const;
  void SetSymbolShape(unsigned int symbolshape);
  bool GetCrosshair() const;
  void SetCrosshair(bool flag);
  void SetActiveColor(const QColor& col);
  void SetPassiveColor(const QColor& col);
private:

  template <class StrategyClass> 
  void ReplaceStrategy();

protected:
  SymbolDrawingStrategyPtr strategy_;
  QColor active_color_;
  QColor passive_color_;
  unsigned int symbolsize_;	
  unsigned int symbolstrength_;
  QMenu* menu_;

  QAction* a_sq_;
  QAction* a_ci_;
  QAction* a_cr_;
  QAction* a_settings_;

  void SetProps(PointlistOverlayBaseSettings* props);
};

}}}  //ns


#endif
