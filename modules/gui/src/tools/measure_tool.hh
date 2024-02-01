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
#ifndef OST_GUI_MEASURE_TOOL_HH
#define OST_GUI_MEASURE_TOOL_HH

/*
  Author: Tobias Schmidt
 */

#include <ost/gui/tools/tool.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/gfx/color.hh>

#include <QList>
#include <QWidget>
namespace ost { namespace gui {

class DLLEXPORT_OST_GUI Measurement {
public:
  typedef enum {
    DIST, ANGL, DIHE
  } Type;

  Measurement(mol::AtomHandleList ahl, gfx::Color col, Real line_width);

  Real GetMeasurement() const;
  Type GetType() const;
  String GetName() const;
  mol::AtomHandleList GetAtoms() const;
  gfx::Color GetColor() const;
  void SetColor(gfx::Color c);
  Real GetLineWidth() const;
  void SetLineWidth(Real w);
  bool IsValid() const;

private:
  mol::AtomHandleList ahl_;
  Type type_;
  String name_;
  Real measurement_;
  bool valid_;
  gfx::Color col_;
  Real line_width_;
};

typedef std::vector<Measurement> MeasurementList;


class DLLEXPORT_OST_GUI MeasureTool : public Tool {
  Q_OBJECT
public:
  MeasureTool();
  virtual ~MeasureTool() {};
  
  int GetMeasureMode();
  virtual void Click(const MouseEvent& event);
  virtual bool CanOperateOn(gfx::NodePtrList nodes);
  virtual QIcon GetIcon() const;

  gfx::Color GetColor() const;
  Real GetLineWidth() const;

  void RenderGL();

public slots:
  void ClearMeasurements();

private:
  int num_clicks_;
  int mode_;
  mol::AtomHandleList sel_atoms_;
  MeasurementList ml_;
};
 
}}

#endif
