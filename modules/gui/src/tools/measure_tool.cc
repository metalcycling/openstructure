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
  Author: Tobias Schmidt
 */
#include "measure_tool.hh"

#include <ost/gfx/gl_helper.hh>
#include <ost/gfx/glext_include.hh>

#include <iostream>
#include <sstream>


#include <ost/platform.hh>
#include <ost/gfx/scene.hh>

#include <ost/geom/vec3.hh>

#include <QDir>
namespace ost { namespace gui {

MeasureTool::MeasureTool()
  : Tool("Measure Tool"), num_clicks_(0)
{
  ToolOptions* opts=this->GetToolOptions();
  ToolOptionEnum* meas_mode(new ToolOptionEnum("meas_mode", "Measure Mode"));
  meas_mode->Add("Distance", 2);
  meas_mode->Add("Angle", 3);  
  meas_mode->Add("Dihedral", 4); 
  meas_mode->SetIndex(0);
  opts->AddOption(meas_mode);
  ToolOptionFloat* col_red(new ToolOptionFloat("col_red", "Red", 1.0, 0.0, 1.0));
  opts->AddOption(col_red);
  ToolOptionFloat* col_green(new ToolOptionFloat("col_green", "Green", 0.5, 0.0, 0.9));
  opts->AddOption(col_green);
  ToolOptionFloat* col_blue(new ToolOptionFloat("col_blue", "Blue", 0.5, 0.0, 1.0));
  opts->AddOption(col_blue);
  ToolOptionFloat* line_width(new ToolOptionFloat("line_width", "Line Width", 4, 0.1, 10.0));
  opts->AddOption(line_width);
  ToolOptionButton* clear_button(new ToolOptionButton("clear_button", "Clear", this, SLOT(ClearMeasurements())));
  opts->AddOption(clear_button);
  mode_=meas_mode->GetValue();
}

int MeasureTool::GetMeasureMode()
{
  ToolOption* opt=this->GetToolOptions()->GetOption("", "meas_mode");
  ToolOptionEnum* meas_mode=dynamic_cast<ToolOptionEnum*>(opt);
  if (meas_mode->GetValue()!=mode_) {
    mode_=meas_mode->GetValue();
    sel_atoms_.clear();
    num_clicks_=0;
  }
  return mode_;
}

gfx::Color MeasureTool::GetColor() const
{
  Real r=dynamic_cast<ToolOptionFloat*>(this->GetToolOptions()->GetOption("", "col_red"))->GetValue();
  Real g=dynamic_cast<ToolOptionFloat*>(this->GetToolOptions()->GetOption("", "col_green"))->GetValue();
  Real b=dynamic_cast<ToolOptionFloat*>(this->GetToolOptions()->GetOption("", "col_blue"))->GetValue();
  return gfx::Color(r,g,b);
}

Real MeasureTool::GetLineWidth() const
{
  return dynamic_cast<ToolOptionFloat*>(this->GetToolOptions()->GetOption("", "line_width"))->GetValue();
}

void MeasureTool::Click(const MouseEvent& event)
{
  int meas_mode=this->GetMeasureMode();
  gfx::Scene& scene=gfx::Scene::Instance();
  scene.SetSelectionMode(1);
  mol::AtomHandle atom_=scene.PickAtom(event.GetPos().x(),scene.GetViewport().height-event.GetPos().y()).second;
  if (atom_.IsValid()) {
    sel_atoms_.push_back(atom_);
    num_clicks_++;
  } else {
    sel_atoms_.clear();
    num_clicks_=0;
  }
  if (num_clicks_>=meas_mode) {
    ml_.push_back(Measurement(sel_atoms_, GetColor(), GetLineWidth()));
    std::stringstream s;
    s << ml_.back().GetType() << " (" << ml_.back().GetName() << ") " << ml_.back().GetMeasurement();
    scene.StatusMessage(s.str());
    num_clicks_=0;
    sel_atoms_.clear();
  }
  scene.RequestRedraw();
}

void MeasureTool::ClearMeasurements()
{
  ml_.clear();
}

void MeasureTool::RenderGL()
{
  glDisable(GL_LIGHTING);
  glPointSize(10.0);
  glBegin(GL_POINTS);
  glColor3f(0.5, 1.0, 0.5);
  mol::AtomHandleList::const_iterator ai;
  for (ai=sel_atoms_.begin();ai!=sel_atoms_.end();++ai) {
    glVertex3(ai->GetPos());
  }
  glEnd();

  glLineStipple(2, 0xEEEE);
  glEnable(GL_LINE_STIPPLE);

  MeasurementList::const_iterator mi;
  for (mi=ml_.begin();mi!=ml_.end();++mi) {
    glLineWidth(mi->GetLineWidth());
    glColor3f(mi->GetColor().Red(), mi->GetColor().Green(),
              mi->GetColor().Blue());
    glBegin(GL_LINE_STRIP);
    mol::AtomHandleList::const_iterator aj;
    mol::AtomHandleList atoms_=mi->GetAtoms();
    for (aj=atoms_.begin(); aj!=atoms_.end();++aj) {
      glVertex3(aj->GetPos());
    }
    glEnd();
  }
  glDisable(GL_LINE_STIPPLE);

  gfx::Scene& scene=gfx::Scene::Instance();
  glDisable(GL_DEPTH_TEST);
  for (mi=ml_.begin();mi!=ml_.end();++mi) {
    std::stringstream ss;
    ss.precision(4);
    ss << mi->GetMeasurement();
    geom::Vec3 pos;
    switch (mi->GetType()) {
      case Measurement::DIST:
        pos=(mi->GetAtoms()[0].GetPos()+mi->GetAtoms()[1].GetPos())/2;
        break;
      case Measurement::ANGL:
        pos=mi->GetAtoms()[1].GetPos();
        break;
      case Measurement::DIHE:
        pos=(mi->GetAtoms()[1].GetPos()+mi->GetAtoms()[2].GetPos())/2;
        break;
    }
    gfx::TextPrim tp(ss.str(), pos, mi->GetColor(), 0.7);
    scene.RenderText(tp);
  }
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glColor3f(1.0f, 1.0f, 1.0f);
  glDisable(GL_TEXTURE_2D);
}


bool MeasureTool::CanOperateOn(gfx::NodePtrList nodes)
{
  return true;
}

QIcon MeasureTool::GetIcon() const
{
  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");
  QIcon icon = QIcon(icon_path.absolutePath()+QDir::separator()+QString("distance_icon.png"));
  return icon;
}

Measurement::Measurement(mol::AtomHandleList ahl, gfx::Color col, Real line_width):
    ahl_(ahl), col_(col), line_width_(line_width)
{
  switch (ahl_.size()) {
    case 2:
    {
      type_=DIST;
      name_=ahl_[0].GetQualifiedName()+" - "+ahl_[1].GetQualifiedName();
      measurement_=geom::Distance(ahl_[0].GetPos(), ahl_[1].GetPos());
      valid_=true;
      break;
    }
    case 3:
    {
      type_=ANGL;
      name_=ahl_[0].GetQualifiedName()+" - "+ahl_[1].GetQualifiedName()+" - "
             +ahl_[2].GetQualifiedName();
      measurement_=57.2958*geom::Angle(ahl_[0].GetPos()-ahl_[1].GetPos(),
                                       ahl_[2].GetPos()-ahl_[1].GetPos());
      valid_=true;
      break;
    }
    case 4:
    {
      type_=DIHE;
      name_=ahl_[0].GetQualifiedName()+" - "+ahl_[1].GetQualifiedName()+" - "
             +ahl_[2].GetQualifiedName()+" - "+ahl_[3].GetQualifiedName();
      geom::Plane pl1 = geom::Plane(ahl_[0].GetPos(),
                                    ahl_[1].GetPos(),
                                    ahl_[2].GetPos());
      geom::Plane pl2 = geom::Plane(ahl_[1].GetPos(),
                                    ahl_[2].GetPos(),
                                    ahl_[3].GetPos());
      measurement_=57.2958*geom::Angle(pl1, pl2);
      geom::Vec3 v = geom::Vec3(ahl_[2].GetPos()-ahl_[3].GetPos());
      if (geom::Dot(v, pl1.GetNormal())>0) {
        measurement_*=-1;
      }
      valid_=true;
      break;
    }
    default:
      std::cout << "WARNING: wrong number of AtomHandles in AtomHandleList" << std::endl;
      valid_=false;
  }
}

bool Measurement::IsValid() const
{
  return valid_;
}

Real Measurement::GetMeasurement() const
{
  return measurement_;
}

String Measurement::GetName() const
{
  return name_;
}

gfx::Color Measurement::GetColor() const
{
  return col_;
}

void Measurement::SetColor(gfx::Color c)
{
  col_=c;
}

Real Measurement::GetLineWidth() const
{
  return line_width_;
}

void Measurement::SetLineWidth(Real w)
{
  line_width_=w;
}

Measurement::Type Measurement::GetType() const
{
  return type_;
}

mol::AtomHandleList Measurement::GetAtoms() const
{
  return ahl_;
}
}}
