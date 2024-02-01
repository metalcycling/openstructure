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
  Author: Marco Biasini
 */

#include "gl_helper.hh"

#include <ost/gfx/gl_include.hh>
#include <ost/gfx/glext_include.hh>
#include <ost/gfx/primitives.hh>


namespace ost { namespace gfx {

Primitive::Primitive(const String& name):
  GfxNode(name), fill_(true), line_(true), outline_color_(1.0,1.0,1.0), 
  fill_color_(0.8,0.8,0.8), line_width_(1.0)
{
  
}

void Primitive::SetFill(const Color& color, bool fill)
{
  fill_color_=color;
  fill_=fill;
}
  
void Primitive::SetOutline(const Color& color, bool line, float width)
{
  outline_color_=color;
  line_=line;
  line_width_=width;
}

const Color& Primitive::GetFillColor() const
{
  return fill_color_;
}
  
const Color& Primitive::GetOutlineColor() const
{
  return outline_color_;
}
  
bool Primitive::HasFill() const
{
  return fill_;
}
bool Primitive::HasOutline() const
{
  return line_;
}

float Primitive::GetLineWidth() const
{
  return line_width_;
}

void Primitive::SetLineWidth(float width)
{
  line_width_=width;
}
  
void Primitive::SetFillColor(const Color& color)
{
  fill_color_=color;
}

void Primitive::SetOutlineColor(const Color& color)
{
  outline_color_=color;
}

Cuboid::Cuboid(const String& name, const geom::Cuboid& cuboid):
  Primitive(name), cuboid_(cuboid)
{ }

void Cuboid::RenderGL(RenderPass pass)
{
  if (!this->IsVisible())
    return;
  glLineWidth(this->GetLineWidth());
  glColor4fv(this->GetOutlineColor());
  geom::Vec3 c=cuboid_.GetCenter();
  geom::Vec3 e=cuboid_.GetHalfExtents();
  geom::Vec3 aa=cuboid_.GetVecA();
  geom::Vec3 ab=cuboid_.GetVecB();
  geom::Vec3 ac=cuboid_.GetVecC();    
  glDisable(GL_LIGHTING);
  if (this->HasOutline()) {
    geom::Vec3 f, t;
    glBegin(GL_LINES);
      f=c-aa*e[0]+ab*e[1]+ac*e[2]; t=c+aa*e[0]+ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);
      f=c-aa*e[0]-ab*e[1]+ac*e[2]; t=c+aa*e[0]-ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);
      f=c-aa*e[0]-ab*e[1]-ac*e[2]; t=c+aa*e[0]-ab*e[1]-ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);
      f=c-aa*e[0]+ab*e[1]-ac*e[2]; t=c+aa*e[0]+ab*e[1]-ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);                  
                                                       
      f=c+aa*e[0]-ab*e[1]+ac*e[2]; t=c+aa*e[0]+ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);
      f=c-aa*e[0]-ab*e[1]+ac*e[2]; t=c-aa*e[0]+ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);      
      f=c+aa*e[0]-ab*e[1]-ac*e[2]; t=c+aa*e[0]+ab*e[1]-ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);      
      f=c-aa*e[0]-ab*e[1]-ac*e[2]; t=c-aa*e[0]+ab*e[1]-ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);                                                                  
      
      f=c+aa*e[0]+ab*e[1]-ac*e[2]; t=c+aa*e[0]+ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);      
      f=c-aa*e[0]+ab*e[1]-ac*e[2]; t=c-aa*e[0]+ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);      
      f=c+aa*e[0]-ab*e[1]-ac*e[2]; t=c+aa*e[0]-ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);      
      f=c-aa*e[0]-ab*e[1]-ac*e[2]; t=c-aa*e[0]-ab*e[1]+ac*e[2];
      glVertex3v(&f[0]); glVertex3v(&t[0]);      
    glEnd();    
  }
  if (this->HasFill()) {
    if (this->GetFillColor().Alpha()<0.9999) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    geom::Vec3 xa, xb, xc, xd;
    glColor4fv(this->GetFillColor());
    glBegin(GL_QUADS);
      xa=c-aa*e[0]+ab*e[1]+ac*e[2]; 
      xb=c+aa*e[0]+ab*e[1]+ac*e[2];
      xc=c+aa*e[0]-ab*e[1]+ac*e[2]; 
      xd=c-aa*e[0]-ab*e[1]+ac*e[2];
      glVertex3v(&xa[0]);
      glVertex3v(&xb[0]);
      glVertex3v(&xc[0]);
      glVertex3v(&xd[0]);
      
      xa=c-aa*e[0]+ab*e[1]-ac*e[2]; 
      xb=c+aa*e[0]+ab*e[1]-ac*e[2];
      xc=c+aa*e[0]-ab*e[1]-ac*e[2]; 
      xd=c-aa*e[0]-ab*e[1]-ac*e[2];
      glVertex3v(&xa[0]);
      glVertex3v(&xb[0]);
      glVertex3v(&xc[0]);
      glVertex3v(&xd[0]);
      
      xa=c-aa*e[0]+ab*e[1]+ac*e[2]; 
      xb=c+aa*e[0]+ab*e[1]+ac*e[2];
      xc=c+aa*e[0]+ab*e[1]-ac*e[2]; 
      xd=c-aa*e[0]+ab*e[1]-ac*e[2];
      glVertex3v(&xa[0]);
      glVertex3v(&xb[0]);
      glVertex3v(&xc[0]);
      glVertex3v(&xd[0]);
      
      xa=c-aa*e[0]-ab*e[1]+ac*e[2]; 
      xb=c+aa*e[0]-ab*e[1]+ac*e[2];
      xc=c+aa*e[0]-ab*e[1]-ac*e[2]; 
      xd=c-aa*e[0]-ab*e[1]-ac*e[2];
      glVertex3v(&xa[0]);
      glVertex3v(&xb[0]);
      glVertex3v(&xc[0]);
      glVertex3v(&xd[0]);
      
      xa=c+aa*e[0]-ab*e[1]+ac*e[2]; 
      xb=c+aa*e[0]+ab*e[1]+ac*e[2];
      xc=c+aa*e[0]+ab*e[1]-ac*e[2]; 
      xd=c+aa*e[0]-ab*e[1]-ac*e[2];
      glVertex3v(&xa[0]);
      glVertex3v(&xb[0]);
      glVertex3v(&xc[0]);
      glVertex3v(&xd[0]);
      
      xa=c-aa*e[0]-ab*e[1]+ac*e[2]; 
      xb=c-aa*e[0]+ab*e[1]+ac*e[2];
      xc=c-aa*e[0]+ab*e[1]-ac*e[2]; 
      xd=c-aa*e[0]-ab*e[1]-ac*e[2];
      glVertex3v(&xa[0]);
      glVertex3v(&xb[0]);
      glVertex3v(&xc[0]);
      glVertex3v(&xd[0]);
      
      xa=c+aa*e[0]-ab*e[1]+ac*e[2]; 
      xb=c+aa*e[0]+ab*e[1]+ac*e[2];
      xc=c+aa*e[0]+ab*e[1]-ac*e[2]; 
      xd=c+aa*e[0]-ab*e[1]-ac*e[2];
      glVertex3v(&xa[0]);
      glVertex3v(&xb[0]);
      glVertex3v(&xc[0]);
      glVertex3v(&xd[0]);
    glEnd();
  }
}


Quad::Quad(const String& name, const geom::Vec3& a, const geom::Vec3& b,
           const geom::Vec3& c, const geom::Vec3& d):
  Primitive(name)
{
  corner_points_[0]=a;
  corner_points_[1]=b;
  corner_points_[2]=c;
  corner_points_[3]=d;      
}

void Quad::RenderGL(RenderPass pass)
{
  if (!this->IsVisible()) {
    return;
  }
  if (!(pass==DEPTH_RENDER_PASS || pass==STANDARD_RENDER_PASS)) {
    return;
  }
    glDisable(GL_LIGHTING);
    if (this->HasFill()) {
      if (this->GetFillColor().Alpha()<0.9999) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }
      glColor4fv(this->GetFillColor());
      glBegin(GL_QUADS);
        for (int i =0; i<4; ++i) {
          glVertex3v(&corner_points_[i][0]);        
        }
      glEnd();
    }  
    if (this->HasOutline()) {
      glEnable(GL_POLYGON_OFFSET_LINE);
      glPolygonOffset(-10, -10.0);
      glLineWidth(this->GetLineWidth());
      glBegin(GL_LINE_LOOP);
        for (int i =0; i<4; ++i) {
          glVertex3v(&corner_points_[i][0]);        
        }
      glEnd();
      glDisable(GL_POLYGON_OFFSET_LINE);    
    } 
}
}}
