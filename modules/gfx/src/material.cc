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
#include "gl_include.hh"
#include "material.hh"

namespace ost { namespace gfx {

Material::Material():
  amb_(0,0,0,1),
  diff_(1,1,1,1),
  spec_(0.3,0.3,0.3,1.0),
  emm_(0,0,0,1),
  shin_(32.0)
{}

Material::Material(const Color& amb, const Color& diff, const Color& spec, float shin, const Color& emm):
  amb_(amb),
  diff_(diff),
  spec_(spec),
  emm_(emm),
  shin_(shin)
{}

Material::Material(float amb, float diff, float spec, float shin, float emm):
  amb_(amb,amb,amb,1.0),
  diff_(diff,diff,diff,1.0),
  spec_(spec,spec,spec,1.0),
  emm_(emm,emm,emm,1.0),
  shin_(shin)
{}

void Material::SetAmb(const Color& c)
{
  amb_=c;
}

void Material::SetAmb(float c)
{
  amb_=Color(c,c,c,1.0);
}

void Material::SetDiff(const Color& c)
{
  diff_=c;
}

void Material::SetDiff(float c)
{
  diff_=Color(c,c,c,1.0);
}

void Material::SetSpec(const Color& c)
{
  spec_=c;
}

void Material::SetSpec(float c)
{
  spec_=Color(c,c,c,1.0);
}

void Material::SetShin(float s)
{
  shin_=s;
}

void Material::SetEmm(const Color& c)
{
  emm_=c;
}

void Material::SetEmm(float c)
{
  emm_=Color(c,c,c,1.0);
}

void Material::RenderGL()
{
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, amb_);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diff_);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_);
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emm_);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shin_);
}


}} // ns
