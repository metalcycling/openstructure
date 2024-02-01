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
  Author: Ansgar Philippsen
*/

#include <ost/log.hh>
#include <ost/gfx/offscreen_buffer.hh>
#include <ost/gfx/scene.hh>

namespace ost { namespace gfx {

static XVisualInfo* get_visinfo(Display* dpy, bool use_accum)
{
  int buf[64];
  unsigned int bufc=0;
  int d[]={16,12,8,6,4,1};
  int c[]={16,12,8,6,5,4,2,1};
  int a[]={8,4};
  XVisualInfo *vis=NULL;

  buf[bufc++]=GLX_RGBA;
  buf[bufc++]=GLX_DEPTH_SIZE;
  unsigned int depthi=bufc++;
  buf[bufc++]=GLX_RED_SIZE;
  unsigned int redi=bufc++;
  buf[bufc++]=GLX_BLUE_SIZE;
  unsigned int bluei=bufc++;
  buf[bufc++]=GLX_GREEN_SIZE;
  unsigned int greeni=bufc++;
  unsigned int accumi=0;
  if(use_accum) {
    buf[bufc++]=GLX_ACCUM_RED_SIZE;
    accumi=bufc++;
    buf[bufc++]=GLX_ACCUM_BLUE_SIZE;
    bufc++;
    buf[bufc++]=GLX_ACCUM_GREEN_SIZE;
    bufc++;
  }
  buf[bufc++]=None;

  for(unsigned int j=0;j<sizeof(c)/sizeof(int);j++) {
    for(unsigned int i=0;i<sizeof(d)/sizeof(int);i++) {
      buf[depthi]=d[i];
      buf[redi]=c[j];
      buf[bluei]=c[j];
      buf[greeni]=c[j];
      if(use_accum) {
        for(unsigned int k=0;k<sizeof(a)/sizeof(int);k++) {
          buf[accumi+0]=a[k];
          buf[accumi+2]=a[k];
          buf[accumi+4]=a[k];
          vis=glXChooseVisual(dpy,DefaultScreen(dpy),buf);
          if(vis) return vis;
        }
      } else {
        vis=glXChooseVisual(dpy,DefaultScreen(dpy),buf);
        if(vis) return vis;
      }
    }
  }
  return NULL;
}

static void create_error(const std::string& m)
{
  LOG_ERROR("error creating offscreen rendering context: " << m);
}

OffscreenBuffer::OffscreenBuffer(unsigned int width, unsigned int height, const OffscreenBufferFormat& f, bool shared):
  width_(width), height_(height), valid_(false), active_(false)
{
  LOG_DEBUG("offscreen buffer: checking for DISPLAY");
  if(getenv("DISPLAY")==NULL) {
    create_error("missing DISPLAY environment variable");
    return;
  }

  LOG_DEBUG("offscreen buffer: XOpenDisplay");
  dpy_ = XOpenDisplay(getenv("DISPLAY"));
  if(dpy_==NULL) {
    create_error("XOpenDisplay() failed");
    return;
  }

#if 0
  std::vector<int> attrib_list;
  attrib_list.push_back(GLX_RENDER_TYPE);
  attrib_list.push_back(GLX_RGBA_BIT);
  attrib_list.push_back(GLX_DRAWABLE_TYPE);
  attrib_list.push_back(GLX_PIXMAP_BIT);
  attrib_list.push_back(GLX_DOUBLEBUFFER);
  attrib_list.push_back(False);
  attrib_list.push_back(GLX_DEPTH_SIZE);
  attrib_list.push_back(f.dbits);
  attrib_list.push_back(GLX_RED_SIZE);
  attrib_list.push_back(f.cbits);
  attrib_list.push_back(GLX_GREEN_SIZE);
  attrib_list.push_back(f.cbits);
  attrib_list.push_back(GLX_BLUE_SIZE);
  attrib_list.push_back(f.cbits);
  attrib_list.push_back(GLX_ALPHA_SIZE);
  attrib_list.push_back(f.abits);
  attrib_list.push_back(0);

  int nelem=0;
  LOG_DEBUG("offscreen buffer: glXChooseFBConfig");
  fbconfig_ =glXChooseFBConfig(dpy_,0,&attrib_list[0],&nelem);
  if(fbconfig_==0 || nelem==0) {
    create_error("glXChooseFBConfig failed");
    return;
  }

#else

  visinfo_ = get_visinfo(dpy_,true);
  if(!visinfo_) {
    create_error("get_visinfo() failed");
    return;
  }

  context_=glXCreateContext(dpy_, visinfo_, 0, False);
  if(!context_) {
    create_error("glxCreateContext() failed");
    return;
  }

  pm_=XCreatePixmap(dpy_,RootWindow(dpy_,DefaultScreen(dpy_)),
                    width_,height_,visinfo_->depth);
  if(!pm_) {
    create_error("XCreatePixmap() failed");
    return;
  }

  glxpm_=glXCreateGLXPixmap(dpy_, visinfo_, pm_);
  if(!glxpm_) {
    create_error("glXCreateGLXPixmap() failed");
    return;
  }
  
#endif

  valid_=true;
}
    
OffscreenBuffer::~OffscreenBuffer()
{
  if(valid_) {
    glXDestroyGLXPixmap(dpy_,glxpm_);
    XFreePixmap(dpy_,pm_);
    glXDestroyContext(dpy_,context_);
    valid_=false;
  }
}

bool OffscreenBuffer::Resize(unsigned int width, unsigned int height)
{
  width_=width;
  height_=height;

  if(pm_) {
    XFreePixmap(dpy_,pm_);
  }
  pm_=XCreatePixmap(dpy_,RootWindow(dpy_,DefaultScreen(dpy_)),
                    width_,height_,visinfo_->depth);
  if(!pm_) {
    create_error("XCreatePixmap() failed");
    return false;
  }

  if(glxpm_) {
    glXDestroyGLXPixmap(dpy_,glxpm_);
  }
  glxpm_=glXCreateGLXPixmap(dpy_, visinfo_, pm_);
  if(!glxpm_) {
    create_error("glXCreateGLXPixmap() failed");
    return false;
  }

  if(active_) {
    MakeActive();
  }

  return true;
}

bool OffscreenBuffer::MakeActive()
{
  if(glXMakeCurrent(dpy_,glxpm_,context_)==False) {
    LOG_ERROR("error activating offscreen rendering: glXMakeCurrent() failed");
    return false;
  }
  active_=true;
  return true;
}

}} // ns
