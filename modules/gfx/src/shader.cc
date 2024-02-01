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

#include <vector>
#include <sstream>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <ost/platform.hh>
#include <ost/log.hh>

#include "impl/scene_fx.hh"
#include "glext_include.hh"

#include "shader.hh"
#include "scene.hh"

namespace ost { namespace gfx {

namespace bf = boost::filesystem;

Shader& Shader::Instance()
{
  static Shader instance;
  return instance;
}

Shader::Shader():
  valid_(false),
  current_program_(0),
  current_name_(""),
  shader_program_map_()
{
  if(!OST_GL_VERSION_2_0) {
    LOG_VERBOSE("OpenGL version less than 2.0, deactivating shader functionality");
  }
}


void Shader::PreGLInit() 
{
#if !defined(__APPLE__)
GLenum err = glewInit();
if (GLEW_OK != err) {
  LOG_ERROR("glew failure: " << glewGetErrorString(err));
  assert(false);
}
#endif
}

bool Shader::Compile(const std::string& shader_name, 
                     const std::string& shader_code,
                     GLenum shader_type,
                     GLuint& shader_id)
{
  if(!OST_GL_VERSION_2_0) return false;
  shader_id = glCreateShader(shader_type);
  const char* tmp_ptr[] = {0};
  tmp_ptr[0]=shader_code.c_str();
  glShaderSource(shader_id, 1, tmp_ptr, NULL);
  glCompileShader(shader_id);

  GLint sh_compiled;
  glGetShaderiv(shader_id, GL_COMPILE_STATUS, &sh_compiled);
  if(sh_compiled==GL_TRUE) {
    LOG_DEBUG("shader [" << shader_name << "] successfully compiled (" << shader_id << ")");
  } else {
    GLint sh_log_length;
    glGetShaderiv(shader_id, GL_INFO_LOG_LENGTH, &sh_log_length);
    std::vector<GLchar> sh_log(sh_log_length+1);
    glGetShaderInfoLog(shader_id, sh_log_length, NULL, &sh_log[0]);
    LOG_WARNING("shader [" << shader_name << "] compilation failed:" << std::endl << String(&sh_log[0]));
    return false;
  }
  return true;
}

void Shader::Link(const std::string& pr_name, const std::string& vs_code, const std::string& fs_code)
{
  if(!OST_GL_VERSION_2_0) return;
  shader_program_map_[pr_name]=0;
  GLuint shader_pr = glCreateProgram();

  GLuint shader_id;
  if(Compile(pr_name+" vs",vs_code,GL_VERTEX_SHADER,shader_id)) {
    LOG_DEBUG("attaching compiled vertex shader id " << shader_id << " to " << shader_pr);
    glAttachShader(shader_pr,shader_id);
  } else {
    LOG_WARNING("skipping [" << pr_name << "] due to error in vertex shader code");
    return;
  }
  if(Compile(pr_name+" fs",fs_code,GL_FRAGMENT_SHADER,shader_id)) {
    LOG_DEBUG("attaching compiled fragment shader id " << shader_id << " to " << shader_pr);
    glAttachShader(shader_pr,shader_id);
  } else {
    LOG_WARNING("skipping [" << pr_name << "] due to error in fragment shader code");
    return;
  }

  glLinkProgram(shader_pr);
  GLint pr_linked;
  glGetProgramiv(shader_pr,GL_LINK_STATUS,&pr_linked);
  if(pr_linked==GL_TRUE) {
    LOG_VERBOSE("shader [" << pr_name << "] sucessfully linked");
  } else {
    GLint pr_log_length;
    glGetProgramiv(shader_pr, GL_INFO_LOG_LENGTH, &pr_log_length);
    std::vector<GLchar> pr_log(pr_log_length+1);
    glGetProgramInfoLog(shader_pr, pr_log_length, NULL, &pr_log[0]);
    LOG_WARNING("shader [" << pr_name << "] linking failed:" << std::endl << String(&pr_log[0]));
    return;
  }
  shader_program_map_[pr_name]=shader_pr;
}

void Shader::Setup() 
{
  if(!OST_GL_VERSION_2_0) return;

  String ost_root = GetSharedDataPath();
  bf::path ost_root_dir(ost_root);
  bf::path shader_dir(ost_root_dir / "shader");

  struct ShaderListEntry {
    const char* file;
    GLenum type;
  } shader_list[] = {
    //////////////////////////////////////////////////////////////////
    // this is the master list of all shader code in shader/

    {"basic_vs.glsl", GL_VERTEX_SHADER},
    {"basic_fs.glsl", GL_FRAGMENT_SHADER},
    {"fraglight_vs.glsl", GL_VERTEX_SHADER},
    {"fraglight_fs.glsl", GL_FRAGMENT_SHADER},
    {"selfx_vs.glsl", GL_VERTEX_SHADER},
    {"selfx_fs.glsl", GL_FRAGMENT_SHADER},
    {"noop_vs.glsl", GL_VERTEX_SHADER},
    {"aaline_fs.glsl", GL_FRAGMENT_SHADER},
    {"iso_vs.glsl", GL_VERTEX_SHADER},
    {"iso_fs.glsl", GL_FRAGMENT_SHADER},
    {"fast_sphere_fs.glsl", GL_FRAGMENT_SHADER},
    {"outline_vs.glsl", GL_VERTEX_SHADER},
    {"dumpnorm_vs.glsl", GL_VERTEX_SHADER},
    {"dumpnorm_fs.glsl", GL_FRAGMENT_SHADER},
    {"quadpp_vs.glsl", GL_VERTEX_SHADER},
    {"convolute1_fs.glsl", GL_FRAGMENT_SHADER},
    {"amboccl_fs.glsl", GL_FRAGMENT_SHADER},
    {"scenefx_vs.glsl", GL_VERTEX_SHADER},
    {"scenefx_fs.glsl", GL_FRAGMENT_SHADER},
    {"beacon_fs.glsl", GL_FRAGMENT_SHADER},
    {"screenblur4_fs.glsl", GL_FRAGMENT_SHADER},
    {"test_tex_fs.glsl", GL_FRAGMENT_SHADER},
    {"material_phong.glsl", GL_FRAGMENT_SHADER},
    {"material_hemi.glsl", GL_FRAGMENT_SHADER},
    {"material_toon1.glsl", GL_FRAGMENT_SHADER},
    {"material_toon2.glsl", GL_FRAGMENT_SHADER},
    {"anaglyph_fs.glsl", GL_FRAGMENT_SHADER}
    //////////////////////////////////////////////////////////////////
  };

  std::map<std::string, std::string> shader_code_map;

  for(unsigned int i=0;i<sizeof(shader_list)/sizeof(ShaderListEntry);++i) {
    String shader_name(shader_list[i].file);
    bf::path shader_file(shader_dir / shader_name);

    if(!bf::exists(shader_file)){
      LOG_WARNING("shader file not found: [" << shader_file.string() << "]");
      shader_code_map[shader_name]=std::string();
      continue;
    }

    bf::ifstream shader_file_stream(shader_file);
    String shader_code;
    String line;
    while (std::getline(shader_file_stream,line)) {
      shader_code += line + "\n";
    }
    shader_code_map[shader_name]=shader_code;
  }

  std::vector<GLuint> shader_program_list;

  Link("basic",
       shader_code_map["basic_vs.glsl"],
       shader_code_map["basic_fs.glsl"]);

  // default frag lighting with phong lighting
  Link("fraglight",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_phong.glsl"]+
       shader_code_map["fraglight_fs.glsl"]);

  // hf lighting model
  Link("hemilight",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_hemi.glsl"]+
       shader_code_map["fraglight_fs.glsl"]);
  // alias
  shader_program_map_["hf"]=shader_program_map_["hemilight"];

  // toon 1 lighting model
  Link("toon1",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_toon1.glsl"]+
       shader_code_map["fraglight_fs.glsl"]);

  // toon 2 lighting model
  Link("toon2",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_toon2.glsl"]+
       shader_code_map["fraglight_fs.glsl"]);

  // raytraced sphere with different materials
  Link("fast_sphere_phong",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_phong.glsl"]+
       shader_code_map["fast_sphere_fs.glsl"]);
  Link("fast_sphere_hf",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_hemi.glsl"]+
       shader_code_map["fast_sphere_fs.glsl"]);
  Link("fast_sphere_toon1",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_toon1.glsl"]+
       shader_code_map["fast_sphere_fs.glsl"]);
  Link("fast_sphere_toon2",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["material_toon2.glsl"]+
       shader_code_map["fast_sphere_fs.glsl"]);

  // selection fx
  Link("selfx",
       shader_code_map["selfx_vs.glsl"],
       shader_code_map["selfx_fs.glsl"]);

  Link("aaline",
       shader_code_map["basic_vs.glsl"],
       shader_code_map["aaline_fs.glsl"]);
  Link("iso",
       shader_code_map["iso_vs.glsl"],
       shader_code_map["iso_fs.glsl"]);

  Link("outline",
       shader_code_map["outline_vs.glsl"],
       shader_code_map["basic_fs.glsl"]);
  Link("dumpnorm",
       shader_code_map["dumpnorm_vs.glsl"],
       shader_code_map["dumpnorm_fs.glsl"]);
  Link("convolute1",
       shader_code_map["quadpp_vs.glsl"],
       shader_code_map["convolute1_fs.glsl"]);
  Link("amboccl",
       shader_code_map["quadpp_vs.glsl"],
       shader_code_map["amboccl_fs.glsl"]);
  Link("beacon",
       shader_code_map["scenefx_vs.glsl"],
       shader_code_map["beacon_fs.glsl"]);
  Link("scenefx",
       shader_code_map["scenefx_vs.glsl"],
       shader_code_map["scenefx_fs.glsl"]);
  Link("screenblur4",
       shader_code_map["quadpp_vs.glsl"],
       shader_code_map["screenblur4_fs.glsl"]);
  Link("test_tex",
       shader_code_map["fraglight_vs.glsl"],
       shader_code_map["test_tex_fs.glsl"]);
  Link("anaglyph",
       shader_code_map["quadpp_vs.glsl"],
       shader_code_map["anaglyph_fs.glsl"]);

  valid_=true;
}

void Shader::Activate(const String& name)
{
  if(!OST_GL_VERSION_2_0) return;

  if(!name.empty()) {
    std::map<String, GLuint>::iterator it = shader_program_map_.find(name);
    if(it!=shader_program_map_.end()) {
      LOG_TRACE("switching to shader [" << name << "]");
      glUseProgram(it->second);
      current_program_=it->second;
      current_name_=name;

      UpdateState();
      return;

    } else {
      LOG_TRACE("shader program [" << name << "] not present");
      return;
    }
  }
  LOG_TRACE("switching to fixed pipeline");
  glUseProgram(0);
  current_program_=0;
  current_name_="";
}

GLuint Shader::GetCurrentProgram() const
{
  return current_program_;
}

String Shader::GetCurrentName() const
{
  return current_name_;
}

bool Shader::IsValid() const
{
  if(!OST_GL_VERSION_2_0) return false;
  return valid_;
}

bool Shader::IsActive() const
{
  return current_program_!=0;
}

void Shader::PushProgram()
{
  if(!OST_GL_VERSION_2_0) return;
  program_stack_.push(current_name_);
}

void Shader::PopProgram()
{
  if(!OST_GL_VERSION_2_0) return;
  if(!program_stack_.empty()) {
    current_name_ = program_stack_.top();
    program_stack_.pop();
    Activate(current_name_);
  }
}

void Shader::UpdateState()
{
  if(!OST_GL_VERSION_2_0) return;
  if(current_program_!=0) {
    // update all settings
    GLint result;
    glGetIntegerv(GL_LIGHTING,&result);
    LOG_TRACE("setting lighting flag to " << result);
    glUniform1i(glGetUniformLocation(current_program_,"lighting_flag"),result);
    GLboolean bresult;
    glGetBooleanv(GL_LIGHT_MODEL_TWO_SIDE,&bresult);
    LOG_TRACE("setting two_sided flag to " << bresult);
    glUniform1i(glGetUniformLocation(current_program_,"two_sided_flag"),bresult);
    if(impl::SceneFX::Instance().WillPostprocess()) {
      result=0;
    } else {
      glGetIntegerv(GL_FOG,&result);
    }
    LOG_TRACE("setting fog flag to " << result);
    glUniform1i(glGetUniformLocation(current_program_,"fog_flag"),result);
    glDisable(GL_COLOR_MATERIAL);

  } else {
    glEnable(GL_COLOR_MATERIAL);
  }
}

}} // ns
