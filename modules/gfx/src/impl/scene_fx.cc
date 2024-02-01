#include <boost/random.hpp>

#include "scene_fx.hh"

#include <ost/log.hh>

#include <ost/gfx/gfx_node.hh>
#include <ost/gfx/shader.hh>
#include <ost/gfx/scene.hh>

#if !GL_VERSION_3_0
#warning using OpenGL 2.0 interface
#endif

namespace ost { namespace gfx { namespace impl {

namespace {
  boost::mt19937 RandomGenerator(time(NULL));
  boost::uniform_01<boost::mt19937> UniformRandom(RandomGenerator);
}

SceneFX& SceneFX::Instance()
{
  static SceneFX inst;
  return inst;
}

SceneFX::SceneFX():
  shadow_flag(false),
  shadow_quality(1),
  shadow_weight(1.0),
  depth_dark_flag(false),
  depth_dark_factor(1.0),
  amb_occl_flag(false),
  amb_occl_factor(0.5),
  amb_occl_mode(1),
  amb_occl_quality(1),
  amb_occl_size(80.0),
  use_beacon(false),
  beacon(),
  scene_tex_id_(),
  depth_tex_id_(),
  shadow_tex_id_(),
  shadow_tex_mat_(),
  occl_tex_id_(),
  dark_tex_id_(),
  norm_tex_id_(),
  kernel_tex_id_(),
  kernel_size_(0),
  kernel2_tex_id_(),
  kernel2_size_(0),
  scene_tex2_id_(),
  norm_tex2_id_(),
  scene_fb_(),
  scene_rb_(),
  depth_rb_(),
  use_fb_(false)
{}

SceneFX::~SceneFX()
{}

void SceneFX::Setup()
{
  if(!OST_GL_VERSION_2_0) return;
  LOG_DEBUG("SceneFX: setup");

  glGenTextures(1,&scene_tex_id_);
  glGenTextures(1,&depth_tex_id_);
  glGenTextures(1,&shadow_tex_id_);
  glGenTextures(1,&occl_tex_id_);
  glGenTextures(1,&dark_tex_id_);
  glGenTextures(1,&norm_tex_id_);
  glGenTextures(1,&kernel_tex_id_);
  glGenTextures(1,&kernel2_tex_id_);

#if GL_VERSION_3_0
  glGenFramebuffers(1,&scene_fb_);
  glGenRenderbuffers(1,&scene_rb_);
  glGenRenderbuffers(1,&depth_rb_);
#else
  glGenFramebuffersEXT(1,&scene_fb_);
  glGenRenderbuffersEXT(1,&scene_rb_);
  glGenRenderbuffersEXT(1,&depth_rb_);
#endif
  glGenTextures(1,&scene_tex2_id_);
  glGenTextures(1,&norm_tex2_id_);

  glBindTexture(GL_TEXTURE_2D, scene_tex_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glBindTexture(GL_TEXTURE_2D, depth_tex_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glBindTexture(GL_TEXTURE_2D, shadow_tex_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glBindTexture(GL_TEXTURE_2D, occl_tex_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glBindTexture(GL_TEXTURE_2D, dark_tex_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glBindTexture(GL_TEXTURE_2D, norm_tex_id_);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_1D, kernel_tex_id_);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);

  std::vector<GLfloat> tmp;
  for(int u=-20;u<=20;++u) {
    for(int v=-20;v<=20;++v) {
      if(UniformRandom()<0.1) {
        // norm to 1
        float x = static_cast<float>(u+20)*0.025;
        float y = static_cast<float>(v+20)*0.025;
        tmp.push_back(x);
        tmp.push_back(y);
        tmp.push_back(0.0);
        tmp.push_back(0.0);
      }
    }
  }

  kernel_size_=tmp.size()/4;
  glTexImage1D(GL_TEXTURE_1D,0,4,kernel_size_,0,GL_RGBA, GL_FLOAT, &tmp[0]);

  glBindTexture(GL_TEXTURE_1D, kernel2_tex_id_);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  
  // depth darkening kernel
  tmp.clear();
  for(int u=-5;u<=5;++u) {
    for(int v=-5;v<=5;++v) {
      if(UniformRandom()<0.5) {
        // norm to -1 ... 1
        float x = static_cast<float>(u)*0.2;
        float y = static_cast<float>(v)*0.2;
        // gaussian kernel with sigma
        float v = exp((-x*x-y*y)*0.5);
        // norml to 0 ... 1
        tmp.push_back(x*0.5+0.5);
        tmp.push_back(y*0.5+0.5);
        tmp.push_back(v);
        tmp.push_back(0.0);
      }
    }
  }

  kernel2_size_=tmp.size()/4;
  glTexImage1D(GL_TEXTURE_1D,0,4,kernel2_size_,0,GL_RGBA, GL_FLOAT, &tmp[0]);
}

void SceneFX::Resize(unsigned int w, unsigned int h)
{
#if 0
  glBindTexture(GL_TEXTURE_2D, scene_tex2_id_);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
  glBindTexture(GL_TEXTURE_2D, norm_tex2_id_);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
  glBindTexture(GL_TEXTURE_2D, 0);

  glBindFramebuffer(GL_FRAMEBUFFER, scene_fb_);

  glBindRenderbuffer(GL_RENDERBUFFER, depth_rb_);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,w,h);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_rb_);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, scene_tex2_id_, 0);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, norm_tex2_id_, 0);

  GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);

  if(status!=GL_FRAMEBUFFER_COMPLETE) {
    use_fb_=false;
    LOG_VERBOSE("SceneFX: framebuffer error code " << status);
  } else {
    use_fb_=true;
  }

  glBindRenderbuffer(GL_RENDERBUFFER, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
#endif
}

void SceneFX::Preprocess() 
{
  if(!OST_GL_VERSION_2_0) return;
  if(use_fb_) {
#if GL_VERSION_3_0
    glBindFramebuffer(GL_FRAMEBUFFER, scene_fb_);
#else
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, scene_fb_);
#endif
  }
}

bool SceneFX::WillPostprocess() const
{
  return shadow_flag || amb_occl_flag || depth_dark_flag || use_beacon;
}

void SceneFX::Postprocess()
{
  if(!OST_GL_VERSION_2_0) return;
  if(use_fb_) {
#if GL_VERSION_3_0
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
#else
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
#endif
  }

  if(!WillPostprocess()) return;

  Viewport vp=Scene::Instance().GetViewport();

  if(!use_fb_) {
    // grab color buffer
    glBindTexture(GL_TEXTURE_2D, scene_tex_id_);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, vp.x, vp.y, vp.width, vp.height, 0);
    // and depth buffer
    glBindTexture(GL_TEXTURE_2D, depth_tex_id_);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, vp.x, vp.y, vp.width, vp.height, 0);

    if(amb_occl_flag) {
      // now for the normal buffer hack if the framebuffer stuff failed
      Shader::Instance().PushProgram();
      Shader::Instance().Activate("dumpnorm");
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      Scene::Instance().GetRootNode()->RenderGL(STANDARD_RENDER_PASS);
      glBindTexture(GL_TEXTURE_2D, norm_tex_id_);
      glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, vp.x, vp.y, vp.width, vp.height, 0);
      Shader::Instance().PopProgram();
    }
  }

  if(shadow_flag) {
    prep_shadow_map();
  }
  if(amb_occl_flag) {
    prep_amb_occlusion();
  }
  if(depth_dark_flag) {
    prep_depth_darkening();
  }

  Shader::Instance().PushProgram();
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_TEXTURE_1D);
  Shader::Instance().Activate("scenefx");
  GLuint cpr=Shader::Instance().GetCurrentProgram();
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D,depth_tex_id_);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D,scene_tex_id_);

  glUniform1i(glGetUniformLocation(cpr,"scene_map"),0);
  glUniform1i(glGetUniformLocation(cpr,"depth_map"),1);

  glUniform2f(glGetUniformLocation(cpr,"i_vp"),
              1.0/static_cast<float>(vp.width),
              1.0/static_cast<float>(vp.height));

  if(shadow_flag) {
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D,shadow_tex_id_);
    glMatrixMode(GL_TEXTURE);
    glPushMatrix();
    glLoadTransposeMatrix(shadow_tex_mat_.Data());
    glMatrixMode(GL_MODELVIEW);
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(glGetUniformLocation(cpr,"shadow_flag"),1);
    glUniform1i(glGetUniformLocation(cpr,"shadow_map"),2);
    glUniform1f(glGetUniformLocation(cpr,"shadow_depth_bias"),0.008);
    glUniform1f(glGetUniformLocation(cpr,"shadow_epsilon"),0.002);
    glUniform1f(glGetUniformLocation(cpr,"shadow_multiplier"),0.4);
    glUniform1f(glGetUniformLocation(cpr,"shadow_weight"),shadow_weight);

  } else {
    glUniform1i(glGetUniformLocation(cpr,"shadow_flag"),0);
  }
  if(amb_occl_flag) {
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D,occl_tex_id_);
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(glGetUniformLocation(cpr,"occl_flag"),1);
    glUniform1i(glGetUniformLocation(cpr,"occl_map"),3);
    glUniform1f(glGetUniformLocation(cpr,"occl_mult"),amb_occl_factor);
  } else {
    glUniform1i(glGetUniformLocation(cpr,"occl_flag"),0);
  }

  if(depth_dark_flag) {
    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_2D,dark_tex_id_);
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(glGetUniformLocation(cpr,"dark_flag"),1);
    glUniform1i(glGetUniformLocation(cpr,"dark_map"),4);
    glUniform1f(glGetUniformLocation(cpr,"dark_mult"),depth_dark_factor);
  } else {
    glUniform1i(glGetUniformLocation(cpr,"dark_flag"),0);
  }

  GLint fog;
  glGetIntegerv(GL_FOG,&fog);
  float znear=Scene::Instance().GetNear();
  float zfar=Scene::Instance().GetFar();
  float fnear=znear+Scene::Instance().GetFogNearOffset();
  float ffar=zfar+Scene::Instance().GetFogFarOffset();
  glUniform1i(glGetUniformLocation(cpr,"fog_flag"),fog);
  glUniform1f(glGetUniformLocation(cpr,"depth_near"),znear);
  glUniform1f(glGetUniformLocation(cpr,"depth_far"),zfar);
  glUniform1f(glGetUniformLocation(cpr,"fog_far"),ffar);
  glUniform1f(glGetUniformLocation(cpr,"fog_scale"),1.0f/(ffar-fnear));

  if(use_beacon) {
    prep_beacon();
  }

  draw_screen_quad(vp.width,vp.height);

  if(shadow_flag) {
    glActiveTexture(GL_TEXTURE2);
    glMatrixMode(GL_TEXTURE);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glActiveTexture(GL_TEXTURE0);
  }

  if(use_beacon) {
    draw_beacon();
  }

  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  Shader::Instance().PopProgram();
}

void SceneFX::DrawTex(unsigned int w, unsigned int h, GLuint texid)
{
  if(!OST_GL_VERSION_2_0) return;
  Shader::Instance().PushProgram();
  Shader::Instance().Activate("");
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texid);

  draw_screen_quad(w,h);

  Shader::Instance().PopProgram();
}

void SceneFX::prep_shadow_map()
{
  GLint smap_size=256 * (1+shadow_quality);

#if GL_VERSION_3_0
  glBindFramebuffer(GL_FRAMEBUFFER, scene_fb_);
  glBindRenderbuffer(GL_RENDERBUFFER, scene_rb_);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA,smap_size,smap_size);
  glBindRenderbuffer(GL_RENDERBUFFER, depth_rb_);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,smap_size,smap_size);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, scene_rb_);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_rb_);
  
  GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
  
  if(status!=GL_FRAMEBUFFER_COMPLETE) {
    LOG_DEBUG("fbo switch for shadow mapping failed, using fallback");
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    smap_size=512;
  }
#else
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, scene_fb_);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, scene_rb_);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_RGBA,smap_size,smap_size);
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, depth_rb_);
  glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT,smap_size,smap_size);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_RENDERBUFFER_EXT, scene_rb_);
  glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, depth_rb_);
  
  GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
  
  if(status!=GL_FRAMEBUFFER_COMPLETE_EXT) {
    LOG_DEBUG("fbo switch for shadow mapping failed, using fallback");
    glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
    smap_size=512;
  }
#endif

  // modelview transform for the lightsource pov
  geom::Transform ltrans(Scene::Instance().GetTransform());
  ltrans.SetRot(Scene::Instance().GetLightRot()*ltrans.GetRot());

  // calculate encompassing box for ortho projection
  geom::AlignedCuboid bb=Scene::Instance().GetBoundingBox(ltrans);
  geom::Vec3 tmin=bb.GetMin();
  geom::Vec3 tmax=bb.GetMax();

  // save overall gl settings
  glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
  // maximize rendering for depth-only information
  glDisable(GL_LIGHTING);
  glDisable(GL_FOG);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_NORMALIZE);
  //glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);

  // render scene with only depth components
  // seen from the light's perspective
  glViewport(0,0,smap_size,smap_size);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(tmin[0],tmax[0],tmin[1],tmax[1],-tmax[2],-tmin[2]);
  //glFrustum(tmin[0],tmax[0],tmin[1],tmax[1],-tmax[2],-tmin[2]);
  float glpmat[16];
  glGetv(GL_PROJECTION_MATRIX, glpmat);
  geom::Mat4 pmat(geom::Transpose(geom::Mat4(glpmat)));

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMultMatrix(ltrans.GetTransposedMatrix().Data());

  // only render non-transparent objects for the shadow map
  Scene::Instance().GetRootNode()->RenderGL(STANDARD_RENDER_PASS);

  // now get the shadow map
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, shadow_tex_id_);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 0,0, smap_size,smap_size, 0);

  // restore settings
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);

  glPopAttrib();
  //glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);

#if GL_VERSION_3_0
  glBindRenderbuffer(GL_RENDERBUFFER, 0);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
#else
  glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
  glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
#endif

  // set up appropriate texture matrix
  geom::Mat4 bias(0.5,0.0,0.0,0.5,
                  0.0,0.5,0.0,0.5,
                  0.0,0.0,0.5,0.5,
                  0.0,0.0,0.0,1.0);
  //shadow_tex_mat_ = bias*pmat*ltrans.GetMatrix();
  glGetv(GL_PROJECTION_MATRIX, glpmat);
  geom::Mat4 pmat2(geom::Transpose(geom::Mat4(glpmat)));
  /*
    given the normalized coordinates in scenefx, the camera projection and modelview transformation
    are first reverted, and then the light modelview and projection are applied, resulting (with the
    bias) in the proper 2D lookup into the shadow map
  */
  try {
    shadow_tex_mat_ = bias*pmat*ltrans.GetMatrix()*geom::Invert(Scene::Instance().GetTransform().GetMatrix())*geom::Invert(pmat2);
  } catch (geom::GeomException& e) {
    LOG_DEBUG("SceneFX: caught matrix inversion exception in prep_shadow_map()");
  }
}

void SceneFX::prep_amb_occlusion()
{
  Viewport vp=Scene::Instance().GetViewport();

  uint qf=1;
  if(amb_occl_quality==0) {qf=4;}
  else if(amb_occl_quality==1) {qf=2;}
  uint width=vp.width/qf;
  uint height=vp.height/qf;

  Shader::Instance().PushProgram();
  Shader::Instance().Activate("amboccl");
  GLuint cpr=Shader::Instance().GetCurrentProgram();
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D,depth_tex_id_);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D,norm_tex_id_);
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_1D,kernel_tex_id_);
  glActiveTexture(GL_TEXTURE0);
  glUniform1i(glGetUniformLocation(cpr,"depth_map"),0);
  glUniform1i(glGetUniformLocation(cpr,"norm_map"),1);
  glUniform1i(glGetUniformLocation(cpr,"kernel"),2);
  glUniform1f(glGetUniformLocation(cpr,"step"),1.0/static_cast<float>(kernel_size_));
  glUniform1f(glGetUniformLocation(cpr,"size"),amb_occl_size);
  glUniform2f(glGetUniformLocation(cpr,"i_vp"),1.0/static_cast<float>(width),1.0/static_cast<float>(height));
  glUniform1i(glGetUniformLocation(cpr,"mode"),amb_occl_mode);
  double pm[16];
  glGetDoublev(GL_PROJECTION_MATRIX,pm);
  glUniform4f(glGetUniformLocation(cpr,"abcd"),pm[0],pm[5],pm[10],pm[14]);

  glMatrixMode(GL_TEXTURE);
  glPushMatrix();
  geom::Mat4 ipm;
  try {
    ipm=(geom::Transpose(geom::Invert(geom::Transpose(geom::Mat4(pm)))));
  } catch (geom::GeomException& e) {
    LOG_DEBUG("SceneFX: caught matrix inversion exception in prep_amb_occlusion()");
  }
  glLoadMatrix(ipm.Data());
  glMatrixMode(GL_MODELVIEW);

  // set up viewport filling quad to run the fragment shader
  draw_screen_quad(width,height);

  glMatrixMode(GL_TEXTURE);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, occl_tex_id_);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0,0, width, height, 0);
  Shader::Instance().PopProgram();
}

void SceneFX::prep_depth_darkening()
{
  Viewport vp=Scene::Instance().GetViewport();

  Shader::Instance().PushProgram();
  Shader::Instance().Activate("convolute1");
  GLuint cpr=Shader::Instance().GetCurrentProgram();
  // assign tex units
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D,depth_tex_id_);
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_1D,kernel2_tex_id_);
  glActiveTexture(GL_TEXTURE0);

  glUniform1i(glGetUniformLocation(cpr,"data"),1);
  glUniform1i(glGetUniformLocation(cpr,"kernel"),2);
  glUniform1f(glGetUniformLocation(cpr,"step"),1.0/static_cast<float>(kernel2_size_));
  glUniform2f(glGetUniformLocation(cpr,"i_vp"),1.0/static_cast<float>(vp.width),1.0/static_cast<float>(vp.height));

  // set up viewport filling quad to run the fragment shader
  draw_screen_quad(vp.width/2,vp.height/2);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, dark_tex_id_);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0,0,vp.width/2, vp.height/2, 0);

  Shader::Instance().PopProgram();
}

void SceneFX::draw_screen_quad(unsigned int w, unsigned int h) 
{
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  
  // setup
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_FOG);
  glDisable(GL_CULL_FACE);
  glDisable(GL_BLEND);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_POINT_SMOOTH);
  glShadeModel(GL_FLAT);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0,1,0,1,-1,1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE, GL_REPLACE);
  // draw
  glColor3f(1.0,0.0,1.0);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0,0.0);
  glVertex2f(0.0,0.0);
  glTexCoord2f(0.0,1.0);
  glVertex2f(0.0,1.0);
  glTexCoord2f(1.0,1.0);
  glVertex2f(1.0,1.0);
  glTexCoord2f(1.0,0.0);
  glVertex2f(1.0,0.0);
  glEnd();
  
  // restore settings
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glPopAttrib();
}

void SceneFX::prep_beacon()
{
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  Scene::Instance().stereo_projection(0);
  GLdouble glpmat[16];
  glGetDoublev(GL_PROJECTION_MATRIX, glpmat);
  geom::Mat4 pmat(geom::Transpose(geom::Mat4(glpmat)));
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  geom::Mat4 mmat=Scene::Instance().GetTransform().GetTransposedMatrix();
  GLdouble glmmat[16];
  for(int i=0;i<16;++i) glmmat[i]=static_cast<GLdouble>(mmat.Data()[i]);
  GLint glvp[4];
  glGetIntegerv(GL_VIEWPORT,glvp);

  double res[3];
  gluUnProject(beacon.wx,beacon.wy,0.0,
	       glmmat,glpmat,glvp,
	       &res[0],&res[1],&res[2]);
  beacon.p0=geom::Vec3(res[0],res[1],res[2]);
  gluUnProject(beacon.wx,beacon.wy,1.0,
	       glmmat,glpmat,glvp,
	       &res[0],&res[1],&res[2]);
  beacon.p1=geom::Vec3(res[0],res[1],res[2]);

  glGetv(GL_PROJECTION_MATRIX, glpmat);
  geom::Mat4 pmat2(geom::Transpose(geom::Mat4(glpmat)));
  try {
    beacon.mat = geom::Invert(Scene::Instance().GetTransform().GetMatrix())*geom::Invert(pmat2);
  } catch (geom::GeomException& e) {
    beacon.mat=geom::Mat4();
  }
}

void SceneFX::draw_beacon()
{
  Viewport vp=Scene::Instance().GetViewport();
  float iw=1.0/static_cast<float>(vp.width);
  float ih=1.0/static_cast<float>(vp.height);
  Shader::Instance().PushProgram();
  Shader::Instance().Activate("beacon");
  uint cpr=Shader::Instance().GetCurrentProgram();

  float rad=0.6;
  float delta = rad/(std::max<float>(1.0,Scene::Instance().znear_) * Scene::Instance().aspect_ratio_ * std::tan(Scene::Instance().fov_*M_PI/360.0))/iw;
  geom::Vec2 q0(beacon.wx-delta,beacon.wy-delta);
  geom::Vec2 q1(beacon.wx-delta,beacon.wy+delta);
  geom::Vec2 q2(beacon.wx+delta,beacon.wy+delta);
  geom::Vec2 q3(beacon.wx+delta,beacon.wy-delta);

  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D,depth_tex_id_);
  glMatrixMode(GL_TEXTURE);
  glPushMatrix();
  glLoadTransposeMatrix(beacon.mat.Data());
  glActiveTexture(GL_TEXTURE0);
  glUniform1i(glGetUniformLocation(cpr,"depth_map"),1);
  glUniform3f(glGetUniformLocation(cpr,"pos"),beacon.p0[0],beacon.p0[1],beacon.p0[2]);
  geom::Vec3 dir=beacon.p1-beacon.p0;
  glUniform3f(glGetUniformLocation(cpr,"dir"),dir[0],dir[1],dir[2]);
  glUniform1f(glGetUniformLocation(cpr,"len"),geom::Length(dir));
  glUniform1f(glGetUniformLocation(cpr,"rad"),rad);

  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0,vp.width,0,vp.height,-1,1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_FOG);
  glDisable(GL_CULL_FACE);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_POINT_SMOOTH);
  glShadeModel(GL_FLAT);
  glViewport(0,0,vp.width,vp.height);
  glEnable(GL_TEXTURE_2D);
  glColor3f(1.0,0.0,1.0);
  glBegin(GL_QUADS);
  glTexCoord2f(iw*q0[0],ih*q0[1]);
  glVertex2f(q0[0],q0[1]);
  glTexCoord2f(iw*q1[0],ih*q1[1]);
  glVertex2f(q1[0],q1[1]);
  glTexCoord2f(iw*q2[0],ih*q2[1]);
  glVertex2f(q2[0],q2[1]);
  glTexCoord2f(iw*q3[0],ih*q3[1]);
  glVertex2f(q3[0],q3[1]);
  glEnd();
  glActiveTexture(GL_TEXTURE1);
  glMatrixMode(GL_TEXTURE);
  glPopMatrix();
  glActiveTexture(GL_TEXTURE0);
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  glPopAttrib();
  Shader::Instance().PopProgram();
}

void SceneFX::screenblur4()
{
  Viewport vp=Scene::Instance().GetViewport();
  glBindTexture(GL_TEXTURE_2D, scene_tex_id_);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, vp.x, vp.y, vp.width, vp.height, 0);
  glBindTexture(GL_TEXTURE_2D, depth_tex_id_);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, vp.x, vp.y, vp.width, vp.height, 0);
  Shader::Instance().PushProgram();
  Shader::Instance().Activate("screenblur4");
  GLuint cpr=Shader::Instance().GetCurrentProgram();
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D,depth_tex_id_);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D,scene_tex_id_);
  glUniform1i(glGetUniformLocation(cpr,"scene_map"),0);
  glUniform1i(glGetUniformLocation(cpr,"depth_map"),1);
  glUniform2f(glGetUniformLocation(cpr,"i_vp"),
              1.0/static_cast<float>(vp.width),
              1.0/static_cast<float>(vp.height));
  glUniform1f(glGetUniformLocation(cpr,"depth_cutoff"),0.5);
  draw_screen_quad(vp.width,vp.height);
  Shader::Instance().PopProgram();
}

}}} // ns
