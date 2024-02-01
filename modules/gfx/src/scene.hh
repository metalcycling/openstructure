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
#ifndef OST_SCENE_HH
#define OST_SCENE_HH

/*
  Author: Ansgar Philippsen
*/

#include <map>
#include <stack>
#include <vector>

#include <boost/shared_array.hpp>

#include <ost/gfx/module_config.hh>
#include <ost/geom/transform.hh>
#include <ost/mol/atom_handle.hh>

#include "gl_include.hh"
#include "color.hh"
#include "gfx_object_fw.hh"
#include "gfx_node_fw.hh"
#include "gfx_node_visitor.hh"
#include "glwin_base.hh"
#include "scene_observer.hh"
#include "gfx_prim.hh"
#include "povray_fw.hh"
#include "exporter_fw.hh"
#include "gradient.hh"
#include "bitmap_io.hh"

namespace ost { namespace gfx {

class InputEvent;

typedef std::vector<SceneObserver*>  SceneObserverList;

struct Viewport {
  int x;
  int y;
  int width;
  int height;
};

namespace impl {class SceneFX;}

/// \brief main class for organization and root for the graphical display
/// 
/// The scene manages graphical objects for rendering. Typical graphical objects 
/// include \ref gfx::Entity "entities", \ref gfx::MapIso "isocontoured maps", 
/// \ref gfx::MapSlab "density slabs", \ref gfx::Surface "molecular surfaces", 
/// or \ref gfx::Primitive "primitives" such as \ref gfx::Cuboid "cuboids", 
/// \ref gfx::Quad "quads" and \ref gfx::PrimList "lines". The nodes are 
/// organized in a tree-like structure whose root 
/// can obtained with GetRootNode().
/// 
/// The center of the eye is controlled with SetCenter() and CenterOn().
/// 
/// By default, the near and far clipping planes are adjusted such that they 
/// contain all visible objects in the scene. This behaviour can be turned of by
/// disabling the AutoAutoslab(). The near and far clipping plane can then be 
/// adjusted manually.
class DLLEXPORT_OST_GFX Scene {
  friend class impl::SceneFX;
 private:
   
  // TODO: this struct may be the seed for a proper
  // refactoring of the scene it into a management
  // and a view part
  struct SceneViewStackEntry {
    geom::Transform transform;
    float fov;
    float znear,zfar;
  };

  typedef std::stack<SceneViewStackEntry> SceneViewStack;

 public:
  static Scene& Instance();

  /// \brief turn fog on or off
  void SetFog(bool f);
  /// \brief check fog status
  bool GetFog() const;
  /// \brief set the fog color
  void SetFogColor(const Color& c);
  /// \brief get the fog color
  Color GetFogColor() const;

  /// \brief turn shadow mapping on and off
  void SetShadow(bool f);
  /// \brief get shadow mapping status
  bool GetShadow() const;
  /// \brief shadow quality from 0 (low) to 3 (high), default=1
  void SetShadowQuality(int q);
  /// \brief get shadow quality
  int GetShadowQuality() const;
  /// \brief multiplier for shadow strength
  void SetShadowWeight(float w);
  /// \brief get shadow strength
  float GetShadowWeight() const;

  /// experimental feature
  void SetDepthDarkening(bool f);
  /// experimental feature
  void SetDepthDarkeningWeight(float f);

  /// experimental feature
  void SetAmbientOcclusion(bool f);
  /// experimental feature
  bool GetAmbientOcclusion() const;
  /// experimental feature
  void SetAmbientOcclusionWeight(float f);
  /// experimental feature
  float GetAmbientOcclusionWeight() const;
  /// experimental feature
  void SetAmbientOcclusionMode(uint m);
  /// experimental feature
  uint GetAmbientOcclusionMode() const;
  /// experimental feature
  void SetAmbientOcclusionQuality(uint q);
  /// experimental feature
  uint GetAmbientOcclusionQuality() const;
  /// experimental feature
  void SetAmbientOcclusionSize(float f);
  /// experimental feature
  float GetAmbientOcclusionSize() const;
 
  void SetHemiParams(const geom::Vec4&);
  geom::Vec4 GetHemiParams() const {return hemi_param_;}

  /// \brief select shading mode
  /// one of fallback, basic, default, hf, toon1, toon2
  void SetShadingMode(const std::string& smode);

  geom::Mat4 GetProjection() const {return pmat_;}
  geom::Mat4 GetInvertedProjection() const {return ipmat_;}

  /// \name clipping planes, fog and field-of-view
  //@{
  /// \brief get near clipping plane
  float GetNear() const;
  
  /// \brief set near clipping plane
  void SetNear(float n);
  
  /// \brief get far clipping plane  
  float GetFar() const;
  
  /// \brief set far clipping plane  
  void SetFar(float f);
  
  /// \brief set near and far clipping plane at once
  void SetNearFar(float n, float f);

  /// \brief set field of view angle
  void SetFOV(float f);

  // \brief get the field of view
  float GetFOV() const;

  float GetAspect() const {return aspect_ratio_;}

  /// \brief offset between near clipping plane and start of fog
  void SetFogNearOffset(float o);

  /// \sa SetFogNearOffset
  float GetFogNearOffset() const;

  /// \brief offset between far clipping plane and end of fog
  void SetFogFarOffset(float o);

  /// \sa SetFogFarOffset
  float GetFogFarOffset() const;

  /// \brief convenciene function to set fog near and far offset
  void SetFogOffsets(float no, float fo);
  
  /// DEPRECATED, use Autoslab() and SetAutoslabMode(int)
  void Autoslab(bool fast);

  /// DEPRECATED, use Autoslab() and SetAutoslabMode(int)
  void Autoslab(bool fast, bool);

  /// DEPRECATED, use SetAutoslabMode(2)
  void AutoslabMax();

  /*!
    \brief adjust near and far clipping plane to fit visible objects

    Use autoslab mode to calculate near and far clipping places; this
    does not need to be called explicitely if AutoAutoslab is active.
    Uses the mode set by \ref SetAutoslabMode
  */
  void Autoslab();

  /*!
    \brief set autoslab mode

    0: fast (default), using only the bounding box
    1: precise, using each graphical element (not implemented)
    2: max, using maximal extent upon rotation
  */
  void SetAutoslabMode(int mode) {
    autoslab_mode_=std::min(2,std::max(0,mode));
  }

  /// \brief return current autoslab mode
  int GetAutoslabMode() const {
    return autoslab_mode_;
  }

  /*!
    \brief turn automatic autoslab'bing on or off for each scene update

    the current autoslab mode is honored \ref SetAutoslabMode(int)
  */
  void AutoAutoslab(bool f);
  
  /// \brief get current state of automatic autoslab'bing
  bool GetAutoAutoslab() const { return auto_autoslab_; }

  //@}

  /// \brief set stereo mode
  /// one of 0 (off), 1 (quad-buffered) 2 (interlaced (for special monitors))
  void SetStereoMode(unsigned int mode);
  int GetStereoMode() const {return stereo_mode_;}

  /// \brief invert stereo eyes for stereo mode=0
  void SetStereoFlip(bool f);
  /// \brief return invert flag for stereo
  bool GetStereoFlip() const {return stereo_inverted_;}
  
  /// \brief stereo view mode
  /// one of 0 (center), -1 (left), 1 (right)
  void SetStereoView(int);
  /// \brief return current stereo view mode
  int GetStereoView() const {return stereo_eye_;}

  /// \brief set stereo eye distance
  void SetStereoIOD(Real);
  /// \brief return current stereo eye distance
  Real GetStereoIOD() const {return stereo_iod_;}

  /// \brief set stereo distance offset from COR
  void SetStereoDistance(Real);
  /// \brief return current stereo distance offset from COR
  Real GetStereoDistance() const {return stereo_distance_;}
  
  /// \brief set stereo algorithm
  /// one of 0 (default) or 1
  void SetStereoAlg(unsigned int);
  /// \brief return current stereo algorithm
  unsigned int GetStereoAlg() const {return stereo_alg_;}
  
  /// \brief set main light direction
  void SetLightDir(const geom::Vec3& dir);
  /// \brief set ambient, diffuse and specular light color
  void SetLightProp(const Color& amb, const Color& diff, const Color& spec);
  /// \brief set ambient, diffuse and specular light intensity
  void SetLightProp(float amb, float diff, float spec);
  /// \brief get main light direction
  geom::Vec3 GetLightDir() const {return light_dir_;}
  /// \brief get main light orientation (internal debugging use)
  geom::Mat3 GetLightRot() const {return light_rot_;}

  /// \brief set the selection mode
  /*
    bad style for now: 0=(reserved), 1=atom, 2=residue, 3=chain, 4=bond, 5=torsion
  */
  void SetSelectionMode(uint m);
  uint GetSelectionMode() const;


  void Export(const String& fname, unsigned int w,
              unsigned int h, bool transparent=false);
  /// \brief export into bitmap, using multisample anti-aliasing
  void Export(const String& fname, unsigned int w,
              unsigned int h, int max_samples, bool transparent=false);

  /// \brief export snapshot of current scene
  void Export(const String& fname, bool transparent=false);

  /// \brief export scene into povray files named fname.pov and fname.inc
  void ExportPov(const std::string& fname, const std::string& wdir=".");

  /// \brief export scene via exporter
  void Export(Exporter* ex) const;

  //@}
  /// \brief entry point for gui events (internal use)
  void OnInput(const InputEvent& e);
  
  /// \brief initialize OpenGL after context has been setup (internal use)
  void InitGL(bool full=true);

  /// \brief handle new viewport size (internal use)
  void Resize(int w, int h);

  /// \brief pick at given mouse coords
  void Pick(int mx, int my, int mask);
  
  float GetDefaultTextSize();

  /// \brief pick atom at given mouse coord
  std::pair<GfxObjP, mol::AtomHandle> PickAtom(int mx, int my);
  
  /// \brief render all gl objects (internal use)
  void RenderGL();

  /// \brief request redraw of gl scene
  void RequestRedraw();
  
  /// \brief send status message to gui
  void StatusMessage(const String& s);

  /// \brief set the viewport; the mapping to the visible window (internal use)
  void SetViewport(int w, int h);

  /// \brief set background color
  void SetBackground(const Color& c);

  /// \brief set background gradient
  void SetBackground(const Gradient& g);

  /// \brief set background image
  void SetBackground(const Bitmap& bm);

  /// \brief get background color
  Color GetBackground() const;

  /// \brief use bg bitmap for stereo mode
  /// this tiles the background bitmap at the far plane
  void SetBackgroundStereoMode(bool);
  bool GetBackgroundStereoMode() const {return bg_stereo_mode_;}

  /// background tile left/right offset
  void SetBackgroundStereoOffset(float);
  float GetBackgroundStereoOffset() const {return bg_stereo_offset_;}

  /// \brief center rotation on the given point
  void SetCenter(const geom::Vec3& cen);
  
  /// \brief retrieve center
  geom::Vec3 GetCenter() const;

  /// \brief center on object of given name
  void CenterOn(const String& s);

  /// \brief center given object
  void CenterOn(const GfxObjP& s);

  /// \brief calculate projection of a point into the scene
  geom::Vec3 Project(const geom::Vec3& v, bool ignore_vp=false) const;
  
  /// \brief calculate unprojected point out of the scene
  geom::Vec3 UnProject(const geom::Vec3& v, bool ignore_vp=false) const;

  /// \brief return bounding box of scene
  /*!
    the sole boolean parameter determines whether or not the scene 
    transformation is applied to calculate the bounding box. Since in
    most cases it should be used, the default value is true.
  */
  geom::AlignedCuboid GetBoundingBox(bool use_tf=true) const;

  /// \brief return bounding box of with a given transform
  geom::AlignedCuboid GetBoundingBox(const geom::Transform& tf) const;

  /// \brief get full underlying transformation
  geom::Transform GetTransform() const;

  /// \brief set transform
  void SetTransform(const geom::Transform& t);

  /// \brief returns a compact, internal representation of the scene orientation
  geom::Mat4 GetRTC() const;

  /// \brief sets a previously retrieved orientation
  void SetRTC(const geom::Mat4& rtc);
  
  /// \brief push the current orientation onto a stack
  void PushView();

  /// \brief retrieve a previously pushed orientation
  void PopView();

  /// brief re-generates the projection matrix (internal use)
  void ResetProjection();

  /// \brief gui glue interface (internal use)
  void Register(GLWinBase* win);
  /// \brief gui glue interface (internal use)
  void Unregister(GLWinBase* win);
  
  /// \name scene graph
  //@{
  /// \brief add graphical object to scene
  void Add(const GfxNodeP& go, bool redraw=true);
  /// \brief remove graphical object from scene
  /// remove graphical object from the scene
  void Remove(const GfxNodeP& go);
  /// remove graphical object from the scene
  void Remove(const String& name);
  
  /// \brief remove all objects from the scene
  void RemoveAll();
  
  /// \brief rename an existing graphical object
  /// defunct for now
  bool Rename(const String& old_name, const String& new_name);

  /// \brief retrieve gfx object by name
  GfxObjP operator[](const String& name);

  /// \brief whether the scene contains a node of the given name
  bool HasNode(const String& name) const;
  
  /// \brief actual event handling for scene (internal use)
  void Apply(const InputEvent& ie, bool request_redraw=true);

  /// \brief apply node visitor to root node
  void Apply(GfxNodeVisitor& v) const;
  
  /// \brief get total number of nodes in scene
  /// 
  /// To obtain the number of top-level nodes, use GfxNode::GetChildCount() of
  /// the root node
  size_t GetNodeCount() const;
  
  /// \brief get root node of scene graph
  GfxNodeP GetRootNode() const;
  //@}
  
  /// \brief observer interface (internal use)
  void AttachObserver(SceneObserver* o);
  /// \brief observer interface (internal use)
  void DetachObserver(SceneObserver* o);

  /// \brief switch into test mode (internal use)
  void SetTestMode(bool t);

  float ElapsedTime() const;

  Viewport GetViewport() const;


  /// \brief show center of rotation of true
  void SetShowCenter(bool f);

  bool GetShowCenter() const {return cor_flag_;}

  /// \brief if true fix center of rotation upon input induced shift
  void SetFixCenter(bool f) {fix_cor_flag_=f;}

  /// \brief return flag
  bool GetFixCenter() const {return fix_cor_flag_;}
  
  /// experimental feature
  void SetBlur(uint n);
  /// experimental feature
  void BlurSnapshot();

  /// internal use
  void RenderText(const TextPrim& t);

  /// experimental feature
  void SetBeacon(int wx, int wy);
  /// experimental feature
  void SetBeaconOff();

  void SetExportAspect(float a);
  float GetExportAspect() const {return export_aspect_;}
  void SetShowExportAspect(bool f);
  bool GetShowExportAspect() const {return show_export_aspect_;}

  bool HasMultisample() const {return ms_flag_;}

  void SetAlphaBias(Real bias);

  void ContextSwitch();
  
protected:
  friend class GfxObj; 
  friend class GfxNode;

  // TODO: this is really a hack and not clean communication
  friend class Entity;
  
  void ObjectChanged(const String& name);
  void SelectionChanged(const String& name, const mol::EntityView& sel);
  void NodeTransformed(const GfxObjP& object);
  void NodeAdded(const GfxNodeP& node);
  void RenderModeChanged(const String& name);

private:  

  template <typename ACTION>
  void NotifyObservers(const ACTION& action) {
    std::for_each(observers_.begin(), observers_.end(), action);
  }
  Scene();
  Scene(const Scene&) {}
  Scene& operator=(const Scene&) {return *this;}
  ~Scene();

  // The GL Context gets typically handled externally, e.g. using QT when
  // calling InitGL(), RenderGL() etc. However, there is no guarantee that the 
  // desired GLContext is active when calling gl related functions in other 
  // occasions. This can be enforced by ActivateGLContext()
  void ActivateGLContext() const;

  GLWinBase*           win_; // target gl window

  mutable GfxNodeP     root_node_; // mutable is slightly hackish
  SceneObserverList    observers_;

  geom::Transform transform_; // overall modelview transformation

  bool gl_init_;

  float fov_; // field of view
  float znear_,zfar_; // near and far clipping plane
  float fnear_,ffar_; // fog near and far offsets

  geom::Mat4 pmat_,ipmat_; // projection and inverted projection matrix
  unsigned int vp_width_,vp_height_; // viewport

  SceneViewStack scene_view_stack_;

  float aspect_ratio_; // aspect ratio for viewport

  Color background_; // background (clear) color

  geom::Vec3 light_dir_; // infinite light source direction
  geom::Mat3 light_rot_; // transform needed for the shadow map
  Color light_amb_;
  Color light_diff_;
  Color light_spec_;
  geom::Vec4 hemi_param_;

  bool cor_flag_;
  bool fix_cor_flag_;
  bool fog_flag_;
  Color fog_color_;
  bool auto_autoslab_; // run autoslab on each scene update
  bool do_autoslab_;   // run autoslab on next scene update
  int autoslab_mode_;  // 0: fast, 1:precise, 2: max

  std::string def_shading_mode_;

  uint selection_mode_;

  bool test_flag_;
  std::vector<unsigned char> tmp_tex_;

  GLuint glyph_tex_id_;
  std::vector<geom::Vec2> glyph_map_;
  float def_text_size_;

  uint blur_count_;
  std::vector<boost::shared_array<unsigned char> > blur_buffer_;

  unsigned int stereo_mode_;
  unsigned int stereo_alg_;
  bool stereo_inverted_;
  int stereo_eye_;
  Real stereo_iod_,stereo_distance_;
  unsigned int scene_left_tex_;
  unsigned int scene_right_tex_;
  unsigned int bg_mode_;
  bool update_bg_;
  bool bg_stereo_mode_;
  float bg_stereo_offset_;
  Gradient bg_grad_;
  Bitmap bg_bm_;
  unsigned int bg_tex_;
  
  bool ms_flag_; // multisample flag

  float export_aspect_; 
  bool show_export_aspect_;

  void set_near(float n);
  void set_far(float f);
  void flag_all_dirty();
  void prep_glyphs();
  void prep_blur();
  void stereo_projection(int view);
  void render_bg();
  void render_export_aspect();
  void render_scene();
  void render_glow();
  void render_stereo();
  void set_bg();
  void do_autoslab();

  bool IsNameAvailable(const String& name) const;

};

}} // ns

#endif
