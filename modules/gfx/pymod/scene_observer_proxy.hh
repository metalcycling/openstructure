#ifndef OST_GFX_SCENE_OBSERVER_PROXY_HH
#define OST_GFX_SCENE_OBSERVER_PROXY_HH

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>

#include <ost/mol/view_type_fw.hh>

#include <ost/gfx/gfx_node.hh>
#include <ost/gfx/scene_observer.hh>

#include <ost/gfx/gfx_node_fw.hh>
#include <ost/gfx/gfx_object_fw.hh>

using namespace ost;
using namespace ost::gfx;


class SceneObserverProxy: public SceneObserver {
public:
  SceneObserverProxy(PyObject *p) :
    self(p) {
  }

  SceneObserverProxy(PyObject *p, const SceneObserver& i) :
    SceneObserver(i), self(p) {
  }

  virtual void NodeAdded(const GfxNodeP& node) {
    call_method<void, GfxNodeP> (self, "NodeAdded", node);
  }

  virtual void NodeRemoved(const GfxNodeP& node) {
    call_method<void, GfxNodeP> (self, "NodeRemoved", node);
  }

  virtual void NodeChanged(const GfxNodeP& node) {
    call_method<void, GfxNodeP> (self, "NodeChanged", node);
  }

  virtual void SelectionModeChanged(const GfxNodeP& node) {
    call_method<void, GfxNodeP> (self, "SelectionModeChanged",
        node);
  }

  virtual void NodeTransformed(const GfxNodeP& node) {
    call_method<void, GfxNodeP> (self, "NodeTransformed", node);
  }

  virtual void RenderModeChanged(const GfxNodeP& node) {
    call_method<void, GfxNodeP> (self, "RenderModeChanged", node);
  }

  void NodeAddedDefault(const GfxNodeP&) {
  }

  void NodeRemovedDefault(const GfxNodeP&) {
  }

  void NodeChangedDefault(const GfxNodeP&) {
  }

  void SelectionModeChangedDefault(const GfxNodeP&) {
  }

  void NodeTransformedDefault(const GfxNodeP&) {
  }

  void RenderModeChangedDefault(const GfxNodeP&) {
  }

private:
  PyObject* self;
};


#endif
