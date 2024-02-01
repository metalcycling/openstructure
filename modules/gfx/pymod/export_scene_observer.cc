#include <boost/python.hpp>
using namespace boost::python;

#include <ost/gfx/scene_observer.hh>
#include <ost/gfx/gfx_node.hh>

#include "scene_observer_proxy.hh"

using namespace ost;
using namespace ost::gfx;

void export_SceneObserver() {
  class_<SceneObserver, boost::shared_ptr<SceneObserverProxy>,
       SceneObserverProxy>("SceneObserver",init<>())
    .def("NodeAdded", &SceneObserverProxy::NodeAddedDefault)
    .def("NodeRemoved", &SceneObserverProxy::NodeRemovedDefault)
    .def("NodeChanged",	&SceneObserverProxy::NodeChangedDefault)
    .def("SelectionModeChanged",	&SceneObserverProxy::SelectionModeChangedDefault)
    .def("NodeTransformed",	&SceneObserverProxy::NodeTransformedDefault)
    .def("RenderModeChanged",	&SceneObserverProxy::RenderModeChangedDefault)
;
}
