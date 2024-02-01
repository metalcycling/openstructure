#ifndef OST_GFX_GLWIN_BASE_PROXY_HH
#define OST_GFX_GLWIN_BASE_PROXY_HH

#include <boost/python.hpp>
#include <boost/python/register_ptr_to_python.hpp>

using namespace boost::python;

#include <ost/gfx/glwin_base.hh>

using namespace ost;
using namespace ost::gfx;

namespace {

class GLWinBaseProxy: public GLWinBase {
public:
  GLWinBaseProxy(PyObject *p) :
    self(p) {
  }

  virtual void DoRefresh() {
    call_method<void>(self, "DoRefresh");
  }
  virtual bool HasStereo() const { 
    return call_method<bool>(self,"HasStereo");
  }
  virtual void StatusMessage(const String& m) {
    call_method<void, const String>(self, "StatusMessage", m);
  }
  virtual bool HasMultisample() const {
    return call_method<bool>(self,"HasMultisample");
  }

private:
  PyObject* self;
};

}

#endif
