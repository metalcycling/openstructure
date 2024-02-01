#include "texture.hh"
#include "bitmap_io.hh"

namespace ost { namespace gfx {

  Texture::Texture(const Bitmap& bm):
    w_(bm.width), h_(bm.height),
    d_()
  {
    if(!bm.data) return;
    d_=boost::shared_array<float>(new float[4*w_*h_]);
    static float f=1.0/255.0;
    for(GLint v=0;v<h_;++v) {
      for(GLint u=0;u<w_;++u) {
        int p=v*w_+u;
        float* c = &d_[p*4];
        if(bm.channels==1) {
          c[0]=f*static_cast<float>(bm.data[p]);
          c[1]=c[0];
          c[2]=c[0];
          c[3]=1.0;
        } else if(bm.channels==2) {
          c[0]=f*static_cast<float>(bm.data[p*2+0]);
          c[1]=c[0];
          c[2]=c[0];
          c[3]=f*static_cast<float>(bm.data[p*2+1]);
        } else if(bm.channels==3) {
          c[0]=f*static_cast<float>(bm.data[p*3+0]);
          c[1]=f*static_cast<float>(bm.data[p*3+1]);
          c[2]=f*static_cast<float>(bm.data[p*3+2]);
          c[3]=1.0;
        } else if(bm.channels==4) {
          c[0]=f*static_cast<float>(bm.data[p*4+0]);
          c[1]=f*static_cast<float>(bm.data[p*4+1]);
          c[2]=f*static_cast<float>(bm.data[p*4+2]);
          c[2]=f*static_cast<float>(bm.data[p*4+3]);
        }
      }
    }
  }

}}
