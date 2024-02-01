uniform sampler2D depth_map;
uniform sampler2D norm_map;
uniform sampler1D kernel;
uniform float step;
uniform float size;
uniform vec2 i_vp;
uniform vec4 abcd;
uniform int mode;

// gl_TexCoord[0] comes from quadpp_vs, i.e. from post processing

vec3 unproject(in vec3 rpos)
{ 
  vec4 coord = gl_TextureMatrix[0]*vec4(rpos.xy*2.0-1.0,rpos.z*2.0-1.0,1.0);
  return coord.xyz/coord.w;
}

float ao(in vec2 tab, in float depth_p, in vec3 norm_p, in vec3 t_pos_p)
{
  vec2 rpos = gl_TexCoord[0].xy+tab.xy*i_vp;
  float depth_q = texture2D(depth_map,rpos).r;
  if(mode==0) {
    if(depth_q>=depth_p) {
      return 0.0;
    } else {
      return 1.0;
    }
  } 

  vec3 t_pos_q = unproject(vec3(rpos.xy,depth_q));
  vec3 diff=t_pos_q-t_pos_p;
  float nd = dot(normalize(norm_p),normalize(diff));

  if(mode==1) {
    return 1.6*max(0.0,nd);
  }

  if(mode==2) {
    float fac=1.0+dot(diff,diff)/100.0;
    return 1.6*max(0.0,nd)/fac;
  }

  if(mode==3) {
    // radius of imposter at distance depth_q to cover a single pixel
    vec4 tmpv=gl_TextureMatrix[0]*vec4(-1.0,0.0,depth_q,1.0);
    float tmpd=tmpv.x/tmpv.w;
    tmpv=gl_TextureMatrix[0]*vec4(1.0,0.0,depth_q,1.0);
    // 10.0 is a fudge factor
    float rad=10.0*i_vp[0]*abs(tmpv.x/tmpv.w-tmpd);
    float s = 6.283*(1.0-sqrt(1.0-min(1.0,rad*rad/dot(diff,diff))));
    return s*max(0.0,nd);
  }

  return 0.0;
}

void main()
{
  float depth_p = texture2D(depth_map,gl_TexCoord[0].xy).r;
  if(depth_p>=1.0) {
    gl_FragColor = vec4(0,0,0,1);
    return;
  }
  vec3 norm_p = (texture2D(norm_map,gl_TexCoord[0].xy).xyz-0.5)*2.0;
  vec3 t_pos_p = unproject(vec3(gl_TexCoord[0].xy,depth_p));

  float i;
  float sum=0.0;
  for(i=0.0;i<1.0;i+=step) {
    vec2 nn=(texture1D(kernel,i).xy-0.5)*size;
    sum+=ao(nn,depth_p,norm_p,t_pos_p);
  }
  sum*=step;

  gl_FragColor.rgb=vec3(sum,sum,sum);
  gl_FragColor.a=1.0;
}
