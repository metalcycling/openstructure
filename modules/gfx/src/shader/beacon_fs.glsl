uniform vec3 pos;
uniform vec3 dir;
uniform float len;
uniform float rad;
uniform sampler2D depth_map;

// gl_TexCoord[0] comes from scenefx_vs, i.e. from post processing

void main()
{
  float depth = texture2D(depth_map,gl_TexCoord[0].xy).r;
  if(depth>=1.0) {
    discard;
  }
  vec4 pcoord = vec4(gl_TexCoord[0].xy*2.0-1.0,depth*2.0-1.0,1.0);
  vec4 coord = gl_TextureMatrix[1]*pcoord;
  coord/=coord.w;

  float d = length(cross(dir,coord.xyz-pos))/len;
  if(d>rad) {
    discard;
  }
  gl_FragColor.rgb=vec3(0,1.0,0);
  gl_FragColor.a=min(0.9,mix(3.0,0.0,d/rad));
}
