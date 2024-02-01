uniform sampler2D scene_map;
uniform sampler2D depth_map;
uniform vec2 i_vp;
uniform float depth_cutoff;

void main()
{
  // reference depth
  float depth = texture2D(depth_map,gl_TexCoord[0].xy).r;

  vec4 scene_color=vec4(0.0,0.0,0.0,0.0);

  int count=0;
  for(float x=-1.5;x<1.51;x+=1.0) {
    for(float y=-1.5;y<1.51;y+=1.0) {
      vec2 p = gl_TexCoord[0].xy+vec2(x*i_vp[0],y*i_vp[1]);
      float dd = abs(depth-texture2D(depth_map,p).r);
      if(dd<depth_cutoff) {
        ++count;
        scene_color+=texture2D(scene_map,p);
      }
    }
  }

  if (count==0) {
    scene_color=texture2D(scene_map,gl_TexCoord[0].xy);
  } else {
    scene_color/=float(count);
  }

  gl_FragColor=scene_color;

}