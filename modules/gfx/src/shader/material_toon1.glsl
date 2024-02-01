void DirectionalLight(in vec3 normal,
                      inout vec4 ambient,
                      inout vec4 diffuse)
{
  float n_vp=max(0.0, dot(normal, vec3(0.0, 0.0, 1.0)));
  if (n_vp>0.6) {
    n_vp=1.0;
  } else if (n_vp>0.5) {
    n_vp=0.8;
  } else if (n_vp>0.1) {
    n_vp=0.1;
  } else {
    n_vp=0.0;
  }
  ambient  += gl_LightSource[0].ambient;
  diffuse  += gl_LightSource[0].diffuse*n_vp;
}

vec4 light(in vec3 normal,
           in vec4 color,
           in bool two_sided)
{
  vec4 amb = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 color2 = vec4(0.0);

  DirectionalLight(normal, amb, diff);

  color2 = (gl_FrontLightModelProduct.sceneColor + 
           amb*color + diff*color);

  if(two_sided) {
    amb=vec4(0.0);
    diff=vec4(0.0);
    DirectionalLight(-normal, amb, diff);
    color2 += (gl_BackLightModelProduct.sceneColor + 
              amb*color + diff*color);
  }

  return clamp(color2,0.0,1.0);
}
