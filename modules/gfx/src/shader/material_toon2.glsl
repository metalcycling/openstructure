bool DirectionalLight(in vec3 normal,
                      in float shin,
                      inout vec4 ambient,
                      inout vec4 diffuse,
                      inout vec4 spec)
{
  float n_vp = max(0.0, dot(normal, normalize(vec3(gl_LightSource[0].position))));
  if (n_vp>0.6) {
    n_vp=0.9;
  } else if (n_vp>0.5) {
    n_vp=0.7;
  } else if (n_vp>0.05) {
    n_vp=0.1;
  } else {
    n_vp=0.0;
  }

  float n_hv = max(0.0, dot(normal, vec3(gl_LightSource[0].halfVector)));
  float pf = shin>0.0 && n_vp>0.0 ? (pow(n_hv,shin)>0.8 ? 1.0 : 0.0) : 0.0;

  ambient  += gl_LightSource[0].ambient;
  diffuse  += gl_LightSource[0].diffuse*n_vp;
  spec     += gl_LightSource[0].specular*pf;
  return true;
}

vec4 light(in vec3 normal,
           in vec4 color,
           in bool two_sided)
{
  vec4 amb = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);
  vec4 color2 = vec4(0.0);

  DirectionalLight(normal, gl_FrontMaterial.shininess, amb, diff, spec);

  color2 = (gl_FrontLightModelProduct.sceneColor + 
           amb*color + diff*color + spec*color);

  if(two_sided) {
    amb=vec4(0.0);
    diff=vec4(0.0);
    spec=vec4(0.0);
    DirectionalLight(-normal, gl_FrontMaterial.shininess, amb, diff, spec);
    color2 += (gl_BackLightModelProduct.sceneColor + 
              amb*color + diff*color + spec*color);
  }

  return clamp(color2,0.0,1.0);
}
