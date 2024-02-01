bool DirectionalLight(in vec3 normal,
                      in int lightn,
                      in float shin,
                      inout vec4 ambient,
                      inout vec4 diffuse,
                      inout vec4 specular,
                      inout bool lflag)
{
  float n_vp = max(0.0, dot(normal, normalize(gl_LightSource[lightn].position.xyz)));

  float pf = 0.0;
  lflag = n_vp>0.0;
  if(n_vp>0.0 && shin>0.0) {
    float n_hv = max(0.0, dot(normal, normalize(gl_LightSource[lightn].halfVector.xyz)));
    pf=pow(n_hv, shin);
  }

  ambient  += gl_LightSource[lightn].ambient;
  diffuse  += gl_LightSource[lightn].diffuse * n_vp;
  specular += gl_LightSource[lightn].specular * pf;

  return true;
}


vec4 light(in vec3 normal,
           in vec4 color,
           in bool two_sided)
{
  vec4 amb = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);
  bool lflag=false;

  DirectionalLight(normal, 0, gl_FrontMaterial.shininess, amb, diff, spec, lflag);

  color  = gl_FrontLightModelProduct.sceneColor  +
           (amb  * gl_FrontMaterial.ambient * color) +
           (diff * gl_FrontMaterial.diffuse * color) +
           (spec * gl_FrontMaterial.specular);

/* TODO
    DirectionalLight(-normal, 0, gl_BackMaterial.shininess, amb, diff, spec, dummy);
    color = gl_BackLightModelProduct.sceneColor  +
            (amb  * gl_BackMaterial.ambient * color) +
            (diff * gl_BackMaterial.diffuse * color) +
            (spec * gl_BackMaterial.specular);
*/

  return color;
}
