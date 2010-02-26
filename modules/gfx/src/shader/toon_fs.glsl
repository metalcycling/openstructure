uniform bool lighting_flag;
uniform bool two_sided_flag;
uniform bool fog_flag;

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

void main()
{
  vec4 amb = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 color = vec4(0.0);
  vec3 normal = normalize(gl_TexCoord[0].stp);

  DirectionalLight(normal, amb, diff);

  color = (gl_FrontLightModelProduct.sceneColor + 
           amb*gl_Color + diff*gl_Color);

  if(two_sided_flag) {
    amb=vec4(0.0);
    diff=vec4(0.0);
    DirectionalLight(-normal, amb, diff);
    color += (gl_BackLightModelProduct.sceneColor + 
              amb*gl_Color + diff*gl_Color);
  }

  gl_FragColor = clamp(color,0.0,1.0);
  gl_FragColor.a = gl_Color.a;
}
