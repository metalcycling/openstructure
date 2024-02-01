uniform bool lighting_flag;
uniform bool two_sided_flag;

void DirectionalLight(in vec3 normal,
                      in float shin,
                      inout vec4 ambient,
                      inout vec4 diffuse,
                      inout vec4 specular)
{
  float n_vp = max(0.0, dot(normal, normalize(gl_LightSource[0].position.xyz)));
  float pf = 0.0;
  if(n_vp>0.0 && shin>0.0) {
    float n_hv = max(0.0, dot(normal, normalize(gl_LightSource[0].halfVector.xyz)));
    pf=pow(n_hv, shin);
  }

  ambient  += gl_LightSource[0].ambient;
  diffuse  += gl_LightSource[0].diffuse * n_vp;
  specular += gl_LightSource[0].specular * pf;
}

void CalcFrontAndBackColor(in vec3 normal)
{
  vec4 amb = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);

  DirectionalLight(normal, gl_FrontMaterial.shininess, amb, diff, spec);

  gl_FrontColor = gl_FrontLightModelProduct.sceneColor + 
                  (amb  * gl_FrontMaterial.ambient * gl_Color) +
                  (diff * gl_FrontMaterial.diffuse * gl_Color) +
                  (spec * gl_FrontMaterial.specular);

  if(two_sided_flag) {
    amb=vec4(0.0);
    diff=vec4(0.0);
    spec=vec4(0.0);
    
    DirectionalLight(-normal, gl_BackMaterial.shininess, amb, diff, spec);
    
    gl_BackColor = gl_BackLightModelProduct.sceneColor + 
                   (amb  * gl_BackMaterial.ambient * gl_Color) +
                   (diff * gl_BackMaterial.diffuse * gl_Color) +
                   (spec * gl_BackMaterial.specular);
  } else {
    gl_BackColor = gl_FrontColor;
  }

}

void main()
{    
  // transformed position
  gl_Position = ftransform();
  vec4 ec_Pos = gl_ModelViewMatrix* gl_Vertex;

  vec3 normal = normalize(gl_NormalMatrix * gl_Normal);
  gl_TexCoord[2].stp = normal;
  
  if(lighting_flag) {
    CalcFrontAndBackColor(normal);
  } else {
    gl_FrontColor = gl_Color;
    gl_BackColor = gl_Color;
  }

  // for some reason, the fog and z coordinate are sign toggled...
  gl_FogFragCoord = -ec_Pos.z;
}

