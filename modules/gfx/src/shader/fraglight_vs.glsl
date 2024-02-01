void main()
{
  // transformed position
  gl_Position = ftransform();

  // eye coordinate position
  vec4 ec_Pos = gl_ModelViewMatrix* gl_Vertex;

  // for some reason, the fog and z coordinate are sign toggled...
  gl_FogFragCoord = -ec_Pos.z;

  // store fragment normal in tex coord 2
  vec3 normal = vec3(0,0,1);
  if(dot(gl_Normal,gl_Normal)>0.001) {
    normal = normalize(gl_NormalMatrix * gl_Normal);
  }
  gl_TexCoord[2].stp=normal;

  // default tex coord
  gl_TexCoord[0] = gl_MultiTexCoord0;

  gl_FrontColor=gl_Color;
  gl_BackColor=gl_Color;
}

