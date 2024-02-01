void main()
{
  gl_Position = ftransform();
  vec3 normal = vec3(0,0,1);
  if(dot(gl_Normal,gl_Normal)>0.001) {
    normal = gl_NormalMatrix * gl_Normal;
  }
  gl_TexCoord[0].stp=normal;
}

