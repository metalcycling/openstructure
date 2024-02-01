uniform float scalef;
uniform vec4 color;

void main()
{    
  vec4 n = vec4(gl_Normal,0.0);
  gl_Position = gl_ModelViewProjectionMatrix*(gl_Vertex+scalef*n);
  gl_FrontColor = color;
  gl_BackColor = color;
  vec4 ec_Pos = gl_ModelViewMatrix* gl_Vertex;
  gl_FogFragCoord = -ec_Pos.z;
}

