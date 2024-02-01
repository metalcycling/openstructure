uniform vec4 sel_color;

void main()
{    
  gl_Position = ftransform();
  vec3 ec_pos = vec3(gl_ModelViewMatrix* gl_Vertex);
  gl_FrontColor = sel_color;
  gl_BackColor = sel_color;
  gl_FogFragCoord = -ec_pos.z;
}

