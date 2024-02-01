uniform bool lighting_flag;

// noop vertex shader

void main()
{    
  gl_Position = gl_Vertex;
  gl_BackColor = gl_Color;
  gl_FrontColor = gl_Color;
}

