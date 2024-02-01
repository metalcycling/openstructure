uniform sampler2D left_scene;
uniform sampler2D right_scene;

void main()
{
  vec4 left_color=texture2D(left_scene,gl_TexCoord[0].xy);
  vec4 right_color=texture2D(right_scene,gl_TexCoord[0].xy);
  gl_FragColor.a = left_color.a+right_color.a;
  gl_FragColor.r=pow(0.7*left_color.g+0.3*left_color.b,1.5);
  gl_FragColor.gb=right_color.gb;
}

