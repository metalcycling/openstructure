void main()
{
  gl_FragColor.rg=fract(gl_TexCoord[0].st);
  gl_FragColor.b=1.0;
  gl_FragColor.a=1.0;
}
