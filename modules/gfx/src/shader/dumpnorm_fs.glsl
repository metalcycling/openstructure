uniform float clip_offset;

void main()
{
  if(gl_FragCoord.z<clip_offset) {
    discard;
  }
  gl_FragColor.rgb=normalize(gl_TexCoord[0]).stp*0.5+0.5;
}

