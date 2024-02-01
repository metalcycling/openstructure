uniform float level;
uniform float tx;
uniform float ty;
uniform sampler2D slab1;
uniform sampler2D slab2;

void main()
{
  float c0=texture2D(slab1,gl_TexCoord[0].xy).s;
  float c1=texture2D(slab1,gl_TexCoord[0].xy+vec2(0.0,ty)).s;
  float c3=texture2D(slab1,gl_TexCoord[0].xy+vec2(tx,0.0)).s;
  float c4=texture2D(slab2,gl_TexCoord[0].xy).s;
  vec3 result=vec3(0.0,0.0,0.0);

  float fr3=0.0;
  float fr0=0.0;
  float fr4=0.0;
  
  if((level>c0 && level<=c3) || level>c3 && level<=c0) {
    fr3=(level-c0)/(c3-c0);
  }

  if((level>c0 && level<=c1) || level>c1 && level<=c0) {
    fr0=(level-c0)/(c1-c0);
  }

  if((level>c0 && level<=c4) || level>c4 && level<=c0) {
    fr4=(level-c0)/(c4-c0);
  }

  float ww=0.0;

  if(c0<level) ww+=1.0;
  if(c1<level) ww+=2.0;
  if(texture2D(slab1,gl_TexCoord[0].xy+vec2(tx,ty)).s<level) ww+=4.0;
  if(c3<level) ww+=8.0;
  if(c4<level) ww+=16.0;
  if(texture2D(slab2,gl_TexCoord[0].xy+vec2(0.0,ty)).s<level) ww+=32.0;
  if(texture2D(slab2,gl_TexCoord[0].xy+vec2(tx,ty)).s<level) ww+=64.0;
  if(texture2D(slab2,gl_TexCoord[0].xy+vec2(tx,0.0)).s<level) ww+=128.0;

  gl_FragColor.rgba=vec4(fr3,fr0,fr4,ww/255.0);
}

