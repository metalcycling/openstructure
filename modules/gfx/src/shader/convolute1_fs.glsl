uniform sampler2D data;
uniform sampler1D kernel;
uniform float step;
uniform vec2 i_vp;

// gl_TexCoord[0] comes from quadpp_vs, i.e. from post processing

void main()
{
  float i;
  float sum=0.0;
  float sum2=0.0;
  for (i=0.0;i<1.0;i=i+step) {
    vec3 kvalue = texture1D(kernel,i).xyz;
    kvalue.xy=(kvalue.xy-0.5)*10.0*i_vp;
    sum+=kvalue.z*texture2D(data,gl_TexCoord[0].xy+kvalue.xy).r;
    sum2+=kvalue.z;
  }
  sum=sum/sum2;
  float delta=5.0*(sum-texture2D(data,gl_TexCoord[0].xy).r);
  if(delta<0.0) {
    gl_FragColor.rgb=vec3(1.0+delta,1.0+delta,1.0+delta);
  } else {
    gl_FragColor.rgb=vec3(1.0,1.0,1.0);
  }
  gl_FragColor.a=1.0;
}
