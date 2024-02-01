/*
  Adapted after Chan and Durand, "Fast Prefiltered Lines",
  in GPU Gems 2
*/
uniform bool fog_flag;
uniform vec3 edge0;
uniform vec3 edge1;
uniform vec3 edge2;
uniform vec3 edge3;
uniform sampler1D table;

void main()
{
  vec3 pos = vec3(gl_FragCoord.x,gl_FragCoord.y,1.0);
  vec4 d = vec4(dot(edge0,pos), dot(edge1,pos),dot(edge2,pos), dot(edge3,pos));
  if(any(lessThan(d,vec4(0.0,0.0,0.0,0.0)))) discard;

  if(fog_flag) {
    float fog = clamp((gl_Fog.end-gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0);
    gl_FragColor.rgb = mix(gl_Fog.color.rgb, gl_Color.rgb, fog);
  } else {
    gl_FragColor.rgb = gl_Color.rgb;
  }
  gl_FragColor.a = texture1D(table,min(d.x,d.z)).r*texture1D(table,min(d.y,d.w)).r;
}
