uniform bool fog_flag;
uniform float clip_offset;

void main()
{
  if(gl_FragCoord.z<clip_offset) {
    discard;
  }

  float fog = fog_flag ? clamp((gl_Fog.end-gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0) : 1.0;
  gl_FragColor.rgb = mix(gl_Fog.color.rgb, gl_Color.rgb, fog);
  gl_FragColor.a = gl_Color.a;
}

