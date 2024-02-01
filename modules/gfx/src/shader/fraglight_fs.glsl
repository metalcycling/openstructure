uniform bool lighting_flag;
uniform bool two_sided_flag;
uniform bool fog_flag;
uniform sampler2D depth_map;
uniform int depth_mode;
uniform bool tex_flag;
uniform sampler2D tex_map;
uniform float clip_offset;

void main()
{
  if(gl_FragCoord.z<clip_offset) {
    discard;
  }
  vec4 color = gl_Color;

  if(tex_flag) {
    color.rgb = texture2D(tex_map,gl_TexCoord[0].st).rgb;
  }

  if(lighting_flag) {
    vec3 normal = normalize(gl_TexCoord[2].stp);
    color=light(normal,color,two_sided_flag);
  }

  gl_FragColor = color;

  float fog = fog_flag ? clamp((gl_Fog.end-gl_FogFragCoord) * gl_Fog.scale, 0.0, 1.0) : 1.0;
  gl_FragColor.rgb = mix(gl_Fog.color.rgb, gl_FragColor.rgb, fog);
  gl_FragColor.a = gl_Color.a;
}
