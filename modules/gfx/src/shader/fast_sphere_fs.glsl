uniform bool lighting_flag;
uniform bool two_sided_flag;
uniform bool fog_flag;
uniform bool write_normals;

// gl_TexCoord[0] is from gl_MultiTexCoord0, which in turn
// is custom crafted in the fast sphere prep routine

void main()
{
  float xx = (gl_TexCoord[0].x-0.5)*2.0;
  float yy = (gl_TexCoord[0].y-0.5)*2.0;
  float xxyy = xx*xx+yy*yy;
  if(xxyy>1.0) discard;
  
  float zz = sqrt(1.0-clamp(xx*xx+yy*yy,0.0,1.0));
  vec3 normal=vec3(xx,yy,zz);

  float z1 = gl_TexCoord[0].z*zz;
  float z2 = gl_TexCoord[0].w*zz;

  gl_FragDepth = gl_FragCoord.z+z2;

  if(write_normals) {
    gl_FragColor.rgb=normal*0.5+0.5;
    return;
  }

  vec4 color = gl_Color;
  if(lighting_flag) {
    color = light(normal,color,two_sided_flag);
  }

  if(fog_flag) {
    float fog = clamp((gl_Fog.end-(gl_FogFragCoord+z1)) * gl_Fog.scale, 0.0, 1.0);
    color.rgb = mix(gl_Fog.color.rgb, color.rgb, fog);
  }
  gl_FragColor.rgb = color.rgb;
  gl_FragColor.a = 1.0;
}
