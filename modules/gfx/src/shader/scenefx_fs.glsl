uniform sampler2D scene_map;
uniform sampler2D depth_map;
uniform vec2 i_vp;
uniform bool shadow_flag;
uniform sampler2D shadow_map;
uniform float shadow_weight;
uniform float shadow_depth_bias;
uniform float shadow_epsilon;
uniform float shadow_multiplier;
uniform bool occl_flag;
uniform sampler2D occl_map;
uniform float occl_mult;
uniform bool dark_flag;
uniform sampler2D dark_map;
uniform float dark_mult;
uniform bool fog_flag;
uniform float depth_near;
uniform float depth_far;
uniform float fog_far;
uniform float fog_scale;

// gl_TexCoord[0] comes from scenefx_vs, i.e. from post processing

float CalcShadowFactor(in vec4 coord, in vec2 o)
{
  // get original depth value of line projected towards light
  float d = texture2D(shadow_map, coord.xy+o*shadow_epsilon).x+shadow_depth_bias;
  return d<=coord.z ? shadow_multiplier : 1.0;
}

void main()
{
  float depth = texture2D(depth_map,gl_TexCoord[0].xy).r;
  if(depth>=1.0) {
    discard;
  }
  vec4 scene_color=texture2D(scene_map,gl_TexCoord[0].xy);
  gl_FragColor.a = scene_color.a;

  float shadow_factor=1.0;
  if(shadow_flag) {
    vec4 pcoord = vec4(gl_TexCoord[0].xy*2.0-1.0,depth*2.0-1.0,1.0);
    vec4 coord = gl_TextureMatrix[2]*pcoord;
    coord/=coord.w;
    shadow_factor = 0.0;
    shadow_factor += 0.18*CalcShadowFactor(coord, vec2(-0.7, -0.7));
    shadow_factor += 0.18*CalcShadowFactor(coord, vec2(0.7, -0.7));
    shadow_factor += 0.18*CalcShadowFactor(coord, vec2(0.7, 0.7));
    shadow_factor += 0.18*CalcShadowFactor(coord, vec2(-0.7, 0.7));
    shadow_factor += 0.28*CalcShadowFactor(coord, vec2(0, 0));
    shadow_factor = mix(1.0, shadow_factor, shadow_weight);
  }

  float occl_factor=1.0;
  if(occl_flag) {
    occl_factor=max(0.0,1.0-occl_mult*texture2D(occl_map,gl_TexCoord[0].xy).r);
  }

  float dark_factor=1.0;
  if(dark_flag) {
    dark_factor=max(0.0,1.0-dark_mult*(1.0-texture2D(dark_map,gl_TexCoord[0].xy).r));
  }

  float fog=1.0;
  if(fog_flag) {
    float z = 2.0*depth_near*depth_far/(-(depth*2.0-1.0)*(depth_far-depth_near)+depth_far+depth_near);
    fog = clamp((fog_far-z) * fog_scale, 0.0, 1.0);
  }

  gl_FragColor.rgb = mix(gl_Fog.color.rgb, 
                         shadow_factor*occl_factor*dark_factor*scene_color.rgb,
                         fog);
}

