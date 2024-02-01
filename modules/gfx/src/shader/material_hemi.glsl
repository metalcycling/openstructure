uniform vec4 hemi_param;

vec4 light(in vec3 normal,
           in vec4 color,
           in bool two_sided)
{
  vec4 ground_color=hemi_param[0]*color;
  vec4 sky_color=hemi_param[1]*color;
  
  vec3 l_dir  = normalize(gl_LightSource[0].position.xyz); // assume directional vector
  // min(1.0,max(0.0,x-u)/(v-u))
  float a=min(1.0,max(0.0,dot(normal,l_dir)-hemi_param[2])/(hemi_param[3]-hemi_param[2]));
  //float a  = 0.5+0.5*dot(normal,l_dir);
  color.rgb = mix(ground_color, sky_color, a).rgb;
  bool backside = false;
  if(two_sided && backside) {
    float a  = 0.5+0.5*dot(-normal,l_dir);
    color.rgb = mix(ground_color, sky_color, a).rgb;
  }
  return color;
}

