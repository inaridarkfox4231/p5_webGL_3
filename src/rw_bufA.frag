// background.

// for noise.
//const vec2 u_10 = vec2(1.0, 0.0);
//const vec2 u_01 = vec2(0.0, 1.0);
//const vec2 u_11 = vec2(1.0, 1.0);
//const vec2 r_vector = vec2(12.9898, 78.233);
const vec3 r_vec_30 = vec3(127.1, 311.7, 251.9);
const vec3 r_vec_31 = vec3(269.5, 183.3, 314.3);
const vec3 r_vec_32 = vec3(419.2, 371.9, 218.4);
const vec3 u_100 = vec3(1.0, 0.0, 0.0);
const vec3 u_010 = vec3(0.0, 1.0, 0.0);
const vec3 u_001 = vec3(0.0, 0.0, 1.0);
const vec3 u_110 = vec3(1.0, 1.0, 0.0);
const vec3 u_101 = vec3(1.0, 0.0, 1.0);
const vec3 u_011 = vec3(0.0, 1.0, 1.0);
const vec3 u_111 = vec3(1.0, 1.0, 1.0);
const float r_coeff = 43758.5453123;
const int octaves = 6;
// hsb to rgb.
vec3 getRGB(float r, float g, float b){
  vec3 c = vec3(r, g, b);
  vec3 rgb = clamp(abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);
    rgb = rgb * rgb * (3.0 - 2.0 * rgb);
    return c.z * mix(vec3(1.0), rgb, c.y);
}
// 3D random vector.(-1.0～1.0)
vec3 random3(vec3 st){
  vec3 v;
  v.x = sin(dot(st, r_vec_30)) * r_coeff;
  v.y = sin(dot(st, r_vec_31)) * r_coeff;
  v.z = sin(dot(st, r_vec_32)) * r_coeff;
  return -1.0 + 2.0 * fract(v); // -1.0～1.0に正規化
}
// simplex noise.
float snoise3(vec3 st){
  vec3 p = st + (st.x + st.y + st.z) / 3.0;
  vec3 f = fract(p);
  vec3 i = floor(p);
  vec3 g0, g1, g2, g3;
  vec4 wt;
  g0 = i;
  g3 = i + u_111;
  if(f.x >= f.y && f.x >= f.z){
    g1 = i + u_100;
    g2 = i + (f.y >= f.z ? u_110 : u_101);
    wt = (f.y >= f.z ? vec4(1.0 - f.x, f.x - f.y, f.y - f.z, f.z) : vec4(1.0 - f.x, f.x - f.z, f.z - f.y, f.y));
  }else if(f.y >= f.x && f.y >= f.z){
    g1 = i + u_010;
    g2 = i + (f.x >= f.z ? u_110 : u_011);
    wt = (f.x >= f.z ? vec4(1.0 - f.y, f.y - f.x, f.x - f.z, f.z) : vec4(1.0 - f.y, f.y - f.z, f.z - f.x, f.x));
  }else{
    g1 = i + u_001;
    g2 = i + (f.x >= f.y ? u_101 : u_011);
    wt = (f.x >= f.y ? vec4(1.0 - f.z, f.z - f.x, f.x - f.y, f.y) : vec4(1.0 - f.z, f.z - f.y, f.y - f.x, f.x));
  }
  float value = 0.0;
  wt = wt * wt * wt * (wt * (wt * 6.0 - 15.0) + 10.0);
  value += wt.x * dot(p - g0, random3(g0));
  value += wt.y * dot(p - g1, random3(g1));
  value += wt.z * dot(p - g2, random3(g2));
  value += wt.w * dot(p - g3, random3(g3));
  return value;
}
// fbm.
float fbm(vec3 st){
  float value = 0.0;
  float amplitude = 0.5;
  for(int i = 0; i < octaves; i++){
    value += amplitude * snoise3(st);
    st *= 2.0;
    amplitude *= 0.5;
  }
  return value;
}
void mainImage(out vec4 fragColor, in vec2 fragCoord){
  float minLength = min(iResolution.x, iResolution.y);
  vec2 p = fragCoord.xy / minLength;
  p = (p + vec2(0.1, 0.18) * iTime) * 3.0;
  float n = 0.5 + 0.5 * fbm(vec3(p, iTime * 0.2)); // ノイズ計算
  vec3 fireColor = getRGB((n - 0.46) * 0.78, 1.0, 1.0);
  vec3 skyColor = vec3(0.0);
  vec3 col = skyColor + (fireColor - skyColor) * smoothstep(0.44, 0.56, n);
  fragColor = vec4(col, 1.0);
}
