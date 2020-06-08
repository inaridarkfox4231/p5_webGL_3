const float pi = 3.14159;
const vec2 r_vector = vec2(12.9898, 78.233);
const float r_coeff = 43758.5453123;
float random(vec2 st){
  return fract(sin(dot(st.xy, r_vector)) * r_coeff);
}
// palette.
const vec3 red = vec3(0.95, 0.3, 0.35);
const vec3 green = vec3(0.3, 0.9, 0.4);
const vec3 blue = vec3(0.2, 0.25, 0.98);
const vec3 gold = vec3(0.85, 0.67, 0.14);
const vec3 white = vec3(1.0);
// hsb to rgb.
vec3 getRGB(float r, float g, float b){
    vec3 c = vec3(r, g, b);
    vec3 rgb = clamp(abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);
    rgb = rgb * rgb * (3.0 - 2.0 * rgb);
    return c.z * mix(vec3(1.0), rgb, c.y);
}
// 菱形タイリング
vec3 diaTiling(in vec2 p, float scale, vec3 tileColor){
  p.x += iTime * 0.5;
  p *= scale;
  vec2 q = (mat2(1.0, 1.0, -sqrt(3.0), sqrt(3.0)) / 3.0) * p;
  vec2 i = floor(q);
  vec2 f = fract(q);
  vec2 e1 = vec2(sqrt(3.0) * 0.5, -0.5);
  vec2 e2 = vec2(sqrt(3.0) * 0.5, 0.5);
  vec3 col = tileColor;
  if(mod(i.x * i.y, 2.0) == 0.0){ col = white; }
  return col;
}
void mainImage(in vec2 fragCoord, out vec4 fragColor){
  vec2 p = (fragCoord.xy * 0.5 - iResolution.xy) / min(iResolution.x, iResolution.y);
  vec3 col_1 = diaTiling(p, 8.0, red);
  float t = pi * 2.0 / 3.0;
  vec3 col_2 = diaTiling((p - vec2(0.0, sqrt(3.0) / 8.0)) * mat2(cos(t), sin(t), -sin(t), cos(t)), 8.0, blue);
  vec3 col_3 = diaTiling((p - vec2(-1.5, 0.5 * sqrt(3.0)) / 8.0) * mat2(cos(t), -sin(t), sin(t), cos(t)), 8.0, green);
  vec3 col = min(col_1, min(col_2, col_3));
  if(col == white){ col = mix(gold, white, (1.4 - length(p)) / 1.4); }
  fragColor = vec4(col, 1.0);
}
