let myShader;

let vs =
"precision mediump float;" +
"attribute vec3 aPosition;" +
"void main(void){" +
"  gl_Position = vec4(aPosition, 1.0);" +
"}";

let fs =
"precision mediump float;" +
"uniform vec2 u_resolution;" +
"uniform vec2 u_mouse;" +
"uniform float u_time;" +
"const vec2 r_vector = vec2(12.9898, 78.233);" +
"const float r_coeff = 43758.5453123;" +
"float random(vec2 st){" +
"  return fract(sin(dot(st.xy, r_vector)) * r_coeff);" +
"}" +
// 横長の菱形
"float dia1(vec2 p){" +
"  return max(abs(p.x + sqrt(3.0) * p.y), abs(p.x - sqrt(3.0) * p.y)) * 2.0 / sqrt(3.0) - 1.0;" +
"}" +
// 右下方向に傾いた菱形
"float dia2(vec2 p){" +
"  return max(abs(4.0 * p.x), abs(2.0 * p.x + 2.0 * sqrt(3.0) * p.y)) / sqrt(3.0) - 1.0;" +
"}" +
// 左下方向に傾いた菱形
"float dia3(vec2 p){" +
"  return max(abs(4.0 * p.x), abs(2.0 * p.x - 2.0 * sqrt(3.0) * p.y)) / sqrt(3.0) - 1.0;" +
"}" +
// hsb to rgb.
"vec3 getRGB(float r, float g, float b){" +
"    vec3 c = vec3(r, g, b);" +
"    vec3 rgb = clamp(abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);" +
"    rgb = rgb * rgb * (3.0 - 2.0 * rgb);" +
"    return c.z * mix(vec3(1.0), rgb, c.y);" +
"}" +
"void main(void){" +
"  vec2 p = (gl_FragCoord.xy - u_resolution) / min(u_resolution.x, u_resolution.y);" +
"  vec3 col = getRGB(0.65, (p.y + 1.0) * 0.5, 1.0);" +
"  gl_FragColor = vec4(col, 1.0);" +
"}";

function setup(){
  createCanvas(400, 400, WEBGL);
  myShader = createShader(vs, fs);
  shader(myShader);
  noLoop();
}

function draw(){
  myShader.setUniform("u_resolution", [width, height]);
  myShader.setUniform("u_mouse", [mouseX, mouseY]);
  myShader.setUniform("u_time", millis() / 1000);
  quad(-1, -1, -1, 1, 1, 1, 1, -1);
}
