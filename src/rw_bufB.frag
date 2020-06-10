// random walk main.
// とりあえず初期配置は関数で。
// たとえば、hueの値を動かすとか。

// constants.
const float pi = 3.14159;
const float gridsize = 10.0;
const float Lfactor = 0.1; // 1.0 / length.
// random.
const vec2 r_vector = vec2(12.9898, 78.233);
const vec3 r3_vector = vec3(124.98, 782.33, 415.69);
const float r_coeff = 43758.5453123;
float random(vec2 st){
  return fract(sin(dot(st.xy, r_vector)) * r_coeff);
}
float random(vec3 st){
  return fract(sin(dot(st, r3_vector)) * r_coeff);
}
// hsb to rgb.
vec3 getRGB(float r, float g, float b){
  vec3 c = vec3(r, g, b);
  vec3 rgb = clamp(abs(mod(c.x * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);
    rgb = rgb * rgb * (3.0 - 2.0 * rgb);
    return c.z * mix(vec3(1.0), rgb, c.y);
}
// マスの値の取得。
// 0.5ずつずらしてセルの真ん中の色を取得する。
vec4 getColor(vec2 i){
  vec2 pos = (i + vec2(0.5, 0.5)) * gridsize * vec2(1.0 / iResolution.x, 1.0 / iResolution.y);
  if(min(pos.x, pos.y) < 0.0 || max(pos.x, pos.y) > 1.0){ return vec4(-1.0); }
  pos.y = pos.y;
  return texture(iChannel0, pos);
}
// 周囲に濃いマスが1つ以下。
bool nonConflict(vec2 i){
  vec2 dir = vec2(1.0, 0.0);
  vec4 neighborColor;
  int check = 0;
  for(int s = 0; s < 4; s++){
    neighborColor = getColor(i + dir);
    dir = dir.yx * vec2(-1.0, 1.0);
    if(neighborColor.a > 1.0 - (0.5 * Lfactor)){ check++; }
    if(check > 1){ return false; }
  }
  return true;
}
// comingWalker.
// 難しい。
vec4 comingWalker(vec2 i){
// まず、周囲に濃いマスが2つ以上ならだめー
  if(!nonConflict(i)){ return vec4(0.0); }
  vec2 dir[4];
  dir[0] = vec2(1.0, 0.0);
  dir[1] = vec2(0.0, 1.0);
  dir[2] = vec2(-1.0, 0.0);
  dir[3] = vec2(0.0, -1.0);
  float k = -1.0;
  vec2 neighbor;
  vec4 neighborColor;
// 周囲に濃いマスがないならblancのまま。
  for(int s = 0; s < 4; s++){
    neighbor = i + dir[s];
    neighborColor = getColor(neighbor);
// 複数の場合はこういうマスが2つ以上無いように気を付ける
    if(neighborColor.a > 1.0 - (0.5 * Lfactor)){ k = float(s); break; }
  }
// 周囲に濃いマスがない時
  if(k < 0.0){ return vec4(0.0); }
  k = mod(k + 2.0, 4.0);
  vec2 around;
  vec4 aroundColor;
  float n = 0.0;
  float tmp[4];
  for(int s = 0; s < 4; s++){ tmp[s] = 4.0; }
  for(int s = 0; s < 4; s++){
    around = neighbor + dir[s];
    aroundColor = getColor(around);
// 複数の場合はブランクなだけでは駄目で、
// なおかつ周囲の濃いマスが1つしかないことが必要。
    if(aroundColor.a > -0.5 && aroundColor.a < 0.5 * Lfactor && nonConflict(around)){ n += 1.0; }
    else{ tmp[s] = float(s); }
  }
  float walkerDir = floor(random(vec3(neighbor, iTime)) * n);
  for(int s = 0; s < 4; s++){
    if(walkerDir >= tmp[s]){ walkerDir += 1.0; }
  }
  if(walkerDir == k){ return neighborColor; }
  return vec4(0.0);
}
// 初期化
vec4 initialize(vec2 p){
  vec2 i = floor(p);
  for(float x = 2.0; x < iResolution.x; x += 6.0){
    if(x * gridsize >= iResolution.x){ break; }
    for(float y = 2.0; y < iResolution.y; y += 6.0){
      if(y * gridsize >= iResolution.y){ break; }
      if(i.x == x && i.y == y){
        float rnd1 = random(vec3(x, y, random(vec2(x, y))));
        float rnd2 = random(vec3(-y, x, random(vec2(-y, x))));
        return vec4(getRGB(0.55 + 0.25 * rnd1, 1.0 - 0.5 * rnd2, 1.0), 1.0);
      }
    }
  }
  return vec4(0.0);
}
// 移動不可能（周囲にブランクマスがない、
// さらにブランクでもその周囲に濃いマスが1つ以下）
bool immovable(vec2 i){
  vec2 dir = vec2(1.0, 0.0);
  vec2 neighbor;
  vec4 neighborColor;
  for(int s = 0; s < 4; s++){
    neighbor = i + dir;
    neighborColor = getColor(neighbor);
    dir = dir.yx * vec2(-1.0, 1.0);
    if(neighborColor.a > -0.5 && neighborColor.a < 0.5 * Lfactor && nonConflict(neighbor)){ return false; }
  }
  return true;
}
void mainImage(out vec4 fragColor, in vec2 fragCoord){
  float minLength = min(iResolution.x, iResolution.y);
  vec2 p = fragCoord.xy / minLength;
  float scale = minLength / gridsize;
  p *= scale;
  vec4 col;
  // 初期化
  if(floor(iTime * 60.0) == 1.0){
    col = initialize(p);
    fragColor = col;
    return;
  }
  vec2 i = floor(p);
  vec4 cur = getColor(i);
  // curはpが属するセルの色である。
  // 変化がなければそのままcurを返す。
// 以下、cur.aの値で処理を分ける。
  if(cur.a > 1.0 - (0.5 * Lfactor)){
     // 移動できないときはremainする
    if(immovable(i)){ col = cur; }else{ col = cur - vec4(vec3(0.0), Lfactor); }
  }else if(cur.a > 1.5 * Lfactor){
    col = cur - vec4(vec3(0.0), Lfactor);
  }else if(cur.a > 0.5 * Lfactor){
    col = vec4(0.0);
  }else{
// iのマスの周囲にこちらにやってくるwalkerがいればその色を設定、
// いなければvec4(0.0)を設定（つまりblanc）
    col = comingWalker(i);
  }
  fragColor = col;
}
