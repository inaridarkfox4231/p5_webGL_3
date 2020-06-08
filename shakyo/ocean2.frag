// また別のやつ。詳しい。
//https://www.shadertoy.com/view/lt3GWj
// うまくいかなさそうだけど・・ていうかみんなAleksanderさんのあれ参考にしてるのね。

// CAUTION!!!: .jsにして読みやすくしてあるけど実行はできないよ。


// A documented, altered, recolored version of "Seascape".
// The famous original at:
// https://www.shadertoy.com/view/Ms2SD1

// "Seascape" by Alexander Alekseev aka TDM - 2014
// Commenting added by bteitler
//  HSV/color adjustments and additional commenting by CaliCoastReplay - 2016

// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

// PI is a mathematical constant relating the ratio of a circle's circumference (distance around
// the edge) to its diameter (distance between two points opposite on the edge).
// Change pi at your own peril, with your own apologies to God.
const float PI	 	= 3.14159265358;

// Can you explain these epsilons to a wide graphics audience?  YOUR comment could go here.
const float EPSILON	= 1e-3;
#define  EPSILON_NRM	(0.5 / iResolution.x)

// Constant indicaing the number of steps taken while marching the light ray.
const int NUM_STEPS = 6;

//Constants relating to the iteration of the heightmap for the wave, another part of the rendering
//process.
const int ITER_GEOMETRY = 2;
const int ITER_FRAGMENT =5;

// Constants that represent physical characteristics of the sea, can and should be changed and
//  played with
const float SEA_HEIGHT = 0.5;
const float SEA_CHOPPY = 3.0;
const float SEA_SPEED = 1.9;
const float SEA_FREQ = 0.24;
const vec3 SEA_BASE = vec3(0.11,0.19,0.22);
const vec3 SEA_WATER_COLOR = vec3(0.55,0.9,0.7);
#define SEA_TIME (iTime * SEA_SPEED)

//Matrix to permute the water surface into a complex, realistic form
mat2 octave_m = mat2(1.7,1.2,-1.2,1.4);

//Space bar key constant
const float KEY_SP    = 32.5/256.0;

//CaliCoastReplay :  These HSV/RGB translation functions are
//from http://gamedev.stackexchange.com/questions/59797/glsl-shader-change-hue-saturation-brightness
//This one converts red-green-blue color to hue-saturation-value color
vec3 rgb2hsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

//CaliCoastReplay :  These HSV/RGB translation functions are
//from http://gamedev.stackexchange.com/questions/59797/glsl-shader-change-hue-saturation-brightness
//This one converts hue-saturation-value color to red-green-blue color
vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

// math
// bteitler: Turn a vector of Euler angles into a rotation matrix
mat3 fromEuler(vec3 ang) {
	vec2 a1 = vec2(sin(ang.x),cos(ang.x));
    vec2 a2 = vec2(sin(ang.y),cos(ang.y));
    vec2 a3 = vec2(sin(ang.z),cos(ang.z));
    mat3 m;
    m[0] = vec3(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x);
	m[1] = vec3(-a2.y*a1.x,a1.y*a2.y,a2.x);
	m[2] = vec3(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y);
	return m;
}

// bteitler: A 2D hash function for use in noise generation that returns range [0 .. 1].  You could
// use any hash function of choice, just needs to deterministic and return
// between 0 and 1, and also behave randomly.  Googling "GLSL hash function" returns almost exactly
// this function: http://stackoverflow.com/questions/4200224/random-noise-functions-for-glsl
// Performance is a real consideration of hash functions since ray-marching is already so heavy.
float hash( vec2 p ) {
    float h = dot(p,vec2(127.1,311.7));
    return fract(sin(h)*83758.5453123);
}

// bteitler: A 2D psuedo-random wave / terrain function.  This is actually a poor name in my opinion,
// since its the "hash" function that is really the noise, and this function is smoothly interpolating
// between noisy points to create a continuous surface.
float noise( in vec2 p ) {
    vec2 i = floor( p );
    vec2 f = fract( p );

    // bteitler: This is equivalent to the "smoothstep" interpolation function.
    // This is a smooth wave function with input between 0 and 1
    // (since it is taking the fractional part of <p>) and gives an output
    // between 0 and 1 that behaves and looks like a wave.  This is far from obvious, but we can graph it to see
    // Wolfram link: http://www.wolframalpha.com/input/?i=plot+x*x*%283.0-2.0*x%29+from+x%3D0+to+1
    // This is used to interpolate between random points.  Any smooth wave function that ramps up from 0 and
    // and hit 1.0 over the domain 0 to 1 would work.  For instance, sin(f * PI / 2.0) gives similar visuals.
    // This function is nice however because it does not require an expensive sine calculation.
    vec2 u = f*f*(3.0-2.0*f);

    // bteitler: This very confusing looking mish-mash is simply pulling deterministic random values (between 0 and 1)
    // for 4 corners of the grid square that <p> is inside, and doing 2D interpolation using the <u> function
    // (remember it looks like a nice wave!)
    // The grid square has points defined at integer boundaries.  For example, if <p> is (4.3, 2.1), we will
    // evaluate at points (4, 2), (5, 2), (4, 3), (5, 3), and then interpolate x using u(.3) and y using u(.1).
    return -1.0+2.0*mix(
                mix( hash( i + vec2(0.0,0.0) ),
                     hash( i + vec2(1.0,0.0) ),
                        u.x),
                mix( hash( i + vec2(0.0,1.0) ),
                     hash( i + vec2(1.0,1.0) ),
                        u.x),
                u.y);
}

// bteitler: diffuse lighting calculation - could be tweaked to taste
// lighting
// diffuseは拡散という意味らしい。
float diffuse(vec3 n,vec3 l,float p) {
    return pow(dot(n,l) * 0.4 + 0.6,p);
}

// bteitler: specular lighting calculation - could be tweaked taste
// specularは鏡面という意味らしい。
float specular(vec3 n,vec3 l,vec3 e,float s) {
    float nrm = (s + 8.0) / (3.1415 * 8.0);
    return pow(max(dot(reflect(e,n),l),0.0),s) * nrm;
}

// bteitler: Generate a smooth sky gradient color based on ray direction's Y value
// sky

// 方向ベクトル（dir）のy成分に応じて空の色を決める感じ。見ての通りy成分以外は一切見てない。
vec3 getSkyColor(vec3 e) {
    e.y = max(e.y,0.0); // yが0.0以下の時は0.0として計算する
    vec3 ret;
    ret.x = pow(1.0-e.y,2.0); // すべての成分は1以下（ノルム1のベクトル）なので(1-y)^2 を取って・・つまり上の方に行くほど濃くなる感じ
    ret.y = 1.0-e.y; // こっちは単純に1-yを取って差をつけている
    ret.z = 0.6+(1.0-e.y)*0.4; // ここで青を表現・・0.6が基本で上の方に行くほど小さくなるようになっている。みたい。
    return ret;
}

// sea
// bteitler: TLDR is that this passes a low frequency random terrain through a 2D symmetric wave function that looks like this:
// http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20
// The <choppy> parameter affects the wave shape.

// 海面の状態をノイズを使って表現する。これはその基本波のようなもので実際は重ね合わせて使う。(fbmのようなもの)
// 跳ねる1-|sin(x)|と|cos(x)|を組み合わせて、2次元にして、そのx成分とy成分を掛けて、0.65乗して滑らかにして、
// それをchoppyしてまた変化させてる（わかりません）。

float sea_octave(vec2 uv, float choppy) {
    // bteitler: Add the smoothed 2D terrain / wave function to the input coordinates
    // which are going to be our X and Z world coordinates.  It may be unclear why we are doing this.
    // This value is about to be passed through a wave function.  So we have a smoothed psuedo random height
    // field being added to our (X, Z) coordinates, and then fed through yet another wav function below.
    uv += noise(uv);
    // Note that you could simply return noise(uv) here and it would take on the characteristics of our
    // noise interpolation function u and would be a reasonable heightmap for terrain.
    // However, that isn't the shape we want in the end for an ocean with waves, so it will be fed through
    // a more wave like function.  Note that although both x and y channels of <uv> have the same value added, there is a
    // symmetry break because <uv>.x and <uv>.y will typically be different values.

    // bteitler: This is a wave function with pointy peaks and curved troughs:
    // http://www.wolframalpha.com/input/?i=1-abs%28cos%28x%29%29%3B
    vec2 wv = 1.0-abs(sin(uv));

    // bteitler: This is a wave function with curved peaks and pointy troughs:
    // http://www.wolframalpha.com/input/?i=abs%28cos%28x%29%29%3B
    vec2 swv = abs(cos(uv));

    // bteitler: Blending both wave functions gets us a new, cooler wave function (output between 0 and 1):
    // http://www.wolframalpha.com/input/?i=abs%28cos%28x%29%29+%2B+abs%28cos%28x%29%29+*+%28%281.0-abs%28sin%28x%29%29%29+-+abs%28cos%28x%29%29%29
    wv = mix(wv,swv,wv);

    // bteitler: Finally, compose both of the wave functions for X and Y channels into a final
    // 1D height value, shaping it a bit along the way.  First, there is the composition (multiplication) of
    // the wave functions: wv.x * wv.y.  Wolfram will give us a cute 2D height graph for this!:
    // http://www.wolframalpha.com/input/?i=%7BAbs%5BCos%5Bx%5D%5D+%2B+Abs%5BCos%5Bx%5D%5D+%28%281.+-+Abs%5BSin%5Bx%5D%5D%29+-+Abs%5BCos%5Bx%5D%5D%29%7D+*+%7BAbs%5BCos%5By%5D%5D+%2B+Abs%5BCos%5By%5D%5D+%28%281.+-+Abs%5BSin%5By%5D%5D%29+-+Abs%5BCos%5By%5D%5D%29%7D
    // Next, we reshape the 2D wave function by exponentiation: (wv.x * wv.y)^0.65.  This slightly rounds the base of the wave:
    // http://www.wolframalpha.com/input/?i=%7B%7BAbs%5BCos%5Bx%5D%5D+%2B+Abs%5BCos%5Bx%5D%5D+%28%281.+-+Abs%5BSin%5Bx%5D%5D%29+-+Abs%5BCos%5Bx%5D%5D%29%7D+*+%7BAbs%5BCos%5By%5D%5D+%2B+Abs%5BCos%5By%5D%5D+%28%281.+-+Abs%5BSin%5By%5D%5D%29+-+Abs%5BCos%5By%5D%5D%29%7D%7D%5E0.65
    // one last final transform (with choppy = 4) results in this which resembles a recognizable ocean wave shape in 2D:
    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5Bx%5D%5D+%2B+Abs%5BCos%5Bx%5D%5D+%28%281.+-+Abs%5BSin%5Bx%5D%5D%29+-+Abs%5BCos%5Bx%5D%5D%29%7D+*+%7BAbs%5BCos%5By%5D%5D+%2B+Abs%5BCos%5By%5D%5D+%28%281.+-+Abs%5BSin%5By%5D%5D%29+-+Abs%5BCos%5By%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4
    // Note that this function is called with a specific frequency multiplier which will stretch out the wave.  Here is the graph
    // with the base frequency used by map and map_detailed (0.16):
    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20
    return pow(1.0-pow(wv.x * wv.y,0.65),choppy);
}

// bteitler: Compute the distance along Y axis of a point to the surface of the ocean
// using a low(er) resolution ocean height composition function (less iterations).

// 海面の高さを計算してる。detailedより短いコードになってますね。
// レイトレで海面上の点を取得するために使ってるのはこっち。
float map(vec3 p) {
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;
    // uvは水平座標で横方向は縮めている（どうしてかは知らない）

    // bteitler: Compose our wave noise generation ("sea_octave") with different frequencies
    // and offsets to achieve a final height map that looks like an ocean.  Likely lots
    // of black magic / trial and error here to get it to look right.  Each sea_octave has this shape:
    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20
    // which should give you an idea of what is going.  You don't need to graph this function because it
    // appears to your left :)
    float d, h = 0.0;
    for(int i = 0; i < ITER_GEOMETRY; i++) {

      // 繰り返し回数は2. detailedでは5に設定されていますね。

        // bteitler: start out with our 2D symmetric wave at the current frequency
    	d = sea_octave((uv+SEA_TIME)*freq,choppy);
        // bteitler: stack wave ontop of itself at an offset that varies over time for more height and wave pattern variance
    	//d += sea_octave((uv-SEA_TIME)*freq,choppy);

      // dが海面の変化でampが海面の高さでそれにより変動を算出してhに加算する感じ

        h += d * amp; // bteitler: Bump our height by the current wave function

        // bteitler: "Twist" our domain input into a different space based on a permutation matrix
        // The scales of the matrix values affect the frequency of the wave at this iteration, but more importantly
        // it is responsible for the realistic assymetry since the domain is shiftly differently.
        // This is likely the most important parameter for wave topology.

      // ツイスト処理。octave_mは行列。要はfbmに似たことをしている。のか？
    	uv *=  octave_m;

      // この辺はfbmぽい。振動数を上げて（2よりはちょっと小さいが）、振幅は減らしている（0.5どころか0.22倍だけど）
        freq *= 1.9; // bteitler: Exponentially increase frequency every iteration (on top of our permutation)
        amp *= 0.22; // bteitler: Lower the amplitude every frequency, since we are adding finer and finer detail
        // bteitler: finally, adjust the choppy parameter which will effect our base 2D sea_octave shape a bit.  This makes
        // the "waves within waves" have different looking shapes, not just frequency and offset
        // choppyも・・0.2倍して0.8を足してる、つまり0.2に近づくように小さくしているみたいね。
        choppy = mix(choppy,1.0,0.2);
    }
    // そんな波を重ね合わせて最終的にhを算出して、それをp.yから引くことで海面からの高さとする。
    return p.y - h;
}

// bteitler: Compute the distance along Y axis of a point to the surface of the ocean
// using a high(er) resolution ocean height composition function (more iterations).

// より詳しく海面の高さを計算するパート。
// この計算はgetNormal内で行っている。
float map_detailed(vec3 p) {
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;

    // bteitler: Compose our wave noise generation ("sea_octave") with different frequencies
    // and offsets to achieve a final height map that looks like an ocean.  Likely lots
    // of black magic / trial and error here to get it to look right.  Each sea_octave has this shape:
    // http://www.wolframalpha.com/input/?i=%7B1-%7B%7B%7BAbs%5BCos%5B0.16x%5D%5D+%2B+Abs%5BCos%5B0.16x%5D%5D+%28%281.+-+Abs%5BSin%5B0.16x%5D%5D%29+-+Abs%5BCos%5B0.16x%5D%5D%29%7D+*+%7BAbs%5BCos%5B0.16y%5D%5D+%2B+Abs%5BCos%5B0.16y%5D%5D+%28%281.+-+Abs%5BSin%5B0.16y%5D%5D%29+-+Abs%5BCos%5B0.16y%5D%5D%29%7D%7D%5E0.65%7D%7D%5E4+from+-20+to+20
    // which should give you an idea of what is going.  You don't need to graph this function because it
    // appears to your left :)
    float d, h = 0.0;
    // こっちは5回であっちより大きく設定されてるね
    for(int i = 0; i < ITER_FRAGMENT; i++) {
        // bteitler: start out with our 2D symmetric wave at the current frequency
    	d = sea_octave((uv+SEA_TIME)*freq,choppy);
        // bteitler: stack wave ontop of itself at an offset that varies over time for more height and wave pattern variance
    	d += sea_octave((uv-SEA_TIME)*freq,choppy);

        h += d * amp; // bteitler: Bump our height by the current wave function

        // bteitler: "Twist" our domain input into a different space based on a permutation matrix
        // The scales of the matrix values affect the frequency of the wave at this iteration, but more importantly
        // it is responsible for the realistic assymetry since the domain is shiftly differently.
        // This is likely the most important parameter for wave topology.
    	uv *= octave_m/1.2;
      // ツイストも弱く設定されてて、なるほど繰り返しが多いからそこら辺考慮してあんま飛び出さないようにしているのね。

        freq *= 1.9; // bteitler: Exponentially increase frequency every iteration (on top of our permutation)
        amp *= 0.22; // bteitler: Lower the amplitude every frequency, since we are adding finer and finer detail
        // bteitler: finally, adjust the choppy parameter which will effect our base 2D sea_octave shape a bit.  This makes
        // the "waves within waves" have different looking shapes, not just frequency and offset
        choppy = mix(choppy,1.0,0.2);
    }
    return p.y - h;
}

// bteitler:
// p: point on ocean surface to get color for // pは海面上の点のようです
// n: normal on ocean surface at <p> // nはpでの海面の法線ベクトルのようですね
// l: light (sun) direction // lは入射光、でいいのかな・・
// eye: ray direction from camera position for this pixel // eyeはカメラ位置からpに向かうベクトル、要はdirのことかと思われる。
// dist: distance from camera to point <p> on ocean surface // distはあっちにも出てくるけどoriからpに向かうベクトルでしょうね。

// 海の色の計算。フレネル反射とかいうのを使ってるんだけど・・
// それ以外にも相当複雑なことやってて手に負えない（（（（

vec3 getSeaColor(vec3 p, vec3 n, vec3 l, vec3 eye, vec3 dist) {
    // bteitler: Fresnel is an exponential that gets bigger when the angle between ocean
    // surface normal and eye ray is smaller
    float fresnel = 1.0 - max(dot(n,-eye),0.0);
    fresnel = pow(fresnel,3.0) * 0.45;

    // bteitler: Bounce eye ray off ocean towards sky, and get the color of the sky
    // 反射
    // getSkyColorとあるから空の色が反射している様子を反映させているのだろうね。空が赤くなったらそれも反映される・・？
    vec3 reflected = getSkyColor(reflect(eye,n))*0.99;

    // bteitler: refraction effect based on angle between light surface normal
    // 屈折
    vec3 refracted = SEA_BASE + diffuse(n,l,80.0) * SEA_WATER_COLOR * 0.27;

    // ここでやってるのはrefracted（屈折）による色とreflected（反射）による色をフレネル値の割合で混ぜ合わせている、と。
    // bteitler: blend the refracted color with the reflected color based on our fresnel term
    vec3 color = mix(refracted,reflected,fresnel);

    // bteitler: Apply a distance based attenuation factor which is stronger
    // at peaks
    // peakにいくほどattenuate（減衰）する効果を付与しているらしい。
    float atten = max(1.0 - dot(dist,dist) * 0.001, 0.0);
    color += SEA_WATER_COLOR * (p.y - SEA_HEIGHT) * 0.15 * atten;

    // bteitler: Apply specular highlight
    // スペキュラーでハイライト（？？）
    // 若干明るくしている（それしかわからん）
    color += vec3(specular(n,l,eye,90.0))*0.5;

    return color;
}

// bteitler: Estimate the normal at a point <p> on the ocean surface using a slight more detailed
// ocean mapping function (using more noise octaves).
// Takes an argument <eps> (stands for epsilon) which is the resolution to use
// for the gradient.  See here for more info on gradients: https://en.wikipedia.org/wiki/Gradient
// tracing

// xzそれぞれの方向に動かして海面の高さを取得してその差でもって法線ベクトルの水平方向成分を
// 近似的に出している。近似なので、その分海面の高さをより詳細に取得する必要があるわけ。
// 垂直方向はepsだけど多分同じことを考えている。。多分。で、正規化。
vec3 getNormal(vec3 p, float eps) {
    // bteitler: Approximate gradient.  An exact gradient would need the "map" / "map_detailed" functions
    // to return x, y, and z, but it only computes height relative to surface along Y axis.  I'm assuming
    // for simplicity and / or optimization reasons we approximate the gradient by the change in ocean
    // height for all axis.
    vec3 n;
    n.y = map_detailed(p); // bteitler: Detailed height relative to surface, temporarily here to save a variable?
    n.x = map_detailed(vec3(p.x+eps,p.y,p.z)) - n.y; // bteitler approximate X gradient as change in height along X axis delta
    n.z = map_detailed(vec3(p.x,p.y,p.z+eps)) - n.y; // bteitler approximate Z gradient as change in height along Z axis delta
    // bteitler: Taking advantage of the fact that we know we won't have really steep waves, we expect
    // the Y normal component to be fairly large always.  Sacrifices yet more accurately to avoid some calculation.
    n.y = eps;
    return normalize(n);

    // bteitler: A more naive and easy to understand version could look like this and
    // produces almost the same visuals and is a little more expensive.
    // vec3 n;
    // float h = map_detailed(p);
    // n.y = map_detailed(vec3(p.x,p.y+eps,p.z)) - h;
    // n.x = map_detailed(vec3(p.x+eps,p.y,p.z)) - h;
    // n.z = map_detailed(vec3(p.x,p.y,p.z+eps)) - h;
    // return normalize(n);
}


//CaliCoastReplay :  Keyboard checking function from the iChannel representing keyboard input
float isKeyPressed(float key)
{
	return texture( iChannel1, vec2(key, 1.0) ).x;
}

// bteitler: Find out where a ray intersects the current ocean
// 光線が海面とどこで衝突するかみたいな。出力は？

// 要するにpからどれだけdirの方向に行ったら海面に達するかを区間縮小で求めていて、
// 最終的にその距離tmidが出力されさらにpはそのときの海面における地点になる。
// 空の場合はtmidに相当する500.0が返されpは特にいじってないね。

//
float heightMapTracing(vec3 ori, vec3 dir, out vec3 p) {
    float tm = 0.0;
    float tx = 500.0; // bteitler: a really far distance, this could likely be tweaked a bit as desired

    // bteitler: At a really far away distance along the ray, what is it's height relative
    // to the ocean in ONLY the Y direction?
    float hx = map(ori + dir * tx);

    // bteitler: A positive height relative to the ocean surface (in Y direction) at a really far distance means
    // this pixel is pure sky.  Quit early and return the far distance constant.
    if(hx > 0.0) return tx;

    // bteitler: hm starts out as the height of the camera position relative to ocean.
    float hm = map(ori + dir * tm);

    // bteitler: This is the main ray marching logic.  This is probably the single most confusing part of the shader
    // since height mapping is not an exact distance field (tells you distance to surface if you drop a line down to ocean
    // surface in the Y direction, but there could have been a peak at a very close point along the x and z
    // directions that is closer).  Therefore, it would be possible/easy to overshoot the surface using the raw height field
    // as the march distance.  The author uses a trick to compensate for this.
    float tmid = 0.0;
    for(int i = 0; i < NUM_STEPS; i++) { // bteitler: Constant number of ray marches per ray that hits the water
        // bteitler: Move forward along ray in such a way that has the following properties:
        // 1. If our current height relative to ocean is higher, move forward more
        // 2. If the height relative to ocean floor very far along the ray is much lower
        //    below the ocean surface, move forward less
        // Idea behind 1. is that if we are far above the ocean floor we can risk jumping
        // forward more without shooting under ocean, because the ocean is mostly level.
        // The idea behind 2. is that if extruding the ray goes farther under the ocean, then
        // you are looking more orthgonal to ocean surface (as opposed to looking towards horizon), and therefore
        // movement along the ray gets closer to ocean faster, so we need to move forward less to reduce risk
        // of overshooting.
        tmid = mix(tm,tx, hm/(hm-hx));
        p = ori + dir * tmid;

    	float hmid = map(p); // bteitler: Re-evaluate height relative to ocean surface in Y axis

        if(hmid < 0.0) { // bteitler: We went through the ocean surface if we are negative relative to surface now
            // bteitler: So instead of actually marching forward to cross the surface, we instead
            // assign our really far distance and height to be where we just evaluated that crossed the surface.
            // Next iteration will attempt to go forward more and is less likely to cross the boundary.
            // A naive implementation might have returned <tmid> immediately here, which
            // results in a much poorer / somewhat indeterministic quality rendering.
            tx = tmid;
            hx = hmid;
        } else {
            // Haven't hit surface yet, easy case, just march forward
            tm = tmid;
            hm = hmid;
        }
    }

    // bteitler: Return the distance, which should be really close to the height map without going under the ocean
    return tmid;
}

// main
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    // bteitler: 2D Pixel location passed in as raw pixel, let's divide by resolution
    // to convert to coordinates between 0 and 1
    vec2 uv = fragCoord.xy / iResolution.xy;

    // -1.0～1.0に正規化したうえでなんかやってるけど、
    // これ要するに縦横の短い方で割ってるのと同じこと。だからまとめて一行で書けるよ。上含めた3行は。
    uv = uv * 2.0 - 1.0; //  bteitler: Shift pixel coordinates from 0 to 1 to between -1 and 1
    uv.x *= iResolution.x / iResolution.y; // bteitler: Aspect ratio correction - if you don't do this your rays will be distorted

    // ここで時間速くしてるね。
    float time = iTime * 2.7; // bteitler: Animation is based on time, but allows you to scrub the animation based on mouse movement

    // ray

    // bteitler: Calculated a vector that smoothly changes over time in a sinusoidal (wave) pattern.
    // This will be used to drive where the user is looking in world space.
   // vec3 ang = vec3(sin(time*3.0)*0.1,sin(time)*0.2+0.3,time);

   // pitchとyawがマウスの位置によって変化しているみたいなんだけどよくわかんね
   // わかりました！！！！（バカみたい）
   // マウスダウン状態でドラッグしないと変化が反映されないようです（単純な話だった）
   // で、マウスダウンで左右に動かしたらぐるぐるした。これがヨーイングか！！！！

   // たとえばrollをなくすと始点が固定されてそのうえでpitchとyawみたいになる？
   // rollがまずあってそのうえでpitchとyawみたいなイメージで捉えるのがよさそうね。続きは帰ってから（予定）

   // まずrollはPIになんかリサージュっぽい感じの変化を
   // まあ基本はPIかな、・・
   // ていうかね。cosとsinの和を括弧でくくったりくくらなかったりしてるのがなんか気持ち悪いのよね。統一して欲しい。
   // そういうのどうでもいいっていう人もいるけど（あれとかこれとか）自分は嫌です。
    float roll = PI + sin(iTime)/14.0 + cos(iTime/2.0)/14.0 ;

   // pitchは・・
   // マウスのyの値が0.8より大きいとpitchが上方修正、0.8より低いと下方修正？？
   // ここら辺は空しか見えなくならないように角度を調整しているっぽいな。
   // でもさぁ、やっぱPI基準で値決めた方が分かりやすいと思うんだよな・・自分で書くときに適当に修正しよ。

   // 大きいほど上を向くようになってるからそういうことなんでしょ。
    float pitch = PI*1.021 + (sin(iTime/2.0)+ cos(iTime))/40.0
        + (iMouse.y/iResolution.y - .8)*PI/3.0  ;

  // yawは・・・まあ、fromEuler見た方が早いかな・・わけわからん。
  // こっちはyと違って完全にxの値に左右されるわけね。0.0から1.0の間のはずだけどね。
  // 0.0からPI * 4.0まで動くから結構広範囲なのね。
    float yaw = iMouse.x/iResolution.x * PI * 4.0;
    vec3 ang = vec3(roll, pitch, yaw);
   // vec3 ang = vec3(roll,pitch,0);

    // bteitler: Calculate the "origin" of the camera in world space based on time.  Camera is located
    // at height 3.5 atx 0 (zero), and flies over the ocean in the z axis over time.
    vec3 ori = vec3(0.0,3.5,time*3.0);

    // bteitler: This is the ray direction we are shooting from the camera location ("ori") that we need to light
    // for this pixel.  The -2.0 indicates we are using a focal length of 2.0 - this is just an artistic choice and
    // results in about a 90 degree field of view.
    //  CaliCoastReplay :  Adjusted slightly to a lower focal length.  Seems to dramatize the scene.
    vec3 dir = normalize(vec3(uv.xy,-1.6));

    // bteitler: Distort the ray a bit for a fish eye effect (if you remove this line, it will remove
    // the fish eye effect and look like a realistic perspective).

    // これをやると凸レンズ？みたいになって水平線がまあるくなる
    // オリジナルがそんな感じなので多分それを真似しているようです（デフォルトではOFFになっている）
   //  dir.z += length(uv) * 0.15;

    // bteitler: Renormalize the ray direction, and then rotate it based on the previously calculated
    // animation angle "ang".  "fromEuler" just calculates a rotation matrix from a vector of angles.
    // if you remove the " * fromEuler(ang)" part, you will disable the camera rotation animation.

    // angから行列を作ってそれによりdirをいじっているんだけどどういうことなんだろう。
    // そもそもどんな行列ができるのかはfromEulerを解析しないとどうしようもない。
    dir = normalize(dir) * fromEuler(ang);

    // tracing

    // bteitler: ray-march to the ocean surface (which can be thought of as a randomly generated height map)
    // and store in p
    vec3 p;
    // pは何も宣言しない場合vec3(0.0)で初期化されるようになっている（他のベクトルも同様）。それは分かるんだけど他のコードでは意図しない挙動をするから
    // こういうのは行儀が悪いよね。自分なら絶対やらない。
    heightMapTracing(ori,dir,p);
    // で、この中でpを決めているんだけど、空に該当する場合pは0のままというわけなんだね。。

    // これpが定まらなかったらどうするんだろ・・ていうか大丈夫なのかこれ
    vec3 dist = p - ori; // bteitler: distance vector to ocean surface for this pixel's ray

    // bteitler: Calculate the normal on the ocean surface where we intersected (p), using
    // different "resolution" (in a sense) based on how far away the ray traveled.  Normals close to
    // the camera should be calculated with high resolution, and normals far from the camera should be calculated with low resolution
    // The reason to do this is that specular effects (or non linear normal based lighting effects) become fairly random at
    // far distances and low resolutions and can cause unpleasant shimmering during motion.
    vec3 n = getNormal(p,
             dot(dist,dist)   // bteitler: Think of this as inverse resolution, so far distances get bigger at an expnential rate
                * EPSILON_NRM // bteitler: Just a resolution constant.. could easily be tweaked to artistic content
           );

    // bteitler: direction of the infinitely far away directional light.  Changing this will change
    // the sunlight direction.
    vec3 light = normalize(vec3(0.0,1.0,0.8));

    // CaliCoastReplay:  Get the sky and sea colors
	vec3 skyColor = getSkyColor(dir);
    vec3 seaColor = getSeaColor(p,n,light,dir,dist);

    //Sea/sky preprocessing

    //CaliCoastReplay:  A distance falloff for the sea color.   Drastically darkens the sea,
    //this will be reversed later based on day/night.
    seaColor /= sqrt(sqrt(length(dist))) ;


    //CaliCoastReplay:  Day/night mode
    bool night;
    if( isKeyPressed(KEY_SP) > 0.0 )    //night mode!
    {
        //Brighten the sea up again, but not too bright at night
    	seaColor *= seaColor * 8.5;

        //Turn down the sky
    	skyColor /= 1.69;

        //Store that it's night mode for later HSV calcc
        night = true;
    }
    else  //day mode!
    {
        //Brighten the sea up again - bright and beautiful blue at day
    	seaColor *= sqrt(sqrt(seaColor)) * 4.0;
        skyColor *= 1.05;
        skyColor -= 0.03;
        night = false;
    }


    //CaliCoastReplay:  A slight "constrasting" for the sky to match the more contrasted ocean
    skyColor *= skyColor;


    //CaliCoastReplay:  A rather hacky manipulation of the high-value regions in the image that seems
    //to add a subtle charm and "sheen" and foamy effect to high value regions through subtle darkening,
    //but it is hacky, and not physically modeled at all.
    vec3 seaHsv = rgb2hsv(seaColor);
    if (seaHsv.z > .75 && length(dist) < 50.0)
       seaHsv.z -= (0.9 - seaHsv.z) * 1.3;
    seaColor = hsv2rgb(seaHsv);

    // bteitler: Mix (linear interpolate) a color calculated for the sky (based solely on ray direction) and a sea color
    // which contains a realistic lighting model.  This is basically doing a fog calculation: weighing more the sky color
    // in the distance in an exponential manner.

    // ----------------- 一番肝心の海の色と空の色を混ぜるところ ------------------ //

    // smoothstepの(0.0, -0.05)とかいう意味不明なことしてる。
    // smoothstep(0.0, 0.05, -dir.y)で書き換えたら同じ挙動を示したのでそういうことなんでしょう。
    // これ以上考えたくないので書き換えちゃいましょう。その方が早い。安定しない挙動はすべて排除！排除！
    vec3 color = mix(
        skyColor,
        seaColor,
    	//pow(smoothstep(0.0,-0.05,dir.y), 0.3) // bteitler: Can be thought of as "fog" that gets thicker in the distance
        pow(smoothstep(0.0, 0.05, -dir.y), 0.3)
    );
    // ということはdir.yが-0.05と0.0の間の所でスムースに変化して0.0より上の場合はskyColorで0.05より下の場合はseaColorで
    // はっきり決まるということね。0.3乗しているのは装飾部分でおそらく排除しても大差ない・・・？
    // と思ったけど排除して1.0にしたら水平線がぼやけちゃった。水平線をくっきりさせるための0.3乗らしい。なるほどねー。
    // 0.1乗だとくっきりしすぎてギザギザ、これは失敗だね。単純なsmoothstepでうまい具合にならなかったのを試行錯誤で調整して見つけ出したっぽいな。

    // Postprocessing （事後処理）

    // 事後処理ばっさりカットしてfragColor = vec4(color, 1.0)したらなんか空が青っぽくなった。
    // 海の色とちがうっぽくなった。
    // でも海面の様子とかは普通だった。どゆこと？？

    // 今見てきたけどオリジナルの方がシンプルだった（つまりここで終わり）。
    // その代わりなんか平均みたいなの取ってたけど。
    // ただでさえ重いこの一連の処理を9回も実行して平均取ってたのよ。そりゃ遅いわ。外したら速くなった。
    // でもこっちの装飾モリモリでも同じ速さだしどうなってるんだろうね。

    // bteitler: Apply an overall image brightness factor as the final color for this pixel.  Can be
    // tweaked artistically.
    fragColor = vec4(pow(color,vec3(0.75)), 1.0);

    // 0.75乗で切ってもやっぱ見た目がなんか違う、でもよくわからんね。
    // まあこの辺はおいおい理解していくことにしてとりあえずskyColorとseaColorの方にフォーカスしますかね。
    // あとクリックで視点が変わるのどうしてなのか知りたい。

    // CaliCoastReplay:  Adjust hue, saturation, and value adjustment for an even more processed look
    // hsv.x is hue, hsv.y is saturation, and hsv.z is value
    vec3 hsv = rgb2hsv(fragColor.xyz);
    //CaliCoastReplay: Increase saturation slightly
    hsv.y += 0.131;
    //CaliCoastReplay:
    //A pseudo-multiplicative adjustment of value, increasing intensity near 1 and decreasing it near
    //0 to achieve a more contrasted, real-world look
    hsv.z *= sqrt(hsv.z) * 1.1;

    if (night)
    {
    ///CaliCoastReplay:
    //Slight value adjustment at night to turn down global intensity
        hsv.z -= 0.045;
        hsv*=0.8;
        hsv.x += 0.12 + hsv.z/100.0;
        //Highly increased saturation at night op, oddly.  Nights appear to be very colorful
        //within their ranges.
        hsv.y *= 2.87;
    }
    else
    {
      //CaliCoastReplay:
        //Add green tinge to the high range
      //Turn down intensity in day in a different way

        hsv.z *= 0.9;

        //CaliCoastReplay:  Hue alteration
        hsv.x -= hsv.z/10.0;
        hsv.x += 0.02 + hsv.z/50.0;
        //Final brightening
        hsv.z *= 1.01;
        //This really "cinemafies" it for the day -
        //puts the saturation on a squared, highly magnified footing.
        //Worth looking into more as to exactly why.
       // hsv.y *= 5.10 * hsv.y * sqrt(hsv.y);
        hsv.y += 0.07;
    }

    //CaliCoastReplay:
    //Replace the final color with the adjusted, translated HSV values
    fragColor.xyz = hsv2rgb(hsv);
}

// ------------- 行列の掛け算や行列の定義について ----------------- //
// math
// bteitler: Turn a vector of Euler angles into a rotation matrix
mat3 fromEuler(vec3 ang) {
	vec2 a1 = vec2(cos(ang.x), sin(ang.x));
  vec2 a2 = vec2(cos(ang.y), sin(ang.y));
  vec2 a3 = vec2(cos(ang.z), sin(ang.z));
  // 画面の横揺れ（roll）
  mat3 m1 = mat3(a1.x, a1.y, 0.0, -a1.y, a1.x, 0.0, 0.0, 0.0, 1.0);
  m1[0] = vec3(a1.x, a1.y, 0.0);
	m1[1] = vec3(-a1.y, a1.x, 0.0);
	m1[2] = vec3(0.0, 0.0, 1.0);
  // 縦揺れ（pitch）
  mat3 m2;
  m2[0] = vec3(1.0, 0.0, 0.0);
	m2[1] = vec3(0.0, a2.x, a2.y);
	m2[2] = vec3(0.0, -a2.y, a2.x);
  // 水平回転（yaw）
  mat3 m3;
  m3[0] = vec3(a3.x, 0.0, a3.y);
	m3[1] = vec3(0.0, 1.0, 0.0);
	m3[2] = vec3(-a3.y, 0.0, a3.x);
  mat3 m;
  m[0] = vec3(a1.x, a2.x * a1.y, a1.y * a2.y);
  m[1] = vec3(-a1.y, a2.x * a1.x, a1.x * a2.y);
  m[2] = vec3(0.0, -a2.y, a2.x);
  // m1, m2, m3の順に適用される
	return m3 * m2 * m1;
}
// まずeuler * dirに直してみた。
// m1を掛ける場合にちゃんとあれしてる・・どうもm1[0]は一列目の縦ベクトル？
// ですね。つまり、上から下、次の列の上から下、次の列の上から下、でOK.
// 次に、掛け算について・・
// OK! 行列の掛け算でm2 * m1を計算してできた結果のmを用意して、
// return のところをm2 * m1とmで比較したら一緒になった。
// つまり縦ベクトルで行列は表現されていて、0とか1でアクセスしているのは
// 縦ベクトルなのだね！そして掛け算もその通りになってるみたい。
// つまりdirを縦ベクトルとみなした場合の計算結果になっている。

// というわけで。
/*
まず、行列の定義の仕方。
mat3 m = mat3(a, b, c, a', b', c', a'', b'', c'');
これでできるのはこれ：
  a  a'  a''
  b  b'  b''
  c  c'  c''
なのです。そして、
m[0] = (a, b, c), m[1] = (a', b', c'), m[2] = (a'', b'', c'')
という感じ。すべてvec3のベクトルですよと。これを使って定義することも可能。
次にベクトルとの掛け算について・・
m * v
とする場合、上記の表現でのmに縦ベクトルとしてのvを掛けてできる縦ベクトルがそのまま演算結果になるっぽい。
なので、m1, m2がある場合、行列の掛け算で上記の記法の下でm1 * m2を計算してできるmに対するmvは、
m * v (= (m1 * m2) * v) = m1 * (m2 * v)
と一緒になるようなのです。だからm1, m2の順に掛けたかったら行列の掛け算としては(m2 * m1)を実行しないといけないみたいね！

テスト完了。これでeuler * dirにしたら意図した通りの挙動になった（ただしrollとpitchからPIを引き算している（当然ね））

*/
