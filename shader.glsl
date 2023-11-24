
uniform float fGlobalTime; // in seconds

uniform float fMidiKnob1;
uniform float fMidiKnob2;
uniform float fMidiKnob3;
uniform float fMidiKnob4;
uniform float fMidiKnob5;
uniform float fMidiKnob6;
uniform float fMidiKnob7;
uniform float fMidiKnob8;
uniform float fMidiPad1;
uniform float fMidiPad2;
uniform float fMidiPad3;
uniform float fMidiPad4;
uniform float fMidiPad5;
uniform float fMidiPad6;
uniform float fMidiPad7;
uniform float fMidiPad8;
uniform vec2 v2Resolution; // viewport resolution (in pixels)

uniform sampler1D texFFT; // towards 0.0 is bass / lower freq, towards 1.0 is higher / treble freq
uniform sampler1D texFFTSmoothed; // this one has longer falloff and less harsh transients
uniform sampler1D texFFTIntegrated; // this is continually increasing

layout(location = 0) out vec4 out_color; // out_color must be written in order to see anything



struct Hit {
  float dist;
  vec3 normal;
  vec3 ambientColor;
  vec3 diffuseColor;
  vec3 specularColor;
  float shininess;
};



/***************************************************************************

        S C E N E S

***************************************************************************/




float hallway( vec3 samplePoint )
{
  float d = sdBox( samplePoint + vec3(0,3,0), vec3(7,0.1,100) );
  d = min( d, sdBox( samplePoint - vec3(0,7,0), vec3(7,0.1,100) ) );
  d = min( d, sdBox( samplePoint - vec3(7,0,0), vec3(0.3,7,100) ) );
  d = min( d, sdBox( samplePoint - vec3(-7,0,0), vec3(0.3,7,100) ) );
  
  d = opSubtraction( sdBox( samplePoint - vec3(-7,2,0), vec3(1.4) ), d );
  return d;
}



float tyyppi( vec3 samplePoint )
{
  //return sphereDist;
  //matrix4x4 m = rotateY( sin(_Time.y+samplePoint.y) );
  //vec3 cubePoint = mul( m, float4(samplePoint, 1.0) ).xyz;


  float f1 = sdRoundBox( samplePoint, vec3(1,0.5,1), 0.9 );
  float m1 = sdVerticalCapsule( samplePoint + vec3(1,2,1), 1.5, 0.2 );
  float m2 = sdVerticalCapsule( samplePoint + vec3(1,2,-1), 1.5, 0.2 );
  float m3 = sdVerticalCapsule( samplePoint + vec3(-1,2,-1), 1.5, 0.2 );
  float m4 = sdVerticalCapsule( samplePoint + vec3(-1,2,1), 1.5, 0.2 );

  float unionRad = 0.9;
  float k = fOpUnionRound( f1, m1, unionRad  );
  k = fOpUnionRound( k, m2, unionRad );
  k = fOpUnionRound( k, m3, unionRad );
  k = fOpUnionRound( k, m4, unionRad );

  return k;
}

float tyyppiScene(vec3 samplePoint) {
  float k = tyyppi( samplePoint - vec3(4,0,0) );

  k = min( hallway(samplePoint), k );
  return k;
}



float sdHexPrismRuuvi_skitso( vec3 p, vec2 h )
{
  vec3 k = vec3(-0.8660254, 0.5, 0.57735);
  pR(k.xy, (cos(p.z+1.57) * cos(p.x) * sin(p.y))* (fGlobalTime) / 100.0 );
  p = abs(p);
  p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
  k = vec3(-0.8660254, 0.5, 0.57735);
  vec2 d = vec2(
       length(p.xy-vec2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
       p.z-h.y );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdRuuvi( vec3 p, vec2 h, float ang )
{
  float rot = atan2(p.x,p.y);
  float d = length(p.xy) - h.x  + sin(p.z*ang+rot)/10.0;
  return max(d, abs(p.z) - h.y);
}
// another direction
float sdRuuviY( vec3 p, vec2 h, float ang )
{
  float rot = atan2(p.x,p.z);
  float l = length(p.xz);
  float d = l - h.x  + cos(p.y*ang+rot)/10.0;
  d = max(d, abs(p.y) - h.y);
  float md = sdCylinder(p - vec3(0,1.1,0), h.x, 0.2);
  return opUnion(md,d);
}



float bolt( vec3 samplePoint, vec3 pos ) {
  float k = sdHexPrism(samplePoint - pos, vec2(1.3,0.6) );
  k = max( k, sdSphere( samplePoint - pos, 1.52 ) );
  k = min( k, sdRuuvi( samplePoint-pos - vec3(0,0,5), vec2(0.8,5), 25 ) );
  return k;
}
float boltNut( vec3 samplePoint, vec3 pos ) {
  float k = sdHexPrism(samplePoint - pos, vec2(1,0.6) );
  k = max( k, sdSphere( samplePoint - pos, 1.2 ) );
  k = opSubtraction( sdRuuvi( samplePoint - pos, vec2(0.8,0.7), 25 ), k );
  return k;
}

float lightBulbKanta( vec3 p, float h, float r )
{
  p.y -= clamp( p.y, 0.0, h );
  return length( p ) - r;
}
float sqr(float a) { return a*a; }
float fOpUnionRound2(float a, float b, float r) {
  //return fOpUnionRound(a, b, r);
  
  float interpolation = clamp(0.5 + 0.5 * (a - b) / r, 0.0, 1.0);
  return mix(a, b, interpolation) - r* interpolation * (1.0 - interpolation);
  
  vec2 u = max(vec2( (r-a), (r-b) ), vec2(0,0));
  //return r-length(u);
  return max(r, min (a, b)) - length(u);
  
  return a*r + b*(1-r);
  //vec2 u = max(vec2(r - a,r - b), vec2(0,0));
  //return max(r, min (a, b)) - length(u);
}
float lightBulb( vec3 samplePoint )
{
  float k = sdSphere( samplePoint, 3 );
  float k2 = sdCylinder( samplePoint - vec3(0,-3,0), 1.75, 1);
  k = fOpUnionRound( k, k2, 1);
  
  float ruuvi = sdRuuviY( samplePoint-vec3(0,-5.3,0), vec2(1.75,1), 10 );
  //k = fOpUnionRound2(k, kanta, 1);
  k = fOpUnion(k, ruuvi);
  //k = kanta;
  //k = fOpUnionRound2(k, sdSphere(samplePoint -vec3(0,-6,0), 1.4), fMidiKnob7 );
  vec3 pos = vec3(0,-6.45,0);
  k = fOpUnionRound2( k, sdEllipsoid( samplePoint-pos, vec3(1.42,0.5,1.42) ), 0.2);
  //k = fOpUnionRound(k, lightBulbKanta(samplePoint -vec3(0,-7,0), 0.15, 1), 1 );
  
  
  
  return k;
}




float sdMouth( vec3 p, float le, float r1, float r2 )
{

  // from opCheapBend https://iquilezles.org/www/articles/distfunctions/distfunctions.htm
  const float ang = 2.0; // or some other amount
  float c = cos(ang*p.x / 18.0);
  float s = sin(ang*p.x / 18.0);
  mat2  m1 = mat2(c,-s,s,c);
  p = vec3(m1*p.xy,p.z);
  //p = vec3(m1*p.xz,p.y).xzy;
  
  //return sdRoundBox(p, vec3(1,0.2,0.5) * 5, 0.1);


  vec3 w = vec3( p.x/3+1, p.y, max(abs(p.z)-le,0.0) );

  float k = length(vec2(length(w.yz)-r1,w.x)) - r2;

  //ristikko
  //if ( sdBox(p,vec3(r2,r1,le)) < 0 )
  if ( sdRoundBox(p,vec3(r2*4.5,r1*0.35,le*1.6),0.2*5) < 0 )
  {
    //horizontal stripes are dy/xy
    float dz = 0.25;
    float dy = 0.3;
    vec3 q = p;
    q.z = abs(mod(p.z+dz, dz*2)-dz);
    q.y = abs(mod(p.y+dy, dy*2)-dy);
    float m = min( length(q.xz)-0.08, length(q.xy) -0.08 );  
    //float m = length(q.xz) -0.02;  
    return min(m,k);
  }

  return k;  
}
float carheadEar( vec3 p )
{
  float le = 0.1 *5;
  float r1 = 0.25 *5;
  float r2 = 0.1 *5;
  vec4 q = vec4( p.x, max(abs(p.y)-le,0.0), p.z, 1 );
  q = q - vec4(0,0.3,0,0) *5;
  q = q * rotateZ( fMidiKnob7 * 2 * PI );
  q = q + vec4(0,0.3,0,0) *5;
  return length(vec2(length(q.xy)-r1,q.z)) - r2;
}

float carhead( vec3 samplePoint ) {

  float k = sdSphere( samplePoint, 5 );
//return min(k,sdMouth( samplePoint - vec3(1.2,-0.8,0), 0.5, 0.12, 0.03 ));
//  float jaw = sdBox( samplePoint - vec3(0.4,-0.6,0), vec3(0.5,0.5,0.4) );
  float jaw = sdRoundBox( samplePoint - vec3(2.5,-3,0), vec3(0.25,0.25,0.2), 2.5 );
  k = fOpUnionSoft( jaw, k, 0.95 );
  // ohimot
  k = opSmoothSubtraction( sdBox(samplePoint - vec3(0,0,6.85), vec3(1.5,1.5,0.5)*5), k, 0.9 );
  k = opSmoothSubtraction( sdBox(samplePoint + vec3(0,0,6.85), vec3(1.5,1.5,0.5)*5), k, 0.9 );
  // nose
  k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(0.96,-0.12,0)*5, 0.5 ), 0.22 );
  k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(1.0,-0.2,0)*5, 0.5 ), 0.52 );
  k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(1.1,-0.35,0)*5, 0.5 ), 0.51 );
  //k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(1.3,-0.4,0.1), 0.01 ), 0.1 );
  // eyes
  k = opSubtraction( sdSphere(samplePoint - vec3(4,0,1.5), 1.5 ), k );
  k = opSubtraction( sdSphere(samplePoint - vec3(4,0,-1.5), 1.5 ), k );
  
  k = opOnion( k, 0.1 );
  //float sdLink( vec3 p, float le, float r1, float r2 )
  vec3 sampleX = samplePoint + vec3(0,0,-4.3);
  sampleX.x = abs(sampleX.x);
  //float ear = sdLink( sampleX, 0.1, 0.25, 0.1 );
  float ear = carheadEar( sampleX );
  k = min(ear, k );

  k = opSubtraction( sdSphere(samplePoint-vec3(4.5,-4.0,0),1.6), k );
  float mouth = sdMouth( samplePoint - vec3(4.0,-4,0), 0.3, 1.2, 0.3 );
  k = fOpUnionSoft( mouth, k, 0.96 );
  //k = min(mouth, k );

  //hajalla?
  /*float size = 100;
  float halfsize = size / 2.0;
  vec2 yz = samplePoint.yz;
  if ( length(samplePoint) < 1 )
  {
    float m = length( mod(samplePoint.yz, 1) ) - 0.1;
    k = min(m,k);
  }*/
  // "windows"
 // k = opSubtraction( sdBox( samplePoint - vec3(0.3,0.4,0), vec3(0.25,0.3,2) ), k );
//  k = opSubtraction( sdBox( samplePoint - vec3(-0.3,0.4,0), vec3(0.25,0.3,2) ), k );

//  k = opUnion( sdSphere(samplePoint - vec3(0.8,0,0.3), 0.4 ), k );
//  k = opUnion( sdSphere(samplePoint - vec3(0.8,0,-0.3), 0.4 ), k );

  return k;
}


float innerSphere( vec3 pos )
{
  return -sdSphere(pos,150);
}



float ruuviScene( vec3 samplePoint ) {

    /*float a = 10.;
    vec3 v3 = vec3(a,a,a);
  vec2 v2 = vec2(a,a);
    float m = length( mod(samplePoint.xy + v2, v2*2) - v2 ) - 0.1;*/
//    return m;

  //return tyyppi(samplePoint);
  //return hallway(samplePoint);
  
//  return carhead(samplePoint);


  //return min(m,k4);

  float k1 = lightBulb( samplePoint );
  //float k2 = sdVerticalCapsule( samplePoint, 5, 1 );

    k1 = min (innerSphere(samplePoint),k1);
  
  vec3 pos = vec3(0,4.5,0);
  float k2 = boltNut( samplePoint, pos );
pos = vec3(-4,0,0);
  float k3 =  bolt( samplePoint, pos );

  return min( k1, min(k2, k3));
}



float sdJalkaSegment(vec3 p, vec3 a, vec3 b) {
  vec3 ab = b - a;
  ab.x += 0.1*sin( length(ab)* fGlobalTime);
  float t = saturate(dot(p - a, ab) / dot(ab, ab));
  return length((ab*t + a) - p);
}

float sdJalkaRuuviSegment( vec3 p, vec2 h, float ang )
{
  float rot = atan2(p.x,p.z);
  float l = length(p.xz);
  float d = l - h.x  + cos(p.y*ang+rot)/10.0;
  d = max(d, abs(p.y) - h.y);
  float md = sdCylinder(p - vec3(0,1.1,0), h.x, 0.2);
  return opUnion(md,d);
}



float jalka( vec3 samplePoint, float t ) {
  
  float ang = 0.35*PI + 0.2 * PI * sin(4*PI*t);
  vec3 knp = vec3(0, -4 * sin(ang), 4 * cos(ang) );
  
  float ang2 = 2.9*PI + t * 2.0 * sin(4*PI*t);
  ang2 = pow( ang2, 2 ) / 40.0;
  vec3 knp2 = vec3(0, -4*sin(ang2), 4*cos(ang2) );
  
  float g1 = sdJalkaSegment( samplePoint, vec3(0,0,0), knp );
  float g2 = sdJalkaSegment( samplePoint, knp, knp + knp2 );
  float g3 = sdJalkaSegment( samplePoint * vec3(1,1,1), knp + knp2, knp+knp2 + vec3(0,0,1) );
  float r = 1;
  float g = fOpUnionSoft(g1,g2,0.4);
  g = min(g,g3);
  //float g = min(min(g1,g2),g3);
  return g - r;
}  

float jalkaScene( vec3 samplePoint ) {
  float t1 = fmod(fGlobalTime / 2.0, 1);
  float t2 = fmod(fGlobalTime / 2.0 + 0.25, 1);
  vec3 p1 = vec3(-2,0,0);
  Hit h;
  return min(
    jalka(samplePoint - p1, t1),
    jalka(samplePoint + p1, t2)
  );
}




const int JALKASCENE = 1;
const int TYYPPISCENE = 2;
const int RUUVISCENE = 3;
const int HEADSCENE = 4;

//const int scene = JALKASCENE;
//const int scene = TYYPPISCENE;
//const int scene = RUUVISCENE;
const int scene = HEADSCENE;




vec3 estimateNormal(vec3 p) {
  //TODO: use larger value further away
  float D = EPSILON;
  switch (scene) {
    case JALKASCENE:
      return normalize(vec3(
          jalkaScene(vec3(p.x + D, p.y, p.z)) - jalkaScene(vec3(p.x - D, p.y, p.z)),
          jalkaScene(vec3(p.x, p.y + D, p.z)) - jalkaScene(vec3(p.x, p.y - D, p.z)),
          jalkaScene(vec3(p.x, p.y, p.z  + D)) - jalkaScene(vec3(p.x, p.y, p.z - D))
      ));
    case RUUVISCENE:
      return normalize(vec3(
          ruuviScene(vec3(p.x + D, p.y, p.z)) - ruuviScene(vec3(p.x - D, p.y, p.z)),
          ruuviScene(vec3(p.x, p.y + D, p.z)) - ruuviScene(vec3(p.x, p.y - D, p.z)),
          ruuviScene(vec3(p.x, p.y, p.z  + D)) - ruuviScene(vec3(p.x, p.y, p.z - D))
      ));
    case TYYPPISCENE :
      return normalize(vec3(
          tyyppiScene(vec3(p.x + D, p.y, p.z)) - tyyppiScene(vec3(p.x - D, p.y, p.z)),
          tyyppiScene(vec3(p.x, p.y + D, p.z)) - tyyppiScene(vec3(p.x, p.y - D, p.z)),
          tyyppiScene(vec3(p.x, p.y, p.z  + D)) - tyyppiScene(vec3(p.x, p.y, p.z - D))
      ));
    case HEADSCENE :
      return normalize(vec3(
          carhead(vec3(p.x + D, p.y, p.z)) - carhead(vec3(p.x - D, p.y, p.z)),
          carhead(vec3(p.x, p.y + D, p.z)) - carhead(vec3(p.x, p.y - D, p.z)),
          carhead(vec3(p.x, p.y, p.z  + D)) - carhead(vec3(p.x, p.y, p.z - D))
      ));
  }
}




Hit sdf( vec3 eye, vec3 viewRayDirection, out float ao )
{
  float start = MIN_DIST;
  float end = MAX_DIST;
  float depth = start;
  float dust = 0.0;
  float aoSum = 0.0;
  float aoMaxSum = 0.0;
  Hit hit;
  hit.dist = start;
  for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
    vec3 pos = eye + depth * viewRayDirection;
//    float dist = scene(pos);
    switch (scene) {
      case JALKASCENE:
        hit.dist = jalkaScene(pos);
        hit.normal = estimateNormal(pos);
        break;
      case TYYPPISCENE:
        hit.dist = tyyppiScene(pos);
        hit.normal = estimateNormal(pos);
        break;
      case RUUVISCENE:
        hit.dist = ruuviScene(pos);
        hit.normal = estimateNormal(pos);
        break;
      case HEADSCENE:
        hit.dist = carhead(pos);
        hit.normal = estimateNormal(pos);
        break;
    }

    //aoSum    += 1. / pow(2., i) * hit.dist;
    //aoMaxSum += 1. / pow(2., i) * (i+1) * hit.dist;
    //ao = aoSum / aoMaxSum;
    ao = 1 - float(i) / (MAX_MARCHING_STEPS-1);

    //hit.dist = legScene( pos );
//    hit.dist = carhead(pos);
    if (hit.dist < EPSILON)
    {
      hit.dist = depth;
      return hit;
    }
    depth += hit.dist;
    dust += 0.001 * hit.dist;
    if (depth >= end)
    {
      hit.dist = depth;
      return hit;
      //return vec2(end,dust);
    }
  }
  return hit;
  //return vec2(end,dust);
}



vec3 texture1( vec2 uv ) {
  return vec3(
    0.7 + 0.2 * sin(uv.x),
    0.6 + 0.2 * cos(uv.y*1.2),
    0.6 + 0.3 * sin(uv.x + cos(uv.y))
  );
}


void main(void)
{
  vec2 uv = vec2(gl_FragCoord.x / v2Resolution.x, gl_FragCoord.y / v2Resolution.y);
  uv -= 0.5;
  uv /= vec2(v2Resolution.y / v2Resolution.x, 1);

  vec2 m;
  m.x = atan(uv.x / uv.y) / 3.14;
  m.y = 1 / length(uv) * .2;
  float d = m.y;

  float f = texture( texFFT, d ).r * 100;
  m.x += sin( fGlobalTime ) * 0.1;
  m.y += fGlobalTime * 0.25;

  vec3 viewDir = rayDirection(75.0, v2Resolution.xy, gl_FragCoord.xy);
  //vec3 eye = vec3(0.0, 4.0 -  cos(fGlobalTime) , 14.0 + sin(fGlobalTime) );

  vec3 tgt = vec3(0.0, 0.0, 0.0);

/*  eye.x += 0*sin(texture( texFFTSmoothed, 0 ).r * 100 );
  eye.y += 0*sin(texture( texFFTSmoothed, 0.1 ).r * 100 );
  eye.z += 0*sin(texture( texFFTSmoothed, 0.3 ).r * 100 );*/
    

float angle = fGlobalTime;
float distance = 5;
angle = cos(angle) + PI / 2 ;
angle = fMidiKnob8 + PI / 2.0;


  vec4 eye = 
    rotateY(fMidiKnob5 * PI * 3 + fGlobalTime/10.0)
    * rotateZ(fMidiKnob1 * PI * 3)
    * vec4( 6+24 * pow(fMidiKnob6+0.7, 3),0,0,1)
    ;




//  mat4 viewToWorld = viewMatrix(eye, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
  mat4 viewToWorld = viewMatrix( eye.xyz, tgt, vec3(0.0, 1.0, 0.0));
  vec3 worldDir = ( ( vec4(viewDir, 0.0) * viewToWorld ) ).xyz;


//worldDir.z *=  sin(v2Resolution.x*100);

/*  eye.xy = uv.xy*5;
  eye.z = -1;
  worldDir = vec3(0,0,1);*/

  float ao;
  Hit hit = sdf( eye.xyz, worldDir, ao );
  float depth = hit.dist;

  vec3 hitPos = eye.xyz + depth * worldDir;

  if (fMidiKnob4 < 0.5 )
  {
    out_color.xyz = hit.normal;
    float dd = out_color.x + out_color.y + out_color.x;
    dd /= 3;
    dd = pow(dd,fMidiKnob3-EPSILON);
    
    dd = clamp(dd,0,1);

    out_color.x = dd * (1-fMidiKnob2) + out_color.x * fMidiKnob2;
    out_color.y = dd * (1-fMidiKnob2) + out_color.y * fMidiKnob2;
    out_color.z = dd * (1-fMidiKnob2) + out_color.z * fMidiKnob2;
    out_color += out_color; // * pow( fMidiPad1, 10.0 );
    out_color.w = 1;
    
  //  depth *= depth;
    //out_color = vec4( c,c,c,1 );
    /*if ( depth > (MAX_DIST - EPSILON) )
    {
        return float4(0,0,1,1);
    }
    vec3 hitPos = eye + depth * worldDir;
    //depth = saturate( depth * 1.0 );
    vec3 norm = estimateNormal(hitPos);

    float d = sdf( )
    vec4 t = plas( m * 3.14, fGlobalTime ) / d;
    t = clamp( t, 0.0, 1.0 );
    out_color = f + t;*/
  }
  else if ( fMidiKnob4 < 1.0 )
  {
    vec3 lightPos = vec3( 10.0 * sin(fGlobalTime), 8.0*cos(fGlobalTime), 4.0*cos(fGlobalTime*0.7));
    vec3 lightDir = hitPos-lightPos;

    float distance = length(lightDir);
    lightDir = normalize(lightDir);
    vec3 ambientColor = vec3( 0.15, 0.15, 0.3 );
    vec3 diffuseColor = vec3( 0.4, 0.3, 0.5 );
  
    vec3 specColor = vec3( 1.0, 1.0, 1.0 );
    float nDotL = dot( hit.normal, lightDir );
    float intensity = saturate(nDotL);

    //float intensity = nDotL;
    vec3 H = normalize(lightDir + viewDir);
    float NdotH = dot( hit.normal, H);
    float spec = pow(saturate(NdotH), 20.1);

    //diffuseColor = hash3( hitPos.xyx / 10000 );

    /*vec3 n = hash3(hitPos);
    vec3 v = normalize(eye.xyz-hitPos);
		float nl = max ( 0.0, dot ( n, lightDir ) );
    vec3 kk = lightDir + v;
    vec3  h  = normalize( kk );
    float hn = max ( 0.0, dot ( h, n ) );
    float sp = pow ( hn, 150.0 );*/


    float vor = voronoi( hitPos.xy );
    float vor2 = voronoi( hitPos.xz );
    ambientColor.x += vor / 10;
    ambientColor.z += vor2 / 10;
    diffuseColor.x += vor;
    diffuseColor.z += vor2;
    
    out_color.xyz = ambientColor +
                     diffuseColor * intensity + // distance +
                     specColor * spec; // / distance;
                     
                   out_color.x *= ao;
                   out_color.y *= ao;
                   out_color.z *= ao;


//    out_color.y = intensity;

  }
}


