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
  const float ang = 10.0; // or some other amount
  float c = cos(ang*p.x);
  float s = sin(ang*p.x);
  mat2  m1 = mat2(c,-s,s,c);
  p = vec3(m1*p.xy,p.z);
  //p = vec3(m1*p.xz,p.y).xzy;
  return sdRoundBox(p, vec3(1,0.2,0.5), 0.1);


  vec3 w = vec3( p.x, p.y, max(abs(p.z)-le,0.0) );

  float k = length(vec2(length(w.yz)-r1,w.x)) - r2;


  //ristikko
  //if ( sdBox(p,vec3(r2,r1,le)) < 0 )
  if ( sdRoundBox(p,vec3(r2*0.9,r1*0.1,le*0.9),0.2) < 0 )
  {
    //horizontal stripes are dy/xy
    float dz = 0.05;
    float dy = 0.06;
    vec3 q = p;
    q.z = abs(mod(p.z+dz, dz*2)-dz);
    q.y = abs(mod(p.y+dy, dy*2)-dy);
    float m = min( length(q.xz)-0.02, length(q.xy) -0.02 );  
    //float m = length(q.xz) -0.02;  
    return min(m,k);
  }

  return k;  
}
float carheadEar( vec3 p )
{
  float le = 0.1;
  float r1 = 0.25;
  float r2 = 0.1;
  vec4 q = vec4( p.x, max(abs(p.y)-le,0.0), p.z, 1 );
  q = q - vec4(0,0.3,0,0);
  q = q * rotateZ( fMidiKnob7 * 2 * PI );
  q = q + vec4(0,0.3,0,0);
  return length(vec2(length(q.xy)-r1,q.z)) - r2;
}

float carhead( vec3 samplePoint ) {

  float k = sdSphere( samplePoint, 1 );

//return min(k,sdMouth( samplePoint - vec3(1.2,-0.8,0), 0.5, 0.12, 0.03 ));

//  float jaw = sdBox( samplePoint - vec3(0.4,-0.6,0), vec3(0.5,0.5,0.4) );
  float jaw = sdRoundBox( samplePoint - vec3(0.6,-0.7,0), vec3(0.05,0.05,0.04), 0.45 );
  k = fOpUnionSoft( jaw, k, 0.25 );

  // ohimot
  k = opSmoothSubtraction( sdBox(samplePoint - vec3(0,0,1.753), vec3(1.5,1.5,0.5)), k, 0.97 );
  k = opSmoothSubtraction( sdBox(samplePoint + vec3(0,0,1.753), vec3(1.5,1.5,0.5)), k, 0.97 );
 

  // nose
  k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(1.15,-0.12,0), 0.0001 ), 0.2 );
  k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(1.20,-0.3,0), 0.02 ), 0.2 );
  k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(1.3,-0.35,0), 0.01 ), 0.1 );
  //k = fOpUnionSoft( k, sdSphere(samplePoint - vec3(1.3,-0.4,0.1), 0.01 ), 0.1 );

  // eyes
  k = opSubtraction( sdSphere(samplePoint - vec3(0.8,0,0.3), 0.5 ), k );
  k = opSubtraction( sdSphere(samplePoint - vec3(0.8,0,-0.3), 0.5 ), k );
  
  k = opOnion( k, 0.1 );

  //float sdLink( vec3 p, float le, float r1, float r2 )

  vec3 sampleX = samplePoint + vec3(0,0,-1);
  sampleX.x = abs(sampleX.x);
  //float ear = sdLink( sampleX, 0.1, 0.25, 0.1 );
  float ear = carheadEar( sampleX );
  k = min(ear, k );


  float mouth = sdMouth( samplePoint - vec3(1.2,-0.8,0), 0.3, 0.2, 0.03 );
  k = min(mouth, k );


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

  vec3 pos = vec3(0,4.5,0);
  float k2 = boltNut( samplePoint, pos );
pos = vec3(-4,0,0);
  float k3 =  bolt( samplePoint, pos );

  return min( k1, min(k2, k3));
}



vec3 estimateNormal(vec3 p) {
  return normalize(vec3(
      ruuviScene(vec3(p.x + EPSILON, p.y, p.z)) - ruuviScene(vec3(p.x - EPSILON, p.y, p.z)),
      ruuviScene(vec3(p.x, p.y + EPSILON, p.z)) - ruuviScene(vec3(p.x, p.y - EPSILON, p.z)),
      ruuviScene(vec3(p.x, p.y, p.z  + EPSILON)) - ruuviScene(vec3(p.x, p.y, p.z - EPSILON))
  ));
  return normalize(vec3(
      tyyppiScene(vec3(p.x + EPSILON, p.y, p.z)) - tyyppiScene(vec3(p.x - EPSILON, p.y, p.z)),
      tyyppiScene(vec3(p.x, p.y + EPSILON, p.z)) - tyyppiScene(vec3(p.x, p.y - EPSILON, p.z)),
      tyyppiScene(vec3(p.x, p.y, p.z  + EPSILON)) - tyyppiScene(vec3(p.x, p.y, p.z - EPSILON))
  ));
}




Hit sdf( vec3 eye, vec3 viewRayDirection )
{
  float start = MIN_DIST;
  float end = MAX_DIST;
  float depth = start;
  float dust = 0.0;
  Hit hit;
  hit.dist = start;
  for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
    vec3 pos = eye + depth * viewRayDirection;
//    float dist = scene(pos);
    hit.dist = ruuviScene(pos);
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
/*  eye.x = distance * sin(angle);
  eye.y = 5 * cos(angle*1.3);
  eye.z = distance * cos(angle);
  eye.x = distance * sin( fMidiKnob1 * PI * 2 ) - cos ( fMidiKnob5 * PI * 2 );
  eye.y = 5  - sin( fMidiKnob5* PI * 2 );
  eye.z = distance * cos( fMidiKnob1 * PI * 2 );*/
  
  vec4 eye = vec4( 0.1+14 * pow(fMidiKnob6+0.5, 7),0,0,1);
  eye = eye * (rotateZ(fMidiKnob1 * PI * 4) * rotateY(fMidiKnob5 * PI * 4));
  




//  mat4 viewToWorld = viewMatrix(eye, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
  mat4 viewToWorld = viewMatrix( eye.xyz, tgt, vec3(0.0, 1.0, 0.0));
  vec3 worldDir = ( ( vec4(viewDir, 0.0) * viewToWorld ) ).xyz;


//worldDir.z *=  sin(v2Resolution.x*100);

/*  eye.xy = uv.xy*5;
  eye.z = -1;
  worldDir = vec3(0,0,1);*/

  Hit hit = sdf( eye.xyz, worldDir );
  float depth = hit.dist;

  vec3 hitPos = eye.xyz + depth * worldDir;
  out_color.xyz = estimateNormal( hitPos );
  
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


