  #version 410 core

uniform float fGlobalTime; // in seconds
uniform vec2 v2Resolution; // viewport resolution (in pixels)

uniform sampler1D texFFT; // towards 0.0 is bass / lower freq, towards 1.0 is higher / treble freq
uniform sampler1D texFFTSmoothed; // this one has longer falloff and less harsh transients
uniform sampler1D texFFTIntegrated; // this is continually increasing

layout(location = 0) out vec4 out_color; // out_color must be written in order to see anything

const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;
const float PI = acos(-1);
const float TAU = (2*PI);
const float PHI = (sqrt(5)*0.5 + 0.5);

#define saturate(x) clamp(x, 0, 1)
#define fmod(x, y) mod(x, y)
#define atan2(x,y) atan(x,y)
float sgn(float x) {
	return (x<0)?-1:1;
}

vec2 sgn(vec2 v) {
	return vec2((v.x<0)?-1:1, (v.y<0)?-1:1);
}


/***************************************************************************

        SIGNED DISTANCE FIELD FUNCTIONS


from:
    MERCURY http://mercury.sexy
    gargaj
    Akleman and Chen


***************************************************************************/


float sdRoundBox( vec3 p, vec3 b, float r )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}

float sdVerticalCapsule( vec3 p, float h, float r )
{
  p.y -= clamp( p.y, 0.0, h );
  return length( p ) - r;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sdSphere( vec3 p, float r )
{
  return length(p) - r;
}

float sdHexPrism( vec3 p, vec2 h )
{
  const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
  p = abs(p);
  p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
  vec2 d = vec2(
       length(p.xy-vec2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
       p.z-h.y );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdHexPrismY( vec3 p, vec2 h )
{
  const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
  p = abs(p);
  p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
  vec2 d = vec2(
       length(p.xy-vec2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
       p.z-h.y );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdLeg( vec3 p, vec3 start, float angle )
{
  float l = 1;
  vec2 f = vec2(cos(angle)*l, sin(angle)*l );
  return length( f - p.xy ) - length( p-start );
//  return mat3(angle * start) - start + p
}


// Plane with normal n (n is normalized) at some distance from the origin
float sdPlane(vec3 p, vec3 n, float distanceFromOrigin) {
  return dot(p, n) + distanceFromOrigin;
}


// infinite box
float sdBox2Cheap(vec2 p, vec2 b) {
  vec2 m = abs(p)-b;
  return max(m.x, m.y);
}

// Endless "corner"
float sdCorner (vec2 p) {
  vec2 m = min(p, vec2(0,0));
  return length(max(p, vec2(0,0))) + max(m.x, m.y);
}


// Blobby ball object. You've probably seen it somewhere. This is not a correct distance bound, beware.
float sdBlob(vec3 p) {
  p = abs(p);
  if (p.x < max(p.y, p.z)) p = p.yzx;
  if (p.x < max(p.y, p.z)) p = p.yzx;
  float b = max(max(max(
    dot(p, normalize(vec3(1, 1, 1))),
    dot(p.xz, normalize(vec2(PHI+1, 1)))),
    dot(p.yx, normalize(vec2(1, PHI)))),
    dot(p.xz, normalize(vec2(1, PHI))));
  float l = length(p);
  return l - 1.5 - 0.2 * (1.5 / 2)* cos(min(sqrt(1.01 - b / l)*(PI / 0.25), PI));
}

// Cylinder standing upright on the xz plane
float sdCylinder(vec3 p, float r, float height) {
  float d = length(p.xz) - r;
  d = max(d, abs(p.y) - height);
  return d;
}

// Capsule: A Cylinder with round caps on both sides
/*float fCapsule(vec3 p, float r, float c) {
  return mix(length(p.xz) - r, length(vec3(p.x, abs(p.y) - c, p.z)) - r, step(c, abs(p.y)));
}*/

// Distance to line segment between <a> and <b>, used for fCapsule() version 2below
float sdLineSegment(vec3 p, vec3 a, vec3 b) {
  vec3 ab = b - a;
  float t = saturate(dot(p - a, ab) / dot(ab, ab));
  return length((ab*t + a) - p);
}

// Capsule version 2: between two end points <a> and <b> with radius r 
float sdCapsule(vec3 p, vec3 a, vec3 b, float r) {
  return sdLineSegment(p, a, b) - r;
}

// Torus in the XZ-plane
float sdTorus(vec3 p, float smallRadius, float largeRadius) {
  return length(vec2(length(p.xz) - largeRadius, p.y)) - smallRadius;
}

// A circle line. Can also be used to make a torus by subtracting the smaller radius of the torus.
float sdCircle(vec3 p, float r) {
  float l = length(p.xz) - r;
  return length(vec2(p.y, l));
}

// A circular disc with no thickness (i.e. a cylinder with no height).
// Subtract some value to make a flat disc with rounded edge.
float sdDisc(vec3 p, float r) {
  float l = length(p.xz) - r;
  return l < 0 ? abs(p.y) : length(vec2(p.y, l));
}


// Hexagonal prism, circumcircle variant
float fHexagonCircumcircle(vec3 p, vec2 h) {
  vec3 q = abs(p);
  return max(q.y - h.y, max(q.x*sqrt(3)*0.5 + q.z*0.5, q.z) - h.x);
  //this is mathematically equivalent to this line, but less efficient:
  //return max(q.y - h.y, max(dot(vec2(cos(PI/3), sin(PI/3)), q.zx), q.z) - h.x);
}

// Hexagonal prism, incircle variant
float fHexagonIncircle(vec3 p, vec2 h) {
  return fHexagonCircumcircle(p, vec2(h.x*sqrt(3)*0.5, h.y));
}

// Cone with correct distances to tip and base circle. Y is up, 0 is in the middle of the base.
float fCone(vec3 p, float radius, float height) {
  vec2 q = vec2(length(p.xz), p.y);
  vec2 tip = q - vec2(0, height);
  vec2 mantleDir = normalize(vec2(height, radius));
  float mantle = dot(tip, mantleDir);
  float d = max(mantle, -q.y);
  float projected = dot(tip, vec2(mantleDir.y, -mantleDir.x));
  
  // distance to tip
  if ((q.y > height) && (projected < 0)) {
    d = max(d, length(tip));
  }
  
  // distance to base ring
  if ((q.x > radius) && (projected > length(vec2(height, radius)))) {
    d = max(d, length(q - vec2(radius, 0)));
  }
  return d;
}


////////////////////////////////////////////////////////////////
//
//                DOMAIN MANIPULATION OPERATORS
//
////////////////////////////////////////////////////////////////
//
// Conventions:
//
// Everything that modifies the domain is named pSomething.
//
// Many operate only on a subset of the three dimensions. For those,
// you must choose the dimensions that you want manipulated
// by supplying e.g. <p.x> or <p.zx>
//
// <inout p> is always the first argument and modified in place.
//
// Many of the operators partition space into cells. An identifier
// or cell index is returned, if possible. This return value is
// intended to be optionally used e.g. as a random seed to change
// parameters of the distance functions inside the cells.
//
// Unless stated otherwise, for cell index 0, <p> is unchanged and cells
// are centered on the origin so objects don't have to be moved to fit.
//
//
////////////////////////////////////////////////////////////////


// rotate 2d space with given angle
void tRotate(inout vec2 p, float angel) {
  vec2 s = sin(vec2(angel, angel + PI * .5));
  p *= mat2(s.y, -s.x, s);
}

// Rotate around a coordinate axis (i.e. in a plane perpendicular to that axis) by angle <a>.
// Read like this: R(p.xz, a) rotates "x towards z".
// This is fast if <a> is a compile-time constant and slower (but still practical) if not.
void pR(inout vec2 p, float a) {
  p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

// Shortcut for 45-degrees rotation
void pR45(inout vec2 p) {
  p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

// Repeat space along one axis. Use like this to repeat along the x axis:
// <float cell = pMod1(p.x,5);> - using the return value is optional.
float pMod1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = fmod(p + halfsize, size) - halfsize;
  return c;
}

// Same, but mirror every second cell so they match at the boundaries
float pModMirror1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = fmod(p + halfsize,size) - halfsize;
  p *= fmod(c, 2.0)*2 - 1;
  return c;
}

// Repeat the domain only in positive direction. Everything in the negative half-space is unchanged.
float pModSingle1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  if (p >= 0)
    p = fmod(p + halfsize, size) - halfsize;
  return c;
}

// Repeat only a few times: from indices <start> to <stop> (similar to above, but more flexible)
float pModInterval1(inout float p, float size, float start, float stop) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = fmod(p+halfsize, size) - halfsize;
  if (c > stop) { //yes, this might not be the best thing numerically.
    p += size*(c - stop);
    c = stop;
  }
  if (c <start) {
    p += size*(c - start);
    c = start;
  }
  return c;
}


// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolar(inout vec2 p, float repetitions) {
  float angle = 2*PI/repetitions;
  float a = atan2(p.y, p.x) + angle/2.;
  float r = length(p);
  float c = floor(a/angle);
  a = fmod(a,angle) - angle/2.;
  p = vec2(cos(a), sin(a))*r;
  // For an odd number of repetitions, fix cell index of the cell in -x direction
  // (cell index would be e.g. -5 and 5 in the two halves of the cell):
  if (abs(c) >= (repetitions/2)) c = abs(c);
  return c;
}

// Repeat in two dimensions
vec2 pMod2(inout vec2 p, vec2 size) {
  vec2 c = floor((p + size*0.5)/size);
  p = fmod(p + size*0.5,size) - size*0.5;
  return c;
}

// Same, but mirror every second cell so all boundaries match
vec2 pModMirror2(inout vec2 p, vec2 size) {
  vec2 halfsize = size*0.5;
  vec2 c = floor((p + halfsize)/size);
  p = fmod(p + halfsize, size) - halfsize;
  p *= fmod(c,vec2(2,2))*2 - vec2(1,1);
  return c;
}

// Same, but mirror every second cell at the diagonal as well
vec2 pModGrid2(inout vec2 p, vec2 size) {
  vec2 c = floor((p + size*0.5)/size);
  p = fmod(p + size*0.5, size) - size*0.5;
  p *= fmod(c,vec2(2,2))*2 - vec2(1,1);
  p -= size/2;
  if (p.x > p.y) p.xy = p.yx;
  return floor(c/2);
}

// Repeat in three dimensions
vec3 pMod3(inout vec3 p, vec3 size) {
  vec3 c = floor((p + size*0.5)/size);
  p = fmod(p + size*0.5, size) - size*0.5;
  return c;
}

// Mirror at an axis-aligned plane which is at a specified distance <dist> from the origin.
float pMirror (inout float p, float dist) {
  float s = sgn(p);
  p = abs(p)-dist;
  return s;
}

// Mirror in both dimensions and at the diagonal, yielding one eighth of the space.
// translate by dist before mirroring.
vec2 pMirrorOctant (inout vec2 p, vec2 dist) {
  vec2 s = sgn(p);
  pMirror(p.x, dist.x);
  pMirror(p.y, dist.y);
  if (p.y > p.x)
    p.xy = p.yx;
  return s;
}

// Reflect space at a plane
float pReflect(inout vec3 p, vec3 planeNormal, float offset) {
  float t = dot(p, planeNormal)+offset;
  if (t < 0) {
    p = p - (2*t)*planeNormal;
  }
  return sgn(t);
}



/****************************************************************************

//             OBJECT COMBINATION OPERATORS

*****************************************************************************/

//
// We usually need the following boolean operators to combine two objects:
// Union: OR(a,b)
// Intersection: AND(a,b)
// Difference: AND(a,!b)
// (a and b being the distances to the objects).
//
// float fTwoBoxes(vec3 p) {
//   float box0 = fBox(p, vec3(1));
//   float box1 = fBox(p-vec3(1), vec3(1));
//   return fOpUnionChamfer(box0, box1, 0.2);
// }
//


// The "Chamfer" flavour makes a 45-degree chamfered edge (the diagonal of a square of size <r>):
float fOpUnionChamfer(float a, float b, float r) {
  return min(min(a, b), (a - r + b)*sqrt(0.5));
}

// Intersection has to deal with what is normally the inside of the resulting object
// when using union, which we normally don't care about too much. Thus, intersection
// implementations sometimes differ from union implementations.
float fOpIntersectionChamfer(float a, float b, float r) {
  return max(max(a, b), (a + r + b)*sqrt(0.5));
}

// Difference can be built from Intersection or Union:
float fOpDifferenceChamfer (float a, float b, float r) {
  return fOpIntersectionChamfer(a, -b, r);
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float fOpUnionRound(float a, float b, float r) {
  vec2 u = max(vec2(r - a,r - b), vec2(0,0));
  return max(r, min (a, b)) - length(u);
}

float fOpIntersectionRound(float a, float b, float r) {
  vec2 u = max(vec2(r + a,r + b), vec2(0,0));
  return min(-r, max (a, b)) + length(u);
}

float fOpDifferenceRound (float a, float b, float r) {
  return fOpIntersectionRound(a, -b, r);
}


// The "Columns" flavour makes n-1 circular columns at a 45 degree angle:
float fOpUnionColumns(float a, float b, float r, float n) {
  if ((a < r) && (b < r)) {
    vec2 p = vec2(a, b);
    float columnradius = r*sqrt(2)/((n-1)*2+sqrt(2));
    pR45(p);
    p.x -= sqrt(2)/2*r;
    p.x += columnradius*sqrt(2);
    if (fmod(n,2) == 1) {
      p.y += columnradius;
    }
    // At this point, we have turned 45 degrees and moved at a point on the
    // diagonal that we want to place the columns on.
    // Now, repeat the domain along this direction and place a circle.
    pMod1(p.y, columnradius*2);
    float result = length(p) - columnradius;
    result = min(result, p.x);
    result = min(result, a);
    return min(result, b);
  } else {
    return min(a, b);
  }
}

float fOpDifferenceColumns(float a, float b, float r, float n) {
  a = -a;
  float m = min(a, b);
  //avoid the expensive computation where not needed (produces discontinuity though)
  if ((a < r) && (b < r)) {
    vec2 p = vec2(a, b);
    float columnradius = r*sqrt(2)/n/2.0;
    columnradius = r*sqrt(2)/((n-1)*2+sqrt(2));

    pR45(p);
    p.y += columnradius;
    p.x -= sqrt(2)/2*r;
    p.x += -columnradius*sqrt(2)/2;

    if (fmod(n,2) == 1) {
      p.y += columnradius;
    }
    pMod1(p.y,columnradius*2);

    float result = -length(p) + columnradius;
    result = max(result, p.x);
    result = min(result, a);
    return -min(result, b);
  } else {
    return -m;
  }
}

float fOpIntersectionColumns(float a, float b, float r, float n) {
  return fOpDifferenceColumns(a,-b,r, n);
}

// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float fOpUnionStairs(float a, float b, float r, float n) {
  float s = r/n;
  float u = b-r;
  return min(min(a,b), 0.5 * (u + a + abs ((fmod (u - a + s, 2 * s)) - s)));
}

// We can just call Union since stairs are symmetric.
float fOpIntersectionStairs(float a, float b, float r, float n) {
  return -fOpUnionStairs(-a, -b, r, n);
}

float fOpDifferenceStairs(float a, float b, float r, float n) {
  return -fOpUnionStairs(-a, b, r, n);
}


// Similar to fOpUnionRound, but more lipschitz-y at acute angles
// (and less so at 90 degrees). Useful when fudging around too much
// by MediaMolecule, from Alex Evans' siggraph slides
float fOpUnionSoft(float a, float b, float r) {
  float e = max(r - abs(a - b), 0);
  return min(a, b) - e*e*0.25/r;
}


// produces a cylindical pipe that runs along the intersection.
// No objects remain, only the pipe. This is not a boolean operator.
float fOpPipe(float a, float b, float r) {
  return length(vec2(a, b)) - r;
}

// first object gets a v-shaped engraving where it intersect the second
float fOpEngrave(float a, float b, float r) {
  return max(a, (a + r - abs(b))*sqrt(0.5));
}

// first object gets a capenter-style groove cut out
float fOpGroove(float a, float b, float ra, float rb) {
  return max(a, min(a + ra, rb - abs(b)));
}

// first object gets a capenter-style tongue attached
float fOpTongue(float a, float b, float ra, float rb) {
  return min(a, max(a - ra, abs(b) - rb));
}



float opOnion( in float sdf, in float thickness )
{
  return abs(sdf)-thickness;
}

float opUnion( float d1, float d2 ) { return min(d1,d2); }

float opSubtraction( float d1, float d2 ) { return max(-d1,d2); }

float opIntersection( float d1, float d2 ) { return max(d1,d2); }


float opRound( in float p, float rad )
{
    return p - rad;
}


/*vec3 opRepLim( in vec3 p, in float c, in vec3 l, in sdf3d primitive )
{
    vec3 q = p-c*clamp(round(p/c),-l,l);
    return primitive( q );
}*/


/***************************************************************************

        M A T H

***************************************************************************/



mat4 rotateY(float theta) {
  float c = cos(theta);
  float s = sin(theta);
  return mat4(
      vec4(c, 0, -s, 0),
      vec4(0, 1, 0, 0),
      vec4(s, 0, c, 0),
      vec4(0, 0, 0, 1)
  );
  return mat4(
      vec4(c, 0, s, 0),
      vec4(0, 1, 0, 0),
      vec4(-s, 0, c, 0),
      vec4(0, 0, 0, 1)
  );
}




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
  float d = length(p.xz) - h.x  + sin(p.y*ang+rot)/10.0;
  return max(d, abs(p.y) - h.y);
}



float bolt( vec3 samplePoint, vec3 pos ) {
  float k = sdHexPrism(samplePoint - pos, vec2(1,0.6) );
  k = max( k, sdSphere( samplePoint - pos, 1.18 ) );
  k = min( k, sdRuuvi( samplePoint-pos - vec3(0,0,2), vec2(0.5,2), 25 ) );
  return k;
}
float boltNut( vec3 samplePoint, vec3 pos ) {
  float k = sdHexPrism(samplePoint - pos, vec2(1,0.6) );
  k = max( k, sdSphere( samplePoint - pos, 1.18 ) );
  k = opSubtraction( sdRuuvi( samplePoint - pos, vec2(0.5,0.7), 25 ), k );
  return k;
}

float lightBulb( vec3 samplePoint )
{
  float k = sdSphere( samplePoint, 3 );
  float k2 = sdCylinder( samplePoint - vec3(0,-3,0), 1.75, 1);
  k = fOpUnionRound( k, k2, 1);
  
  float kanta = sdRuuviY( samplePoint-vec3(0,-5,0), vec2(1.75,1), 10 );
  return fOpUnionRound(k, kanta, 1);
  return k;
}

float ruuviScene( vec3 samplePoint ) {

  return lightBulb( samplePoint );

  vec3 pos = vec3(0,0,0);
  return boltNut( samplePoint, pos );
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




/***************************************************************************

        R E N D E R E R

***************************************************************************/




vec2 sdf( vec3 eye, vec3 viewRayDirection )
{
    float start = MIN_DIST;
    float end = MAX_DIST;
    float depth = start;
    float dust = 0.0;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        vec3 pos = eye + depth * viewRayDirection;
//        float dist = scene(pos);
	    float dist = ruuviScene(pos);
        if (dist < EPSILON)
            return vec2(depth,dust);
        depth += dist;
        dust += 0.001 * dist;
        if (depth >= end)
          return vec2(end,dust);
    }
    return vec2(end,dust);
}

vec3 rayDirection(float fieldOfView, vec2 size, vec2 fragCoord) {
  vec2 xy = fragCoord - size / 2.0;
  float z = size.y / tan(radians(fieldOfView) / 2.0);
  return normalize(vec3(xy, -z));
}


mat4 viewMatrix(vec3 eye, vec3 center, vec3 up) {
  vec3 f = normalize(center - eye);
  vec3 s = normalize(cross(f, up));
  vec3 u = cross(s, f);
  return mat4(
    vec4(s.x,u.x,-f.x, 0.0),
    vec4(s.y,u.y,-f.y, 0.0),
    vec4(s.z,u.z,-f.z, 0.0),
    vec4(0.0, 0.0, 0.0, 1) );
  return mat4(
      vec4(s, 0.0),
      vec4(u, 0.0),
      vec4(-f, 0.0),
      vec4(0.0, 0.0, 0.0, 1)
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
  vec3 eye = vec3(0.0, 4.0 -  cos(fGlobalTime) , 14.0 + sin(fGlobalTime) );
  vec3 tgt = vec3(0.0, 0.0, 0.0);

    
  eye.x += 0*sin(texture( texFFTSmoothed, 0 ).r * 100 );
  eye.y += 0*sin(texture( texFFTSmoothed, 0.1 ).r * 100 );
  eye.z += 0*sin(texture( texFFTSmoothed, 0.3 ).r * 100 );

  eye.x = 34.0 * sin(fGlobalTime);
  eye.y = 0;
  eye.z = 34.0 * cos(fGlobalTime);




//  mat4 viewToWorld = viewMatrix(eye, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
  mat4 viewToWorld = viewMatrix( eye, tgt, vec3(0.0, 1.0, 0.0));
  vec3 worldDir = ( ( vec4(viewDir, 0.0) * viewToWorld ) ).xyz;


//worldDir.z *=  sin(v2Resolution.x*100);

/*  eye.xy = uv.xy*5;
  eye.z = -1;
  worldDir = vec3(0,0,1);*/

  vec2 depth = sdf( eye, worldDir );

  vec3 hitPos = eye + depth.x * worldDir;
  out_color.xyz = estimateNormal( hitPos );
  out_color.w = 1;
  
//  depth *= depth;
  //out_color = vec4( c,c,c,1 );
  /*if ( depth.x > (MAX_DIST - EPSILON) )
  {
      return float4(0,0,1,1);
  }
  vec3 hitPos = eye + depth.x * worldDir;
  //depth.x = saturate( depth.x * 1.0 );
  vec3 norm = estimateNormal(hitPos);

  float d = sdf( )
	vec4 t = plas( m * 3.14, fGlobalTime ) / d;
	t = clamp( t, 0.0, 1.0 );
	out_color = f + t;*/
}


