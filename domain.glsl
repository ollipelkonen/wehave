#version 410 core

struct Hit {
  float dist;
  vec3 ambientColor;
  vec3 diffuseColor;
  vec3 specularColor;
  float shininess;
};

const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 200.0;
const float EPSILON = 0.0001;
const float PI = acos(-1.0);
const float TAU = (2*PI);
const float PHI = (sqrt(5)*0.5 + 0.5);

#define saturate(x) clamp(x, 0, 1)
#define fmod(x, y) mod(x, y)
#define atan2(x,y) atan(x,y)
float sgn(float x) {
       return (x<0)?-1.0:1.0;
}

vec2 sgn(vec2 v) {
       return vec2((v.x<0)?-1:1, (v.y<0)?-1:1);
}

vec3 hash3( vec3 n )
{
    return fract(sin(n)*vec3(158.5453123,278.1459123,341.3490423));
}
vec2 hash2( vec2 n )
{
    return fract(sin(n)*vec2(158.5453123,278.1459123));
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









/***************************************************************************

        M A T H

***************************************************************************/



mat4 rotateY(float theta) {
  float c = cos(theta);
  float s = sin(theta);
  return mat4(
      vec4( c,  0, -s,  0),
      vec4( 0,  1,  0,  0),
      vec4( s,  0,  c,  0),
      vec4( 0,  0,  0,  1)
  );
}
mat4 rotateX(float theta) {
  float c = cos(theta);
  float s = sin(theta);
  return mat4(
      vec4( 1,  0,  0,  0),
      vec4( 0,  c,  s,  0),
      vec4( 0, -s,  c,  0),
      vec4( 0,  0,  0,  1)
  );
}
mat4 rotateZ(float theta) {
  float c = cos(theta);
  float s = sin(theta);
  return mat4(
      vec4( c, -s,  0,  0),
      vec4( s,  c,  0,  0),
      vec4( 0,  0,  1,  0),
      vec4( 0,  0,  0,  1)
  );
}






/***************************************************************************

        R E N D E R E R

***************************************************************************/


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



