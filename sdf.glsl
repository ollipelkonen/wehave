
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

float sdEllipsoid( in vec3 p, in vec3 r )
{
    return (length(p/r) - 1.0) * min(min(r.x,r.y),r.z);
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

float sdLink( vec3 p, float le, float r1, float r2 )
{
  vec3 q = vec3( p.x, max(abs(p.y)-le,0.0), p.z );
  return length(vec2(length(q.xy)-r1,q.z)) - r2;
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

