
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
  //float interpolation = clamp(0.5 + 0.5 * (a - b) / r, 0.0, 1.0);
  //return mix(a, b, interpolation) - r* interpolation * (1.0 - interpolation);
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

#define fOpUnion min

// Similar to fOpUnionRound, but more lipschitz-y at acute angles
// (and less so at 90 degrees). Useful when fudging around too much
// by MediaMolecule, from Alex Evans' siggraph slides
float fOpUnionSoft(float a, float b, float r) {
float e = max(r - abs(a - b), 0);
return min(a, b) - e*e*0.25/r;
}

float opSmoothUnion( float d1, float d2, float k ) {
  float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
  return mix( d2, d1, h ) - k*h*(1.0-h); }

float opSmoothSubtraction( float d1, float d2, float k ) {
  float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
  return mix( d2, -d1, h ) + k*h*(1.0-h); }

float opSmoothIntersection( float d1, float d2, float k ) {
  float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
  return mix( d2, d1, h ) + k*h*(1.0-h); }


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


float fOpInterpolate(float shape1, float shape2, float amount){
  return mix(shape1, shape2, amount);
}

float opOnion( in float sdf, in float thickness )
{
return abs(sdf)-thickness;
}

float opUnion( float d1, float d2 ) { return min(d1,d2); }

float opSubtraction( float d1, float d2 ) { return max(-d1,d2); }

float opIntersection( float d1, float d2 ) { return max(d1,d2); }

float opRound( in float p, float rad ) { return p - rad; }


//TODO:
/*float opDisplace( in sdf3d primitive, in vec3 p )
{
float d1 = primitive(p);
float d2 = displacement(p);
return d1+d2;
}
float opTwist( in sdf3d primitive, in vec3 p )
{
const float k = 10.0; // or some other amount
float c = cos(k*p.y);
float s = sin(k*p.y);
mat2  m = mat2(c,-s,s,c);
vec3  q = vec3(m*p.xz,p.y);
return primitive(q);
}
float opCheapBend( in sdf3d primitive, in vec3 p )
{
const float k = 10.0; // or some other amount
float c = cos(k*p.x);
float s = sin(k*p.x);
mat2  m = mat2(c,-s,s,c);
vec3  q = vec3(m*p.xy,p.z);
return primitive(q);
}*/


/*vec3 opRepLim( in vec3 p, in float c, in vec3 l, in sdf3d primitive )
{
  vec3 q = p-c*clamp(round(p/c),-l,l);
  return primitive( q );
}*/
// polynomial smooth min (k = 0.1);
float sminPoly( float a, float b, float k )
{
  float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
  return mix( b, a, h ) - k*h*(1.0-h);
}

// polynomial smooth min (k = 0.1);
float sminPoly_MediaMolecule( float a, float b, float k )
{
  float h = max( k-abs(a-b), 0.0 )/k;
  return min( a, b ) - h*h*k*(1.0/4.0);
}

// return a factor for blending material too.
float sminQuadratic ( float a, float b, float k )
{
  float h = max( k-abs(a-b), 0.0 )/k;
  float m = h*h*0.5;
  float s = m*k*(1.0/2.0);
  return (a<b) ? vec2(a-s,m).x : vec2(b-s,1.0-m).x;
}

//--------------------------------------------------------

// polynomial smooth min (k = 0.1);
float sminCubic( float a, float b, float k )
{
  float h = max( k-abs(a-b), 0.0 )/k;
  return min( a, b ) - h*h*h*k*(1.0/6.0);
}

// return a factor for blending material too.
float sminCubic2( float a, float b, float k )
{
  float h = max( k-abs(a-b), 0.0 )/k;
  float m = h*h*h*0.5;
  float s = m*k*(1.0/3.0); 
  return (a<b) ? vec2(a-s,m).x : vec2(b-s,1.0-m).x;
}

// return a factor for blending material too.
float sminN( float a, float b, float k )
{
  const float n = 3.; //3.
  float h = max( k-abs(a-b), 0.0 )/k;
  float m = pow(h, n)*0.5;
  float s = m*k/n; 
  return (a<b) ? vec2(a-s,m).x : vec2(b-s,1.0-m).x;
}

// return a factor for blending material too.
float sminExp( in float a, in float b, in float k )
{
  float f1 = exp2( -k*a );
  float f2 = exp2( -k*b );
  return vec2(-log2(f1+f2)/k,f2).x;
}
//--------------------------------------------------------

// power smooth min (k = 8);
float sminPow( float a, float b, float k )
{
  k/=2.;
  a = pow( a, k ); b = pow( b, k );
  return pow( (a*b)/(a+b), 1.0/k );
}
