Shader "Hidden/gnabiili"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        // No culling or depth
        Cull Off ZWrite Off ZTest Always



        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"
            #include "kaikk.cginc"
            #include "matrix.cginc"
            #include "snoise.cginc"

            static const int MAX_MARCHING_STEPS = 255;
            static const float MIN_DIST = 0.0;
            static const float MAX_DIST = 100.0;
            static const float EPSILON = 0.0001;

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
            };

            float intersectSDF(float distA, float distB) {
                return max(distA, distB);
            }

            float unionSDF(float distA, float distB) {
                return min(distA, distB);
            }

            float differenceSDF(float distA, float distB) {
                return max(distA, -distB);
            }


            float huimaSDF( float3 p )
            {
                //p.xyz = fmod( p.xyz, 3 );
                
                p = mul( rotateY(p.y*sin(_Time.x)), p / 1.5 );
                p.xyz = pMod3( p.xyz, float3(1,1,1) );
                float3 d = abs(p) - (float3(1,1,1) / 2.0);
                
                float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
                float outsideDistance = length(max(d, 0.0));

                float p2 = mul( rotateY(p.y*sin(_Time.x)), p / 1.5 );

                return fOpDifferenceRound(insideDistance + outsideDistance, length(p)-1.0, 0.2 + 0.2*sin(_Time.w));
            }

            float sphereSDF(float3 p) {
                //p.xyz = fmod(p.xyz,1); 
                p = mul( rotateY(p.y*sin(_Time.x)), p );
                //p.xyz = pMod3( p.xyz, float3(1,1,1) );
                //p = mul( rotateY(p.y*sin(_Time.x)), p );
                return length(p) - 1.0;
            }

            float planeSDF( float3 p, float y )
            {
                float ip = 1;
                float d = fmod(abs(p.z),ip);
                return p.y-y+step(d,0.2)+sin(p.x)*0.1;
            }            
            
            float cubeSDF(float3 p, float3 size) {
                //p.xyz = fmod( p.xyz, 3 );
                //p = mul( rotateY(p.y*sin(_Time.x)), p / 2 );
                p = mul( rotateY(p.y*sin(_Time.x)), p );
                //p.xyz = pMod3( p.xyz, float3(1,1,1) );
                float3 d = abs(p) - (size / 2.0);
                
                /*float3 nois = p;
                noise(nois);
                d += 0.2 * nois;*/
                // Assuming p is inside the cube, how far is it from the surface?
                // Result will be negative or zero.
                float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
                
                // Assuming p is outside the cube, how far is it from the surface?
                // Result will be positive or zero.
                float outsideDistance = length(max(d, 0.0));
                
                return insideDistance + outsideDistance;
            }

            

            float fractalSDF( float3 p )
            {
                /*float e = 0.1;
                return sinh(p) / 2 * pow(e) * abs((')lG'(z) / logG(z)*/
                float3 w = p;

                for ( int i=0; i<20; i++ )
                {
                    float wr = sqrt(dot(w,w));
                    float wo = acos(p.y/wr);
                    float wi = atan2(p.x,p.z);

                    // scale and rotate the point
                    wr = pow( wr, 8.0 );
                    wo = wo * 8.0;
                    wi = wi * 8.0;

                    // convert back to cartesian coordinates
                    w.x = wr * sin(wo)*sin(wi);
                    w.y = wr * cos(wo);
                    w.z = wr * sin(wo)*cos(wi);
                    w += float3(1,0,0) * 0.1;
                }
                //return length(w);
                return distance(w,p);


                float r;
                int n = 0;
                float3 z = p;
                float Scale = 2.0;
                float3 Offset = float3(1,1,1);
                float Iterations = 30;
                while (n < Iterations) {
                    //if(z.x+z.y<0) z.xy = -z.yx; // fold 1
                    if(z.x+z.y<0) z.xy = abs(z.yx-Offset)-Offset*(Scale-1.); 
                    if(z.x+z.z<0) z.xz = -z.zx;
                    if(z.y+z.z<0) z.zy = -z.yz;
                    z = z*Scale - Offset*(Scale-1.0);
                    //float3 n1 = float3(1,0,0);
                    //z-=2.0 * min(0.0, dot(z, n1)) * n1;
                    n++;
                }
                //return length(z);
                return (length(z) ) * pow(Scale, -float(n));
            }

            float4x4 viewMatrix(float3 eye, float3 center, float3 up) {
                float3 f = normalize(center - eye);
                float3 s = normalize(cross(f, up));
                float3 u = cross(s, f);
                return float4x4(
                    float4(s, 0.0),
                    float4(u, 0.0),
                    float4(-f, 0.0),
                    float4(0.0, 0.0, 0.0, 1)
                );
            }

float DE(float3 pos) {
	float3 z = pos;
	float dr = 1.0;
	float r = 0.0;
    float Iterations = 30;
    float Bailout = 1.15;
    float Power = 8;

	for (int i = 0; i < Iterations ; i++) {
		r = length(z);
		if (r>Bailout) break;
		
		// convert to polar coordinates
		float theta = acos(z.z/r);
		float phi = atan2(z.y,z.x);
		dr =  pow( r, Power-1.0)*Power*dr + 1.0;
		
		// scale and rotate the point
		float zr = pow( r,Power);
		theta = theta*Power;
		phi = phi*Power;
		
		// convert back to cartesian coordinates
		z = zr*float3(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
		z+=pos;
	}
	return 0.5*log(r)*r/dr;
}


            float sceneSDF(float3 samplePoint) {
                //return sphereDist;
                float4x4 m = rotateY( sin(_Time.y+samplePoint.y) );
                float3 cubePoint = mul( m, float4(samplePoint, 1.0) ).xyz;


                //return huimaSDF( cubePoint );
                //return fractalSDF( cubePoint );
                float sphereDist = sphereSDF(samplePoint / float3(1.1,2,1.1)) * float3(1.1,2,1.1);
                float cubeDist = cubeSDF(cubePoint,float3(2,4,2));
                //float cubeDist = cubeSDF(cubePoint,float3(1,1,1));
                //return unionSDF(cubeDist, sphereDist);

                float planeDist = planeSDF(samplePoint,-1.5);
                return min( planeDist, fOpDifferenceRound(cubeDist, sphereDist, 0.2 + 0.2*sin(_Time.w)) );
                return fOpDifferenceRound(cubeDist, sphereDist, 0.2 + 0.2*sin(_Time.w));
                return fOpIntersectionColumns(cubeDist,sphereDist, 0.7 + 0.5*sin(_Time.w), 4);
                return fOpDifferenceStairs(cubeDist,sphereDist, 0.3 + 0.2*sin(_Time.w), 6);
                return fOpDifferenceChamfer(cubeDist, sphereDist, 0.3 + 0.2*sin(_Time.y));
                return differenceSDF(cubeDist, sphereDist);
            }

            float3 estimateNormal(float3 p) {
                return normalize(float3(
                    sceneSDF(float3(p.x + EPSILON, p.y, p.z)) - sceneSDF(float3(p.x - EPSILON, p.y, p.z)),
                    sceneSDF(float3(p.x, p.y + EPSILON, p.z)) - sceneSDF(float3(p.x, p.y - EPSILON, p.z)),
                    sceneSDF(float3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(float3(p.x, p.y, p.z - EPSILON))
                ));
            }


            float2 sdf( float3 eye, float3 viewRayDirection )
            {
                float start = MIN_DIST;
                float end = MAX_DIST;
                float depth = start;
                float dust = 0;
                for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
                    float3 pos = eye + depth * viewRayDirection;
                    float dist = sceneSDF(pos);
                    if (dist < EPSILON)
                        return float2(depth,dust);
                    depth += dist;
                    dust += snoise(pos) * dist;
                    if (depth >= end)
                        return float2(end,dust);
                }
                return end; 
            }

            float3 rayDirection(float fieldOfView, float2 size, float2 fragCoord) {
                float2 xy = fragCoord - size / 2.0;
                float z = size.y / tan(radians(fieldOfView) / 2.0);
                return normalize(float3(xy, -z));
            }


            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = v.vertex;
                o.vertex.y *= 4;
                o.vertex.y -= 1;
                o.vertex.x *= 2;
                o.uv = v.uv;
                return o;
            }

            sampler2D _MainTex;


            float shadow( float3 pos, float3 light )
            {
                float depth = 0;
                float maxDepth = distance( pos, light );
                float minDepth = MAX_MARCHING_STEPS;
                float3 viewRayDirection = normalize( light-pos );
                pos += viewRayDirection * EPSILON*2;
                for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
                    float dist = sceneSDF(pos + depth * viewRayDirection);
                    if (i>1 && dist < EPSILON) {
                        return dist/EPSILON;
                    }
                    minDepth = min(dist,minDepth);
                    depth += dist;
                    /*if (depth >= maxDepth) */
                    if ( depth >= (MAX_DIST - EPSILON) )
                        return 1.0;
                    
                }
                //return 1.0/min(1,minDepth);
                float s = min(1,minDepth);
                return s*s*(3.0-2.0*s);
            }

            float4 f( fixed4 f ) 
            {
                return f;
            }




            fixed4 frag (v2f i) : SV_Target
            {
                fixed4 tex = tex2D(_MainTex, i.uv);
                float4 col = float4(tex);
                col = float4(1,1,1,1);
                //col.rgb = sdf().xxx;
                float3 lightPos = float3(5,12,2);
                //float3 lightPos = float3(12*sin(_Time.w),12*cos(_Time.w),5+4*sin(_Time.y));

                float3 viewDir = rayDirection(45.0, _ScreenParams.xy, i.vertex.xy);
                float3 eye = float3(4.0, 2.0, 5.0);
                float4x4 viewToWorld = viewMatrix(eye, float3(0.0, 0.0, 0.0), float3(0.0, 1.0, 0.0));
                float3 worldDir = ( mul( float4(viewDir, 0.0), viewToWorld ) ).xyz;

                float2 depth = sdf( eye, worldDir );
                if ( depth.x > (MAX_DIST - EPSILON) )
                {
                    return float4(0,0,1,1);
                }
                float3 hitPos = eye + depth.x * worldDir;
                //depth.x = saturate( depth.x * 1.0 );
                float3 norm = estimateNormal(hitPos);

                col = float4( norm, 1 );
                float4 lighting = lit( dot(norm,lightPos), dot(norm,worldDir), 0 );
                //col = float4( col.xyz * 0.5 + float3(0.2,0.1,0.3) * lighting.y + float3(0.6,0.6,0.6) * lighting.z, 1 );
                col.xyz = 0.5 + float3(0.2,0.1,0.3) * lighting.y + float3(0.6,0.6,0.6) * lighting.z;
                col.xyz *= shadow( hitPos, lightPos );
                //col = float4( depth.x.xxx, 1 );

                col.xyz = saturate(depth.yyy);

                //col += f(tex);

                return fixed4(col.xyz,1);
            }
            ENDCG
        }
    }
}
