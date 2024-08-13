
// Include Deliot & Belcour implementation
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Glints2023.hlsl"

float _OverrideDMax;
//float _LogSinSunAngle;
float _SunSolidAngle;
float _ZeroIfPgt1;
int _GlintReferenceLog2SampleCount;
int _GlintNDFIntegrationMode;

float ComputeTotalNDF(float roughness)
{
#if 1
    return 1.0 + pow(roughness, 0.5*2.821827587526826);
#else
    // https://www.wolframalpha.com/input?i=integrate+2*+a+%2F+%28%28x*a+-+x%29*x+%2B+1%29**2+dx+from+0+to+1
    // https://proofwiki.org/wiki/Arctangent_of_Imaginary_Number
    float a2 = Sq(roughness);
    float sqrt_term = sqrt(1-a2);
    // return 1 + 0.5*a2*log((1+sqrt_term)/(1-sqrt_term))/sqrt_term;
    return 1 + a2 * log((1 + sqrt_term) / roughness) / sqrt_term;
#endif
}

#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/GlintsSubdivision.hlsl"



float SampleGlints2024NDF(float3 halfwayTS, float LdotH, float roughness, float2 uv, float2 duvdx, float2 duvdy)
{
    float Dtarget = D_GGX(halfwayTS.z, roughness);
    float Dmax = D_GGX(1, roughness);

    // Approx. \int_hemisphere D(h) dh
    float R = _MicrofacetRoughness;
    float SD = ComputeTotalNDF(roughness);
    // Sun:
    // omega = 6.8e-5 sr.
    // sin(gamma) = 4.6e-3 == gamma
    //float sinGammaSq = exp(2*_LogSinSunAngle);
    //float cosGamma = sqrt(1 - sinGammaSq);
    //float Al = 2*PI * (1.0 - cosGamma);
    float Al = _SunSolidAngle;
    float Ah = Al / abs(4*LdotH);
    float newDmax = R * SD / Ah;
    Dmax = lerp(Dmax, newDmax, _OverrideDMax);
    if (Dtarget > Dmax)
        Dtarget *= 1-_ZeroIfPgt1;
    Dmax = max(Dmax, Dtarget);
    //Dtarget = min(Dtarget, Dmax);

    // Explicitly simulate smooth surface if _LogMicrofacetDensity < 0
    if (_LogMicrofacetDensity <= 0)
        return 1;

    float D = SampleGlints2023NDF(halfwayTS, Dtarget, Dmax, uv, duvdx, duvdy);
    return D;
}

float SampleGlints2024NDF_2023Wrapped(float3 halfwayTS, float LdotH, float roughness, float2 uv, float2 duvdx, float2 duvdy)
{
    // Wrapper for previous method with our interface...
    float Dtarget = D_GGX(halfwayTS.z, roughness);
    float Dmax = D_GGX(1, roughness);
	return SampleGlints2023NDF(halfwayTS, Dtarget, Dmax, uv, duvdx, duvdy);
}

float SampleGlints2024NDF_Optimized(float3 halfwayTS, float LdotH, float roughness, float2 uv, float2 duvdx, float2 duvdy)
{
    float Dtarget = D_GGX(halfwayTS.z, roughness);

    // Approx. \int_hemisphere D(h) dh
    float SD = ComputeTotalNDF(roughness);
    // Sun:
    // omega = 6.8e-5 sr.
    // sin(gamma) = 4.6e-3 == gamma
    //float sinGammaSq = exp(2*_LogSinSunAngle);
    //float cosGamma = sqrt(1 - sinGammaSq);
    //float Al = 2*PI * (1.0 - cosGamma);
    float Al = _SunSolidAngle;
    float Ah = Al / abs(4*LdotH);
    float Dmax = SD / Ah;
    Dmax = max(Dmax, Dtarget);

	float successProb = Dtarget / Dmax;
	float result = SampleGlints2023NDF_Internal(halfwayTS, successProb, uv, duvdx, duvdy);
	float D = result * Dmax;
    return D;
}



void IntegrateDGGXOnly_AreaRef(float NdotV, float4x3 lightVerts, float roughness, float3 fresnel0, out float ndfValue, out float3 bsdfValue, uint sampleCount = 512)
{
    float accNDF_LightSampling = 0;
    float accNDF_LightSampling2 = 0; // second moment
    float accNDF_BSDFSampling = 0;
    float accNDF_BSDFSampling2 = 0; // second moment

    float3 accBSDF_LightSampling = float3(0,0,0);
    float3 accBSDF_LightSampling2 = float3(0,0,0); // second moment
    float3 accBSDF_BSDFSampling = float3(0,0,0);
    float3 accBSDF_BSDFSampling2 = float3(0,0,0); // second moment

    // lightVerts are in local coordinate system!
    // see GetOrthoBasisViewNormal
    // X axis ~ view direction
    // Z axis ~ normal

    float3 V   = real3(sqrt(1 - NdotV * NdotV), 0, NdotV);
    float3x3 localToWorld = k_identity3x3;

    // coneNormals point inwards
    float4x3 coneNormals;
    coneNormals[0] = cross(lightVerts[0], lightVerts[1]);
    coneNormals[1] = cross(lightVerts[1], lightVerts[2]);
    coneNormals[2] = cross(lightVerts[2], lightVerts[3]);
    coneNormals[3] = cross(lightVerts[3], lightVerts[0]);


    // vxy
    float3 v00 = normalize(lightVerts[0]); //--
    float3 v01 = normalize(lightVerts[1]); //-+
    float3 v11 = normalize(lightVerts[2]); //++
    float3 v10 = normalize(lightVerts[3]); //+-

    // "Lower" triangle
    float3 e0x = v10 - v00;
    float3 e0y = v01 - v00;
    float3 n0 = normalize(cross(e0y, e0x));//normalize(v00 + v10 + v01);
    float A20 = length(cross(e0y, e0x));

    // "Upper" triangle
    float3 e1x = v01 - v11;
    float3 e1y = v10 - v11;
    float3 n1 = normalize(cross(e1y, e1x));//normalize(v11 + v10 + v01);
    float A21 = length(cross(e1y, e1x));


    // Cache for V term...
    float partLambdaV = GetSmithJointGGXPartLambdaV(NdotV, roughness);

    [loop]
    for (uint i = 0; i < sampleCount; ++i)
    {
        float2 u = Hammersley2d(i, sampleCount);

        // BSDF importance sampling
        {
            float3 L;
            float NdotL;
            float NdotH;
            float VdotH;
            SampleGGXDir(u, V, localToWorld, roughness, L, NdotL, NdotH, VdotH);
            if (NdotL > 0.0)
            {
                // step(a, b) == a <= b ? 1 : 0
                float ndfWeightOverPdf = rcp(NdotH);

                float V = V_SmithJointGGX(NdotL, NdotV, roughness, partLambdaV);
                //float V = V_SmithJointGGX(NdotL, NdotV, roughness);
                float3 F = F_Schlick(fresnel0, VdotH);
                float3 bsdfWeightOverPdf = F * (4.0 * V * NdotL * VdotH / NdotH);

                // coneNormals point inwards
                float mask = 1.0;
                mask *= (dot(L, coneNormals[0]) > 0) ? 1 : 0;
                mask *= (dot(L, coneNormals[1]) > 0) ? 1 : 0;
                mask *= (dot(L, coneNormals[2]) > 0) ? 1 : 0;
                mask *= (dot(L, coneNormals[3]) > 0) ? 1 : 0;
                ndfWeightOverPdf *= mask;
                bsdfWeightOverPdf *= mask;

                accNDF_BSDFSampling += ndfWeightOverPdf;
                accNDF_BSDFSampling2 += ndfWeightOverPdf*ndfWeightOverPdf;

                accBSDF_BSDFSampling += bsdfWeightOverPdf;
                accBSDF_BSDFSampling2 += bsdfWeightOverPdf*bsdfWeightOverPdf;
            }
        }

        // Light source importances sampling
        {
            float3 v_sampled;
            float3 v_normal;
            float v_A2;
            if (u.x + u.y < 1)
            {
                // v00, v01, v10
                v_sampled = v00 + u.x * (v01-v00) + u.y * (v10-v00);
                v_normal = n0;
                v_A2 = A20;
            }
            else
            {
                // v11, v10, v01
                u = 1-u;
                v_sampled = v11 + u.x * (v01-v11) + u.y * (v10-v11);
                v_normal = n1;
                v_A2 = A21;
            }
            // \int_\Omega \cos\theta_i \dx{\omega} = \int_A \cos\theta_i \cos\theta_o / \|v-x\|^2 \dx{v}
            //v_weight = area * v_cosTheta / dot(v_sampled, v_sampled);
            float v_weight = v_A2 / dot(v_sampled, v_sampled);
            v_sampled = normalize(v_sampled);
            v_weight *= abs(dot(v_sampled, v_normal));

            float3 L = v_sampled;
            float3 H = normalize(V+L);

            float VdotH = dot(V, H);
            float NdotH = H.z;
            float NdotL = L.z;

            float D = D_GGX(NdotH, roughness);
            float V = V_SmithJointGGX(NdotL, NdotV, roughness, partLambdaV);
            //float V = V_SmithJointGGX(NdotL, NdotV, roughness);
            float3 F = F_Schlick(fresnel0, VdotH);

            float ndfWeightOverPdf = abs(v_weight) * D / abs(4*VdotH);
            float3 brdfWeightOverPdf = F * (abs(v_weight) * D * V * max(0, NdotL));

            accNDF_LightSampling += ndfWeightOverPdf;
            accNDF_LightSampling2 += ndfWeightOverPdf*ndfWeightOverPdf;

            accBSDF_LightSampling += brdfWeightOverPdf;
            accBSDF_LightSampling2 += brdfWeightOverPdf*brdfWeightOverPdf;
        }
    }
    accNDF_LightSampling /= sampleCount;
    accNDF_BSDFSampling /= sampleCount;

    accBSDF_LightSampling /= sampleCount;
    accBSDF_BSDFSampling /= sampleCount;

    // roughness = 0 -> bsdf sampling
    // roughness = 1 -> light sampling
    // also for smaller roughnesses -> light sampling

    ndfValue = lerp(accNDF_BSDFSampling, accNDF_LightSampling, pow(roughness, 0.1));
    bsdfValue = lerp(accBSDF_BSDFSampling, accBSDF_LightSampling, pow(roughness, 0.1));
}

float SampleGlints2024NDF_Area(float3 halfwayTS, float LdotH, float roughness, float integratedNDF, float4x3 lightVerts, float2 uv, float2 duvdx, float2 duvdy)
{
    // integratedNDF = \int_light D(H)/(4*LdotH) dL

    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);

    if (integratedNDF > totalNDF)
        integratedNDF *= 1-_ZeroIfPgt1;
    float p = saturate(integratedNDF / totalNDF); // = R*Dtarget/Dmax

    // Skip computation below if probability is zero!
    if (p == 0)
        return 0;

    // Explicitly simulate smooth surface if _LogMicrofacetDensity < 0
    if (_LogMicrofacetDensity <= 0)
        return 1;

    // Division by microfacet count handled internally.
    float D = SampleGlints2023NDF_Internal(halfwayTS, p, uv, duvdx, duvdy) / p;
    return D;
}

float SampleGlints2024NDF_Area_2023Wrapped(float3 halfwayTS, float LdotH, float roughness, float integratedNDF, float4x3 lightVerts, float2 uv, float2 duvdx, float2 duvdy)
{
    // Wrapper for previous method with our interface...
    float Dtarget = D_GGX(halfwayTS.z, roughness);
    float Dmax = D_GGX(1, roughness);
	return SampleGlints2023NDF(halfwayTS, Dtarget, Dmax, uv, duvdx, duvdy);
}

float SampleGlints2024NDF_Area_Optimized(float3 halfwayTS, float LdotH, float roughness, float integratedNDF, float4x3 lightVerts, float2 uv, float2 duvdx, float2 duvdy)
{
    // integratedNDF = \int_light D(H)/(4*LdotH) dL

    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);

    float p = saturate(integratedNDF / totalNDF); // = R*Dtarget/Dmax

    // Skip computation below if probability is zero!
    if (p == 0)
        return 0;

    // Division by microfacet count handled internally.
    float D = SampleGlints2023NDF_Internal(halfwayTS, p, uv, duvdx, duvdy) / p;
    return D;
}
