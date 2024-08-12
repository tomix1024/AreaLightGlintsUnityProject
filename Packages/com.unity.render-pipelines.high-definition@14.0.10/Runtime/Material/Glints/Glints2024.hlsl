
// Include Deliot & Belcour implementation
#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/Glints/Glints2023.hlsl"

float _OverrideDMax;
//float _LogSinSunAngle;
float _SunSolidAngle;
float _ZeroIfPgt1;

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
