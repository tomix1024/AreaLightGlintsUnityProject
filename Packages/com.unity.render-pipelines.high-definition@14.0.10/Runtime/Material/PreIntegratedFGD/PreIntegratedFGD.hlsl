#include "Packages/com.unity.render-pipelines.high-definition/Runtime/Material/PreIntegratedFGD/PreIntegratedFGD.cs.hlsl"

TEXTURE2D(_PreIntegratedFGD_GGXDisneyDiffuse);

// For image based lighting, a part of the BSDF is pre-integrated.
// This is done both for specular GGX height-correlated and DisneyDiffuse
// reflectivity is  Integral{(BSDF_GGX / F) - use for multiscattering
void GetPreIntegratedFGDGGXAndDisneyDiffuse(float NdotV, float perceptualRoughness, float3 fresnel0, out float3 specularFGD, out float diffuseFGD, out float reflectivity)
{
    // We want the LUT to contain the entire [0, 1] range, without losing half a texel at each side.
    float2 coordLUT = Remap01ToHalfTexelCoord(float2(sqrt(NdotV), perceptualRoughness), FGDTEXTURE_RESOLUTION);

    float3 preFGD = SAMPLE_TEXTURE2D_LOD(_PreIntegratedFGD_GGXDisneyDiffuse, s_linear_clamp_sampler, coordLUT, 0).xyz;

    // Pre-integrate GGX FGD
    // Integral{BSDF * <N,L> dw} =
    // Integral{(F0 + (1 - F0) * (1 - <V,H>)^5) * (BSDF / F) * <N,L> dw} =
    // (1 - F0) * Integral{(1 - <V,H>)^5 * (BSDF / F) * <N,L> dw} + F0 * Integral{(BSDF / F) * <N,L> dw}=
    // (1 - F0) * x + F0 * y = lerp(x, y, F0)
    specularFGD = lerp(preFGD.xxx, preFGD.yyy, fresnel0);

    // Pre integrate DisneyDiffuse FGD:
    // z = DisneyDiffuse
    // Remap from the [0, 1] to the [0.5, 1.5] range.
    diffuseFGD = preFGD.z + 0.5;

    reflectivity = preFGD.y;
}

void GetPreIntegratedFGDGGXAndLambert(float NdotV, float perceptualRoughness, float3 fresnel0, out float3 specularFGD, out float diffuseFGD, out float reflectivity)
{
    GetPreIntegratedFGDGGXAndDisneyDiffuse(NdotV, perceptualRoughness, fresnel0, specularFGD, diffuseFGD, reflectivity);
    diffuseFGD = 1.0;
}

TEXTURE2D(_PreIntegratedFGD_CharlieAndFabric);

void GetPreIntegratedFGDCharlieAndFabricLambert(float NdotV, float perceptualRoughness, float3 fresnel0, out float3 specularFGD, out float diffuseFGD, out float reflectivity)
{
    // Read the texture
    float3 preFGD = SAMPLE_TEXTURE2D_LOD(_PreIntegratedFGD_CharlieAndFabric, s_linear_clamp_sampler, float2(NdotV, perceptualRoughness), 0).xyz;

    specularFGD = lerp(preFGD.xxx, preFGD.yyy, fresnel0) * 2.0f * PI;

    // z = FabricLambert
    diffuseFGD = preFGD.z;

    reflectivity = preFGD.y;
}

TEXTURE2D(_PreIntegratedFGD_DGGXOnly);

void GetPreIntegratedFGDDGGXOnly(float NdotV, float perceptualRoughness, out float FGD)
{
    // Read the texture
    // Texture stores values in [0, 1], but represents values in [0, 2].
    float3 preFGD = 2*SAMPLE_TEXTURE2D_LOD(_PreIntegratedFGD_DGGXOnly, s_linear_clamp_sampler, float2(sqrt(NdotV), perceptualRoughness), 0).xyz;

    // Only single channel, no fresnelF0 interpolation!
    FGD = preFGD.x; // preFGD.y;
}
