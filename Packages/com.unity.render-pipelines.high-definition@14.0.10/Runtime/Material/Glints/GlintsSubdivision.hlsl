#ifndef GLINTS_SUBDIVISION_HLSL
#define GLINTS_SUBDIVISION_HLSL

#include "Packages/com.unity.render-pipelines.core/ShaderLibrary/AreaLighting.hlsl"

#define GLINT_SUBDIVISION_COUNT_X 4
#define GLINT_SUBDIVISION_COUNT_Y 4
#define GLINT_SUBDIVISION_COUNT (GLINT_SUBDIVISION_COUNT_X + GLINT_SUBDIVISION_COUNT_Y)
#define GLINT_SUBDIVISION_CELL_COUNT (1<<(GLINT_SUBDIVISION_COUNT))
#define GLINT_SUBDIVISION_CELL_COUNT_X (1<<(GLINT_SUBDIVISION_COUNT_X))
#define GLINT_SUBDIVISION_CELL_COUNT_Y (1<<(GLINT_SUBDIVISION_COUNT_Y))
#define GLINT_SUBDIVISION_MIPMAP_CELL_COUNT ((1<<(GLINT_SUBDIVISION_COUNT+1)) - 1) // sum_i=0^GLINT_SUBDIVISION_COUNT 2**i
#define GLINT_SUBDIVISION_MIPMAP_SIZE(mip) (1<<(mip))
#define GLINT_SUBDIVISION_MIPMAP_OFFSET(mip) ((1<<(mip)) - 1)

float UintToNormalizedFloat(uint x)
{
    return asfloat(0x3f800000 | (x >> 9)) - 1.0f;
}

float RandomFloat01(inout uint state)
{
    return UintToNormalizedFloat(XorShift(state));
}

void PolygonIrradiance_Subdivided(float4x3 lightVerts, float3x3 ltcTransform, out float ltcValue_subdivided[GLINT_SUBDIVISION_CELL_COUNT])
{
    // float4x3 LS = mul(lightVerts, ltcTransform);
    // TODO APPLY SUBDIVISION IN WORLD SPACE, BEFORE LTC TRANSFORM!

    // vxy
    float3 v00 = normalize(lightVerts[0]); //--
    float3 v01 = normalize(lightVerts[1]); //-+
    float3 v11 = normalize(lightVerts[2]); //++
    float3 v10 = normalize(lightVerts[3]); //+-
    float cosTheta0 = dot(v00, v01);
    float cosTheta1 = dot(v10, v11);
    float sinTheta0 = sqrt(1 - cosTheta0*cosTheta0);
    float sinTheta1 = sqrt(1 - cosTheta1*cosTheta1);
    float theta0 = acos(cosTheta0);
    float theta1 = acos(cosTheta1);
    for (int y = 0; y < GLINT_SUBDIVISION_CELL_COUNT_Y; ++y)
    {
        // Slerp along y axis
        float tya = float(y) / float(GLINT_SUBDIVISION_CELL_COUNT_Y);
        float tyb = float(y+1) / float(GLINT_SUBDIVISION_CELL_COUNT_Y);
        float3 v0a = normalize(sin((1-tya)*theta0)/sinTheta0 * v00 + sin(tya*theta0)/sinTheta0 * v01);
        float3 v0b = normalize(sin((1-tyb)*theta0)/sinTheta0 * v00 + sin(tyb*theta0)/sinTheta0 * v01);
        float3 v1a = normalize(sin((1-tya)*theta1)/sinTheta1 * v10 + sin(tya*theta1)/sinTheta1 * v11);
        float3 v1b = normalize(sin((1-tyb)*theta1)/sinTheta1 * v10 + sin(tyb*theta1)/sinTheta1 * v11);

        float cosTheta_a = dot(v0a, v1a);
        float cosTheta_b = dot(v0b, v1b);
        float sinTheta_a = sqrt(1 - cosTheta_a*cosTheta_a);
        float sinTheta_b = sqrt(1 - cosTheta_b*cosTheta_b);
        float theta_a = acos(cosTheta_a);
        float theta_b = acos(cosTheta_b);
        for (int x = 0; x < GLINT_SUBDIVISION_CELL_COUNT_X; ++x)
        {
            // Slerp along x axis
            float txa = float(x) / float(GLINT_SUBDIVISION_CELL_COUNT_X);
            float txb = float(x+1) / float(GLINT_SUBDIVISION_CELL_COUNT_X);
            float3 vaa = normalize(sin((1-txa)*theta_a)/sinTheta_a * v0a + sin(txa*theta_a)/sinTheta_a * v1a);
            float3 vba = normalize(sin((1-txb)*theta_a)/sinTheta_a * v0a + sin(txb*theta_a)/sinTheta_a * v1a);
            float3 vab = normalize(sin((1-txa)*theta_b)/sinTheta_b * v0b + sin(txa*theta_b)/sinTheta_b * v1b);
            float3 vbb = normalize(sin((1-txb)*theta_b)/sinTheta_b * v0b + sin(txb*theta_b)/sinTheta_b * v1b);

            float4x3 lightVerts_cell;
            lightVerts_cell[0] = vaa;
            lightVerts_cell[1] = vab;
            lightVerts_cell[2] = vbb;
            lightVerts_cell[3] = vba;

            float4x3 LS_cell = mul(lightVerts_cell, ltcTransform);
            float ltcValue_cell = PolygonIrradiance(LS_cell);
            ltcValue_subdivided[y*GLINT_SUBDIVISION_CELL_COUNT_X + x] = ltcValue_cell;
        }
    }
}

float GenerateBinomialValue_DualGated(inout uint rng, float leftProb, float microfacetCount)
{
#if 1
    float result = 0;
    for (int i = 0; i < microfacetCount; ++i)
    {
        float rand = RandomFloat01(rng);
        result += rand < leftProb;
    }
    return result;
#else
    // Apparently this can happen and it leads to bad visual artefacts :/
    //if (microfacetCount < 1)
    //    return leftProb * microfacetCount;
    float rand0 = RandomFloat01(rng); // uniform in [0, 1]
    float rand1 = RandomFloat01(rng); // uniform in [0, 1]
    float rand2 = RandomFloat01(rng); // uniform in [0, 1]

    // microfacetCount "should" be integer-valued
    if (microfacetCount < 5)
    {
        float result = 0;//leftProb * frac(microfacetCount);
        if (microfacetCount >= 1)
            result += (rand0 < leftProb);
        if (microfacetCount >= 2)
            result += (rand1 < leftProb);
        if (microfacetCount >= 3)
            result += (rand2 < leftProb);
        if (microfacetCount >= 4)
        {
            float rand3 = RandomFloat01(rng); // uniform in [0, 1]
            result += (rand3 < leftProb);
        }
        return result;
    }
	// Compute binomial properties
	float probOneHitLeft = 1 - pow(1 - leftProb, microfacetCount); // Probability of hitting at least one microfacet in left node.
    float probOneHitRight = 1 - pow(leftProb, microfacetCount-1); // Probability of hitting at leat one microfacet in the right node (given that we have at least one hit in the left node).
	float leftCountMean = 1 + (microfacetCount - 2) * leftProb; // Expected value of hits in the left node given that we have at least one hit in the left, AND one hit in the right node.
	float leftCountVar = (microfacetCount - 2.0) * leftProb * (1.0 - leftProb); // Standard deviation of number of hits in the footprint given already one hit
    float leftCountStdDev = sqrt(max(leftCountVar, 0));
    float binomialSmoothWidthLeft = 0;  //0.1 * clamp(probOneHitLeft * 10, 0.0, 1.0) * clamp((1.0 - probOneHitLeft) * 10, 0.0, 1.0);
    float binomialSmoothWidthRight = 0; //0.1 * clamp(probOneHitRight * 10, 0.0, 1.0) * clamp((1.0 - probOneHitRight) * 10, 0.0, 1.0);

    // Compute gating on both sides!
    float randGatingLeft = rand0; // uniform in [0, 1]
    float randGatingRight = rand1; // uniform in [0, 1]
    float gatingLeft = 1.0; // at least one left.
    float gatingRight = 1.0; // at least one right.
    if (binomialSmoothWidthLeft > 0.0000001)
    {
        gatingLeft = saturate(RemapTo01(randGatingLeft, probOneHitLeft + binomialSmoothWidthLeft, probOneHitLeft - binomialSmoothWidthLeft));
    }
    else
    {
        gatingLeft = randGatingLeft < probOneHitLeft;
    }
    if (binomialSmoothWidthRight > 0.0000001)
    {
        gatingRight = saturate(RemapTo01(randGatingRight, probOneHitRight + binomialSmoothWidthRight, probOneHitRight - binomialSmoothWidthRight));
    }
    else
    {
        gatingRight = randGatingRight > probOneHitRight;
    }

    // Transform to gaussian with mean + stddev.
    float leftCountGauss = sampleNormalDistribution(rand2, leftCountMean, leftCountStdDev);
    // NOTE: We have adjusted the mean such that we don't need to adjust the clamp limits!
    // YES WE DO NEED TO ADJUST CLAMP LIMITS! at least one microfacet is left and at least one microfacet is right!
    leftCountGauss = clamp(leftCountGauss, 1, microfacetCount-1);
    // Combine the three gating cases...
    float result = gatingLeft * lerp(microfacetCount, leftCountGauss, gatingRight);
    return round(result);
#endif
}

void ApplySubdivisionOnLevel(int level, inout float result_subdivided[GLINT_SUBDIVISION_CELL_COUNT], in float targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_CELL_COUNT], inout uint rngState)
{
    int coarse_level = level;
    int fine_level = level+1;
    // Count backwards! -> inplace computation!
    for (int i = GLINT_SUBDIVISION_MIPMAP_SIZE(coarse_level)-1; i >= 0; --i)
    {
        float parent_N = result_subdivided[/*GLINT_SUBDIVISION_MIPMAP_OFFSET(coarse_level) +*/ i];

        // Randomly swap left and right children in `GenerateBinomialValue_DualGated`.
        // For small microfacetCounts(<1) there seems to be some kind of bias...
        int swapLR = (RandomFloat01(rngState) < 0.5) ? 1 : 0;

        float pdfCombined = targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_OFFSET(coarse_level) + i];
        float pdfLeft = targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) + 2*i + swapLR];
        float pdfRight = targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) + 2*i+1 - swapLR];

        float relativeProbLeft = pdfLeft / pdfCombined; // = pdfLeft / (pdfLeft + pdfRight);

        float left_N = GenerateBinomialValue_DualGated(rngState, relativeProbLeft, parent_N);
        float right_N = parent_N - left_N;

        result_subdivided[/*GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) +*/ 2*i + swapLR] = left_N;
        result_subdivided[/*GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) +*/ 2*i+1 - swapLR] = right_N;
    }
}

void ApplySubdivisionOnLevel(int level, inout float4 result_subdivided[GLINT_SUBDIVISION_CELL_COUNT], in float targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_CELL_COUNT], inout uint rngState)
{
    int coarse_level = level;
    int fine_level = level + 1;
    // Count backwards! -> inplace computation!
    for (int i = GLINT_SUBDIVISION_MIPMAP_SIZE(coarse_level) - 1; i >= 0; --i)
    {
        float4 parent_N = result_subdivided[/*GLINT_SUBDIVISION_MIPMAP_OFFSET(coarse_level) +*/ i];

        // Randomly swap left and right children in `GenerateBinomialValue_DualGated`.
        // For small microfacetCounts(<1) there seems to be some kind of bias...
        int swapLR = (RandomFloat01(rngState) < 0.5) ? 1 : 0;

        float pdfCombined = targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_OFFSET(coarse_level) + i];
        float pdfLeft = targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) + 2 * i + swapLR];
        float pdfRight = targetProb_mipmap[GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) + 2 * i + 1 - swapLR];

        float relativeProbLeft = pdfLeft / pdfCombined; // = pdfLeft / (pdfLeft + pdfRight);

        float4 left_N = float4(
            GenerateBinomialValue_DualGated(rngState, relativeProbLeft, parent_N.x),
            GenerateBinomialValue_DualGated(rngState, relativeProbLeft, parent_N.y),
            GenerateBinomialValue_DualGated(rngState, relativeProbLeft, parent_N.z),
            GenerateBinomialValue_DualGated(rngState, relativeProbLeft, parent_N.w)
        );
        float4 right_N = parent_N - left_N;

        result_subdivided[/*GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) +*/ 2 * i + swapLR] = left_N;
        result_subdivided[/*GLINT_SUBDIVISION_MIPMAP_OFFSET(fine_level) +*/ 2 * i + 1 - swapLR] = right_N;
    }
}

void GenerateAngularBinomialValueForSurfaceCell_Subdivided(uint rngStateSubdiv, float4 randB, float4 randG, float2 slopeLerp, float footprintOneHitProba, float binomialSmoothWidth, float footprintMean, float footprintSTD, float microfacetCount, float targetNDF_subdivided[GLINT_SUBDIVISION_MIPMAP_CELL_COUNT], out float result[GLINT_SUBDIVISION_CELL_COUNT])
{
    float4 gating;
    if (binomialSmoothWidth > 0.0000001)
        gating = saturate(RemapTo01(randB, footprintOneHitProba + binomialSmoothWidth, footprintOneHitProba - binomialSmoothWidth));
    else
        gating = randB < footprintOneHitProba;

    float4 gauss = randG * footprintSTD + footprintMean;
    // The microfacet count should be reduced by 1 after the gating!
    float fixedMicrofacetCount = lerp(microfacetCount, max(0, microfacetCount - 1), _FixSampledMicrofacetCount);
    gauss = clamp(gauss, 0, fixedMicrofacetCount);
    float4 results = gating * (1.0 + gauss);
    // For subdivision experiment we need to round to integer here!
    // In practice rounding leads to artifacts in some situations!
    results = round(results);

    // Subdivide!
    float4 results_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    results_subdivided[0] = results;
    // Subdivide...
    for (int level = 0; level < GLINT_SUBDIVISION_COUNT; ++level)
    {
        ApplySubdivisionOnLevel(level, results_subdivided, targetNDF_subdivided, rngStateSubdiv); // TODO separate RNG state for each thing!
    }

    for (int i = 0; i < GLINT_SUBDIVISION_CELL_COUNT; ++i)
    {
        result[i] = BilinearLerp(results_subdivided[i], slopeLerp);
    }
}


void SampleGlintGridSimplex_Subdivided(float2 uv, uint gridSeed, float2 slope, float footprintArea, float targetNDF_subdivided[GLINT_SUBDIVISION_MIPMAP_CELL_COUNT], float gridWeight, out float result[GLINT_SUBDIVISION_CELL_COUNT])
{
    // Get surface space glint simplex grid cell
    const float2x2 gridToSkewedGrid = float2x2(1.0, -0.57735027, 0.0, 1.15470054);
    float2 skewedCoord = mul(gridToSkewedGrid, uv);
    int2 baseId = int2(floor(skewedCoord));
    float3 temp = float3(frac(skewedCoord), 0.0);
    temp.z = 1.0 - temp.x - temp.y;
    float s = step(0.0, -temp.z);
    float s2 = 2.0 * s - 1.0;
    int2 glint0 = baseId + int2(s, s);
    int2 glint1 = baseId + int2(s, 1.0 - s);
    int2 glint2 = baseId + int2(1.0 - s, s);
    float3 barycentrics = float3(-temp.z * s2, s - temp.y * s2, s - temp.x * s2);

    // Generate per surface cell random numbers
    float3 rand0 = pcg3dFloat(uint3(glint0 + 2147483648, gridSeed)); // TODO : optimize away manual seeds
    float3 rand1 = pcg3dFloat(uint3(glint1 + 2147483648, gridSeed));
    float3 rand2 = pcg3dFloat(uint3(glint2 + 2147483648, gridSeed));

    // rng state during subdivision...
    uint rngState0 = JenkinsHash(uint3(glint0 + 2147483648, gridSeed));
    uint rngState1 = JenkinsHash(uint3(glint1 + 2147483648, gridSeed));
    uint rngState2 = JenkinsHash(uint3(glint2 + 2147483648, gridSeed));

    // Get per surface cell per slope cell random numbers
    float4 rand0SlopesB, rand1SlopesB, rand2SlopesB, rand0SlopesG, rand1SlopesG, rand2SlopesG;
    float2 slopeLerp0, slopeLerp1, slopeLerp2;
    CustomRand4Texture(slope, rand0.yz, rand0SlopesB, rand0SlopesG, slopeLerp0);
    CustomRand4Texture(slope, rand1.yz, rand1SlopesB, rand1SlopesG, slopeLerp1);
    CustomRand4Texture(slope, rand2.yz, rand2SlopesB, rand2SlopesG, slopeLerp2);

    // Compute microfacet count with randomization
    float3 logDensityRand = clamp(sampleNormalDistribution(float3(rand0.x, rand1.x, rand2.x), _LogMicrofacetDensity, _DensityRandomization), 0.0, 50.0); // TODO : optimize sampleNormalDist
    float3 microfacetCount = max(0.0.rrr, footprintArea.rrr * exp(logDensityRand));
    float3 microfacetCountBlended = microfacetCount * gridWeight;

    // Compute binomial properties
    float hitProba = targetNDF_subdivided[0]; // probability of hitting desired half vector in NDF distribution for ANY of the subdivision cells!
    float3 footprintOneHitProba = (1.0 - pow(1.0 - hitProba.rrr, microfacetCountBlended)); // probability of hitting at least one microfacet in footprint
    float3 footprintMean = (microfacetCountBlended - 1.0) * hitProba.rrr; // Expected value of number of hits in the footprint given already one hit
    float3 footprintSTD = sqrt((microfacetCountBlended - 1.0) * hitProba.rrr * (1.0 - hitProba.rrr)); // Standard deviation of number of hits in the footprint given already one hit
    float3 binomialSmoothWidth = 0.1 * clamp(footprintOneHitProba * 10, 0.0, 1.0) * clamp((1.0 - footprintOneHitProba) * 10, 0.0, 1.0);

    // This does the trick! It's even continuous around N=1
    footprintOneHitProba = lerp(footprintOneHitProba, hitProba * microfacetCountBlended, microfacetCountBlended < 1);

    // Generate numbers of reflecting microfacets
    float result0_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    float result1_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    float result2_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
#if 1
    // Apply subdivision to result{0,1,2} separately here!
    GenerateAngularBinomialValueForSurfaceCell_Subdivided(rngState0, rand0SlopesB, rand0SlopesG, slopeLerp0, footprintOneHitProba.x, binomialSmoothWidth.x, footprintMean.x, footprintSTD.x, microfacetCountBlended.x, targetNDF_subdivided, result0_subdivided);
    GenerateAngularBinomialValueForSurfaceCell_Subdivided(rngState1, rand1SlopesB, rand1SlopesG, slopeLerp1, footprintOneHitProba.y, binomialSmoothWidth.y, footprintMean.y, footprintSTD.y, microfacetCountBlended.y, targetNDF_subdivided, result1_subdivided);
    GenerateAngularBinomialValueForSurfaceCell_Subdivided(rngState2, rand2SlopesB, rand2SlopesG, slopeLerp2, footprintOneHitProba.z, binomialSmoothWidth.z, footprintMean.z, footprintSTD.z, microfacetCountBlended.z, targetNDF_subdivided, result2_subdivided);
#else
    // Applying subdivision after any blending increases error...
    float result0, result1, result2;
    result0 = GenerateAngularBinomialValueForSurfaceCell(rand0SlopesB, rand0SlopesG, slopeLerp0, footprintOneHitProba.x, binomialSmoothWidth.x, footprintMean.x, footprintSTD.x, microfacetCountBlended.x);
    result1 = GenerateAngularBinomialValueForSurfaceCell(rand1SlopesB, rand1SlopesG, slopeLerp1, footprintOneHitProba.y, binomialSmoothWidth.y, footprintMean.y, footprintSTD.y, microfacetCountBlended.y);
    result2 = GenerateAngularBinomialValueForSurfaceCell(rand2SlopesB, rand2SlopesG, slopeLerp2, footprintOneHitProba.z, binomialSmoothWidth.z, footprintMean.z, footprintSTD.z, microfacetCountBlended.z);

    // Initialize with "full" area light.
    result0_subdivided[0] = result0;
    result1_subdivided[0] = result1;
    result2_subdivided[0] = result2;

    // Subdivide...
    for (int level = 0; level < GLINT_SUBDIVISION_COUNT; ++level)
    {
        ApplySubdivisionOnLevel(level, result0_subdivided, targetNDF_subdivided, rngState0);
        ApplySubdivisionOnLevel(level, result1_subdivided, targetNDF_subdivided, rngState1);
        ApplySubdivisionOnLevel(level, result2_subdivided, targetNDF_subdivided, rngState2);
    }
#endif

    for (int i = 0; i < GLINT_SUBDIVISION_CELL_COUNT; ++i)
    {
        // Interpolate result for glint grid cell
        float3 results = float3(result0_subdivided[i], result1_subdivided[i], result2_subdivided[i]) / microfacetCount.xyz;
        result[i] = dot(results, barycentrics);
    }
}

void SampleGlints2023NDF_Subdivided(float3 localHalfVector, float successProb_subdivided[GLINT_SUBDIVISION_MIPMAP_CELL_COUNT], float2 uv, float2 duvdx, float2 duvdy, out float D_subdivided[GLINT_SUBDIVISION_CELL_COUNT])
{
    // ACCURATE PIXEL FOOTPRINT ELLIPSE
    float2 ellipseMajorAxis, ellipseMinorAxis;
    float ellipseMajorLength, ellipseMinorLength;
    GetGradientEllipseNew(duvdx, duvdy, ellipseMajorAxis, ellipseMinorAxis, ellipseMajorLength, ellipseMinorLength);
    float ellipseRatio = ellipseMajorLength / ellipseMinorLength;

    // SHARED GLINT NDF VALUES
    //float halfScreenSpaceScaler = _ScreenSpaceScale * 0.5;
    //float footprintArea = ellipseMajorLength * halfScreenSpaceScaler * ellipseMinorLength * halfScreenSpaceScaler * 4.0;
    float2 slope = localHalfVector.xy; // Orthogrtaphic slope projected grid

    // MANUAL LOD COMPENSATION
    float lod = log2(ellipseMinorLength * _ScreenSpaceScale);
    float lod0 = floor(lod); //lod >= 0.0 ? (int)(lod) : (int)(lod - 1.0);
    float lod1 = lod0 + 1;
    float divLod0 = exp2(lod0);
    float divLod1 = exp2(lod1);
    float lodLerp = frac(lod);
    float footprintAreaLOD0 = divLod0*divLod0;// = pow(exp2(lod0), 2.0);
    float footprintAreaLOD1 = divLod1*divLod1;// = pow(exp2(lod1), 2.0);

    // MANUAL ANISOTROPY RATIO COMPENSATION
    float ratio0 = max(pow(2.0, (int)log2(ellipseRatio)), 1.0);
    float ratio1 = ratio0 * 2.0;
    float ratioLerp = clamp(Remap(ellipseRatio, ratio0, ratio1, 0.0, 1.0), 0.0, 1.0);

    // MANUAL ANISOTROPY ROTATION COMPENSATION
    float2 v1 = float2(0.0, 1.0);
    float2 v2 = ellipseMajorAxis;
    float theta = atan2(v1.x * v2.y - v1.y * v2.x, v1.x * v2.x + v1.y * v2.y) * RAD2DEG;
    float thetaGrid = 90.0 / max(ratio0, 2.0);
    float thetaBin = (int)(theta / thetaGrid) * thetaGrid;
    thetaBin = thetaBin + (thetaGrid / 2.0);
    float thetaBin0 = theta < thetaBin ? thetaBin - thetaGrid / 2.0 : thetaBin;
    float thetaBinH = thetaBin0 + thetaGrid / 4.0;
    float thetaBin1 = thetaBin0 + thetaGrid / 2.0;
    float thetaBinLerp = Remap(theta, thetaBin0, thetaBin1, 0.0, 1.0);
    thetaBin0 = thetaBin0 <= 0.0 ? 180.0 + thetaBin0 : thetaBin0;

    // TETRAHEDRONIZATION OF ROTATION + RATIO + LOD GRID
    bool centerSpecialCase = (ratio0 == 1.0);
    float2 divLods = float2(divLod0, divLod1);
    float2 footprintAreas = float2(footprintAreaLOD0, footprintAreaLOD1);
    float2 ratios = float2(ratio0, ratio1);
    float4 thetaBins = float4(thetaBin0, thetaBinH, thetaBin1, 0.0); // added 0.0 for center singularity case
    float3 tetraA, tetraB, tetraC, tetraD;
    GetAnisoCorrectingGridTetrahedron(centerSpecialCase, thetaBinLerp, ratioLerp, lodLerp, tetraA, tetraB, tetraC, tetraD);
    if (centerSpecialCase == true) // Account for center singularity in barycentric computation
        thetaBinLerp = Remap01To(thetaBinLerp, 0.0, ratioLerp);
    float4 tetraBarycentricWeights = GetBarycentricWeightsTetrahedron(float3(thetaBinLerp, ratioLerp, lodLerp), tetraA, tetraB, tetraC, tetraD); // Compute barycentric coordinates within chosen tetrahedron

    // PREPARE NEEDED ROTATIONS
    tetraA.x *= 2; tetraB.x *= 2; tetraC.x *= 2; tetraD.x *= 2;
    if (centerSpecialCase == true) // Account for center singularity (if center vertex => no rotation)
    {
        tetraA.x = (tetraA.y == 0) ? 3 : tetraA.x;
        tetraB.x = (tetraB.y == 0) ? 3 : tetraB.x;
        tetraC.x = (tetraC.y == 0) ? 3 : tetraC.x;
        tetraD.x = (tetraD.y == 0) ? 3 : tetraD.x;
    }
    float2 uvRotA = RotateUV(uv, thetaBins[tetraA.x] * DEG2RAD, 0.0.rr);
    float2 uvRotB = RotateUV(uv, thetaBins[tetraB.x] * DEG2RAD, 0.0.rr);
    float2 uvRotC = RotateUV(uv, thetaBins[tetraC.x] * DEG2RAD, 0.0.rr);
    float2 uvRotD = RotateUV(uv, thetaBins[tetraD.x] * DEG2RAD, 0.0.rr);

    // SAMPLE GLINT GRIDS
    uint gridSeedA = HashWithoutSine13(float3(log2(divLods[tetraA.z]), thetaBins[tetraA.x] % 360, ratios[tetraA.y])) * 4294967296.0;
    uint gridSeedB = HashWithoutSine13(float3(log2(divLods[tetraB.z]), thetaBins[tetraB.x] % 360, ratios[tetraB.y])) * 4294967296.0;
    uint gridSeedC = HashWithoutSine13(float3(log2(divLods[tetraC.z]), thetaBins[tetraC.x] % 360, ratios[tetraC.y])) * 4294967296.0;
    uint gridSeedD = HashWithoutSine13(float3(log2(divLods[tetraD.z]), thetaBins[tetraD.x] % 360, ratios[tetraD.y])) * 4294967296.0;
    float sampleA_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    SampleGlintGridSimplex_Subdivided(uvRotA / divLods[tetraA.z] / float2(1.0, ratios[tetraA.y]), gridSeedA, slope, ratios[tetraA.y] * footprintAreas[tetraA.z], successProb_subdivided, tetraBarycentricWeights.x, sampleA_subdivided);
    float sampleB_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    SampleGlintGridSimplex_Subdivided(uvRotB / divLods[tetraB.z] / float2(1.0, ratios[tetraB.y]), gridSeedB, slope, ratios[tetraB.y] * footprintAreas[tetraB.z], successProb_subdivided, tetraBarycentricWeights.y, sampleB_subdivided);
    float sampleC_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    SampleGlintGridSimplex_Subdivided(uvRotC / divLods[tetraC.z] / float2(1.0, ratios[tetraC.y]), gridSeedC, slope, ratios[tetraC.y] * footprintAreas[tetraC.z], successProb_subdivided, tetraBarycentricWeights.z, sampleC_subdivided);
    float sampleD_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    SampleGlintGridSimplex_Subdivided(uvRotD / divLods[tetraD.z] / float2(1.0, ratios[tetraD.y]), gridSeedD, slope, ratios[tetraD.y] * footprintAreas[tetraD.z], successProb_subdivided, tetraBarycentricWeights.w, sampleD_subdivided);
    for (int i = 0; i < GLINT_SUBDIVISION_CELL_COUNT; ++i)
    {
        D_subdivided[i] = (sampleA_subdivided[i] + sampleB_subdivided[i] + sampleC_subdivided[i] + sampleD_subdivided[i]);
    }
}

void BuildMipChain(in float input[GLINT_SUBDIVISION_CELL_COUNT], out float output[GLINT_SUBDIVISION_MIPMAP_CELL_COUNT])
{
    for (int i = 0; i < GLINT_SUBDIVISION_CELL_COUNT; ++i)
        output[GLINT_SUBDIVISION_MIPMAP_OFFSET(GLINT_SUBDIVISION_COUNT) + i] = input[i];
    for (int level = GLINT_SUBDIVISION_COUNT; level > 0; --level)
    {
        for (int i = 0; i < GLINT_SUBDIVISION_MIPMAP_SIZE(level-1); ++i)
        {
            float left_value = output[GLINT_SUBDIVISION_MIPMAP_OFFSET(level) + 2*i];
            float right_value = output[GLINT_SUBDIVISION_MIPMAP_OFFSET(level) + 2*i+1];
            float value = left_value + right_value;
            output[GLINT_SUBDIVISION_MIPMAP_OFFSET(level-1) + i] = value;
        }
    }
}

float3 ComputeGlints2024_Area_Subdivided(float3 halfwayTS, float LdotH, float roughness, float3 ltcValue_subdivided[GLINT_SUBDIVISION_CELL_COUNT], float integratedNDF_subdivided[GLINT_SUBDIVISION_CELL_COUNT], float2 uv, float2 duvdx, float2 duvdy)
{
    // Approx. \int_hemisphere D(h) dh
    float totalNDF = ComputeTotalNDF(roughness);

    float R = _MicrofacetRoughness;

    float p_subdivided[GLINT_SUBDIVISION_CELL_COUNT];
    for (int i = 0; i < GLINT_SUBDIVISION_CELL_COUNT; ++i)
        p_subdivided[i] = saturate(integratedNDF_subdivided[i] / totalNDF); // = A = R*Dtarget/Dmax

    // Compute mipmap hierarchy on p
    float p_mipmap[GLINT_SUBDIVISION_MIPMAP_CELL_COUNT];
    BuildMipChain(p_subdivided, p_mipmap);

    // Skip computation below if probability is zero!
    if (p_mipmap[0] == 0)
        return 0;

    // Division by microfacet count handled internally.
    // Explicitly simulate smooth surface if _LogMicrofacetDensity < 0
    float D_subdivided[GLINT_SUBDIVISION_CELL_COUNT] = p_subdivided;
    if (_LogMicrofacetDensity > 0)
        SampleGlints2023NDF_Subdivided(halfwayTS, p_mipmap, uv, duvdx, duvdy, D_subdivided);
    float3 ltcValue = float3(0, 0, 0);
    for (int i = 0; i < GLINT_SUBDIVISION_CELL_COUNT; ++i)
    {
        if (p_subdivided[i] > 0)
            ltcValue += ltcValue_subdivided[i] * D_subdivided[i] / p_subdivided[i];
    }
    ltcValue = float3(SanitizePositiveFinite(ltcValue.x), SanitizePositiveFinite(ltcValue.y), SanitizePositiveFinite(ltcValue.z));
    return ltcValue;
}

#endif // GLINTS_SUBDIVISION_HLSL
