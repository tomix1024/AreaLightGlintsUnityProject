using UnityEngine;
using UnityEditor;
using System.Collections.Generic;

[ExecuteInEditMode]
class TakeScreenshotMaterials : TakeScreenshotBase
{
    public Camera[] cameras;
    public string[] cameraNames;

    public string outputPath;

    public Material[] materials;
    public string[] materialNames;
    public MeshRenderer renderer;

    public float[] smoothnesses;

    public float[] microfacetDensities;

    public GameObject[] enableObjects;
    public GameObject[] disableObjects;

    public bool updateTextureOffset;

    private Dictionary<string, string> pathReplacement = new Dictionary<string, string>();

    private IEnumerable<Material> ForMaterials()
    {
        if (materials == null || materials.Length == 0)
        {
            Material material = renderer.material;
            yield return material;
        }
        else
        {
            for (int i = 0; i < materials.Length; ++i)
            {
                Material material = materials[i];
                if (material == null)
                    continue;
                if (materialNames != null && i < materialNames.Length)
                    pathReplacement["{mat}"] = materialNames[i];
                else
                    pathReplacement.Remove("{mat}");
                renderer.material = material;
                yield return material;
            }
        }
        pathReplacement.Remove("{mat}");
    }

    private IEnumerable<float> ForSmoothness(Material mat)
    {
        if (smoothnesses == null || smoothnesses.Length == 0)
        {
            float smoothness = mat.GetFloat("_Smoothness");
            pathReplacement["{smoothness}"] = $"{smoothness}";
            yield return smoothness;
        }
        else
        {
            foreach (var smoothness in smoothnesses)
            {
                pathReplacement["{smoothness}"] = $"{smoothness}";
                mat.SetFloat("_Smoothness", smoothness);
                yield return smoothness;
            }
        }
        pathReplacement.Remove("{smoothness}");
    }

    private IEnumerable<float> ForMicrofacetDensities(Material mat)
    {
        if (microfacetDensities == null || microfacetDensities.Length == 0)
        {
            float density = mat.GetFloat("_LogMicrofacetDensity");
            pathReplacement["{density}"] = $"{density}";
            yield return density;
        }
        else
        {
            foreach (var density in microfacetDensities)
            {
                pathReplacement["{density}"] = $"{density}";
                mat.SetFloat("_LogMicrofacetDensity", density);
                yield return density;
            }
        }
    }

    private IEnumerable<Camera> ForCameras()
    {
        for (int i = 0; i < cameras.Length; ++i)
        {
            Camera cam = cameras[i];
            string camName = cameraNames != null && i < cameraNames.Length ? cameraNames[i] :  $"cam{i:02d}";
            pathReplacement["{cam}"] = camName;
            yield return cam;
        }
        pathReplacement.Remove("{cam}");
    }

    public override void Run()
    {
        pathReplacement.Clear();
        CreateTextures();

        foreach (var obj in enableObjects)
        {
            obj.SetActive(true);
        }
        foreach (var obj in disableObjects)
        {
            obj.SetActive(false);
        }

        foreach (var material in ForMaterials())
        {
            int textureOffset = 0;
            foreach (var smoothness in ForSmoothness(material))
            {
                foreach (var density in ForMicrofacetDensities(material))
                {
                    textureOffset += 10;
                    if (updateTextureOffset)
                        material.SetTextureOffset("_BaseColorMap", new Vector2(textureOffset, 0));
                    foreach (var cam in ForCameras())
                    {
                        // Render current camera
                        RenderCameraToTexture(cam);

                        // Customize path
                        string tempPath = outputPath;
                        foreach (var entry in pathReplacement)
                        {
                            Debug.Log($"replace {entry.Key} with {entry.Value}");
                            tempPath = tempPath.Replace(entry.Key, entry.Value);
                        }

                        Debug.Log($"output path {tempPath}");
                        // Save to custom path
                        SaveToPathNoExt(tempPath);
                    }
                }
            }
        }
    }

}
