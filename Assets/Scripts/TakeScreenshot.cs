using UnityEngine;
using UnityEditor;

[ExecuteInEditMode]
class TakeScreenshot : TakeScreenshotBase
{
    public Camera camera;

    public string outputPath;

    public override void Run()
    {
        if (camera == null)
            camera = GetComponentInChildren<Camera>();

        CreateTextures();

        RenderCameraToTexture(camera);

        SaveToPathNoExt(outputPath);
    }

}
