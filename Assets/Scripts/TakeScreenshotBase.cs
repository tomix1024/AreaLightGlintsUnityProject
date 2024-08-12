using UnityEngine;
using UnityEditor;

[ExecuteInEditMode]
class TakeScreenshotBase : MonoBehaviour
{
    public bool hdr = false;
    public int width = 1024;
    public int height = 1024;

    public bool runOnce = false;

    public void Update()
    {
        if (runOnce)
        {
            runOnce = false;
            Run();
        }
    }

    private RenderTexture renderTexture;
    private Texture2D texture;

    protected void CreateTextures()
    {
        RenderTextureFormat renderTextureFormat = hdr ? RenderTextureFormat.DefaultHDR : RenderTextureFormat.Default;
        TextureFormat textureFormat = hdr ? TextureFormat.RGBAHalf : TextureFormat.RGB24;
        renderTexture = new RenderTexture(width, height , 16, renderTextureFormat);
        texture = new Texture2D(width, height, textureFormat, mipChain: false, linear: hdr);
    }

    protected void RenderCameraToTexture(Camera camera)
    {
        camera.targetTexture = renderTexture;
        camera.Render();
        camera.targetTexture = null;

        RenderTexture.active = renderTexture;
        texture.ReadPixels(new Rect(0, 0, width, height), 0, 0);
        RenderTexture.active = null;
    }

    protected void SaveToPathNoExt(string path)
    {
        byte[] bytes = hdr ? texture.EncodeToEXR() : texture.EncodeToPNG();
        string ext = hdr ? "exr" : "png";
        System.IO.File.WriteAllBytes($"{path}.{ext}", bytes);
    }

    public virtual void Run() {}
}
