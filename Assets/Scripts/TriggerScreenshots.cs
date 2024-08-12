using UnityEngine;
using UnityEditor;

[ExecuteInEditMode]
class TriggerScreenshots : MonoBehaviour
{
    public TakeScreenshotBase[] screenshots;

    public bool runOnce = false;

    public void Update()
    {
        if (runOnce)
        {
            runOnce = false;
            Run();
        }
    }

    public void Run()
    {
        foreach (var screenshot in screenshots)
        {
            if (screenshot != null)
                screenshot.runOnce = true;
        }
    }
}
