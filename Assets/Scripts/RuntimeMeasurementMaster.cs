using UnityEngine;
using System.Collections;

public class RuntimeMeasurementMaster : MonoBehaviour
{
    public MaterialController materialController;
    public LightController lightController;
    public MeasureRuntime measureRuntime;

    public void Start()
    {
        if (materialController == null)
            materialController = GetComponent<MaterialController>();
        if (lightController == null)
            lightController = GetComponent<LightController>();
        if (measureRuntime == null)
            measureRuntime = GetComponent<MeasureRuntime>();

        StartCoroutine(MeasurementLoop());
    }

    IEnumerator MeasurementLoop()
    {
        yield return null; // wait for next frame

        // warmup GPU!
        yield return new WaitForSeconds(60);

        for (int lightIndex = 0; lightIndex < lightController.lights.Length; ++lightIndex)
        {
            lightController.currentIndex = lightIndex;
            string lightName = lightController.lightNames[lightIndex];
            lightController.AssignMaterial();

            for (int materialIndex = 0; materialIndex < materialController.materials.Length; ++materialIndex)
            {
                materialController.currentIndex = materialIndex;
                string matName = materialController.materialNames[materialIndex];

                bool[] enableNormalMapOptions = { false, true };
                foreach (bool enableNormalMap in enableNormalMapOptions)
                {
                    string flatrandom = enableNormalMap ? "rand" : "flat";
                    materialController.enableNormalMap = enableNormalMap;
                    materialController.AssignMaterial();

                    // Wait before starting measurement
                    yield return new WaitForSeconds(10);

                    measureRuntime.filename = $"runtime_{lightName}_{matName}_{flatrandom}.csv";
                    measureRuntime.StartMeas();
                    while (measureRuntime.MeasInProgress())
                    {
                        yield return null;
                    }
                }
            }
        }
        Application.Quit();
    }
}
