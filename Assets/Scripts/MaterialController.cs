using UnityEngine;
using UnityEngine.UI;

public class MaterialController : MonoBehaviour
{
    public Material[] materials;
    public string[] materialNames;
    public int currentIndex = 0;

    public Renderer[] affectedRenderers;

    public KeyCode key = KeyCode.Space;
    public KeyCode switchNormalKey = KeyCode.N;

    public bool enableNormalMap = false;

    public void Awake()
    {
        Debug.Log("Setup!");

        AssignMaterial();
    }

    public void Update()
    {
        if (Input.GetKeyDown(key))
        {
            currentIndex = (currentIndex + 1) % materials.Length;
            AssignMaterial();
        }
        if (Input.GetKeyDown(switchNormalKey))
        {
            enableNormalMap = !enableNormalMap;
            AssignMaterial();
        }

        KeyCode[] numberCodes = { KeyCode.Alpha1, KeyCode.Alpha2, KeyCode.Alpha3, KeyCode.Alpha4, KeyCode.Alpha5, KeyCode.Alpha6, KeyCode.Alpha7, KeyCode.Alpha8, KeyCode.Alpha9, KeyCode.Alpha0, KeyCode.Minus, KeyCode.Equals };

        for (int i = 0; i < numberCodes.Length; ++i)
        {
            if (i >= materials.Length)
                break;

            if (Input.GetKeyDown(numberCodes[i]))
            {
                currentIndex = i;
                AssignMaterial();
                break;
            }
        }
    }

    public string CurrentMaterialName()
    {
        if (materialNames != null && materialNames.Length > currentIndex && materialNames[currentIndex] != "")
        {
            return materialNames[currentIndex];
        }
        else
        {
            return materials[currentIndex].name;
        }
    }

    public void AssignMaterial()
    {
        var material = materials[currentIndex];
        if (enableNormalMap)
        {
            material.SetFloat("_NormalScale", 1.0f);
            // material.SetFloat("_UseNormalMap", 1.0f);
            // material.EnableKeyword("USE_NORMALMAP");
        }
        else
        {
            material.SetFloat("_NormalScale", 0.0f);
            // material.SetFloat("_UseNormalMap", 0.0f);
            // material.DisableKeyword("USE_NORMALMAP");
        }
        foreach (var renderer in affectedRenderers)
        {
            renderer.material = material;
        }
    }
}
