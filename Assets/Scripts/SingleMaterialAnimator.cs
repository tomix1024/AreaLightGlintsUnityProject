using UnityEngine;
using System.Collections;

[ExecuteInEditMode]
public class SingleMaterialAnimator : MonoBehaviour
{
    // Workaround because of animated properties!
    public string key = "_Smoothness";
    public float value = 0;

    public MeshRenderer renderer;
    public int index;

    void Update()
    {
        if (renderer == null)
            return;
        var material = renderer.sharedMaterials[index];
        if (material == null)
            return;

        Debug.Log($"set {material.name} {key}={value}");
        material.SetFloat(key, value);
    }
}
