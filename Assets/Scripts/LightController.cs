using UnityEngine;
using UnityEngine.UI;

public class LightController : MonoBehaviour
{
    public GameObject[] lights; // Actually any type supported...
    public string[] lightNames;
    public int currentIndex = 0;

    public KeyCode key = KeyCode.L;

    public void Awake()
    {
        Debug.Log("Setup LightController!");

        AssignMaterial();
    }

    public void Update()
    {
        if (Input.GetKeyDown(key))
        {
            currentIndex = (currentIndex + 1) % lights.Length;
            AssignMaterial();
        }
    }

    public string CurrentLightName()
    {
        if (lightNames != null && lightNames.Length > currentIndex && lightNames[currentIndex] != "")
        {
            return lightNames[currentIndex];
        }
        else
        {
            return lights[currentIndex].name;
        }
    }

    public void AssignMaterial()
    {
        for (int i = 0; i < lights.Length; ++i)
        {
            lights[i].SetActive(i == currentIndex);
        }
    }
}
