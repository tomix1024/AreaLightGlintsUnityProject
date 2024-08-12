using UnityEngine;
using System.Collections;
using UnityEngine.UI;
using TMPro;

[ExecuteInEditMode]
public class DisplayValue : MonoBehaviour
{
    // Workaround because of animated properties!
    public float value = 0;

    public TextMeshProUGUI uitext = null;

    void Update()
    {
        if (uitext == null)
            return;
        string text = value.ToString("0.00");
        uitext.text = text;
    }
}
