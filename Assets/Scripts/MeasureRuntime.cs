using UnityEngine;
using System.Collections;
using UnityEngine.UI;
using System.IO;

public class MeasureRuntime : MonoBehaviour
{
    // In seconds:
    public KeyCode startKey = KeyCode.E;
    public float recordDuration = 1;

    public string filename = "runtime.csv";

    public int bufferSize = 1024;

    private float startTime = float.NegativeInfinity;

    private int bufferIndex = 0;
    private float[] buffer = null;

    private bool hasData = false;

    private bool triggerStart = false;
    private bool measurementInProgress = false;

    public void StartMeas()
    {
        triggerStart = true;
        measurementInProgress = true;
    }

    public bool MeasInProgress()
    {
        return measurementInProgress;
    }

    void Awake()
    {
        Reset();
    }

    void Update()
    {
        float timePoint = Time.realtimeSinceStartup;
        if (triggerStart || Input.GetKeyDown(startKey))
        {
            triggerStart = false;
            startTime = timePoint;
        }

        if (timePoint >= startTime && timePoint < startTime + recordDuration)
        {
            hasData = true;
            if (bufferIndex < buffer.Length)
                buffer[bufferIndex++] = Time.unscaledDeltaTime;
        }
        else if (hasData)
        {
            SaveData();
            Reset();
        }
    }


    void SaveData()
    {
        if (bufferIndex > 0)
        {
            using(var file = new StreamWriter(filename))
            {
                file.WriteLine("time deltaTime");

                float time = 0;
                for (int i = 0; i < bufferIndex; ++i)
                {
                    float deltaTime = buffer[i];
                    file.WriteLine(time + " " + deltaTime);
                    time += deltaTime;
                }
            }
        }
    }

    void Reset()
    {
        measurementInProgress = false;
        hasData = false;
        if (buffer == null || buffer.Length != bufferSize)
            buffer = new float[bufferSize];
        bufferIndex = 0;
    }
}
