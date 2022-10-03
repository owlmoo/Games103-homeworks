using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Test : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        Vector4 vector4 = new Vector4(1, 2, 3, 4);
        Vector3 vector3 = vector4;
        print(vector3);
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
