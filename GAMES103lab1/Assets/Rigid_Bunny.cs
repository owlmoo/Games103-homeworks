using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	
	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.99f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.5f;                 // for collision

	bool firstDt = true;
	Vector3 a = new Vector3(0, -9.8f, 0);   //加速度
	Vector3[] verticer;   // 每个顶点在初始状态下相较于局部原点的位置
	// Use this for initialization
	void Start () 
	{		
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		verticer = vertices;
		float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
		Time.fixedDeltaTime = dt;
	}
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}

	// In this function, update v and w by the impulse due to the collision with
	//a plane <P, N>
	void Collision_Impulse(Vector3 P, Vector3 N)
	{
		int colliderCount = 0;
		Vector3 averager = Vector3.zero, averagev = Vector3.zero;
		Matrix4x4 matrix4 = Matrix4x4.Rotate(transform.rotation);

		for (int i = 0;i < verticer.Length; i++)
        {
			Vector4 r = new Vector4(verticer[i].x, verticer[i].y, verticer[i].z, 1.0f);
			Vector3 relax = matrix4 * r;
			Vector3 xi = transform.position + relax;

			if(DetectX(xi, P, N))
            {
				Vector3 relav = Get_Cross_Matrix(w) * matrix4 * r;
				Vector3 vi = v + relav;

				if(DetectV(vi, N))
                {
					colliderCount += 1;
					averager += verticer[i];
					averagev += vi;
				}
            }
        }

		if(colliderCount > 0)
        {
			averager /= colliderCount;
			averagev /= colliderCount;

			// 计算希望得到的新的速度
			Vector3 vn = Vector3.Dot(averagev, N) * N;
			Vector3 vt = averagev - vn;

			float alpha = Mathf.Max(0.0f, 1 - 0.2f * (vn.magnitude / vt.magnitude));
			vn *= -0.9f;
			vt *= 0.8f;

			Vector3 vnew = vn + vt;

			// 根据希望得到的速度差，反向计算冲量j
			Matrix4x4 I = matrix4 * I_ref * matrix4.transpose;

			Vector4 ri = new Vector4(averager.x, averager.y, averager.z, 1.0f);
			Matrix4x4 Rri = Get_Cross_Matrix(matrix4 * ri);
			Matrix4x4 K = SubMatrix(DivideMatrix(Matrix4x4.identity, mass), Rri * I.inverse * Rri);

			Vector3 vbetween = vnew - averagev;
			Vector3 j = K.inverse * new Vector4(vbetween.x, vbetween.y, vbetween.z, 1.0f);

			// 更新v和w
			v += j / mass;
			
			Vector3 dw = I.inverse * Rri * new Vector4(j.x, j.y, j.z, 1.0f);
			w += dw;
        }
	}

	bool DetectX(Vector3 xi, Vector3 P, Vector3 N)
    {
		Vector3 X = xi - P;
		float direction = Vector3.Dot(X, N);
		return direction < 0;
    }

	bool DetectV(Vector3 vi, Vector3 N)
    {
		return Vector3.Dot(vi, N) < 0;
    }

	Matrix4x4 DivideMatrix(Matrix4x4 matrix, float mass)
    {
		matrix.m00 *= 1.0f / mass;
		matrix.m01 *= 1.0f / mass;
		matrix.m02 *= 1.0f / mass;
		matrix.m03 *= 1.0f / mass;
		matrix.m10 *= 1.0f / mass;
		matrix.m11 *= 1.0f / mass;
		matrix.m12 *= 1.0f / mass;
		matrix.m13 *= 1.0f / mass;
		matrix.m20 *= 1.0f / mass;
		matrix.m21 *= 1.0f / mass;
		matrix.m22 *= 1.0f / mass;
		matrix.m23 *= 1.0f / mass;
		matrix.m30 *= 1.0f / mass;
		matrix.m31 *= 1.0f / mass;
		matrix.m32 *= 1.0f / mass;
		matrix.m33 *= 1.0f / mass;
		return matrix;
    }

	Matrix4x4 SubMatrix(Matrix4x4 matrix1, Matrix4x4 matrix2)
    {
		Matrix4x4 matrix3 = new Matrix4x4();
		matrix3.m00 = matrix1.m00 - matrix2.m00;
		matrix3.m01 = matrix1.m01 - matrix2.m01;
		matrix3.m02 = matrix1.m02 - matrix2.m02;
		matrix3.m03 = matrix1.m03 - matrix2.m03;
		matrix3.m10 = matrix1.m10 - matrix2.m10;
		matrix3.m11 = matrix1.m11 - matrix2.m11;
		matrix3.m12 = matrix1.m12 - matrix2.m12;
		matrix3.m13 = matrix1.m13 - matrix2.m13;
		matrix3.m20 = matrix1.m20 - matrix2.m20;
		matrix3.m21 = matrix1.m21 - matrix2.m21;
		matrix3.m22 = matrix1.m22 - matrix2.m22;
		matrix3.m23 = matrix1.m23 - matrix2.m23;
		matrix3.m30 = matrix1.m30 - matrix2.m30;
		matrix3.m31 = matrix1.m31 - matrix2.m31;
		matrix3.m32 = matrix1.m32 - matrix2.m32;
		matrix3.m33 = matrix1.m33 - matrix2.m33;
		return matrix3;
    }

	void FixedUpdate() 
	{
		//Game Control
		if(Input.GetKey(KeyCode.R))
		{
			transform.position = new Vector3 (0, 0.6f, 0);
			restitution = 0.5f;
			launched=false;
		}
		if(Input.GetKey(KeyCode.S))
		{
			v = new Vector3 (5, 2, 0);
			launched = true;
			firstDt = true;
		}

        // Part I: Update velocities
        if (launched)
        {
			if (firstDt)
            {
                v -= a * (Time.fixedDeltaTime / 2);


				firstDt = false;
            }
            v += a * Time.fixedDeltaTime;
        }

        // Part II: Collision Impulse

        if (launched)
        {
			Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
			Collision_Impulse(new Vector3(1.99f, 0, 0), new Vector3(-1, 0, 0));
		}
		// Part III: Update position & orientation

		Vector3 x = transform.position;
		Quaternion q = transform.rotation;

		if (launched)
        {
			v *= linear_decay;
			w *= angular_decay;

			x += v * Time.fixedDeltaTime;

			Quaternion dq = new Quaternion(w.x * Time.fixedDeltaTime, w.y * Time.fixedDeltaTime, w.z * Time.fixedDeltaTime, 0.0f) * q;

			q = new Quaternion(q.x + dq.x, q.y + dq.y, q.z + dq.z, q.w + dq.w);
			q = q.normalized;
        }

		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}
