using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class FVM : MonoBehaviour
{
	float dt 			= 0.003f;
    float mass 			= 1;
	float stiffness_0	= 20000.0f;
    float stiffness_1 	= 5000.0f;
    float damp			= 0.999f;

	int[] 		Tet;
	int tet_number;			//The number of tetrahedra

	Vector3[] 	Force;
	Vector3[] 	V;
	Vector3[] 	X;
	int number;				//The number of vertices

	Matrix4x4[] inv_Dm;

	SVD svd = new SVD();

    float g = -9.8f;        //The gravity acceleration
    float lambda = 0.2f;    //The laplacian blender parameters

    List<HashSet<int>> neighbors;  //The neighbors of every vertex

    void Start()
    {
    	// FILO IO: Read the house model from files.
    	// The model is from Jonathan Schewchuk's Stellar lib.
    	{
    		string fileContent = File.ReadAllText("Assets/house2.ele");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		
    		tet_number=int.Parse(Strings[0]);
        	Tet = new int[tet_number*4];

    		for(int tet=0; tet<tet_number; tet++)
    		{
				Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;
				Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
				Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
				Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
			}
    	}
    	{
			string fileContent = File.ReadAllText("Assets/house2.node");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		number = int.Parse(Strings[0]);
    		X = new Vector3[number];
       		for(int i=0; i<number; i++)
       		{
       			X[i].x=float.Parse(Strings[i*5+5])*0.4f;
       			X[i].y=float.Parse(Strings[i*5+6])*0.4f;
       			X[i].z=float.Parse(Strings[i*5+7])*0.4f;
       		}
    		//Centralize the model.
	    	Vector3 center=Vector3.zero;
	    	for(int i=0; i<number; i++)		center+=X[i];
	    	center=center/number;
	    	for(int i=0; i<number; i++)
	    	{
	    		X[i]-=center;
	    		float temp=X[i].y;
	    		X[i].y=X[i].z;
	    		X[i].z=temp;
	    	}
		}

        //Create triangle mesh.
       	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];

        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }

        int[] triangles = new int[tet_number*12];
        for(int t=0; t<tet_number*4; t++)
        {
        	triangles[t*3+0]=t*3+0;
        	triangles[t*3+1]=t*3+1;
        	triangles[t*3+2]=t*3+2;
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.triangles = triangles;
		mesh.RecalculateNormals ();


		V 	  = new Vector3[number];
        Force = new Vector3[number];

        inv_Dm = new Matrix4x4[tet_number];

        for(int tet = 0; tet < tet_number; tet++)
        {
            inv_Dm[tet] = Build_Edge_Matrix(tet).inverse;
        }

        // Add neighbors
        neighbors = new List<HashSet<int>>();
        for (int i = 0; i < number; i++) neighbors.Add(new HashSet<int>());

        for(int tet = 0; tet < tet_number; tet++)
        {
            for(int i = 0; i < 4; i++)
            {
                for(int j = 0; j < 4; j++)
                {
                    neighbors[Tet[tet * 4 + i]].Add(Tet[tet * 4 + j]);
                }
            }
        }
    }

    Matrix4x4 Build_Edge_Matrix(int tet)
    {
        Matrix4x4 ret = Matrix4x4.zero;
        Vector3[] XX = new Vector3[3];

        XX[0] = X[Tet[tet * 4 + 0]] - X[Tet[tet * 4 + 1]];
        XX[1] = X[Tet[tet * 4 + 0]] - X[Tet[tet * 4 + 2]];
        XX[2] = X[Tet[tet * 4 + 0]] - X[Tet[tet * 4 + 3]];

        for(int i = 0; i < 3; i++)  // 按列
        {
            for (int j = 0; j < 3; j++)  // 按行
            {
                ret[j, i] = XX[i][j];
            }
        }
        ret[3, 3] = 1.0f;
		return ret;
    }

    Matrix4x4 Minus_Matrix(Matrix4x4 matrix1, Matrix4x4 matrix2)
    {
        Matrix4x4 matrix3 = new Matrix4x4();
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                matrix3[i, j] = matrix1[i, j] - matrix2[i, j];
            }
        }
        return matrix3;
    }

    Matrix4x4 Add_Matrix(Matrix4x4 matrix1, Matrix4x4 matrix2)
    {
        Matrix4x4 matrix3 = new Matrix4x4();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                matrix3[i, j] = matrix1[i, j] + matrix2[i, j];
            }
        }
        return matrix3;
    }

    Matrix4x4 Multiply_Float(Matrix4x4 matrix1, float a)
    {
        Matrix4x4 matrix2 = new Matrix4x4();
        for (int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                matrix2[i, j] = a * matrix1[i, j];
            }
        }
        return matrix2;
    }

    float Matrix_Trace(Matrix4x4 matrix)
    {
        float result = 0.0f;
        for(int i = 0; i < 3; i++)  // 忽略最后一位的迹
        {
            result += matrix[i, i];
        }
        return result;
    }


    void _Update()
    {
    	// Jump up.
		if(Input.GetKeyDown(KeyCode.Space))
    	{
    		for(int i=0; i<number; i++)
    			V[i].y+=0.2f;
    	}

    	for(int i=0 ;i<number; i++)
    	{
            Force[i] = new Vector3(0, mass * g, 0);
    	}

    	for(int tet=0; tet<tet_number; tet++)
    	{
            Matrix4x4 tempM = inv_Dm[tet].inverse;
            Matrix4x4 XX = Build_Edge_Matrix(tet);
            Matrix4x4 F = XX * inv_Dm[tet];

            Matrix4x4 G = Multiply_Float(Minus_Matrix(F.transpose * F, Matrix4x4.identity), 0.5f);

            Matrix4x4 S = Add_Matrix(Multiply_Float(G, 2 * stiffness_1), Multiply_Float(Matrix4x4.identity, stiffness_0 * Matrix_Trace(G)));

            Matrix4x4 P = F * S;
            Matrix4x4 result = Multiply_Float(P * inv_Dm[tet].transpose, -1.0f / (6 * Matrix4x4.Determinant(inv_Dm[tet])));

            // 接下来就是从result中提取各个force
            Vector3 force1 = new Vector3(result[0, 0], result[1, 0], result[2, 0]);
            Vector3 force2 = new Vector3(result[0, 1], result[1, 1], result[2, 1]);
            Vector3 force3 = new Vector3(result[0, 2], result[1, 2], result[2, 2]);

            Vector3 force0 = -(force1 + force2 + force3);

            Force[Tet[tet * 4 + 0]] += force0;
            Force[Tet[tet * 4 + 1]] += force1;
            Force[Tet[tet * 4 + 2]] += force2;
            Force[Tet[tet * 4 + 3]] += force3;

        }

        Vector3[] tmpV = new Vector3[number];

    	for(int i=0; i<number; i++)
    	{
            Vector3 dV = (Force[i] / mass) * dt;
            
            V[i] += dV;

            X[i] += V[i] * dt;

            if (X[i].y <= -3.0f)
            {
                V[i] = -V[i] * damp;
                X[i] = new Vector3(X[i].x, -3.0f, X[i].z);
            }

            tmpV[i] = Vector3.zero;
            foreach(int neighbor in neighbors[i])
            {
                tmpV[i] += V[neighbor];
            }
            tmpV[i] /= neighbors[i].Count;
        }

        for (int i = 0; i < number; i++) V[i] = lambda * tmpV[i] + (1 - lambda) * V[i];
    }

    void Update()
    {
    	for(int l=0; l<10; l++)
    		 _Update();

    	// Dump the vertex array for rendering.
    	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.RecalculateNormals ();
    }
}


// 遇到的问题：
// 1、公式写错了
// 2、laplacian 平滑的时候直接用的替代，没有blender