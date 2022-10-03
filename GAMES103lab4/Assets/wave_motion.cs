using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.004f;
	float damping 	= 0.996f;
	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag=true;

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;

	float liftlevel = 1.0f;      // 随机升起的高度水平 (0.2 ~ 2) * liftlevel
	Vector3[] globalX;
	GameObject cubeObj;
	GameObject blockObj;

	int[,] neighbor = new int[4, 2]
	{
		{-1, 0 },
		{0, -1 },
		{1, 0 },
		{0, 1 }
	};

	// Use this for initialization
	void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		Vector3[] X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}
		globalX = new Vector3[size * size];
		for(int i = 0; i < size; i++)
        {
			globalX[i] = X[i];
        }

		cubeObj = GameObject.Find("Cube");
		blockObj = GameObject.Find("Block");
	}

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}

	void Shallow_Wave(float[,] old_h, float[,] h, float [,] new_h)
	{
		//Step 1:
		//TODO: Compute new_h based on the shallow wave model.
		for(int i = 0; i < size; i++)
        {
			for(int j = 0; j < size; j++)
            {
                new_h[i, j] = h[i, j] + (h[i,j] - old_h[i,j]) * damping;

				for(int k = 0; k < 4; k++)
                {
					int tmpi = i + neighbor[k, 0];
					int tmpj = j + neighbor[k, 1];

					if(tmpi >= 0 && tmpi < size && tmpj >= 0 && tmpj < size)
                    {
						new_h[i, j] += (h[tmpi, tmpj] - h[i, j]) * rate;
                    }
				}
            }
        }

        //Step 2: Block->Water coupling
        //TODO: for block 1, calculate low_h.
        //TODO: then set up b and cg_mask for conjugate gradient.
        //TODO: Solve the Poisson equation to obtain vh (virtual height).


        // 先找到contact有哪些点，做y=0的一个截面，看正方体的截面是什么
        // 为简便起见，截面恒为正方形，下降为0.5f
        cube_v = cubeObj.transform.position;
        Vector3[] X = GetComponent<MeshFilter>().mesh.vertices;

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                if ((Mathf.Abs(X[i * size + j].x - cube_v.x) <= 0.5f) && Mathf.Abs(X[i * size + j].z - cube_v.z) <= 0.5f)
                {
                    cg_mask[i, j] = true;
                    low_h[i, j] = -0.5f;
                    b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
                }
                else
                {
					cg_mask[i, j] = false;
					vh[i, j] = 0;
                }
            }
        }

        Conjugate_Gradient(cg_mask, b, vh, 0, size - 1, 0, size - 1);



		//TODO: for block 2, calculate low_h.
		//TODO: then set up b and cg_mask for conjugate gradient.
		//TODO: Solve the Poisson equation to obtain vh (virtual height).
		cube_w = blockObj.transform.position;
        
		for(int i = 0; i < size; i++)
        {
			for(int j = 0; j < size; j++)
            {
                if ((Mathf.Abs(X[i * size + j].x - cube_w.x) <= 0.5f) && Mathf.Abs(X[i * size + j].z - cube_w.z) <= 0.5f)
                {
					cg_mask[i, j] = true;
					low_h[i, j] = -0.5f;
					b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
                }
            }
        }

		Conjugate_Gradient(cg_mask, b, vh, 0, size - 1, 0, size - 1);
        //TODO: Diminish vh.

		for(int i = 0; i < size; i++)
        {
			for(int j = 0; j < size; j++)
            {
				vh[i, j] *= gamma;
            }
        }

        //TODO: Update new_h by vh.

        for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					int tmpi = i + neighbor[k, 0];
					int tmpj = j + neighbor[k, 1];

					if (tmpi >= 0 && tmpi < size && tmpj >= 0 && tmpj < size)
					{
						new_h[i, j] += (vh[tmpi, tmpj] - vh[i, j]) * rate; 
					}
				}
			}
		}

		//Step 3
		//TODO: old_h <- h; h <- new_h;
		for (int i = 0; i < size; i++)
        {
			for (int j = 0; j < size; j++)
			{
				old_h[i, j] = h[i, j];
				h[i, j] = new_h[i, j];
			}
		}
		
		//Step 4: Water->Block coupling.
		//More TODO here.
	}
	

	// Update is called once per frame
	void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

        //TODO: Load X.y into h.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                h[i, j] = X[i * size + j].y;
            }
        }

        if (Input.GetKeyDown ("r")) 
		{
			//TODO: Add random water.
			int i = (int)Random.Range(0, size);
			int j = (int)Random.Range(0, size);

			float r = Random.Range(0.2f, 2.0f) * liftlevel;

			h[i, j] += r;

			//将高度提升平均分配给邻居
			List<Vector2> neighbors = new List<Vector2>();

			for(int k = 0; k < 4; k++)
            {
				int newi = i + neighbor[k, 0];
				int newj = j + neighbor[k, 1];

				if(newi >= 0 && newi < size && newj >= 0 && newj < size)
                {
					neighbors.Add(new Vector2(newi, newj));
                }
            }
			float drop = r / neighbors.Count;
			for(int k = 0; k < neighbors.Count; k++)
            {
				h[(int)neighbors[k].x, (int)neighbors[k].y] -= drop;
            }
		}
	
		for(int l=0; l<8; l++)
		{
			Shallow_Wave(old_h, h, new_h);
		}

        //TODO: Store h back into X.y and recalculate normal.

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                X[i * size + j].y = h[i, j];
            }
        }
		mesh.vertices = X;
		mesh.RecalculateNormals();
    }
}
