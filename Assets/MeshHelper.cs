// C#
using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class MeshHelper {
	static List<Vector3> vertices;
	static List<Vector3> normals;
	// [... all other vertex data arrays you need]

	static List<int> indices;
	static Dictionary<uint, int> newVectices;

	static int GetNewVertex(int i1, int i2) {
		// We have to test both directions since the edge
		// could be reversed in another triangle
		uint t1 = ((uint) i1 << 16) | (uint) i2;
		uint t2 = ((uint) i2 << 16) | (uint) i1;
		if(newVectices.ContainsKey(t2))
			return newVectices[t2];
		if(newVectices.ContainsKey(t1))
			return newVectices[t1];
		// generate vertex:
		int newIndex = vertices.Count;
		newVectices.Add(t1, newIndex);

		// calculate new vertex
		vertices.Add((vertices[i1] + vertices[i2]) * 0.5f);
		normals.Add((normals[i1] + normals[i2]).normalized);
		// [... all other vertex data arrays]

		return newIndex;
	}


	public static void Subdivide(Mesh mesh) {
		newVectices = new Dictionary<uint, int>();

		vertices = new List<Vector3>(mesh.vertices);
		normals = new List<Vector3>(mesh.normals);
		// [... all other vertex data arrays]
		indices = new List<int>();

		int[] triangles = mesh.triangles;
		for(int i = 0; i < triangles.Length; i += 3) {
			int i1 = triangles[i + 0];
			int i2 = triangles[i + 1];
			int i3 = triangles[i + 2];

			int a = GetNewVertex(i1, i2);
			int b = GetNewVertex(i2, i3);
			int c = GetNewVertex(i3, i1);
			indices.Add(i1);
			indices.Add(a);
			indices.Add(c);
			indices.Add(i2);
			indices.Add(b);
			indices.Add(a);
			indices.Add(i3);
			indices.Add(c);
			indices.Add(b);
			indices.Add(a);
			indices.Add(b);
			indices.Add(c); // center triangle
		}
		mesh.vertices = vertices.ToArray();
		mesh.normals = normals.ToArray();
		// [... all other vertex data arrays]
		mesh.triangles = indices.ToArray();

		// since this is a static function and it uses static variables
		// we should erase the arrays to free them:
		newVectices = null;
		vertices = null;
		normals = null;
		// [... all other vertex data arrays]

		indices = null;
	}

	/*
	 * Creates a square mesh from (0, 0, 0) to (width, height, 0) with two triangles to which numSubdivisions subdivisions are applied.
	 */
	public static Mesh createSquareMesh(float width, float height, int numSubdivisions) {
		Mesh mesh = new Mesh();
		mesh.vertices = new Vector3[] {
			new Vector3(0, 0, 0),
			new Vector3(width, 0, 0),
			new Vector3(0, height, 0),
			new Vector3(width, height, 0)
		};
		mesh.triangles = new int[] {
			0, 2, 1,
			2, 3, 1
		};
		mesh.uv = new Vector2[] {
			new Vector2(0, 0),
			new Vector2(1, 0),
			new Vector2(0, 1),
			new Vector2(1, 1)
		};
		mesh.RecalculateNormals();

		for(int i = 0; i < numSubdivisions; i++) {
			MeshHelper.Subdivide(mesh);
		}
		return mesh;
	}

	/*
	 * Creates a triangle mesh with corners (0, 0, 0), (width, 0, 0) and (0, height, 0) to which numSubdivisions subdivisions are applied.
	 */
	public static Mesh createTriangleMesh(float width, float height, int numSubdivisions) {
		Mesh mesh = new Mesh();
		mesh.vertices = new Vector3[] {
			new Vector3(0, 0, 0),
			new Vector3(width, 0, 0),
			new Vector3(0, height, 0)
		};
		mesh.triangles = new int[] {
			0, 2, 1
		};
		mesh.uv = new Vector2[] {
			new Vector2(0, 0),
			new Vector2(1, 0),
			new Vector2(0, 1)
		};
		mesh.RecalculateNormals();
		
		for(int i = 0; i < numSubdivisions; i++) {
			MeshHelper.Subdivide(mesh);
		}
		return mesh;
	}

	/**
	 * Creates an array of vertex constraints constraining all outer edges of the given mesh.
	 * The returned array contains one element per vertex, where true means constrained and false unconstrained.
	 */
	public static bool[] createOuterEdgeVertexContraints(Vec3D[] vertices, List<Edge> edges) {
		bool[] verticesMovementConstraints = new bool[vertices.Length];
		for(int i = 0; i < verticesMovementConstraints.Length; i++) {
			verticesMovementConstraints[i] = false;
		}
		foreach(Edge edge in edges) {
			if(!edge.hasSideFlaps()) {
				verticesMovementConstraints[edge.ve1] = true;
				verticesMovementConstraints[edge.ve2] = true;
			}
		}
		return verticesMovementConstraints;
	}
}
