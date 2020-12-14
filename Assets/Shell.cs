using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Shell : MonoBehaviour {

    public float kLength;
    public float kArea;
    public float kBend;

    private GameObject shellObj;
    private List<int>[] vertexTriangles;
    private Vector3[] originalVertices; // Vertices in undeformed state.
    private Vector3[] verticesVelocity;
    private Vector3[] verticesAcceleration;

    private bool doUpdate = true;

    // Start is called before the first frame update.
    void Start() {

        // Initialize fields. Doing this overwrites the values set in Unity's inspector.
        this.kLength = 0f;
        this.kArea = 0f;
        this.kBend = 10f;
        
        // Create a new object in the scene.
        this.shellObj = new GameObject();
        MeshRenderer meshRenderer = this.shellObj.AddComponent<MeshRenderer>();
        meshRenderer.sharedMaterial = new Material(Shader.Find("Custom/StandardTwoSides"));
        MeshFilter meshFilter = this.shellObj.AddComponent<MeshFilter>();
        Mesh mesh = this.createMesh();
        //Mesh mesh = this.createSingleTriangleMesh();

        // Store undeformed vertices.
        this.originalVertices = (Vector3[]) mesh.vertices.Clone();

        // Perform mesh deformation and translation (this should be replaced with some meaningful force/deformation).
        Vector3[] verts = mesh.vertices;

        //verts[0].x += 2.5f;
        //verts[0].y += 2.5f;
        //verts[0].z -= 5f;

        verts[1].z += 1f;
        verts[2].z += 1f;

        //verts[1].z += 1f;

        mesh.vertices = verts;
        mesh.RecalculateNormals();
        meshFilter.mesh = mesh;
        this.shellObj.transform.position = new Vector3(0, 1, 0);

        // Create triangle list per vertex.
        this.vertexTriangles = new List<int>[mesh.vertexCount];
        for(int i = 0; i < this.vertexTriangles.Length; i++) {
            this.vertexTriangles[i] = new List<int>();
        }
        for(int i = 0; i < mesh.triangles.Length; i += 3) {
            int triangleId = i / 3;
            this.vertexTriangles[mesh.triangles[i]].Add(triangleId);
            this.vertexTriangles[mesh.triangles[i + 1]].Add(triangleId);
            this.vertexTriangles[mesh.triangles[i + 2]].Add(triangleId);
        }

        // Sort triangle list per vertex in clockwise order.
        for(int vertexInd = 0; vertexInd < this.vertexTriangles.Length; vertexInd++) {
            List<int> triangleList = this.vertexTriangles[vertexInd];
            if(triangleList.Count > 0) {
                List<int> sortedTriangleList = new List<int>(triangleList.Count);
                sortedTriangleList.Add(triangleList[0]);
                for(int i = 1; i < triangleList.Count; i++) {

                    // Get triangle vertices.
                    int triangleBaseInd1 = triangleList[i] * 3;
                    int v11 = mesh.triangles[triangleBaseInd1];
                    int v12 = mesh.triangles[triangleBaseInd1 + 1];
                    int v13 = mesh.triangles[triangleBaseInd1 + 2];

                    // Get the vertex indices of the other vertices that are connected to this vertex's edges.
                    int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
                    int otherVertexAntiClockwiseInd1 = (vertexInd == v11 ? v12 : (vertexInd == v12 ? v13 : v11));

                    // Insert the triangle into the sorted triangles list.
                    bool inserted = false;
                    for(int j = 0; j < sortedTriangleList.Count; j++) {

                        // Get triangle vertices.
                        int triangleBaseInd2 = sortedTriangleList[j] * 3;
                        int v21 = mesh.triangles[triangleBaseInd2];
                        int v22 = mesh.triangles[triangleBaseInd2 + 1];
                        int v23 = mesh.triangles[triangleBaseInd2 + 2];

                        // Get the vertex indices of the other vertices that are connected to this vertex's edges.
                        int otherVertexClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));
                        int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));

                        // Insert the triangle before the current triangle if it is placed next to it in anti-clockwise rotation.
                        if(otherVertexClockwiseInd1 == otherVertexAntiClockwiseInd2) {
                            sortedTriangleList.Insert(j, triangleList[i]);
                            inserted = true;
                            break;
                        }

                        // Insert the triangle after the current triangle if it is placed next to it in clockwise rotation.
                        if(otherVertexAntiClockwiseInd1 == otherVertexClockwiseInd2) {
                            sortedTriangleList.Insert(j + 1, triangleList[i]);
                            inserted = true;
                            break;
                        }
                    }
                    if(!inserted) {
                        sortedTriangleList.Add(triangleList[i]);
                    }
                }
            }
        }

        // Initialize vertex velocity and acceleration.
        int numVertices = mesh.vertices.Length;
        this.verticesVelocity = new Vector3[numVertices];
        this.verticesAcceleration = new Vector3[numVertices];
        for(int i = 0; i < numVertices; i++) {
            this.verticesVelocity[i] = new Vector3(0, 0, 0);
            this.verticesAcceleration[i] = new Vector3(0, 0, 0);
        }
    }

    private Mesh createMesh() {
        Mesh mesh = new Mesh();
        float width = 5;
        float height = 5;
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

        //for(int i = 0; i < 2; i++) {
        //    MeshHelper.Subdivide(mesh);
        //}
        return mesh;
    }

    private Mesh createSingleTriangleMesh() {
        Mesh mesh = new Mesh();
        float width = 5;
        float height = 5;
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

        return mesh;
    }

    // Update is called once per frame.
    void Update() {

        // Don't update if no noticable time has passed, or when the game is paused.
        if(Time.deltaTime == 0 || Time.timeScale == 0
                || (!this.doUpdate && !Input.GetKey(KeyCode.F) && !Input.GetKeyDown(KeyCode.G))) {
            return;
        }

        // Get the mesh.
        Mesh mesh = this.getMesh();

        // Compute vertex energy gradient array.
        Vector3[] vertexEnergyGradient = new Vector3[mesh.vertexCount];
        for(int i = 0; i < vertexEnergyGradient.Length; i++) {
            vertexEnergyGradient[i] = new Vector3(0, 0, 0);
        }
        for(int vertexInd = 0; vertexInd < this.vertexTriangles.Length; vertexInd++) {
            
            for(int i = 0; i < this.vertexTriangles[vertexInd].Count; i++) {

                // Get two possibly adjacent triangles in clockwise direction.
                int triangleId = this.vertexTriangles[vertexInd][i];
                int nextTriangleId = this.vertexTriangles[vertexInd][(i + 1) % this.vertexTriangles[vertexInd].Count];

                // Get triangle vertices.
                int triangleBaseInd1 = triangleId * 3;
                int triangleBaseInd2 = nextTriangleId * 3;
                int v11 = mesh.triangles[triangleBaseInd1];
                int v12 = mesh.triangles[triangleBaseInd1 + 1];
                int v13 = mesh.triangles[triangleBaseInd1 + 2];
                int v21 = mesh.triangles[triangleBaseInd2];
                int v22 = mesh.triangles[triangleBaseInd2 + 1];
                int v23 = mesh.triangles[triangleBaseInd2 + 2];

                // Get the vertex indices of the other vertices that are connected to the possibly shared triangle edge.
                int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
                int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));

                // Handle the edge, or both edges if they are not the same.
                bool edgeSharedByTriangles = (otherVertexClockwiseInd1 == otherVertexAntiClockwiseInd2);

                // Calculate spring energy gradient in the edge.
                vertexEnergyGradient[vertexInd] += kLength * this.getEdgeLengthEnergyGradient(vertexInd, otherVertexClockwiseInd1);

                if(edgeSharedByTriangles) {

                    // Calculate bending energy gradient.
                    vertexEnergyGradient[vertexInd] += kBend
                            * this.getBendingEnergyGradient(v11, v12, v13, v21, v22, v23, vertexInd, otherVertexClockwiseInd1);
                } else {

                    // Calculate spring energy gradient in the second edge.
                    vertexEnergyGradient[vertexInd] += kLength * this.getEdgeLengthEnergyGradient(vertexInd, otherVertexAntiClockwiseInd2);
                }

                // Calculate the area energy gradient in the triangle.
                if(vertexInd == v11) {
                    vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(v11, v12, v13);
                } else if(vertexInd == v12) {
                    vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(v12, v13, v11);
                } else {
                    vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(v13, v11, v12);
                }
            }
        }

        // Perform Newmark Time Stepping (ODE integration).
        float mass = 0.1f; // Some random constant mass. TODO - Replace this with some area-dependent mass.
        float gamma = 1f; // 0.5f;
        float beta = 0f; // 0.25f;
        Vector3[] vertices = mesh.vertices;
        for(int i = 0; i < vertices.Length; i++) {

            // Calculate acceleration.
            Vector3 newAcceleration = mass * -vertexEnergyGradient[i];

            // Update position.
            vertices[i] += Time.deltaTime * this.verticesVelocity[i]
                    + Time.deltaTime * Time.deltaTime * ((1f / 2f - beta) * this.verticesAcceleration[i] + beta * newAcceleration);

            // Update velocity.
            this.verticesVelocity[i] += Time.deltaTime * ((1f - gamma) * this.verticesAcceleration[i] + gamma * newAcceleration);

            // Update acceleration.
            this.verticesAcceleration[i] = newAcceleration;
        }
        mesh.vertices = vertices;
        mesh.RecalculateNormals();
    }

    private Mesh getMesh() {
        return this.shellObj.GetComponent<MeshFilter>().mesh;
    }

    private float getEdgeLength(int v1, int v2) {
        Vector3[] vertices = this.getMesh().vertices;
        Vector3 diff = vertices[v1] - vertices[v2];
        return diff.magnitude;
    }

    private float getEdgeLengthEnergy(int v1, int v2) {
        Vector3[] vertices = this.getMesh().vertices;
        Vector3 edge = vertices[v2] - vertices[v1]; // Vector from v1 to v2.
        Vector3 undeformedEdge = this.originalVertices[v2] - this.originalVertices[v1];
        float edgeLength = edge.magnitude;
        float undeformedEdgeLength = undeformedEdge.magnitude;
        return Mathf.Pow(1f - edgeLength / undeformedEdgeLength, 2) * undeformedEdgeLength;
    }

    /**
     * Computes the edge length energy gradient of vertex v1, for the edge between vertices v1 and v2.
     */
    private Vector3 getEdgeLengthEnergyGradient(int v1, int v2) {
        Vector3[] vertices = this.getMesh().vertices;
        Vector3 edge = vertices[v2] - vertices[v1]; // Vector from v1 to v2.
        Vector3 undeformedEdge = this.originalVertices[v2] - this.originalVertices[v1];
        float edgeLength = edge.magnitude;
        float undeformedEdgeLength = undeformedEdge.magnitude;
        Vector3 dEdgeLength = (vertices[v1] - vertices[v2]) / edgeLength;
        return (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength;
    }

    /**
     * Computes the triangle area energy gradient of vertex v1, for the triangle defined by vertices v1, v2 and v3.
     * Note that 1/3 of the area energy is used for v1 (the other two parts will be used for v2 and v3).
     */
    private Vector3 getTriangleAreaEnergyGradient(int v1, int v2, int v3) {
        Vector3[] vertices = this.getMesh().vertices;

        // Get two edges, where only one is dependent on vertex v1.
        Vector3 edge21 = vertices[v1] - vertices[v2]; // dEdge21 / dv1 = {1, 1, 1}.
        Vector3 edge23 = vertices[v3] - vertices[v2]; // dEdge23 / dv1 = {0, 0, 0}.

        // Calculate the triangle area gradient.
        float crossProdLength = Vector3.Cross(edge21, edge23).magnitude;
        Vector3 dCrossProdLength = 1f / crossProdLength * new Vector3(
                (edge21.x * edge23.z - edge23.x * edge21.z) * edge23.z + (edge21.x * edge23.y - edge23.x * edge21.y) * edge23.y,
                (edge21.y * edge23.z - edge23.y * edge21.z) * edge23.z + (edge21.x * edge23.y - edge23.x * edge21.y) * -edge23.x,
                (edge21.y * edge23.z - edge23.y * edge21.z) * -edge23.y + (edge21.x * edge23.z - edge23.x * edge21.z) * -edge23.x);
        float triangleArea = crossProdLength / 6f; // Area of triangle is half the cross product length, and we only look at a third.
        Vector3 dTriangleArea = dCrossProdLength / 6f; // Area of triangle is half the cross product length, and we only look at a third.

        // Calculate undeformed triangle area.
        // TODO - Determine beforehand and store for performance?
        float undeformedTriangleArea = Vector3.Cross(
                this.originalVertices[v1] - this.originalVertices[v2],
                this.originalVertices[v3] - this.originalVertices[v2]
            ).magnitude / 6f; // Area of triangle is half the cross product length, and we only look at a third.

        // Calculate the area energy gradient.
        return (2 * triangleArea / undeformedTriangleArea - 2) * dTriangleArea;
    }

    /**
     * Computes the bending energy gradient of the adjacent triangles defined by vertices v1_ and v2_.
     * Vertices ve1 and ve2 are the vertices that define the edge that is shared between the two triangles.
     * Edge ve1 -> ve2 belongs to vertex v1_ and edge ve2 -> ve1 belongs to vertex v2_ (this matters for the direction of the normals).
     */
    private Vector3 getBendingEnergyGradient(int v11, int v12, int v13, int v21, int v22, int v23, int ve1, int ve2) {
        Vector3[] vertices = this.getMesh().vertices;

        // Shared edge, clockwise for triangle 1, anti-clockwise for triangle 2.
        Vector3 e1 = vertices[ve2] - vertices[ve1];
        Vector3 e1_undeformed = this.originalVertices[ve2] - this.originalVertices[ve1];

        // Edge in triangle 1 that does not include ve1.
        Vector3 e2 = (v11 == ve1 ? vertices[v13] - vertices[v12] : (v12 == ve1 ? vertices[v11] - vertices[v13] : vertices[v12] - vertices[v11]));
        Vector3 e2_undeformed = (v11 == ve1 ? this.originalVertices[v13] - this.originalVertices[v12]
                : (v12 == ve1 ? this.originalVertices[v11] - this.originalVertices[v13] : this.originalVertices[v12] - this.originalVertices[v11]));

        // Edge in triangle 2 that does not include ve1.
        Vector3 e3 = (v21 == ve1 ? vertices[v23] - vertices[v22] : (v22 == ve1 ? vertices[v21] - vertices[v23] : vertices[v22] - vertices[v21]));
        Vector3 e3_undeformed = (v21 == ve1 ? this.originalVertices[v23] - this.originalVertices[v22]
                : (v22 == ve1 ? this.originalVertices[v21] - this.originalVertices[v23] : this.originalVertices[v22] - this.originalVertices[v21]));

        // Triangle normals.
        Vector3 cross_e1_e2 = Vector3.Cross(e1, e2);
        Vector3 cross_e1_e3 = Vector3.Cross(e1, e3);
        Vector3 n1 = cross_e1_e2.normalized;
        Vector3 n2 = cross_e1_e3.normalized;

        // Calculate d_teta_d_ve1, based on rewriting d_teta_d_x1 in paper: http://ddg.math.uni-goettingen.de/pub/bendingCAGD.pdf
        float dot_e1_norm_e2_norm = Vector3.Dot(e1.normalized, e2.normalized);
        float dot_e1_norm_e3_norm = Vector3.Dot(e1.normalized, e3.normalized);
        if(Mathf.Abs(dot_e1_norm_e2_norm) == 1f || dot_e1_norm_e3_norm == 1f) {
            return new Vector3(0f, 0f, 0f); // Triangle vertices are on a single line, gradient is 0 here.
        }
        Vector3 d_teta_d_ve1 = -1 / e1.magnitude * (dot_e1_norm_e2_norm / Mathf.Sqrt(1 - (dot_e1_norm_e2_norm * dot_e1_norm_e2_norm)) * n1
                - dot_e1_norm_e3_norm / Mathf.Sqrt(1 - (dot_e1_norm_e3_norm * dot_e1_norm_e3_norm)) * n2);

        // Undeformed triangle normals. These don't have to be using the same edges, as long as they are clockwise as well (or have a minus sign).
        Vector3 n1_undeformed = Vector3.Cross(this.originalVertices[v12] - this.originalVertices[v11],
                this.originalVertices[v13] - this.originalVertices[v12]).normalized;
        Vector3 n2_undeformed = Vector3.Cross(this.originalVertices[v22] - this.originalVertices[v21],
                this.originalVertices[v23] - this.originalVertices[v22]).normalized;

        // Angle between triangle normals.
        Vector3 v_triangle2_unshared = (v21 == ve1 ? vertices[v23] : (v22 == ve1 ? vertices[v21] : vertices[v22]));
        float teta_e_sign = Mathf.Sign(Vector3.Dot(n1, v_triangle2_unshared - vertices[ve1])); // 1 if teta_e positive, -1 if negative.
        float teta_e = Mathf.Acos(Vector3.Dot(n1, n2)) * teta_e_sign;
        Vector3 v_triangle2_undeformed_unshared = (v21 == ve1 ? this.originalVertices[v23] : (v22 == ve1 ? this.originalVertices[v21] : this.originalVertices[v22]));
        float teta_e_undeformed_sign = Mathf.Sign(Vector3.Dot(n1_undeformed, v_triangle2_undeformed_unshared - this.originalVertices[ve1])); // 1 if teta_e_undeformed positive, -1 if negative.
        float teta_e_undeformed = Mathf.Acos(Vector3.Dot(n1_undeformed, n2_undeformed)) * teta_e_undeformed_sign;

        // bending energy gradient.
        float h_e_undeformed = (cross_e1_e2.magnitude + cross_e1_e3.magnitude) / e1.magnitude / 6f;
        float d_W_bending_energy_edge_d_teta_e = 2 * (teta_e - teta_e_undeformed) * e1_undeformed.magnitude / h_e_undeformed;

        // Return the result.
        return d_W_bending_energy_edge_d_teta_e * d_teta_d_ve1;
    }
}
