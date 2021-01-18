using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Shell : MonoBehaviour {

    // Material property representing constants.
    public float kLength;
    public float kArea;
    public float kBend;
    public float shellThickness = 0.002f; // [m].
    public float shellMaterialDensity = 40f; // [kg/m^3].
    public bool useFlatUndeformedBendState = true; // When true, assumes the bending energy to be 0 when the object is flat.

    // Game object to simulate.
    public GameObject shellObj = null;

    // Vertex related fields.
    private List<int>[] vertexTriangles;
    private Vector3[] originalVertices; // Vertices in undeformed state.
    private Vector3[] verticesVelocity;
    private Vector3[] verticesAcceleration;
    private bool[] verticesMovementConstraints; // When true, movement for the corresponding vertex is prohibited.

    // Simulation update loop settings.
    private bool doUpdate = false;
    private bool doGradientDescent = true;
    public float kGradientDescent;
    public float maxGradientDescentStep;
    public float timeScale = 1f;
    public Vector3 windPressure; // [N/m^2]. TODO - Could also apply scalar pressure in triangle normal directions.

    // Cached vertex/triangle properties.
    private Vector3[] triangleNormals;
    private Vector3[] undeformedTriangleNormals;
    private float[] triangleAreas;
    private float[] undeformedTriangleAreas;

    // Start is called before the first frame update.
    void Start() {

        // Initialize fields. Doing this overwrites the values set in Unity's inspector.
        this.kLength = 10f;
        this.kArea = 10f;
        this.kBend = 0.1f;
        this.kGradientDescent = 0.1f;
        this.maxGradientDescentStep = 0.01f;
        this.windPressure = new Vector3(0f, 0f, 100f);

        // Create the new object in the scene.
        Mesh mesh;
        if(this.shellObj == null) {
            this.shellObj = new GameObject();
            this.shellObj.AddComponent<MeshRenderer>();
            this.shellObj.AddComponent<MeshFilter>();
            mesh = MeshHelper.createSquareMesh(2);
            //Mesh mesh = MeshHelper.createSingleTriangleMesh();
        } else {
            mesh = this.shellObj.GetComponent<MeshFilter>().mesh;
        }
        MeshRenderer meshRenderer = this.shellObj.GetComponent<MeshRenderer>();
        meshRenderer.sharedMaterial = new Material(Shader.Find("Custom/StandardTwoSides"));
        MeshFilter meshFilter = this.shellObj.GetComponent<MeshFilter>();

        // Store undeformed vertices.
        this.originalVertices = (Vector3[]) mesh.vertices.Clone();

        // Perform mesh deformation and translation (this should be replaced with some meaningful force/deformation).
        Vector3[] verts = mesh.vertices;

        //verts[0].x += 2.5f;
        //verts[0].y += 2.5f;
        //verts[0].z -= 5f;

        mesh.vertices = verts;
        mesh.RecalculateNormals();

        // Add the mesh to the mesh filter.
        meshFilter.mesh = mesh;

        // Set shell position.
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
                        int otherVertexClockwiseInd2 = (vertexInd == v21 ? v23 : (vertexInd == v23 ? v22 : v21));
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

                // Replace vertex triangles list with sorted vertex triangles list.
                this.vertexTriangles[vertexInd] = sortedTriangleList;
            }
        }

        // Initialize additional vertex and triangle data.
        int numVertices = mesh.vertices.Length;
        this.verticesVelocity = new Vector3[numVertices];
        this.verticesAcceleration = new Vector3[numVertices];
        this.verticesMovementConstraints = new bool[numVertices];
        for(int i = 0; i < numVertices; i++) {
            this.verticesVelocity[i] = new Vector3(0, 0, 0);
            this.verticesAcceleration[i] = new Vector3(0, 0, 0);
            this.verticesMovementConstraints[i] = false;
        }
        int numTriangles = mesh.triangles.Length;
        this.triangleNormals = new Vector3[numTriangles];
        this.triangleAreas = new float[numTriangles];
        this.undeformedTriangleNormals = new Vector3[numTriangles];
        this.undeformedTriangleAreas = new float[numTriangles];
        int[] triangles = mesh.triangles;
        Vector3[] vertices = mesh.vertices;
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];
            Vector3 crossProd = Vector3.Cross(vertices[v2] - vertices[v1], vertices[v3] - vertices[v1]);
            float crossProdMag = crossProd.magnitude;

            // Store the triangle normal and area if they exist (i.e. if the triangle has a normal and therefore an area).
            if(float.IsNaN(crossProdMag)) {
                this.undeformedTriangleNormals[triangleId] = Vector3.zero;
                this.undeformedTriangleAreas[triangleId] = 0f;
            } else {
                this.undeformedTriangleNormals[triangleId] = crossProd / crossProdMag;
                this.undeformedTriangleAreas[triangleId] = crossProdMag / 2f; // Triangle area is half of the cross product of any two of its edges.
            }
        }

        // Add movement constraints on mesh edge vertices.
        for(int vertexInd = 0; vertexInd < this.vertexTriangles.Length; vertexInd++) {
            List<int> triangleList = this.vertexTriangles[vertexInd];
            for(int i = 0; i < triangleList.Count; i++) {

                // Get triangle vertices.
                int triangleBaseInd1 = triangleList[i] * 3;
                int v11 = mesh.triangles[triangleBaseInd1];
                int v12 = mesh.triangles[triangleBaseInd1 + 1];
                int v13 = mesh.triangles[triangleBaseInd1 + 2];

                // Get next triangle vertices.
                int triangleBaseInd2 = triangleList[(i + 1) % triangleList.Count] * 3;
                int v21 = mesh.triangles[triangleBaseInd2];
                int v22 = mesh.triangles[triangleBaseInd2 + 1];
                int v23 = mesh.triangles[triangleBaseInd2 + 2];

                // Get the vertex indices of the other vertices that are connected to this vertex's edges.
                int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
                int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));

                // Constrain the vertex if a gap between triangles has been found.
                if(otherVertexClockwiseInd1 != otherVertexAntiClockwiseInd2) {
                    this.verticesMovementConstraints[vertexInd] = true;
                }
            }
        }
        //mesh.vertices = vertices;
    }

    private void reset() {
        Mesh mesh = this.getMesh();
        Vector3[] vertices = mesh.vertices;
        for(int i = 0; i < mesh.vertices.Length; i++) {
            this.verticesVelocity[i] = new Vector3(0, 0, 0);
            this.verticesAcceleration[i] = new Vector3(0, 0, 0);
            vertices[i] = new Vector3(this.originalVertices[i].x, this.originalVertices[i].y, this.originalVertices[i].z);
        }
        mesh.vertices = vertices;
        mesh.RecalculateNormals();
    }

    // Update is called once per frame.
    private void Update() {

        // Handle reset hotkey.
        if(Input.GetKeyDown(KeyCode.R)) {
            this.reset();
            print("Mesh reset.");
        }

        // Handle simulation pause/resume hotkey.
        if(Input.GetKeyDown(KeyCode.P)) {
            this.doUpdate = !this.doUpdate;
            print((this.doUpdate ? "Simulation resumed." : "Simulation paused."));
        }

        // Handle single step hotkey.
        if(Input.GetKeyDown(KeyCode.G)) {
            this.simulationStep(Time.fixedDeltaTime);
            print("Performed single step (" + Time.fixedDeltaTime + "s).");
        }
    }

    // FixedUpdate is called every fixed interval (Edit -> Project Settings -> Time -> Fixed Timestep).
    void FixedUpdate() {

        // Don't update if no noticable time has passed or when the simulation has been paused.
        if(Time.deltaTime == 0 || Time.timeScale == 0 || (!this.doUpdate && !Input.GetKey(KeyCode.F))) {
            return;
        }

        // Perform a simulation step.
        this.simulationStep(Time.deltaTime);
    }

    private void simulationStep(float deltaTime) {
        deltaTime *= this.timeScale;

        // Get the mesh.
        Mesh mesh = this.getMesh();
        int[] triangles = mesh.triangles;
        Vector3[] vertices = mesh.vertices;

        // Compute triangle normals and areas.
        this.recalcTriangleNormalsAndAreas(triangles, vertices);

        // Calculate vertex energy gradient array and vertex wind force array.
        Vector3[] vertexEnergyGradient = this.calcVertexEnergyGradient(triangles, vertices);
        Vector3[] vertexWindForce = this.calcVertexWindForce(triangles, vertices);

        // Update the vertices using discrete integration or gradient descent.
        if(this.doGradientDescent) {

            // Perform grafient descent.
            for(int i = 0; i < vertices.Length; i++) {

                // Skip vertex if it is constrained.
                if(this.verticesMovementConstraints[i]) {
                    continue;
                }

                // Calculate vertex area (a third of the area of triangles that this vertex is part of).
                float vertexArea = 0f;
                foreach(int triangleId in this.vertexTriangles[i]) {
                    vertexArea += this.triangleAreas[triangleId];
                }
                vertexArea /= 3f;

                // Update vertex position.
                Vector3 step = this.kGradientDescent * (vertexWindForce[i] - vertexEnergyGradient[i]);
                if(step.magnitude > maxGradientDescentStep) {
                    step = step.normalized * maxGradientDescentStep;
                }
                vertices[i] += step;
            }

        } else {

            // Perform Newmark Time Stepping (ODE integration).
            float gamma = 0.5f;
            float beta = 0.25f;
            for(int i = 0; i < vertices.Length; i++) {

                // Skip vertex if it is constrained.
                if(this.verticesMovementConstraints[i]) {
                    this.verticesVelocity[i] = new Vector3(0f, 0f, 0f);
                    this.verticesAcceleration[i] = new Vector3(0f, 0f, 0f);
                    continue;
                }

                // Calculate lumped vertex mass (a third of the area of triangles that this vertex is part of).
                float vertexArea = 0f;
                foreach(int triangleId in this.vertexTriangles[i]) {
                    vertexArea += this.triangleAreas[triangleId];
                }
                vertexArea /= 3f;
                float mass = vertexArea * this.shellThickness * this.shellMaterialDensity;

                // Calculate vertex acceleration.
                Vector3 newAcceleration = (vertexWindForce[i] - vertexEnergyGradient[i]) / mass;

                // Update vertex position.
                vertices[i] += deltaTime * this.verticesVelocity[i]
                        + deltaTime * deltaTime * ((1f / 2f - beta) * this.verticesAcceleration[i] + beta * newAcceleration);

                // Update vertex velocity.
                this.verticesVelocity[i] += deltaTime * ((1f - gamma) * this.verticesAcceleration[i] + gamma * newAcceleration);
                this.verticesVelocity[i] *= 0.9f; // TODO - Replace this constant damping with something more realistic friction-based damping.

                // Update vertex acceleration.
                this.verticesAcceleration[i] = newAcceleration;
            }
        }
    }

    private Mesh getMesh() {
        return this.shellObj.GetComponent<MeshFilter>().mesh;
    }

    private void recalcTriangleNormalsAndAreas(int[] triangles, Vector3[] vertices) {
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];
            Vector3 crossProd = Vector3.Cross(vertices[v2] - vertices[v1], vertices[v3] - vertices[v1]);
            float crossProdMag = crossProd.magnitude;

            // Store the triangle normal and area if they exist (i.e. if the triangle has a normal and therefore an area).
            if(float.IsNaN(crossProdMag)) {
                this.triangleNormals[triangleId] = Vector3.zero;
                this.triangleAreas[triangleId] = 0f;
            } else {
                this.triangleNormals[triangleId] = crossProd / crossProdMag;
                this.triangleAreas[triangleId] = crossProdMag / 2f; // Triangle area is half of the cross product of any two of its edges.
            }
        }
    }

    private Vector3[] calcVertexEnergyGradient(int[] triangles, Vector3[] vertices) {

        // Initialize vertex energy gradient array.
        Vector3[] vertexEnergyGradient = new Vector3[vertices.Length];
        for(int i = 0; i < vertexEnergyGradient.Length; i++) {
            vertexEnergyGradient[i] = Vector3.zero;
        }

        // Compute vertex energy gradient array.
        for(int vertexInd = 0; vertexInd < this.vertexTriangles.Length; vertexInd++) {
            for(int i = 0; i < this.vertexTriangles[vertexInd].Count; i++) {

                // Get two possibly adjacent triangles in clockwise direction.
                int triangleId = this.vertexTriangles[vertexInd][i];
                int nextTriangleId = this.vertexTriangles[vertexInd][(i + 1) % this.vertexTriangles[vertexInd].Count];

                // Get triangle vertices.
                int triangleBaseInd1 = triangleId * 3;
                int triangleBaseInd2 = nextTriangleId * 3;
                int v11 = triangles[triangleBaseInd1];
                int v12 = triangles[triangleBaseInd1 + 1];
                int v13 = triangles[triangleBaseInd1 + 2];
                int v21 = triangles[triangleBaseInd2];
                int v22 = triangles[triangleBaseInd2 + 1];
                int v23 = triangles[triangleBaseInd2 + 2];

                // Get the vertex indices of the other vertices that are connected to the possibly shared triangle edge.
                int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
                int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));

                // Handle the edge, or both edges if they are not the same.
                bool edgeSharedByTriangles = (otherVertexClockwiseInd1 == otherVertexAntiClockwiseInd2);

                // Calculate spring energy gradient in the edge.
                vertexEnergyGradient[vertexInd] += kLength * this.getEdgeLengthEnergyGradient(vertices, vertexInd, otherVertexClockwiseInd1);

                if(edgeSharedByTriangles) {

                    // Calculate bending energy gradient.
                    vertexEnergyGradient[vertexInd] += kBend * this.getBendingEnergyGradient(
                        vertices, triangleId, nextTriangleId, v11, v12, v13, v21, v22, v23, vertexInd, otherVertexClockwiseInd1);
                } else {

                    // Calculate spring energy gradient in the second edge.
                    vertexEnergyGradient[vertexInd] += kLength * this.getEdgeLengthEnergyGradient(vertices, vertexInd, otherVertexAntiClockwiseInd2);
                }

                // Calculate the area energy gradient in the triangle.
                if(vertexInd == v11) {
                    vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v11, v12, v13);
                } else if(vertexInd == v12) {
                    vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v12, v13, v11);
                } else {
                    vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v13, v11, v12);
                }
            }
        }
        return vertexEnergyGradient;
    }

    private Vector3[] calcVertexWindForce(int[] triangles, Vector3[] vertices) {

        // Initialize vertex wind force array.
        Vector3[] vertexWindForce = new Vector3[vertices.Length];
        for(int i = 0; i < vertexWindForce.Length; i++) {
            vertexWindForce[i] = Vector3.zero;
        }

        // Compute vertex wind force array.
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];

            // Compute projection of the wind pressure vector on the triangle normal.
            Vector3 triangleNormal = this.triangleNormals[triangleId];
            float triangleArea = this.triangleAreas[triangleId];
            Vector3 totalTriangleWindForce = Vector3.Dot(this.windPressure, triangleNormal) * triangleArea * triangleNormal;

            // Add a third of the total triangle wind force to each of its vertices.
            vertexWindForce[v1] += totalTriangleWindForce / 3f;
            vertexWindForce[v2] += totalTriangleWindForce / 3f;
            vertexWindForce[v3] += totalTriangleWindForce / 3f;
        }
        return vertexWindForce;
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
    private Vector3 getEdgeLengthEnergyGradient(Vector3[] vertices, int v1, int v2) {
        Vector3 edge = vertices[v2] - vertices[v1]; // Vector from v1 to v2.
        Vector3 undeformedEdge = this.originalVertices[v2] - this.originalVertices[v1];
        float edgeLength = edge.magnitude;
        if(float.IsNaN(edgeLength)) {
            return Vector3.zero; // Edge is zero-length, so the gradient is 0.
        }
        float undeformedEdgeLength = undeformedEdge.magnitude;
        Vector3 dEdgeLength = (vertices[v1] - vertices[v2]) / edgeLength;
        return (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength;
    }

    /**
     * Computes the triangle area energy gradient of vertex v1, for the triangle defined by vertices v1, v2 and v3.
     * Note that 1/3 of the area energy is used for v1 (the other two parts will be used for v2 and v3).
     */
    private Vector3 getTriangleAreaEnergyGradient(Vector3[] vertices, int triangleId, int v1, int v2, int v3) {

        // Get two edges, where only one is dependent on vertex v1.
        Vector3 edge21 = vertices[v1] - vertices[v2]; // dEdge21 / dv1 = {1, 1, 1}.
        Vector3 edge23 = vertices[v3] - vertices[v2]; // dEdge23 / dv1 = {0, 0, 0}.

        // Calculate the triangle area gradient.
        if(this.triangleAreas[triangleId] == 0f) {
            return Vector3.zero; // Area is 0 m^2, so the gradient is 0.
        }
        float crossProdLength = this.triangleAreas[triangleId] * 2f;
        Vector3 dCrossProdLength = 1f / crossProdLength * new Vector3(
                (edge21.x * edge23.z - edge23.x * edge21.z) * edge23.z + (edge21.x * edge23.y - edge23.x * edge21.y) * edge23.y,
                (edge21.y * edge23.z - edge23.y * edge21.z) * edge23.z + (edge21.x * edge23.y - edge23.x * edge21.y) * -edge23.x,
                (edge21.y * edge23.z - edge23.y * edge21.z) * -edge23.y + (edge21.x * edge23.z - edge23.x * edge21.z) * -edge23.x);
        Vector3 dTriangleArea = dCrossProdLength / 6f; // Area of triangle is half the cross product length, and we only look at a third.

        // Calculate the area energy gradient.
        return (2 * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2) * dTriangleArea;
    }

    /**
     * Computes the bending energy gradient of the adjacent triangles defined by vertices v1_ and v2_.
     * Vertices ve1 and ve2 are the vertices that define the edge that is shared between the two triangles.
     * Edge ve1 -> ve2 belongs to vertex v1_ and edge ve2 -> ve1 belongs to vertex v2_ (this matters for the direction of the normals).
     */
    private Vector3 getBendingEnergyGradient(Vector3[] vertices, int triangleId1, int triangleId2,
            int v11, int v12, int v13, int v21, int v22, int v23, int ve1, int ve2) {

        // Shared edge, clockwise for triangle 1, anti-clockwise for triangle 2.
        Vector3 e1 = vertices[ve2] - vertices[ve1];
        Vector3 e1_undeformed = this.originalVertices[ve2] - this.originalVertices[ve1];

        // Edge in triangle 1 that does not include ve1.
        Vector3 e2 = (v11 == ve1 ? vertices[v13] - vertices[v12] : (v12 == ve1 ? vertices[v11] - vertices[v13] : vertices[v12] - vertices[v11]));

        // Edge in triangle 2 that does not include ve1.
        Vector3 e3 = (v21 == ve1 ? vertices[v23] - vertices[v22] : (v22 == ve1 ? vertices[v21] - vertices[v23] : vertices[v22] - vertices[v21]));

        // Return if any of the vertices overlap. This means that no triangle normals exist.
        if(e1.magnitude == 0f || e2.magnitude == 0f || e3.magnitude == 0f) {
            return Vector3.zero;
        }

        // Triangle normals.
        Vector3 n1 = this.triangleNormals[triangleId1];
        Vector3 n2 = this.triangleNormals[triangleId2];

        // Calculate d_teta_d_ve1, based on rewriting d_teta_d_x1 in paper: http://ddg.math.uni-goettingen.de/pub/bendingCAGD.pdf
        float dot_e1_norm_e2_norm = Vector3.Dot(e1.normalized, e2.normalized);
        float dot_e1_norm_e3_norm = Vector3.Dot(e1.normalized, e3.normalized);
        float dot_e1_norm_e2_norm_square = dot_e1_norm_e2_norm * dot_e1_norm_e2_norm;
        float dot_e1_norm_e3_norm_square = dot_e1_norm_e3_norm * dot_e1_norm_e3_norm;
        Vector3 d_teta_d_ve1 = -1 / e1.magnitude * (dot_e1_norm_e2_norm / Mathf.Sqrt(1f - dot_e1_norm_e2_norm_square) * n1
                - dot_e1_norm_e3_norm / Mathf.Sqrt(1f - dot_e1_norm_e3_norm_square) * n2);
        if(float.IsNaN(d_teta_d_ve1.x) || float.IsNaN(d_teta_d_ve1.y) || float.IsNaN(d_teta_d_ve1.z)
                || float.IsInfinity(d_teta_d_ve1.x) || float.IsInfinity(d_teta_d_ve1.y) || float.IsInfinity(d_teta_d_ve1.z)) {
            return Vector3.zero; // Triangle vertices are on a single line. Gradient is 0 here.
        }

        // Angle between triangle normals.
        Vector3 v_triangle2_unshared = (v21 == ve1 ? vertices[v23] : (v22 == ve1 ? vertices[v21] : vertices[v22]));
        float teta_e_sign = Mathf.Sign(Vector3.Dot(n1, v_triangle2_unshared - vertices[ve1])); // 1 if teta_e positive, -1 if negative.
        float teta_e = Mathf.Acos(Vector3.Dot(n1, n2)) * teta_e_sign;
        if(float.IsNaN(teta_e)) {
            teta_e = Mathf.PI; // The triangles are on top of each other, which is both 180 and -180 degrees.
        }

        // Angle between undeformed triangle normals, or 0 when assuming a flat rest state for bending.
        float teta_e_undeformed;
        if(!this.useFlatUndeformedBendState) {

            // Undeformed triangle normals. These don't have to be using the same edges, as long as they are clockwise as well (or have a minus sign).
            Vector3 n1_undeformed = this.undeformedTriangleNormals[triangleId1];
            Vector3 n2_undeformed = this.undeformedTriangleNormals[triangleId2];

            // Angle between undeformed triangle normals.
            Vector3 v_triangle2_undeformed_unshared = (v21 == ve1 ? this.originalVertices[v23] : (v22 == ve1 ? this.originalVertices[v21] : this.originalVertices[v22]));
            float teta_e_undeformed_sign = Mathf.Sign(Vector3.Dot(n1_undeformed, v_triangle2_undeformed_unshared - this.originalVertices[ve1])); // 1 if teta_e_undeformed positive, -1 if negative.
            teta_e_undeformed = Mathf.Acos(Vector3.Dot(n1_undeformed, n2_undeformed)) * teta_e_undeformed_sign;
        } else {
            teta_e_undeformed = 0f;
        }

        // bending energy gradient.
        // h_e_undeformed is a third of the average triangle height, where the height is twice the triangle area divided by the triangle width.
        float h_e_undeformed = (this.undeformedTriangleAreas[triangleId1] + this.undeformedTriangleAreas[triangleId2]) / e1_undeformed.magnitude / 3f;
        float d_W_bending_energy_edge_d_teta_e = 2f * (teta_e - teta_e_undeformed) * e1_undeformed.magnitude / h_e_undeformed;

        // Return the result.
        return d_W_bending_energy_edge_d_teta_e * d_teta_d_ve1;
    }
}
