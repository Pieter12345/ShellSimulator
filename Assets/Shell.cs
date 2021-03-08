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
    private List<int>[] sortedVertexTriangles;
    private List<Edge> edges;
    private Vector3[] originalVertices; // Vertices in undeformed state.
    private bool[] verticesMovementConstraints; // When true, movement for the corresponding vertex is prohibited.

    private Vec3D[] vertexPositions; // Format: {{x1, y1, z1}, {x2, y2, z2}, ...}.
    private VecD vertexVelocities; // Format: {vx1, vy1, vz1, vx2, v2y, v2z, ...}.
    private VecD vertexAccelerations; // Format: {ax1, ay1, az1, ax2, a2y, a2z, ...}.

    // Simulation update loop settings.
    public TimeSteppingMethod timeSteppingMethod = TimeSteppingMethod.GRADIENT_DESCENT;
    private bool doUpdate = false;
    public double kGradientDescent;
    public double maxGradientDescentStep;
    public double timeScale = 1f;
    public float dampingFactor = 0.99f; // [F * s / m] = [kg / s] ? Factor applied to vertex velocity per time step. TODO - Replace with proper energy dissipation.
    public Vector3 windPressure; // [N/m^2]. TODO - Could also apply scalar pressure in triangle normal directions.

    // Cached vertex/triangle properties.
    private Vec3D[] triangleNormals;
    private Vector3[] triangleNormals_old;
    private Vector3[][] dTriangleNormals_dv1;
    private Vector3[][] dTriangleNormals_dv2;
    private Vector3[][] dTriangleNormals_dv3;
    private Vector3[] undeformedTriangleNormals;
    private float[] triangleAreas;
    private Vector3[] dTriangleAreas_dv1;
    private Vector3[] dTriangleAreas_dv2;
    private Vector3[] dTriangleAreas_dv3;
    private float[] undeformedTriangleAreas;

    void Awake() {
        QualitySettings.vSyncCount = 0; // Disable V-sync.
        Application.targetFrameRate = 100; // Set max framerate.
        System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US"); // To print decimal points instead of commas.
    }

    // Start is called before the first frame update.
    void Start() {

        // Initialize fields. Doing this overwrites the values set in Unity's inspector.
        this.kLength = 0f;
        this.kArea = 10f;
        this.kBend = 0f;
        this.kGradientDescent = 1f;
        this.maxGradientDescentStep = 0.01f;
        this.windPressure = new Vector3(0f, 0f, 10f);

        // Create the new object in the scene.
        Mesh mesh;
        if(this.shellObj == null) {
            this.shellObj = new GameObject();
            this.shellObj.AddComponent<MeshRenderer>();
            this.shellObj.AddComponent<MeshFilter>();
            //mesh = MeshHelper.createTriangleMesh(5, 5, 0);
            //mesh = MeshHelper.createSquareMesh(5, 5, 1);
            mesh = MeshHelper.createTriangleMesh(5, 5, 3); // 5 subdivisions leads to 561 vertices and 3072 triangles.
        } else {
            mesh = this.shellObj.GetComponent<MeshFilter>().mesh;
        }
        MeshRenderer meshRenderer = this.shellObj.GetComponent<MeshRenderer>();
        meshRenderer.sharedMaterial = new Material(Shader.Find("Custom/StandardTwoSides"));
        MeshFilter meshFilter = this.shellObj.GetComponent<MeshFilter>();
        print("Mesh loaded with " + mesh.vertices.Length + " vertices and " + mesh.triangles.Length + " triangles.");

        // Store undeformed vertices.
        this.originalVertices = (Vector3[]) mesh.vertices.Clone();

        // Perform mesh deformation and translation (this should be replaced with some meaningful force/deformation).
        Vector3[] verts = mesh.vertices;

        //verts[0].x += 2.5f;
        //verts[0].y += 2.5f;
        //verts[0].z -= 5f;
        //verts[1].y += 1f;
        //verts[1].x += 1f;
        //verts[1].z += 1f;

        mesh.vertices = verts;
        mesh.RecalculateNormals();

        // Add the mesh to the mesh filter.
        meshFilter.mesh = mesh;

        // Get mesh variables for convenience.
        int[] triangles = mesh.triangles;
        Vector3[] vertices = mesh.vertices;
        int numVertices = mesh.vertexCount;

        // Set shell position.
        this.shellObj.transform.position = new Vector3(0, 1, 0);

        // Create clockwise sorted triangle list per vertex.
        this.sortedVertexTriangles = MeshUtils.getSortedVertexTriangles(triangles, numVertices);

        // Cache mesh edges.
        this.edges = MeshUtils.getEdges(this.sortedVertexTriangles, triangles);

        // Initialize additional vertex and triangle data.
        this.vertexPositions = new Vec3D[numVertices];
        this.verticesMovementConstraints = new bool[numVertices];
        for(int i = 0; i < numVertices; i++) {
            Vector3 vertex = vertices[i];
            this.vertexPositions[i] = new Vec3D(vertex.x, vertex.y, vertex.z);
            this.verticesMovementConstraints[i] = false;
        }
        this.vertexVelocities = new VecD(numVertices * 3);
        this.vertexAccelerations = new VecD(numVertices * 3);
        for(int i = 0; i < numVertices * 3; i++) {
            this.vertexVelocities[i] = 0;
            this.vertexAccelerations[i] = 0;
        }
        int numTriangles = triangles.Length;
        this.triangleNormals_old = new Vector3[numTriangles];
        this.triangleNormals = new Vec3D[numTriangles];
        this.dTriangleNormals_dv1 = new Vector3[numTriangles][];
        this.dTriangleNormals_dv2 = new Vector3[numTriangles][];
        this.dTriangleNormals_dv3 = new Vector3[numTriangles][];
        this.triangleAreas = new float[numTriangles];
        this.dTriangleAreas_dv1 = new Vector3[numTriangles];
        this.dTriangleAreas_dv2 = new Vector3[numTriangles];
        this.dTriangleAreas_dv3 = new Vector3[numTriangles];
        this.undeformedTriangleNormals = new Vector3[numTriangles];
        this.undeformedTriangleAreas = new float[numTriangles];
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];
            Vector3 crossProd = Vector3.Cross(this.originalVertices[v2] - this.originalVertices[v1], this.originalVertices[v3] - this.originalVertices[v1]);
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
        for(int vertexInd = 0; vertexInd < this.sortedVertexTriangles.Length; vertexInd++) {
            List<int> triangleList = this.sortedVertexTriangles[vertexInd];
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

        //// Add sail mesh specific movement constraints on the sail corner vertices.
        //Vector3[] vertices = mesh.vertices;
        //int vertMaxX = 0;
        //int vertMinX = 0;
        //int vertMaxY = 0;
        //int vertMinY = 0;
        //int vertMaxZ = 0;
        //int vertMinZ = 0;
        //for(int i = 0; i < numVertices; i++) {
        //    Vector3 vertex = vertices[i];
        //    if(vertex.x > vertices[vertMaxX].x) {
        //        vertMaxX = i;
        //    } else if(vertex.x < vertices[vertMinX].x) {
        //        vertMinX = i;
        //    }
        //    if(vertex.y > vertices[vertMaxY].y) {
        //        vertMaxY = i;
        //    } else if(vertex.y < vertices[vertMinY].y) {
        //        vertMinY = i;
        //    }
        //    if(vertex.z > vertices[vertMaxZ].z) {
        //        vertMaxZ = i;
        //    } else if(vertex.z < vertices[vertMinZ].z) {
        //        vertMinZ = i;
        //    }
        //}
        //this.verticesMovementConstraints[vertMinY] = true; // Mast/boom intersection.
        //this.verticesMovementConstraints[vertMaxX] = true; // Boom end.
        //this.verticesMovementConstraints[vertMinX] = true; // Mast top.
        ////vertices[vertMinY].z += 1f;
        ////vertices[vertMaxX].z += 2f;
        ////vertices[vertMinX].z += 3f;
        ////mesh.vertices = vertices;

        //// Visualize vertex constraints.
        //Vector3[] vertices = mesh.vertices;
        //for(int i = 0; i < this.verticesMovementConstraints.Length; i++) {
        //    if(this.verticesMovementConstraints[i]) {
        //        vertices[i].z += 2;
        //    }
        //}
        //mesh.vertices = vertices;
        //mesh.RecalculateNormals();
    }

    private void reset() {
        Mesh mesh = this.getMesh();
        Vector3[] vertices = mesh.vertices;
        for(int i = 0; i < mesh.vertices.Length; i++) {
            vertices[i] = new Vector3(this.originalVertices[i].x, this.originalVertices[i].y, this.originalVertices[i].z);

            this.vertexPositions[i] = new Vec3D(this.originalVertices[i].x, this.originalVertices[i].y, this.originalVertices[i].z);
            for(int coord = 0; coord < 3; coord++) {
                this.vertexVelocities[3 * i + coord] = 0;
            }
        }
        //this.vertexPositions[5][2] += 5; // TODO - Remove test.
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

    private void simulationStep(double deltaTime) {
        deltaTime *= this.timeScale;

        // Get the mesh.
        Mesh mesh = this.getMesh();
        int[] triangles = mesh.triangles;
        Vector3[] vertices = mesh.vertices;

        // Compute triangle normals and areas.
        this.recalcTriangleNormalsAndAreas(triangles, vertices);

        // Update the vertices using the chosen time stepping method.
        switch(this.timeSteppingMethod) {
            case TimeSteppingMethod.GRADIENT_DESCENT: {
                this.doGradientDescentStep();
                break;
            }
            case TimeSteppingMethod.EXPLICIT_NEWMARK: {
                this.doExplititIntegrationNewmarkStep(deltaTime);
                break;
            }
            case TimeSteppingMethod.IMPLICIT: {
                this.doImplicitIntegrationStep(deltaTime);
                break;
            }
        }
        

        // TODO - This implementation is wrong. Use things like the mass calculation if useful. Remove after having used all useful parts.
        /*
         * Implicit integration.
         * Solve: x(n+1) - deltaTime * x'(n) - deltaTime^2 * (windForce(x(n+1)) - energyGradient(x(n+1))) / mass(x(n+1)) - x(n) = f(x(n+1)) = 0
         * Parameters:
         *     x(n+1): Mesh vertex positions after this step (initialize with current positions).
         *     x'(n): Mesh vertex velocities before this step.
         *     (windForce(x(n+1)) - energyGradient(x(n+1))) / mass(x(n+1)) = f(x(n+1)): Acceleration calculated using vertex positions after this step.
         * Solving:
         *     Repeat x(n+1) -= f(x(n+1)) / f'(x(n+1)) until f(x(n+1)) is close enough to 0.
         *  Determine velocity afterwards (fill in): x'(n+1) = x'(n) + deltaTime * (windForce(x(n+1)) - energyGradient(x(n+1))) / mass(x(n+1))
         * /
        
        // Calculate vertex energy gradient array and vertex wind force array.
        Vector3[] vertexEnergyGradient = this.calcVertexEnergyGradient(triangles, vertices);
        Vector3[] vertexWindForce = this.calcVertexWindForce(triangles, vertices);

        // As a first guess for implicit integration, we take the current vertex positions.
        Vector3[] newVertices = (Vector3[]) vertices.Clone();
        Vector3[] newVertexWindForce = vertexWindForce;
        Vector3[] newVertexEnergyGradient = vertexEnergyGradient;

        float absError = 0f;
        float maxNormalizedError = 0.1f;
        int iterations = 0;
        int maxIterations = 20;
        do {
            if(iterations++ >= maxIterations) {
                print("Terminating Newton's method since " + maxIterations + " iterations have been reached. Current normalized error: " + (absError / vertices.Length));
                newVertices = vertices; // Do not make any changes.
                break;
            }

            // Calculate energy Hessian and wind force Hessian for the approximation of the new position.
            Vector3[][] newEnergyHess = this.calcVertexEnergyHessian(triangles, newVertices);
            Vector3[][] newWindForceHess = this.calcVertexWindForceHessian(triangles, newVertices);

            // Make the Hessians positive definite. This affects the path that Newton's method takes, but not its target.
            // TODO - Test whether this works and potentially make changes to the method to make the Hessian positive definite.
            for(int i = 0; i < vertices.Length; i++) {
                makeHessPositiveDefinite(newEnergyHess[i]);
                makeHessPositiveDefinite(newWindForceHess[i]);
            }

            absError = 0f;
            for(int i = 0; i < vertices.Length; i++) {

                // Skip vertex if it is constrained.
                if(this.verticesMovementConstraints[i]) {
                    this.verticesVelocity[i] = new Vector3(0f, 0f, 0f);
                    this.verticesAcceleration[i] = new Vector3(0f, 0f, 0f);
                    continue;
                }

                // Calculate lumped vertex mass (a third of the area of triangles that this vertex is part of).
                float newVertexArea = 0f;
                float dTriangleArea_dx = 0f;
                float dTriangleArea_dy = 0f;
                float dTriangleArea_dz = 0f;
                foreach(int triangleId in this.sortedVertexTriangles[i]) {
                    newVertexArea += this.triangleAreas[triangleId];
                    if(i == triangles[triangleId * 3]) {
                        dTriangleArea_dx += this.dTriangleAreas_dv1[triangleId].x;
                        dTriangleArea_dy += this.dTriangleAreas_dv1[triangleId].y;
                        dTriangleArea_dz += this.dTriangleAreas_dv1[triangleId].z;
                    } else if(i == triangles[triangleId * 3 + 1]) {
                        dTriangleArea_dx += this.dTriangleAreas_dv2[triangleId].x;
                        dTriangleArea_dy += this.dTriangleAreas_dv2[triangleId].y;
                        dTriangleArea_dz += this.dTriangleAreas_dv2[triangleId].z;
                    } else if(i == triangles[triangleId * 3 + 2]) {
                        dTriangleArea_dx += this.dTriangleAreas_dv3[triangleId].x;
                        dTriangleArea_dy += this.dTriangleAreas_dv3[triangleId].y;
                        dTriangleArea_dz += this.dTriangleAreas_dv3[triangleId].z;
                    }
                }
                newVertexArea /= 3f;
                float dNewMass_dx = dTriangleArea_dx * this.shellThickness * this.shellMaterialDensity;
                float dNewMass_dy = dTriangleArea_dy * this.shellThickness * this.shellMaterialDensity;
                float dNewMass_dz = dTriangleArea_dz * this.shellThickness * this.shellMaterialDensity;
                float newMass = newVertexArea * this.shellThickness * this.shellMaterialDensity;

                // Calculate vertex acceleration.
                Vector3 newAcceleration = (newVertexWindForce[i] - newVertexEnergyGradient[i]) / newMass;

                // Calculate differential function f for the guessed/approximated vertex positions.
                Vector3 f = newVertices[i] - (float) deltaTime * this.verticesVelocity[i] - (float) (deltaTime * deltaTime) * newAcceleration - vertices[i];

                // Update absolute error.
                absError += f.magnitude;

                // Calculate differential function f gradient.
                Vector3[] dNewAcceleration = new Vector3[] {
                        (newWindForceHess[i][0] - newEnergyHess[i][0]) / newMass - (newVertexWindForce[i] - newVertexEnergyGradient[i]) * dNewMass_dx / (newMass * newMass),
                        (newWindForceHess[i][1] - newEnergyHess[i][1]) / newMass - (newVertexWindForce[i] - newVertexEnergyGradient[i]) * dNewMass_dy / (newMass * newMass),
                        (newWindForceHess[i][2] - newEnergyHess[i][2]) / newMass - (newVertexWindForce[i] - newVertexEnergyGradient[i]) * dNewMass_dz / (newMass * newMass)
                };
                Vector3[] df = new Vector3[] {
                        // d_f(u)_dux = {d_f_x(u)_dux, d_f_y(u)_dux, d_f_z(u)_dux}.
                        new Vector3(1, 0, 0) - (float) (deltaTime * deltaTime) * dNewAcceleration[0],

                        // d_f(u)_duy = {d_f_x(u)_duy, d_f_y(u)_duy, d_f_z(u)_duy}.
                        new Vector3(0, 1, 0) - (float) (deltaTime * deltaTime) * dNewAcceleration[1],

                        // d_f(u)_duz = {d_f_x(u)_duz, d_f_y(u)_duz, d_f_z(u)_duz}.
                        new Vector3(0, 0, 1) - (float) (deltaTime * deltaTime) * dNewAcceleration[2]
                };

                // Update newVertices following Newton's method: numVertices[i] -= f(u) / f'(u).
                // TODO - Validate that f'(u) is symmetrical, and otherwise use division if necessary.
                newVertices[i] -= f.x * df[0] + f.y * df[1] + f.z * df[2];
                //print("Stepping: " + (f.x * df[0] + f.y * df[1] + f.z * df[2]));

                // Update triangle normals and areas.
                this.recalcTriangleNormalsAndAreas(triangles, newVertices);

                // Update vertex energy gradient array and vertex wind force array.
                newVertexEnergyGradient = this.calcVertexEnergyGradient(triangles, newVertices);
                newVertexWindForce = this.calcVertexWindForce(triangles, newVertices);
            }
        } while(absError / vertices.Length > maxNormalizedError);

        // Update velocity: x'(n+1) = x'(n) + deltaTime * (windForce(x(n+1)) - energyGradient(x(n+1))) / mass(x(n+1))
        for(int i = 0; i < vertices.Length; i++) {

            // Skip vertex if it is constrained.
            if(this.verticesMovementConstraints[i]) {
                continue;
            }

            // Calculate lumped vertex mass (a third of the area of triangles that this vertex is part of).
            float newVertexArea = 0f;
            foreach(int triangleId in this.sortedVertexTriangles[i]) {
                newVertexArea += this.triangleAreas[triangleId];
            }
            newVertexArea /= 3f;
            float newMass = newVertexArea * this.shellThickness * this.shellMaterialDensity;

            // Update the velocity.
            this.verticesVelocity[i] += (float) deltaTime * (newVertexWindForce[i] - newVertexEnergyGradient[i]) / newMass;
        }

        // Update the vertices.
        mesh.vertices = newVertices;
        mesh.RecalculateNormals();
        */
    }

    private void doGradientDescentStep() {
        int[] triangles = this.getMesh().triangles;
        VecD vertexEnergyGradient = this.calcVertexEnergyGradient(triangles, this.vertexPositions);
        VecD vertexWindForce = this.calcVertexWindForce(triangles, this.vertexPositions);
        VecD step = this.kGradientDescent * (vertexWindForce - vertexEnergyGradient);
        for(int vertexInd = 0; vertexInd < this.vertexPositions.Length; vertexInd++) {
            if(!this.verticesMovementConstraints[vertexInd]) {
                VecD stepVec = new VecD(step[3 * vertexInd], step[3 * vertexInd + 1], step[3 * vertexInd + 2]);
                double stepVecMag = stepVec.magnitude;

                // Disallow NaN and infinite steps.
                if(stepVec.containsNaN()) {
                    print("NaN step: " + stepVec);
                    continue;
                }
                if(stepVec.containsInf()) {
                    print("Infinite step: " + step);
                    continue;
                }

                // Limit step size.
                if(stepVecMag > this.maxGradientDescentStep) {
                    stepVec = stepVec / stepVecMag * this.maxGradientDescentStep;
                }

                // Apply step.
                for(int coord = 0; coord < 3; coord++) {
                    this.vertexPositions[vertexInd][coord] += stepVec[coord];
                }
            }
        }
        this.updateMesh();
    }

    private void doExplititIntegrationNewmarkStep(double deltaTime) {
        
        // Calculate vertex energy gradient array and vertex wind force array.
        int[] triangles = this.getMesh().triangles;
        VecD vertexEnergyGradient = this.calcVertexEnergyGradient(triangles, this.vertexPositions);
        VecD vertexWindForce = this.calcVertexWindForce(triangles, this.vertexPositions);

        // Perform Newmark Time Stepping (ODE integration).
        double gamma = 0.5d;
        double beta = 0.25d;
        for(int i = 0; i < this.vertexPositions.Length; i++) {

            // Skip vertex if it is constrained.
            if(this.verticesMovementConstraints[i]) {
                this.vertexVelocities[3 * i] = 0;
                this.vertexVelocities[3 * i + 1] = 0;
                this.vertexVelocities[3 * i + 2] = 0;
                this.vertexAccelerations[3 * i] = 0;
                this.vertexAccelerations[3 * i + 1] = 0;
                this.vertexAccelerations[3 * i + 2] = 0;
                continue;
            }

            // Calculate lumped vertex mass (a third of the area of triangles that this vertex is part of).
            double vertexArea = 0f;
            foreach(int triangleId in this.sortedVertexTriangles[i]) {
                vertexArea += this.triangleAreas[triangleId];
            }
            vertexArea /= 3d;
            double mass = vertexArea * this.shellThickness * this.shellMaterialDensity;

            // Calculate vertex acceleration.
            Vec3D newAcceleration = new Vec3D();
            for(int coord = 0; coord < 3; coord++) {
                newAcceleration[coord] = (vertexWindForce[3 * i + coord] - vertexEnergyGradient[3 * i + coord]) / mass;
            }

            // Update vertex position.
            for(int coord = 0; coord < 3; coord++) {
                this.vertexPositions[i][coord] += deltaTime * this.vertexVelocities[3 * i + coord]
                        + deltaTime * deltaTime * ((1d / 2d - beta) * this.vertexAccelerations[3 * i + coord] + beta * newAcceleration[coord]);
            }

            // Update vertex velocity.
            for(int coord = 0; coord < 3; coord++) {
                this.vertexVelocities[3 * i + coord] += deltaTime * ((1d - gamma) * this.vertexAccelerations[3 * i + coord] + gamma * newAcceleration[coord]);
                this.vertexVelocities[3 * i + coord] *= dampingFactor; // TODO - Replace this constant damping with something more realistic friction-based damping.
            }

            // Update vertex acceleration.
            for(int coord = 0; coord < 3; coord++) {
                this.vertexAccelerations[3 * i + coord] = newAcceleration[coord];
            }
        }
        this.updateMesh();
    }

    private void doImplicitIntegrationStep(double deltaTime) {
        
        /*
         * Calculate all requires mesh things here, following the paper:
         * https://studios.disneyresearch.com/wp-content/uploads/2019/03/Discrete-Bending-Forces-and-Their-Jacobians-Paper.pdf
         * This should eventually be moved to methods and fields where convenient.
         */
        int[] triangles = this.getMesh().triangles;
        int numVertices = this.vertexPositions.Length;

        // TODO - Create a system-wide length energy Hessian.
        MatD forceHess = new MatD(numVertices * 3, numVertices * 3);
        foreach(Edge edge in this.edges) {

            // The length energy Hessian consists of 4 (3x3) parts that have to be inserted into the matrix.
            MatD lengthEnergyHess = this.getEdgeLengthEnergyHess(this.vertexPositions, edge.ve1, edge.ve2);
            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < 3; j++) {

                    // ddLengthEnergy_dv1_dv1.
                    forceHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] = lengthEnergyHess[i, j];

                    // ddLengthEnergy_dv2_dv2.
                    forceHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] = lengthEnergyHess[i + 3, j + 3];

                    // ddLengthEnergy_dv1_dv2.
                    forceHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] = lengthEnergyHess[i, j + 3];

                    // ddLengthEnergy_dv2_dv1.
                    forceHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] = lengthEnergyHess[i + 3, j];
                }
            }
        }

        // Calculate lumped vertex mass (a third of the area of triangles that each vertex is part of).
        // TODO - Look into adding constraints by making 1/mass === 0 for constrained variables (can be done per x, y, z of every vertex).
        MatD inverseMassMatrix = new MatD(numVertices * 3, numVertices * 3); // Diagonal matrix with pairs of 3 equal 1/mass_i (for each x, y and z).
        for(int i = 0; i < numVertices; i++) {
            double vertexArea = 0d;
            foreach(int triangleId in this.sortedVertexTriangles[i]) {
                vertexArea += this.triangleAreas[triangleId];
            }
            vertexArea /= 3d;
            double mass = vertexArea * this.shellThickness * this.shellMaterialDensity;
            double inverseMass = 1d / mass;
            inverseMassMatrix[3 * i, 3 * i] = inverseMass;
            inverseMassMatrix[3 * i + 1, 3 * i + 1] = inverseMass;
            inverseMassMatrix[3 * i + 2, 3 * i + 2] = inverseMass;
        }




        // Implicit differentiation following paper: https://www.cs.cmu.edu/~baraff/papers/sig98.pdf
        //VecD vertexForces = new VecD(vertices.Length * 3); // Format: {fx1, fy1, fz1, fx2, ...}.
        // TODO - vertexForces = vertexWindForce - vertexEnergyGradient; // TODO - THIS IS ESSENTIAL TO HAVE FORCES.
        VecD vertexEnergyGradient = this.calcVertexEnergyGradient(triangles, this.vertexPositions);
        VecD vertexForces = -vertexEnergyGradient;

        /*
         * Linear equation: (I - h * M_inverse * d_f_d_v - h^2 * M_inverse * d_f_d_x) * delta_v = h * M_inverse * (f(t0) + h * d_f_d_x * v(t0))
         * Format: mat * delta_v = vec
         * mat = I - h * M_inverse * d_f_d_v - h^2 * M_inverse * d_f_d_x
         * vec = h * M_inverse * (f(t0) + h * d_f_d_x * v(t0))
         */
        MatD d_vertexForces_d_velocity = new MatD(numVertices * 3, numVertices * 3); // This is a zero matrix until some velocity based damping force will be implemented.
        MatD d_vertexForces_d_pos = -forceHess; // This is the energy Hessian (switch sign from energy to force).
        MatD identMat = new MatD(inverseMassMatrix.numRows, inverseMassMatrix.numColumns).addDiag(1.0);
        MatD mat = identMat - deltaTime * inverseMassMatrix * d_vertexForces_d_velocity - deltaTime * deltaTime * inverseMassMatrix * d_vertexForces_d_pos;
        VecD vec = deltaTime * inverseMassMatrix * (vertexForces + deltaTime * d_vertexForces_d_pos * this.vertexVelocities);
        VecD deltaVelocity = this.sparseLinearSolve(mat, vec);
        //print("mat: " + mat);
        //print("vec: " + vec);
        //print("deltaVelocity: " + deltaVelocity);

        // Update vertex velocity.
        this.vertexVelocities += deltaVelocity;

        // Update vertex positions: pos(i + 1) = pos(i) + deltaTime * velocity(i + 1).
        for(int i = 0; i < this.vertexPositions.Length; i++) {
            for(int coord = 0; coord < 3; coord++) {
                this.vertexPositions[i][coord] += deltaTime * this.vertexVelocities[3 * i + coord];
            }
        }

        this.updateMesh();
    }

    /*
     * Updates the current mesh with the current vertices positions.
     * After setting the new vertices, the mesh normals are recalculated.
     */
    private void updateMesh() {
        Mesh mesh = this.getMesh();
        Vector3[] meshVerts = mesh.vertices;
        for(int i = 0; i < meshVerts.Length; i++) {
            for(int coord = 0; coord < 3; coord++) {
                meshVerts[i][coord] = (float) this.vertexPositions[i][coord];
            }
        }
        mesh.vertices = meshVerts;
        mesh.RecalculateNormals();
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
            Vector3 e12 = vertices[v2] - vertices[v1];
            Vector3 e13 = vertices[v3] - vertices[v1];
            Vector3 cross_e12_e13 = Vector3.Cross(e12, e13);
            float crossprod_length = cross_e12_e13.magnitude; // Length is the same, regardless of which edges are used.

            // Store the triangle normal and area if they exist (i.e. if the triangle has a normal and therefore an area).
            if(float.IsNaN(crossprod_length)) {
                print("Encountered zero-area triangle.");

                // The triangle area is 0 and the triangle has infinitely many normals in a circle (two edges parallel) or sphere (all vertices on the same point).
                this.triangleNormals_old[triangleId] = Vector3.zero;
                this.triangleNormals[triangleId] = new Vec3D(0, 0, 0);
                this.triangleAreas[triangleId] = 0f;

                // Technically, there are infinitely many options for the normal to change, and the area will almust always grow.
                // This simplification might introduce a small error in one timestep, but only in this exact zero-area case.
                this.dTriangleNormals_dv1[triangleId] = new Vector3[] {Vector3.zero, Vector3.zero, Vector3.zero};
                this.dTriangleNormals_dv2[triangleId] = new Vector3[] {Vector3.zero, Vector3.zero, Vector3.zero};
                this.dTriangleNormals_dv3[triangleId] = new Vector3[] {Vector3.zero, Vector3.zero, Vector3.zero};
                this.dTriangleAreas_dv1[triangleId] = Vector3.zero;
                this.dTriangleAreas_dv2[triangleId] = Vector3.zero;
                this.dTriangleAreas_dv3[triangleId] = Vector3.zero;
                continue;
            }
            this.triangleNormals_old[triangleId] = cross_e12_e13 / crossprod_length;
            this.triangleNormals[triangleId] = new Vec3D(this.triangleNormals_old[triangleId][0], this.triangleNormals_old[triangleId][1], this.triangleNormals_old[triangleId][2]);
            this.triangleAreas[triangleId] = crossprod_length / 2f; // Triangle area is half of the cross product of any two of its edges.

            /*
             * Triangle normal partial derivatives:
             * e1 === e12 and e2 === e13 and crossprod_length === cross_e1_e2_length.
             * n = cross_e1_e2 / cross_e1_e2_length = {e1y * e2z - e2y * e1z, e1x * e2z - e2x * e1z, e1x * e2y - e2x * e1y} / cross_e1_e2_length // Vector3.
             * d_n_d_e1 = {d_n_d_e1x, d_n_d_e1y, d_n_d_e1z} // 3x3 matrix.
             * d_n_d_e1x = (cross_e1_e2_length * d_cross_e1_e2_d_e1x - cross_e1_e2 * d_cross_e1_e2_length_d_e1x) / cross_e1_e2_length^2 // Vector3.
             * d_n_d_e1y = (cross_e1_e2_length * d_cross_e1_e2_d_e1y - cross_e1_e2 * d_cross_e1_e2_length_d_e1y) / cross_e1_e2_length^2 // Vector3.
             * d_n_d_e1z = (cross_e1_e2_length * d_cross_e1_e2_d_e1z - cross_e1_e2 * d_cross_e1_e2_length_d_e1z) / cross_e1_e2_length^2 // Vector3.
             * d_cross_e1_e2_d_e1x = {0, e2z, e2y} // Vector3.
             * d_cross_e1_e2_d_e1y = {e2z, 0, -e2x} // Vector3.
             * d_cross_e1_e2_d_e1z = {-e2y, -e2x, 0} // Vector3.
             */
            Vector3 d_cross_e12_e13_length_d_e12 = 1f / crossprod_length * new Vector3(
                    (e12.x * e13.z - e13.x * e12.z) * e13.z + (e12.x * e13.y - e13.x * e12.y) * e13.y,
                    (e12.y * e13.z - e13.y * e12.z) * e13.z + (e12.x * e13.y - e13.x * e12.y) * -e13.x,
                    (e12.y * e13.z - e13.y * e12.z) * -e13.y + (e12.x * e13.z - e13.x * e12.z) * -e13.x);
            Vector3 d_cross_e12_e13_length_d_e13 = 1f / crossprod_length * new Vector3(
                    (e12.y * e13.z - e13.y * e12.z) * 0      + (e12.x * e13.z - e13.x * e12.z) * -e12.z + (e12.x * e13.y - e13.x * e12.y) * -e12.y,
                    (e12.y * e13.z - e13.y * e12.z) * -e12.z + (e12.x * e13.z - e13.x * e12.z) * 0      + (e12.x * e13.y - e13.x * e12.y) * e12.x,
                    (e12.y * e13.z - e13.y * e12.z) * e12.y  + (e12.x * e13.z - e13.x * e12.z) * e12.x  + (e12.x * e13.y - e13.x * e12.y) * 0);
            Vector3 d_cross_e12_e13_d_e12x = new Vector3(0f, e13.z, e13.y);
            Vector3 d_cross_e12_e13_d_e12y = new Vector3(e13.z, 0f, -e13.x);
            Vector3 d_cross_e12_e13_d_e12z = new Vector3(-e13.y, -e13.x, 0f);
            Vector3 d_cross_e12_e13_d_e13x = new Vector3(0f, -e12.z, -e12.y);
            Vector3 d_cross_e12_e13_d_e13y = new Vector3(-e12.z, 0f, e12.x);
            Vector3 d_cross_e12_e13_d_e13z = new Vector3(e12.y, e12.x, 0f);
            Vector3 d_n_d_e12x = (crossprod_length * d_cross_e12_e13_d_e12x - cross_e12_e13 * d_cross_e12_e13_length_d_e12.x) / (crossprod_length * crossprod_length);
            Vector3 d_n_d_e12y = (crossprod_length * d_cross_e12_e13_d_e12y - cross_e12_e13 * d_cross_e12_e13_length_d_e12.y) / (crossprod_length * crossprod_length);
            Vector3 d_n_d_e12z = (crossprod_length * d_cross_e12_e13_d_e12z - cross_e12_e13 * d_cross_e12_e13_length_d_e12.z) / (crossprod_length * crossprod_length);
            Vector3 d_n_d_e13x = (crossprod_length * d_cross_e12_e13_d_e13x - cross_e12_e13 * d_cross_e12_e13_length_d_e13.x) / (crossprod_length * crossprod_length);
            Vector3 d_n_d_e13y = (crossprod_length * d_cross_e12_e13_d_e13y - cross_e12_e13 * d_cross_e12_e13_length_d_e13.y) / (crossprod_length * crossprod_length);
            Vector3 d_n_d_e13z = (crossprod_length * d_cross_e12_e13_d_e13z - cross_e12_e13 * d_cross_e12_e13_length_d_e13.z) / (crossprod_length * crossprod_length);
            Vector3[] d_n_d_e12 = new Vector3[] {d_n_d_e12x, d_n_d_e12y, d_n_d_e12z};
            Vector3[] d_n_d_e13 = new Vector3[] {d_n_d_e13x, d_n_d_e13y, d_n_d_e13z};
            Vector3 d_e12_d_v1 = new Vector3(-1f, -1f, -1f);
            Vector3 d_e12_d_v2 = new Vector3(1f, 1f, 1f);
            Vector3 d_e13_d_v3 = new Vector3(1f, 1f, 1f);
            this.dTriangleNormals_dv1[triangleId] = new Vector3[] {d_n_d_e12[0] * d_e12_d_v1.x, d_n_d_e12[1] * d_e12_d_v1.y, d_n_d_e12[2] * d_e12_d_v1.z};
            this.dTriangleNormals_dv2[triangleId] = new Vector3[] {d_n_d_e12[0] * d_e12_d_v2.x, d_n_d_e12[1] * d_e12_d_v2.y, d_n_d_e12[2] * d_e12_d_v2.z};
            this.dTriangleNormals_dv3[triangleId] = new Vector3[] {d_n_d_e13[0] * d_e13_d_v3.x, d_n_d_e13[1] * d_e13_d_v3.y, d_n_d_e13[2] * d_e13_d_v3.z};

            /*
             * Triangle area partial derivatives:
             * e1 === e12 and e2 === e13 and crossprod_length === cross_e1_e2_length.
             * triangleArea = crossprod_length / 2f
             * dTriangleArea_dv1 = {d_cross_e12_e13_length_d_e12x / 2f, d_cross_e12_e13_length_d_e12x / 2f, d_cross_e12_e13_length_d_e12z / 2f} .* d_e12_d_v1
             * dTriangleArea_dv2 = {d_cross_e12_e13_length_d_e12x / 2f, d_cross_e12_e13_length_d_e12y / 2f, d_cross_e12_e13_length_d_e12z / 2f} .* d_e12_d_v2
             * dTriangleArea_dv3 = {d_cross_e12_e13_length_d_e13x / 2f, d_cross_e12_e13_length_d_e13y / 2f, d_cross_e12_e13_length_d_e13z / 2f} .* d_e13_d_v3
             */
            this.dTriangleAreas_dv1[triangleId] = new Vector3(
                    d_cross_e12_e13_length_d_e12.x * d_e12_d_v1.x,
                    d_cross_e12_e13_length_d_e12.y * d_e12_d_v1.y,
                    d_cross_e12_e13_length_d_e12.z * d_e12_d_v1.z) / 2f;
            this.dTriangleAreas_dv2[triangleId] = new Vector3(
                    d_cross_e12_e13_length_d_e12.x * d_e12_d_v2.x,
                    d_cross_e12_e13_length_d_e12.y * d_e12_d_v2.y,
                    d_cross_e12_e13_length_d_e12.z * d_e12_d_v2.z) / 2f;
            this.dTriangleAreas_dv3[triangleId] = new Vector3(
                    d_cross_e12_e13_length_d_e13.x * d_e13_d_v3.x,
                    d_cross_e12_e13_length_d_e13.y * d_e13_d_v3.y,
                    d_cross_e12_e13_length_d_e13.z * d_e13_d_v3.z) / 2f;
        }
    }

    private Vector3[] calcVertexEnergyGradient(int[] triangles, Vector3[] vertices) {

        // Initialize vertex energy gradient array.
        Vector3[] vertexEnergyGradient = new Vector3[vertices.Length];
        for(int i = 0; i < vertexEnergyGradient.Length; i++) {
            vertexEnergyGradient[i] = Vector3.zero;
        }

        // Compute vertex energy gradient array.
        for(int vertexInd = 0; vertexInd < this.sortedVertexTriangles.Length; vertexInd++) {
            for(int i = 0; i < this.sortedVertexTriangles[vertexInd].Count; i++) {

                // Get two possibly adjacent triangles in clockwise direction.
                int triangleId = this.sortedVertexTriangles[vertexInd][i];
                int nextTriangleId = this.sortedVertexTriangles[vertexInd][(i + 1) % this.sortedVertexTriangles[vertexInd].Count];

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
                    if(kBend != 0f) {
                        vertexEnergyGradient[vertexInd] += kBend * this.getBendingEnergyGradient(
                            vertices, triangleId, nextTriangleId, v11, v12, v13, v21, v22, v23, vertexInd, otherVertexClockwiseInd1);
                    }
                } else {

                    // Calculate spring energy gradient in the second edge.
                    if(kLength != 0f) {
                        vertexEnergyGradient[vertexInd] += kLength * this.getEdgeLengthEnergyGradient(vertices, vertexInd, otherVertexAntiClockwiseInd2);
                    }
                }

                // Calculate the area energy gradient in the triangle.
                if(kArea != 0f) {
                    if(vertexInd == v11) {
                        vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v11, v12, v13);
                    } else if(vertexInd == v12) {
                        vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v12, v13, v11);
                    } else {
                        vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v13, v11, v12);
                    }
                }
            }
        }
        return vertexEnergyGradient;
    }

    private VecD calcVertexEnergyGradient(int[] triangles, Vec3D[] vertices) {

        // Initialize vertex energy gradient array.
        VecD vertexEnergyGradient = new VecD(vertices.Length * 3);

        // Compute edge length and bending energy gradient.
        foreach(Edge edge in this.edges) {

            // Compute edge length energy gradient.
            if(this.kLength != 0f) {
                VecD edgeLengthEnergyGrad = this.kLength * this.getEdgeLengthEnergyGradient(vertices, edge.ve1, edge.ve2);
                vertexEnergyGradient[3 * edge.ve1] += edgeLengthEnergyGrad[0];
                vertexEnergyGradient[3 * edge.ve1 + 1] += edgeLengthEnergyGrad[1];
                vertexEnergyGradient[3 * edge.ve1 + 2] += edgeLengthEnergyGrad[2];
                vertexEnergyGradient[3 * edge.ve2] += edgeLengthEnergyGrad[3];
                vertexEnergyGradient[3 * edge.ve2 + 1] += edgeLengthEnergyGrad[4];
                vertexEnergyGradient[3 * edge.ve2 + 2] += edgeLengthEnergyGrad[5];
            }

            // Compute edge bending energy gradient.
            if(edge.hasSideFlaps() && this.kBend != 0f) {
                VecD edgeBendEnergyGrad = this.kBend * this.getBendingEnergyGradient(vertices, edge);
                vertexEnergyGradient[3 * edge.ve1] += edgeBendEnergyGrad[0];
                vertexEnergyGradient[3 * edge.ve1 + 1] += edgeBendEnergyGrad[1];
                vertexEnergyGradient[3 * edge.ve1 + 2] += edgeBendEnergyGrad[2];
                vertexEnergyGradient[3 * edge.ve2] += edgeBendEnergyGrad[3];
                vertexEnergyGradient[3 * edge.ve2 + 1] += edgeBendEnergyGrad[4];
                vertexEnergyGradient[3 * edge.ve2 + 2] += edgeBendEnergyGrad[5];
                vertexEnergyGradient[3 * edge.vf1] += edgeBendEnergyGrad[6];
                vertexEnergyGradient[3 * edge.vf1 + 1] += edgeBendEnergyGrad[7];
                vertexEnergyGradient[3 * edge.vf1 + 2] += edgeBendEnergyGrad[8];
                vertexEnergyGradient[3 * edge.vf2] += edgeBendEnergyGrad[9];
                vertexEnergyGradient[3 * edge.vf2 + 1] += edgeBendEnergyGrad[10];
                vertexEnergyGradient[3 * edge.vf2 + 2] += edgeBendEnergyGrad[11];
            }
        }

        // Compute triangle area energy gradient.
        if(this.kArea != 0f) {
            for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
                int v1 = triangles[3 * triangleId];
                int v2 = triangles[3 * triangleId + 1];
                int v3 = triangles[3 * triangleId + 2];
                VecD triangleAreaEnergyGrad = this.kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v1, v2, v3);
                for(int coord = 0; coord < 3; coord++) {
                    vertexEnergyGradient[3 * v1 + coord] += triangleAreaEnergyGrad[coord];
                    vertexEnergyGradient[3 * v2 + coord] += triangleAreaEnergyGrad[3 + coord];
                    vertexEnergyGradient[3 * v3 + coord] += triangleAreaEnergyGrad[6 + coord];
                }
            }
        }

        // TODO - Add area and bending energy gradients. Strip out the logic that's no longer needed (since the edge list has been added).
        // Compute vertex energy gradient array.
        //for(int vertexInd = 0; vertexInd < this.vertexTriangles.Length; vertexInd++) {
        //    for(int i = 0; i < this.vertexTriangles[vertexInd].Count; i++) {

        //        // Get two possibly adjacent triangles in clockwise direction.
        //        int triangleId = this.vertexTriangles[vertexInd][i];
        //        int nextTriangleId = this.vertexTriangles[vertexInd][(i + 1) % this.vertexTriangles[vertexInd].Count];

        //        // Get triangle vertices.
        //        int triangleBaseInd1 = triangleId * 3;
        //        int triangleBaseInd2 = nextTriangleId * 3;
        //        int v11 = triangles[triangleBaseInd1];
        //        int v12 = triangles[triangleBaseInd1 + 1];
        //        int v13 = triangles[triangleBaseInd1 + 2];
        //        int v21 = triangles[triangleBaseInd2];
        //        int v22 = triangles[triangleBaseInd2 + 1];
        //        int v23 = triangles[triangleBaseInd2 + 2];

        //        // Get the vertex indices of the other vertices that are connected to the possibly shared triangle edge.
        //        int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
        //        int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));

        //        // Handle the edge, or both edges if they are not the same.
        //        bool edgeSharedByTriangles = (otherVertexClockwiseInd1 == otherVertexAntiClockwiseInd2);

        //        if(edgeSharedByTriangles) {

        //            // Calculate bending energy gradient.
        //            if(kBend != 0f) {
        //                vertexEnergyGradient[vertexInd] += kBend * this.getBendingEnergyGradient(
        //                    vertices, triangleId, nextTriangleId, v11, v12, v13, v21, v22, v23, vertexInd, otherVertexClockwiseInd1);
        //            }
        //        }

        //        // Calculate the area energy gradient in the triangle.
        //        if(kArea != 0f) {
        //            if(vertexInd == v11) {
        //                vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v11, v12, v13);
        //            } else if(vertexInd == v12) {
        //                vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v12, v13, v11);
        //            } else {
        //                vertexEnergyGradient[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v13, v11, v12);
        //            }
        //        }
        //    }
        //}
        return vertexEnergyGradient;
    }

    private Vector3[][] calcVertexEnergyHessian(int[] triangles, Vector3[] vertices) {

        // Initialize vertex energy Hessian.
        Vector3[][] vertexEnergyHess = new Vector3[vertices.Length][];
        for(int i = 0; i < vertexEnergyHess.Length; i++) {
            vertexEnergyHess[i] = new Vector3[] {Vector3.zero, Vector3.zero, Vector3.zero};
        }

        // Compute vertex energy Hessian.
        for(int vertexInd = 0; vertexInd < this.sortedVertexTriangles.Length; vertexInd++) {
            for(int i = 0; i < this.sortedVertexTriangles[vertexInd].Count; i++) {

                // Get two possibly adjacent triangles in clockwise direction.
                int triangleId = this.sortedVertexTriangles[vertexInd][i];
                int nextTriangleId = this.sortedVertexTriangles[vertexInd][(i + 1) % this.sortedVertexTriangles[vertexInd].Count];

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

                // Calculate spring energy Hessian in the edge.
                Vector3[] lengthHess = this.getEdgeLengthEnergyHess(vertices, vertexInd, otherVertexClockwiseInd1);
                vertexEnergyHess[vertexInd][0] += kLength * lengthHess[0];
                vertexEnergyHess[vertexInd][1] += kLength * lengthHess[1];
                vertexEnergyHess[vertexInd][2] += kLength * lengthHess[2];

                // TODO - Implement area and bending Hessian.
                if(edgeSharedByTriangles) {

                    // Calculate bending energy Hessian.
                    //vertexEnergyHess[vertexInd] += kBend * this.getBendingEnergyGradient(
                    //    vertices, triangleId, nextTriangleId, v11, v12, v13, v21, v22, v23, vertexInd, otherVertexClockwiseInd1);
                } else {

                    // Calculate spring energy Hessian in the second edge.
                    //vertexEnergyHess[vertexInd] += kLength * this.getEdgeLengthEnergyGradient(vertices, vertexInd, otherVertexAntiClockwiseInd2);
                }

                // Calculate the area energy Hessian in the triangle.
                if(vertexInd == v11) {
                    //vertexEnergyHess[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v11, v12, v13);
                } else if(vertexInd == v12) {
                    //vertexEnergyHess[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v12, v13, v11);
                } else {
                    //vertexEnergyHess[vertexInd] += kArea * this.getTriangleAreaEnergyGradient(vertices, triangleId, v13, v11, v12);
                }
            }
        }
        return vertexEnergyHess;
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
            Vector3 triangleNormal = this.triangleNormals_old[triangleId];
            float triangleArea = this.triangleAreas[triangleId];
            Vector3 totalTriangleWindForce = Vector3.Dot(this.windPressure, triangleNormal) * triangleArea * triangleNormal;

            // Add a third of the total triangle wind force to each of its vertices.
            vertexWindForce[v1] += totalTriangleWindForce / 3f;
            vertexWindForce[v2] += totalTriangleWindForce / 3f;
            vertexWindForce[v3] += totalTriangleWindForce / 3f;
        }
        return vertexWindForce;
    }

    /*
     * Calculates the wind force acting on each vertex.
     * Returns the wind force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
     */
    private VecD calcVertexWindForce(int[] triangles, VecD[] vertices) {

        // Initialize vertex wind force array.
        VecD vertexWindForce = new VecD(vertices.Length * 3);

        // Compute vertex wind force.
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];

            // Compute projection of the wind pressure vector on the triangle normal.
            VecD triangleNormal = this.triangleNormals[triangleId];
            if(triangleNormal == null) {
                continue; // Triangle has a zero-area and no normal, so the projected wind force is zero as well.
            }
            float triangleArea = this.triangleAreas[triangleId];
            VecD totalTriangleWindForce = VecD.dot(new VecD(this.windPressure), triangleNormal) * triangleArea * triangleNormal;

            // Add a third of the total triangle wind force to each of its vertices.
            VecD totalTriangleWindForcePart = totalTriangleWindForce / 3d;
            for(int coord = 0; coord < 3; coord++) {
                vertexWindForce[3 * v1 + coord] += totalTriangleWindForcePart[coord];
                vertexWindForce[3 * v2 + coord] += totalTriangleWindForcePart[coord];
                vertexWindForce[3 * v3 + coord] += totalTriangleWindForcePart[coord];
            }
        }
        return vertexWindForce;
    }

    private Vector3[][] calcVertexWindForceHessian(int[] triangles, Vector3[] vertices) {

        // Initialize vertex wind force Hessian.
        Vector3[][] vertexWindForceHess = new Vector3[vertices.Length][];
        for(int i = 0; i < vertexWindForceHess.Length; i++) {
            vertexWindForceHess[i] = new Vector3[] { Vector3.zero, Vector3.zero, Vector3.zero };
        }

        // Compute vertex wind force Hessian.
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];

            // Compute projection of the wind pressure vector on the triangle normal.
            Vector3 triangleNormal = this.triangleNormals_old[triangleId];
            Vector3[] dTriangleNormal_dv1 = this.dTriangleNormals_dv1[triangleId];
            Vector3[] dTriangleNormal_dv2 = this.dTriangleNormals_dv2[triangleId];
            Vector3[] dTriangleNormal_dv3 = this.dTriangleNormals_dv3[triangleId];
            float triangleArea = this.triangleAreas[triangleId];
            Vector3 dTriangleArea_dv1 = this.dTriangleAreas_dv1[triangleId];
            Vector3 dTriangleArea_dv2 = this.dTriangleAreas_dv2[triangleId];
            Vector3 dTriangleArea_dv3 = this.dTriangleAreas_dv3[triangleId];
            Vector3 totalTriangleWindForce = Vector3.Dot(this.windPressure, triangleNormal) * triangleArea * triangleNormal;

            /*
             * Vector3 totalTriangleWindForce = Vector3.Dot(this.windPressure, triangleNormal) * triangleArea * triangleNormal
             * = (windPressure.x * triangleNormal.x + windPressure.y * triangleNormal.y + windPressure.z * triangleNormal.z) * triangleArea * triangleNormal
             * 
             * a = windPressure.x * triangleNormal.x + windPressure.y * triangleNormal.y + windPressure.z * triangleNormal.z // float.
             * da_dv1 = {da_dv1x, da_dv1y, da_dv1z} // Vector3.
             * da_dv1x = windPressure.x * dTriangleNormal_x_dv1x + windPressure.y * dTriangleNormal_y_dv1x + windPressure.z * dTriangleNormal_z_dv1x // float.
             * da_dv1y = windPressure.x * dTriangleNormal_x_dv1y + windPressure.y * dTriangleNormal_y_dv1y + windPressure.z * dTriangleNormal_z_dv1y // float.
             * da_dv1z = windPressure.x * dTriangleNormal_x_dv1z + windPressure.y * dTriangleNormal_y_dv1z + windPressure.z * dTriangleNormal_z_dv1z // float.
             * 
             * dTotalTriangleWindForce_dv1 = da_dv1 * triangleArea * triangleNormal + a * dTriangleArea_dv1 * triangleNormal + a * triangleArea * d_triangleNormal_d_v1
             * dTotalTriangleWindForce_dv1x = da_dv1x * triangleArea * triangleNormal + a * dTriangleArea_dv1x * triangleNormal + a * triangleArea * d_triangleNormal_dv1x
             * dTotalTriangleWindForce_dv1y = da_dv1y * triangleArea * triangleNormal + a * dTriangleArea_dv1y * triangleNormal + a * triangleArea * d_triangleNormal_dv1y
             * dTotalTriangleWindForce_dv1z = da_dv1z * triangleArea * triangleNormal + a * dTriangleArea_dv1z * triangleNormal + a * triangleArea * d_triangleNormal_dv1z
             * 
             */
            float a = Vector3.Dot(this.windPressure, triangleNormal);
            // TODO - Validate that dTriangleNormal_dv1[0] is dTriangleNormal_x_dv1, or otherwise use first column.
            Vector3 da_dv1 = windPressure.x * dTriangleNormal_dv1[0] + windPressure.y * dTriangleNormal_dv1[1] + windPressure.z * dTriangleNormal_dv1[2];
            Vector3 da_dv2 = windPressure.x * dTriangleNormal_dv2[0] + windPressure.y * dTriangleNormal_dv2[1] + windPressure.z * dTriangleNormal_dv2[2];
            Vector3 da_dv3 = windPressure.x * dTriangleNormal_dv3[0] + windPressure.y * dTriangleNormal_dv3[1] + windPressure.z * dTriangleNormal_dv3[2];

            Vector3[] dTriangleVertexWindForce_dv1 = new Vector3[] {
                (da_dv1.x * triangleArea * triangleNormal + a * dTriangleArea_dv1.x * triangleNormal + a * triangleArea * dTriangleNormal_dv1[0]) / 3f,
                (da_dv1.y * triangleArea * triangleNormal + a * dTriangleArea_dv1.y * triangleNormal + a * triangleArea * dTriangleNormal_dv1[1]) / 3f,
                (da_dv1.z * triangleArea * triangleNormal + a * dTriangleArea_dv1.z * triangleNormal + a * triangleArea * dTriangleNormal_dv1[2]) / 3f
            };
            Vector3[] dTriangleVertexWindForce_dv2 = new Vector3[] {
                (da_dv2.x * triangleArea * triangleNormal + a * dTriangleArea_dv2.x * triangleNormal + a * triangleArea * dTriangleNormal_dv2[0]) / 3f,
                (da_dv2.y * triangleArea * triangleNormal + a * dTriangleArea_dv2.y * triangleNormal + a * triangleArea * dTriangleNormal_dv2[1]) / 3f,
                (da_dv2.z * triangleArea * triangleNormal + a * dTriangleArea_dv2.z * triangleNormal + a * triangleArea * dTriangleNormal_dv2[2]) / 3f
            };
            Vector3[] dTriangleVertexWindForce_dv3 = new Vector3[] {
                (da_dv3.x * triangleArea * triangleNormal + a * dTriangleArea_dv3.x * triangleNormal + a * triangleArea * dTriangleNormal_dv3[0]) / 3f,
                (da_dv3.y * triangleArea * triangleNormal + a * dTriangleArea_dv3.y * triangleNormal + a * triangleArea * dTriangleNormal_dv3[1]) / 3f,
                (da_dv3.z * triangleArea * triangleNormal + a * dTriangleArea_dv3.z * triangleNormal + a * triangleArea * dTriangleNormal_dv3[2]) / 3f
            };

            // Add a third of the total triangle wind force to each of its vertices.
            for(int i = 0; i < 3; i++) {
                vertexWindForceHess[v1][i] += dTriangleVertexWindForce_dv1[i];
                vertexWindForceHess[v2][i] += dTriangleVertexWindForce_dv2[i];
                vertexWindForceHess[v3][i] += dTriangleVertexWindForce_dv3[i];
            }
        }
        return vertexWindForceHess;
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
        Vector3 result = (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength;
        if(float.IsNaN(result.x) || float.IsNaN(result.y) || float.IsNaN(result.z)) {
            print("NaN length gradient: " + result + " undeformedEdgeLength: " + undeformedEdgeLength
                    + " edgeLength: " + edgeLength);
            return Vector3.zero;
        }
        if(float.IsInfinity(result.x) || float.IsInfinity(result.y) || float.IsInfinity(result.z)) {
            print("Infinite length gradient: " + result + " undeformedEdgeLength: " + undeformedEdgeLength
                    + " edgeLength: " + edgeLength);
            return Vector3.zero;
        }
        return result;
    }

    /**
     * Computes the edge length energy gradient of the edge between vertices v1 and v2, towards {v1x, v1y, v1z, v2x, v2y, v2z}.
     * The result is a vector of length 6.
     */
    private VecD getEdgeLengthEnergyGradient(Vec3D[] vertexPositions, int v1, int v2) {
        Vec3D edge = vertexPositions[v2] - vertexPositions[v1]; // Vector from v1 to v2.
        double edgeLength = edge.magnitude;
        if(double.IsNaN(edgeLength)) {
            return new VecD(0, 0, 0, 0, 0, 0); // Edge is zero-length, so the gradient is 0.
        }
        double undeformedEdgeLength = (this.originalVertices[v2] - this.originalVertices[v1]).magnitude;
        Vec3D dEdgeLength_dv1 = (vertexPositions[v1] - vertexPositions[v2]) / edgeLength;
        Vec3D dEdgeLength_dv2 = -dEdgeLength_dv1;
        VecD dEdgeLength_dv1v2 = new VecD(dEdgeLength_dv1, dEdgeLength_dv2); // Partial derivative towards {v1x, v1y, v1z, v2x, v2y, v2z}.

        VecD result = (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_dv1v2;
        for(int i = 0; i < result.length; i++) {
            if(double.IsNaN(result[i])) {
                print("NaN length gradient: " + result + " undeformedEdgeLength: " + undeformedEdgeLength
                        + " edgeLength: " + edgeLength);
                return new VecD(0, 0, 0, 0, 0, 0);
            }
            if(double.IsInfinity(result[i])) {
                print("Infinite length gradient: " + result + " undeformedEdgeLength: " + undeformedEdgeLength
                        + " edgeLength: " + edgeLength);
                return new VecD(0, 0, 0, 0, 0, 0);
            }
        }
        return result;
    }

    /**
     * Computes the edge length energy Hessian of the edge between vertices v1 and v2, towards {v1x, v1y, v1z, v2x, v2y, v2z}.
     * The result is a 6x6 matrix containing all combinations of double partial derivatives towards {v1x, v1y, v1z, v2x, v2y, v2z}.
     */
    private MatD getEdgeLengthEnergyHess(VecD[] vertices, int v1, int v2) {
        /*
         * edgeLength (float):
         *     sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
         * 
         * dEdgeLength_dv1x (float):
         *     (v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2) = (v1x - v2x) / edgeLength
         * 
         * dEdgeLength_dv1 (Vector3):
         *     (v1 - v2) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2) = (v1 - v2) / edgeLength
         * 
         * dEdgeLength_dv2 (Vector3):
         *     -dEdgeLength_dv1
         * 
         * ddEdgeLength_dv1x_dv1y (float):
         *    ((v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
         *    (edgeLength * (v1x - v2x)' - (v1x - v2x) * edgeLength') / edgeLength^2
         *    (edgeLength * 0 - (v1x - v2x) * dEdgeLength.y) / edgeLength^2
         *    (v2x - v1x) * dEdgeLength.y / edgeLength^2
         *
         * ddEdgeLength_dv1y_dv1x (float):
         *    ((v1y - v2y) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
         *    (edgeLength * (v1y - v2y)' - (v1y - v2y) * edgeLength') / edgeLength^2
         *    (edgeLength * 0 - (v1y - v2y) * dEdgeLength.x) / edgeLength^2
         *    (v2y - v1y) * dEdgeLength.x / edgeLength^2 = (v2y - v1y) * (v1x - v2x) / edgeLength / edgeLength^2 = (v2y - v1y) * (v1x - v2x) / edgeLength^3
         *    ddEdgeLength_dv1x_dv1y // The resulting Hessian is symmetrical.
         * 
         * ddEdgeLength_dv1x_dv1x (float):
         *     ((v1x - v2x) / edgeLength)'
         *     (edgeLength * (v1x - v2x)' - (v1x - v2x) * dEdgeLength.x) / edgeLength^2
         *     (edgeLength - (v1x - v2x) * dEdgeLength.x) / ((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
         *         = (edgeLength - (v1x - v2x) * dEdgeLength.x) / edgeLength^2
         *         = (edgeLength - (v1x - v2x) * (v1x - v2x) / edgeLength) / edgeLength^2
         *         = (edgeLength^2 - (v1x - v2x)^2) / edgeLength^3
         * 
         * ddEdgeLength_dv1y_dv1y (float):
         *     ((v1y - v2y) / edgeLength)'
         *     (edgeLength * (v1y - v2y)' - (v1y - v2y) * dEdgeLength.y) / edgeLength^2
         *     (edgeLength - (v1y - v2y) * dEdgeLength.y) / ((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
         *     (edgeLength - (v1y - v2y) * dEdgeLength.y) / edgeLength^2
         *     (edgeLength - (v1y - v2y) * (v1y - v2y) / edgeLength) / edgeLength^2
         *     (edgeLength^2 - (v1y - v2y)^2) / edgeLength^3
         * 
         * ddEdgeLength_dv2x_dv2y (float):
         *    (-(v1x - v2x) / edgeLength)'
         *    (edgeLength * -(v1x - v2x)' - -(v1x - v2x) * edgeLength') / edgeLength^2
         *    (edgeLength * 0 + (v1x - v2x) * dEdgeLength_dv2y) / edgeLength^2
         *    (v1x - v2x) * dEdgeLength_dv2y / edgeLength^2
         *    (v1x - v2x) * (-(v1y - v2y) / edgeLength) / edgeLength^2
         *    (v1x - v2x) * (v2y - v1y) / edgeLength^3
         *    -ddEdgeLength_dv1x_dv1y
         * 
         * ddEdgeLength_dv1x_dv2x (float):
         *    ((v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
         *    (edgeLength * (v1x - v2x)' - (v1x - v2x) * edgeLength') / edgeLength^2
         *    (edgeLength * -1 - (v1x - v2x) * dEdgeLength_dv2x) / edgeLength^2
         *    (-edgeLength - (v1x - v2x) * dEdgeLength_dv2x) / edgeLength^2
         *    (-edgeLength - (v1x - v2x) * -(v1x - v2x) / edgeLength) / edgeLength^2
         *    (-edgeLength^2 + (v1x - v2x)^2) / edgeLength^3
         * 
         * ddEdgeLength_dv1x_dv2y (float):
         *    ((v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
         *    (edgeLength * (v1x - v2x)' - (v1x - v2x) * edgeLength') / edgeLength^2
         *    (edgeLength * 0 - (v1x - v2x) * dEdgeLength_dv2y) / edgeLength^2
         *    (v2x - v1x) * dEdgeLength_dv2y / edgeLength^2
         *    (v2x - v1x) * ((v1y - v2y) / edgeLength) / edgeLength^2
         *    (v2x - v1x) * (v1y - v2y) / edgeLength^3
         *    ddEdgeLength_dv1x_dv1y // x and y are just independent?
         * 
         * ddEdgeLength_dv2x_dv1y (float):
         *    (-(v1x - v2x) / edgeLength)'
         *    (edgeLength * -(v1x - v2x)' - -(v1x - v2x) * edgeLength') / edgeLength^2
         *    (edgeLength * 0 + (v1x - v2x) * dEdgeLength.y) / edgeLength^2
         *    (v1x - v2x) * dEdgeLength_dv1y / edgeLength^2
         *    ddEdgeLength_dv1y_dv2x // The resulting Hessian is symmetrical.
         * 
         * ddEdgeLength_dv2x_dv2x (float):
         *     (-(v1x - v2x) / edgeLength)'
         *     (edgeLength * -(v1x - v2x)' - -(v1x - v2x) * dEdgeLength_dv2x) / edgeLength^2
         *     (edgeLength + (v1x - v2x) * dEdgeLength_dv2x) / ((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
         *         = (edgeLength + (v1x - v2x) * dEdgeLength_dv2x) / edgeLength^2
         *         = (edgeLength + (v1x - v2x) * -(v1x - v2x) / edgeLength) / edgeLength^2
         *         = (edgeLength^2 - (v1x - v2x)^2) / edgeLength^3
         *         = ddEdgeLength_dv1x_dv1x
         * 
         * ddEdgeLength_dv2i_dv2i (Vector3):
         *     ddEdgeLength_dv1i_dv1i // The 3 diagonals are equal, together these form the diagonals of the 6x6 length Hessian.
         * 
         * ddEdgeLength_dv2i_dv2j where i != j (3x3 matrix excluding diagonals):
         *     -ddEdgeLength_dv1i_dv1j
         * 
         * ddEdgeLength_dv1_dv1 (symmetric 3x3 matrix):
         *     | ddEdgeLength_dv1x_dv1x, ddEdgeLength_dv1y_dv1x, ddEdgeLength_dv1z_dv1x |
         *     | ddEdgeLength_dv1x_dv1y, ddEdgeLength_dv1y_dv1y, ddEdgeLength_dv1z_dv1y |
         *     | ddEdgeLength_dv1x_dv1z, ddEdgeLength_dv1y_dv1z, ddEdgeLength_dv1z_dv1z |
         *     =
         *     | edgeLength^2 - (v1x - v2x)^2, (v2y - v1y) * (v1x - v2x)   , (v2z - v1z) * (v1x - v2x) |
         *     | (v2x - v1x) * (v1y - v2y)   , edgeLength^2 - (v2y - v1y)^2, (v2z - v1z) * (v1x - v2x) | / edgeLength^3
         *     | (v2x - v1x) * (v1z - v2z)   , (v2y - v1y) * (v1z - v2z), edgeLength^2 - (v2z - v1z)^2 |
         * 
         * ddEdgeLength_dv2_dv2 (symmetric 3x3 matrix):
         *     | ddEdgeLength_dv2x_dv2x, ddEdgeLength_dv2y_dv2x, ddEdgeLength_dv2z_dv2x |
         *     | ddEdgeLength_dv2x_dv2y, ddEdgeLength_dv2y_dv2y, ddEdgeLength_dv2z_dv2y |
         *     | ddEdgeLength_dv2x_dv2z, ddEdgeLength_dv2y_dv2z, ddEdgeLength_dv2z_dv2z |
         *     =
         *     |  ddEdgeLength_dv1x_dv1x, -ddEdgeLength_dv1y_dv1x, -ddEdgeLength_dv1z_dv1x |
         *     | -ddEdgeLength_dv1x_dv1y,  ddEdgeLength_dv1y_dv1y, -ddEdgeLength_dv1z_dv1y |
         *     | -ddEdgeLength_dv1x_dv1z, -ddEdgeLength_dv1y_dv1z,  ddEdgeLength_dv1z_dv1z |
         * 
         * ddEdgeLength_dv1_dv2 (symmetric 3x3 matrix):
         *     | ddEdgeLength_dv1x_dv2x, ddEdgeLength_dv1y_dv2x, ddEdgeLength_dv1z_dv2x |
         *     | ddEdgeLength_dv1x_dv2y, ddEdgeLength_dv1y_dv2y, ddEdgeLength_dv1z_dv2y |
         *     | ddEdgeLength_dv1x_dv2z, ddEdgeLength_dv1y_dv2z, ddEdgeLength_dv1z_dv2z |
         *     =
         *     | -ddEdgeLength_dv1x_dv1x,  ddEdgeLength_dv1y_dv1x,  ddEdgeLength_dv1z_dv1x |
         *     |  ddEdgeLength_dv1x_dv1y, -ddEdgeLength_dv1y_dv1y,  ddEdgeLength_dv1z_dv1y |
         *     |  ddEdgeLength_dv1x_dv1z,  ddEdgeLength_dv1y_dv1z, -ddEdgeLength_dv1z_dv1z |
         * 
         * ddEdgeLength_dv2_dv1 (symmetric 3x3 matrix):
         *     ddEdgeLength_dv1_dv2 // Symmetrical.
         * 
         * edgeLengthHess (ddEdgeLength_dv1v2_dv1v2) (6x6 matrix, composed of 4 3x3 matrices):
         *     | ddEdgeLength_dv1_dv1, ddEdgeLength_dv1_dv2 |
         *     | ddEdgeLength_dv2_dv1, ddEdgeLength_dv2_dv2 |
         * 
         * Length energy gradient (dLengthEnergy_dv1v2) (Vector with length 6):
         *     Substitute: dEdgeLength_dv1v2 = [dEdgeLength_dv1, dEdgeLength_dv2]
         *     (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_dv1v2
         *     dLengthEnergy_di = (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_di // Double, for i = {v1x, v1y, v1z, v2x, v2y, v2z}.
         * 
         * Length energy Hessian (ddLengthEnergy_dv1v2_dv1v2) (6x6 matrix):
         *     
         *     Consider per element i in length energy gradient, towards each element j (elements are {v1x, v1y, v1z, v2x, v2y, v2z}):
         *         ((2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_di)'
         *         (2 * edgeLength / undeformedEdgeLength - 2)' * dEdgeLength_di + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_di_dj
         *         2 / undeformedEdgeLength * dEdgeLength_dj * dEdgeLength_di + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_di_dj
         *             // For all 6 coordinates, this becomes symmetrical due to ddEdgeLength_di_dj == ddEdgeLength_dj_di and dEdgeLength_dj * dEdgeLength_di == dEdgeLength_di * dEdgeLength_dj
         *     
         *     ddLengthEnergy_di_dj = 2 / undeformedEdgeLength * dEdgeLength_dj * dEdgeLength_di + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_di_dj
         *         // For i and j = {v1x, v1y, v1z, v2x, v2y, v2z}.
         */

        // TODO - This is a copy from the energy gradient code. Combine this in a way to prevent double calculations.
        VecD edge = vertexPositions[v2] - vertexPositions[v1]; // Vector from v1 to v2.
        double edgeLength = edge.magnitude;
        if(double.IsNaN(edgeLength)) {
            return new MatD(new double[,] {
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0}
            });
        }
        double undeformedEdgeLength = (this.originalVertices[v2] - this.originalVertices[v1]).magnitude;
        VecD dEdgeLength_dv1 = (vertexPositions[v1] - vertexPositions[v2]) / edgeLength;
        VecD dEdgeLength_dv2 = -dEdgeLength_dv1;
        VecD dEdgeLength_dv1v2 = new VecD(dEdgeLength_dv1, dEdgeLength_dv2); // Partial derivative towards {v1x, v1y, v1z, v2x, v2y, v2z}.

        VecD dEdgeEnergy_dv1v2 = (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_dv1v2;
        // TODO - Copied gradient code ends here (See TODO above).

        // Calculate edge length Hessian.
        double edgeLengthSquare = edgeLength * edgeLength;
        double edgeLengthCube = edgeLengthSquare * edgeLength;
        MatD ddEdgeLength_dv1_dv1 = new MatD(new double[,] {
            {edgeLengthSquare - edge[0] * edge[0],                  - edge[1] * edge[0],                  - edge[2] * edge[0]},
            {                 - edge[0] * edge[1], edgeLengthSquare - edge[1] * edge[1],                  - edge[2] * edge[1]},
            {                 - edge[0] * edge[2],                  - edge[1] * edge[2], edgeLengthSquare - edge[2] * edge[2]}
        }) / edgeLengthCube;
        
        MatD ddEdgeLength_dv1_dv2 = ddEdgeLength_dv1_dv1.Clone();
        ddEdgeLength_dv1_dv2[0, 0] *= -1;
        ddEdgeLength_dv1_dv2[1, 1] *= -1;
        ddEdgeLength_dv1_dv2[2, 2] *= -1;
        MatD ddEdgeLength_dv2_dv1 = ddEdgeLength_dv1_dv2;
        MatD ddEdgeLength_dv2_dv2 = -ddEdgeLength_dv1_dv2;

        MatD edgeLengthHess = new MatD(6, 6);
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                edgeLengthHess[i, j] = ddEdgeLength_dv1_dv1[i, j];
                edgeLengthHess[i + 3, j] = ddEdgeLength_dv1_dv2[i, j];
                edgeLengthHess[i, j + 3] = ddEdgeLength_dv2_dv1[i, j];
                edgeLengthHess[i + 3, j + 3] = ddEdgeLength_dv2_dv2[i, j];
            }
        }

        // Calculate length energy Hessian.
        // ddLengthEnergy_di_dj = 2 / undeformedEdgeLength * dEdgeLength_dj * dEdgeLength_di + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_di_dj
        MatD lengthEnergyHess = new MatD(6, 6);
        double a = (2 * edgeLength / undeformedEdgeLength - 2);
        for(int i = 0; i < 6; i++) {
            for(int j = 0; j < 6; j++) {
                lengthEnergyHess[i, j] = 2 / undeformedEdgeLength * dEdgeLength_dv1v2[i] * dEdgeLength_dv1v2[j] + a * edgeLengthHess[i, j];
            }
        }
        return lengthEnergyHess;
    }

    /**
     * Computes the edge length energy Hessian of vertex v1, for the edge between vertices v1 and v2.
     */
    private Vector3[] getEdgeLengthEnergyHess(Vector3[] vertices, int v1, int v2) {
        /*
         * EdgeLength (float):
         *     sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
         * 
         * dEdgeLength_dv1x (float):
         *     (v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2) = (v1x - v2x) / edgeLength
         * 
         * dEdgeLength_dv1 (Vector3):
         *     (v1 - v2) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2) = (v1 - v2) / edgeLength
         * 
         * ddEdgeLength_dv1x_dv1y (float):
         *    ((v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
         *    (edgeLength * (v1x - v2x)' - (v1x - v2x) * edgeLength') / edgeLength^2
         *    (edgeLength * 0 - (v1x - v2x) * dEdgeLength.y) / edgeLength^2
         *    (v2x - v1x) * dEdgeLength.y / edgeLength^2
         *
         * ddEdgeLength_dv1y_dv1x (float):
         *    ((v1y - v2y) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
         *    (edgeLength * (v1y - v2y)' - (v1y - v2y) * edgeLength') / edgeLength^2
         *    (edgeLength * 0 - (v1y - v2y) * dEdgeLength.x) / edgeLength^2
         *    (v2y - v1y) * dEdgeLength.x / edgeLength^2 = (v2y - v1y) * (v1x - v2x) / edgeLength / edgeLength^2 = (v2y - v1y) * (v1x - v2x) / edgeLength^3
         * 
         * ddEdgeLength_dv1x_dv1x (float):
         *     ((v1x - v2x) / edgeLength)'
         *     (edgeLength * (v1x - v2x)' - (v1x - v2x) * dEdgeLength.x) / edgeLength^2
         *     (edgeLength - (v1x - v2x) * dEdgeLength.x) / ((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
         *         = (edgeLength - (v1x - v2x) * dEdgeLength.x) / edgeLength^2
         *         = (edgeLength - (v1x - v2x) * (v1x - v2x) / edgeLength) / edgeLength^2
         *         = (edgeLength^2 - (v1x - v2x)^2) / edgeLength^3
         * 
         * ddEdgeLength_dv1y_dv1y (float):
         *     ((v1y - v2y) / edgeLength)'
         *     (edgeLength * (v1y - v2y)' - (v1y - v2y) * dEdgeLength.y) / edgeLength^2
         *     (edgeLength - (v1y - v2y) * dEdgeLength.y) / ((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
         *     (edgeLength - (v1y - v2y) * dEdgeLength.y) / edgeLength^2
         *     (edgeLength - (v1y - v2y) * (v1 - v2) / edgeLength) / edgeLength^2
         *     (edgeLength^2 - (v1y - v2y)^2) / edgeLength^3
         * 
         * ddEdgeLength_dv1_dv1 (symmetric 3x3 matrix, Hessian):
         *     | ddEdgeLength_dv1x_dv1x, ddEdgeLength_dv1y_dv1x, ddEdgeLength_dv1z_dv1x |
         *     | ddEdgeLength_dv1x_dv1y, ddEdgeLength_dv1y_dv1y, ddEdgeLength_dv1z_dv1y |
         *     | ddEdgeLength_dv1x_dv1z, ddEdgeLength_dv1y_dv1z, ddEdgeLength_dv1z_dv1z |
         *     =
         *     | edgeLength^2 - (v1x - v2x)^2, (v2y - v1y) * (v1x - v2x)   , (v2z - v1z) * (v1x - v2x) |
         *     | (v2x - v1x) * (v1y - v2y)   , edgeLength^2 - (v2y - v1y)^2, (v2z - v1z) * (v1x - v2x) | / edgeLength^3
         *     | (v2x - v1x) * (v1z - v2z)   , (v2y - v1y) * (v1z - v2z), edgeLength^2 - (v2z - v1z)^2 |
         * 
         * Length energy gradient (Vector3):
         *     (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_dv1
         * 
         * Length energy Hessian (3x3 matrix):
         *     ((2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_dv1)'
         *     (2 * edgeLength / undeformedEdgeLength - 2)' * dEdgeLength_dv1 + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1
         *     2 / undeformedEdgeLength * dEdgeLength_dv1 + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1
         *     =
         *     | 2 / undeformedEdgeLength * dEdgeLength_dv1x + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1x |
         *     | 2 / undeformedEdgeLength * dEdgeLength_dv1y + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1y |
         *     | 2 / undeformedEdgeLength * dEdgeLength_dv1z + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1z |
         */

        Vector3 e_v1_v2 = vertices[v2] - vertices[v1];
        Vector3 undeformedEdge = this.originalVertices[v2] - this.originalVertices[v1];
        float edgeLength = e_v1_v2.magnitude;
        if(float.IsNaN(edgeLength)) {
            return new Vector3[] {Vector3.zero, Vector3.zero, Vector3.zero};
        }
        float undeformedEdgeLength = undeformedEdge.magnitude;
        Vector3 dEdgeLength_dv1 = (vertices[v1] - vertices[v2]) / edgeLength;

        float edgeLengthSquare = edgeLength * edgeLength;
        float edgeLengthCube = edgeLengthSquare * edgeLength;
        Vector3[] ddEdgeLength_dv1_dv1 = new Vector3[] {
            new Vector3(edgeLengthSquare - e_v1_v2.x * e_v1_v2.x,                  - e_v1_v2.y * e_v1_v2.x,                  - e_v1_v2.z * e_v1_v2.x) / edgeLengthCube,
            new Vector3(                 - e_v1_v2.x * e_v1_v2.y, edgeLengthSquare - e_v1_v2.y * e_v1_v2.y,                  - e_v1_v2.z * e_v1_v2.y) / edgeLengthCube,
            new Vector3(                 - e_v1_v2.x * e_v1_v2.z,                  - e_v1_v2.y * e_v1_v2.z, edgeLengthSquare - e_v1_v2.z * e_v1_v2.z) / edgeLengthCube
        };

        // Calculate length energy Hessian. This matrix is symmetrical, so transposing it has no effect.
        Vector3 a = 2 / undeformedEdgeLength * dEdgeLength_dv1;
        float b = 2 * edgeLength / undeformedEdgeLength - 2;
        Vector3[] lengthEnergyHess = new Vector3[] {
            a + b * ddEdgeLength_dv1_dv1[0],
            a + b * ddEdgeLength_dv1_dv1[1],
            a + b * ddEdgeLength_dv1_dv1[2]
        };
        return lengthEnergyHess;
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
        //Vector3 dTriangleArea = dCrossProdLength / 6f; // Area of triangle is half the cross product length, and we only look at a third.
        Vector3 dTriangleArea = dCrossProdLength / 2f; // Area of triangle is half the cross product length.

        // Calculate the area energy gradient.
        Vector3 result = (2 * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2) * dTriangleArea;
        if(float.IsNaN(result.x) || float.IsNaN(result.y) || float.IsNaN(result.z)) {
            print("NaN area gradient: " + result + " triangleArea: " + this.triangleAreas[triangleId]
                    + " dTriangleArea: " + dTriangleArea + " crossProdLength: " + crossProdLength);
            return Vector3.zero;
        }
        if(float.IsInfinity(result.x) || float.IsInfinity(result.y) || float.IsInfinity(result.z)) {
            print("Infinite area gradient: " + result + " triangleArea: " + this.triangleAreas[triangleId]
                    + " dTriangleArea: " + dTriangleArea + " crossProdLength: " + crossProdLength);
            return Vector3.zero;
        }
        return result;
    }

    /**
     * Computes the triangle area energy gradient for the triangle defined by vertices v1, v2 and v3.
     * Returns the gradient towards {v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z}.
     */
    private VecD getTriangleAreaEnergyGradient(Vec3D[] vertices, int triangleId, int v1Ind, int v2Ind, int v3Ind) {

        // Return if the area is zero.
        if(this.triangleAreas[triangleId] == 0d) {
            return new VecD(0, 0, 0, 0, 0, 0, 0, 0, 0); // Area is 0 m^2, so the gradient is 0.
        }
        
        // Get triangle vertices.
        Vec3D v1 = vertices[v1Ind];
        Vec3D v2 = vertices[v2Ind];
        Vec3D v3 = vertices[v3Ind];

        // Compute triangle energy gradient. See thesis notes for derivation.
        double crossProdLength = this.triangleAreas[triangleId] * 2d;
        VecD a = new VecD(0          , v3.z - v2.z, v2.y - v3.y, 0                        , v2.z - v3.z + v1.z - v2.z, v2.y - v1.y + v3.y - v2.y, 0          , v2.z - v1.z, v1.y - v2.y);
        VecD b = new VecD(v2.z - v3.z, 0          , v3.x - v2.x, v2.z - v1.z + v3.z - v2.z, 0                        , v2.x - v3.x + v1.x - v2.x, v1.z - v2.z, 0          , v2.x - v1.x);
        VecD c = new VecD(v3.y - v2.y, v2.x - v3.x, 0          , v2.y - v3.y + v1.y - v2.y, v2.x - v1.x + v3.x - v2.x, 0                        , v2.y - v1.y, v1.x - v2.x, 0          );
        double d = ((v1.y - v2.y) * (v3.z - v2.z) - (v3.y - v2.y) * (v1.z - v2.z));
        double e = ((v2.x - v1.x) * (v3.z - v2.z) + (v3.x - v2.x) * (v1.z - v2.z));
        double f = ((v1.x - v2.x) * (v3.y - v2.y) - (v3.x - v2.x) * (v1.y - v2.y));
        VecD d_crossProdLength_dv1v2v3 = (d * a + e * b + f * c) / crossProdLength;
        VecD d_triangleArea_dv1v2v3 = d_crossProdLength_dv1v2v3 / 2d;
        VecD d_triangleEnergy_dv1v2v3 = (2d * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2d) * d_triangleArea_dv1v2v3;
        return d_triangleEnergy_dv1v2v3;
    }

    /**
     * Computes the triangle area energy Hessian of vertex v1, for the triangle defined by vertices v1, v2 and v3.
     * Note that 1/3 of the area energy is used for v1 (the other two parts will be used for v2 and v3).
     */
    private Vector3[] getTriangleAreaEnergyHessian(Vector3[] vertices, int triangleId, int v1, int v2, int v3) {

        /*
         * Triangle energy gradient in v1 (Vector3):
         *     (2 * triangleArea / undeformedTriangleArea - 2) * dTriangleArea_dv1
         *     = 2 * triangleArea / undeformedTriangleArea * dTriangleArea_dv1 - 2 * dTriangleArea_dv1
         *     = 2 / undeformedTriangleArea * triangleArea * dTriangleArea_dv1 - 2 * dTriangleArea_dv1
         * 
         * Assumption: The triangle Hessian in v1 only depends on the gradient in v1, and not on the gradients in v2 and v3.
         * 
         * Triangle energy Hessian (3x3 matrix):
         *     | ddTriangleAreaEnergy_dv1x_dv1x, ddTriangleAreaEnergy_dv1x_dv1y, ddTriangleAreaEnergy_dv1x_dv1z |
         *     | ddTriangleAreaEnergy_dv1y_dv1x, ddTriangleAreaEnergy_dv1y_dv1y, ddTriangleAreaEnergy_dv1y_dv1z |
         *     | ddTriangleAreaEnergy_dv1z_dv1x, ddTriangleAreaEnergy_dv1z_dv1y, ddTriangleAreaEnergy_dv1z_dv1z |
         * 
         * ddTriangleAreaEnergy_dv1x_dv1x = d((2 * triangleArea / undeformedTriangleArea - 2) * dTriangleArea_dv1x) / dv1x
         *     = 2 * dTriangleArea_dv1x / undeformedTriangleArea * dTriangleArea_dv1x + (2 * triangleArea / undeformedTriangleArea - 2) * ddTriangleArea_dv1x_dv1x
         * ddTriangleAreaEnergy_dv1x_dv1y = d((2 * triangleArea / undeformedTriangleArea - 2) * dTriangleArea_dv1x) / dv1y
         *     = 2 * dTriangleArea_dv1y / undeformedTriangleArea * dTriangleArea_dv1x + (2 * triangleArea / undeformedTriangleArea - 2) * ddTriangleArea_dv1x_dv1y
         * 
         * Cache ddTriangleArea_dv1_dv1 if this is used in bending energy as well? Doesn't look like it. Or if in the caching code more variables are available already.
         * 
         * 
         */
        // TODO - Implement.
        /*
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
        Vector3 result = (2 * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2) * dTriangleArea;
        if(float.IsNaN(result.x) || float.IsNaN(result.y) || float.IsNaN(result.z)) {
            print("NaN area gradient: " + result + " triangleArea: " + this.triangleAreas[triangleId]
                    + " dTriangleArea: " + dTriangleArea + " crossProdLength: " + crossProdLength);
            return Vector3.zero;
        }
        if(float.IsInfinity(result.x) || float.IsInfinity(result.y) || float.IsInfinity(result.z)) {
            print("Infinite area gradient: " + result + " triangleArea: " + this.triangleAreas[triangleId]
                    + " dTriangleArea: " + dTriangleArea + " crossProdLength: " + crossProdLength);
            return Vector3.zero;
        }
        return result;
        */
        return null;
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
        Vector3 n1 = this.triangleNormals_old[triangleId1];
        Vector3 n2 = this.triangleNormals_old[triangleId2];

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
        Vector3 result = d_W_bending_energy_edge_d_teta_e * d_teta_d_ve1;
        if(float.IsNaN(result.x) || float.IsNaN(result.y) || float.IsNaN(result.z)) {
            print("NaN bending gradient: " + result + " d_W_bending_energy_edge_d_teta_e: " + d_W_bending_energy_edge_d_teta_e
                    + " d_teta_d_ve1: " + d_teta_d_ve1);
            return Vector3.zero;
        }
        if(float.IsInfinity(result.x) || float.IsInfinity(result.y) || float.IsInfinity(result.z)) {
            print("Infinite bending gradient: " + result + " d_W_bending_energy_edge_d_teta_e: " + d_W_bending_energy_edge_d_teta_e
                    + " d_teta_d_ve1: " + d_teta_d_ve1);
            return Vector3.zero;
        }
        return result;
    }

    /**
     * Computes the bending energy gradient of the given edge. This gradient is taken towards the edge-defining vertices ve1 and ve2, as well
     * as the vertices vf1 and vf2, completing triangles 1 and 2 connected to the edge respectively.
     * The return values are the bending energy gradient towards: {ve1x. ve1y, ve1z, ve2x. ve2y, ve2z, vf1x. vf1y, vf1z, vf2x. vf2y, vf2z}.
     */
    private VecD getBendingEnergyGradient(Vec3D[] vertices, Edge edge) {

        // Define required constants.
        double undeformedEdgeLength = (this.originalVertices[edge.ve2] - this.originalVertices[edge.ve1]).magnitude;
        // h_e_undeformed is a third of the average triangle height, where the height is twice the triangle area divided by the triangle width.
        double h_e_undeformed = (this.undeformedTriangleAreas[edge.triangleId1] + this.undeformedTriangleAreas[edge.triangleId2]) / undeformedEdgeLength / 3d;
        double teta_e_undeformed = 0; // TODO - Add option to use a non-flat rest state depending on the useFlatUndeformedBendState field.

        // Get triangle normals.
        Vec3D n1 = this.triangleNormals[edge.triangleId1];
        Vec3D n2 = this.triangleNormals[edge.triangleId2];

        // Return if at least one of the triangles has a zero area and no normal.
        if(n1 == null || n2 == null) {
            return new VecD(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }
        
        // Calculate hinge angle.
        double teta_e_sign = Mathf.Sign((float) VecD.dot(n1, vertices[edge.vf2] - vertices[edge.ve1])); // 1 if teta_e positive, -1 if negative.
        double dot_n1_n2 = VecD.dot(n1, n2);
        double teta_e = -Mathf.Acos((float) (dot_n1_n2 > 1d ? 1d : (dot_n1_n2 < -1d ? -1d : dot_n1_n2))) * teta_e_sign; // Limit the dot product of the normals at 1 (fix precision errors).

        // TODO - Remove this check if it can't happen. Check what happens for non-existing normals though.
        if(double.IsNaN(teta_e)) {
            print("NaN teta_e. sign: " + teta_e_sign + ", dot(n1, n2): " + VecD.dot(n1, n2));
            // teta_e = Mathf.PI; // The triangles are on top of each other, which is both 180 and -180 degrees.
            return new VecD(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }
        
        // Calculate energy derivative towards hinge angle teta.
        double d_fi_d_teta = 2d * (teta_e - teta_e_undeformed) * undeformedEdgeLength / h_e_undeformed;

        // Return if the energy gradient is zero.
        if(d_fi_d_teta == 0d) {
            return new VecD(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }

        // Define edges, angles and edge heights (following paper https://studios.disneyresearch.com/wp-content/uploads/2019/03/Discrete-Bending-Forces-and-Their-Jacobians-Paper.pdf).
        Vec3D e0 = vertices[edge.ve2] - vertices[edge.ve1]; // Middle horizontal edge.
        Vec3D t1_e1 = vertices[edge.vf1] - vertices[edge.ve2]; // Top right edge.
        Vec3D t1_e2 = vertices[edge.vf1] - vertices[edge.ve1]; // Top left edge.
        Vec3D t2_e1 = vertices[edge.vf2] - vertices[edge.ve2]; // Bottom right edge.
        Vec3D t2_e2 = vertices[edge.vf2] - vertices[edge.ve1]; // Bottom left edge.
        float t1_alpha1 = Mathf.Acos((float) (VecD.dot(e0, t1_e2) / (e0.magnitude * t1_e2.magnitude)));
        float t1_alpha2 = Mathf.Acos((float) (VecD.dot(-e0, t1_e1) / (e0.magnitude * t1_e1.magnitude)));
        float t2_alpha1 = Mathf.Acos((float) (VecD.dot(e0, t2_e2) / (e0.magnitude * t2_e2.magnitude)));
        float t2_alpha2 = Mathf.Acos((float) (VecD.dot(-e0, t2_e1) / (e0.magnitude * t2_e1.magnitude)));
        double t1_h0 = 2d * this.triangleAreas[edge.triangleId1] / e0.magnitude;
        double t1_h1 = 2d * this.triangleAreas[edge.triangleId1] / t1_e1.magnitude;
        double t1_h2 = 2d * this.triangleAreas[edge.triangleId1] / t1_e2.magnitude;
        double t2_h0 = 2d * this.triangleAreas[edge.triangleId2] / e0.magnitude;
        double t2_h1 = 2d * this.triangleAreas[edge.triangleId2] / t2_e1.magnitude;
        double t2_h2 = 2d * this.triangleAreas[edge.triangleId2] / t2_e2.magnitude;

        // Calculate derivatives of teta towards all 4 vertices.
        VecD d_teta_d_ve1 = Mathf.Cos(t1_alpha2) / t1_h1 * n1 + Mathf.Cos(t2_alpha2) / t2_h1 * n2;
        VecD d_teta_d_ve2 = Mathf.Cos(t1_alpha1) / t1_h2 * n1 + Mathf.Cos(t2_alpha1) / t2_h2 * n2;
        VecD d_teta_d_vf1 = -n1 / t1_h0;
        VecD d_teta_d_vf2 = -n2 / t2_h0;
        
        // Assemble hinge angle gradient.
        VecD d_teta_d_ve1ve2vf1vf2 = new VecD(d_teta_d_ve1, d_teta_d_ve2, d_teta_d_vf1, d_teta_d_vf2);

        // Return the bending energy gradient.
        return d_fi_d_teta * d_teta_d_ve1ve2vf1vf2;
    }

    /**
     * Makes the given Hessian positive definite by adding the identity matrix to it until it is positive definite.
     */
    private static void makeHessPositiveDefinite(Vector3[] hess) {
        if(hess[0].x <= 0) {
            float amount = 0.01f - hess[0].x; // 0.01f to ensure a positive non-zero value.
            hess[0].x += amount;
            hess[1].y += amount;
            hess[2].z += amount;
        }
        float subRowOneFromTwoAmount = hess[1].x / hess[0].x;
        float subRowOneFromThreeAmount = hess[2].x / hess[0].x;
        float hess11 = hess[1].y - subRowOneFromTwoAmount * hess[0].y;
        if(hess11 <= 0) {
            float amount = 0.01f - hess11; // 0.01f to ensure a positive non-zero value.
            hess[0].x += amount;
            hess[1].y += amount;
            hess[2].z += amount;
        }
        float subRowTwoFromThreeAmount = hess[2].y / hess[1].y;
        float hess22 = hess[2].z - subRowOneFromThreeAmount * hess[0].z - subRowTwoFromThreeAmount * hess[1].z;
        if(hess22 <= 0) {
            float amount = 0.01f - hess22; // 0.01f to ensure a positive non-zero value.
            hess[0].x += amount;
            hess[1].y += amount;
            hess[2].z += amount;
        }
    }

    /*
     * Solves the given linear equation using a sparse linear solver.
     * Equation: mat * x = vec
     * Returns vector x.
     */
    private VecD sparseLinearSolve(MatD mat, VecD vec) {
        // See example code: https://www.alglib.net/translator/man/manual.csharp.html#example_linlsqr_d_1
        alglib.sparsematrix algMat;
        alglib.sparsecreate(mat.numRows, mat.numColumns, out algMat);
        for(int row = 0; row < mat.numRows; row++) {
            for(int col = 0; col < mat.numColumns; col++) {
                if(mat[row, col] != 0) {
                    alglib.sparseset(algMat, row, col, mat[row, col]);
                }
            }
        }
        alglib.sparseconverttocrs(algMat);
        alglib.linlsqrstate solverObj;
        alglib.linlsqrreport report;
        double[] x;
        alglib.linlsqrcreate(mat.numRows, mat.numColumns, out solverObj);
        alglib.linlsqrsolvesparse(solverObj, algMat, vec.asDoubleArray());
        alglib.linlsqrresults(solverObj, out x, out report);
        return new VecD(x);
    }
}
