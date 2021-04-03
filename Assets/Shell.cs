﻿using System;
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
    private double undeformedEdgeLengthFactor = 1d; // Factor to multiple all undeformed edge lengths with to make artificial edge-constrained sails non-flat.
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
    public double gravityConstant = 9.81d;

    // Cached vertex/triangle properties.
    private Vec3D[] triangleNormals;
    private Vec3D[][] dTriangleNormals_dv1;
    private Vec3D[][] dTriangleNormals_dv2;
    private Vec3D[][] dTriangleNormals_dv3;
    private Vec3D[] undeformedTriangleNormals;
    private double[] triangleAreas;
    private Vec3D[] dTriangleAreas_dv1;
    private Vec3D[] dTriangleAreas_dv2;
    private Vec3D[] dTriangleAreas_dv3;
    private double[] undeformedTriangleAreas;

    void Awake() {
        QualitySettings.vSyncCount = 0; // Disable V-sync.
        Application.targetFrameRate = 100; // Set max framerate.
        System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US"); // To print decimal points instead of commas.
    }

    // Start is called before the first frame update.
    void Start() {

        // Initialize fields. Doing this overwrites the values set in Unity's inspector.
        this.kLength = 10f;
        this.kArea = 0f;
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
            this.undeformedEdgeLengthFactor = 1.1d;
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
        this.triangleNormals = new Vec3D[numTriangles];
        this.dTriangleNormals_dv1 = new Vec3D[numTriangles][];
        this.dTriangleNormals_dv2 = new Vec3D[numTriangles][];
        this.dTriangleNormals_dv3 = new Vec3D[numTriangles][];
        this.triangleAreas = new double[numTriangles];
        this.dTriangleAreas_dv1 = new Vec3D[numTriangles];
        this.dTriangleAreas_dv2 = new Vec3D[numTriangles];
        this.dTriangleAreas_dv3 = new Vec3D[numTriangles];
        this.undeformedTriangleNormals = new Vec3D[numTriangles];
        this.undeformedTriangleAreas = new double[numTriangles];
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];
            Vec3D crossProd = Vec3D.cross(
                    (new Vec3D(this.originalVertices[v2]) - new Vec3D(this.originalVertices[v1])) * this.undeformedEdgeLengthFactor,
                    (new Vec3D(this.originalVertices[v3]) - new Vec3D(this.originalVertices[v1])) * this.undeformedEdgeLengthFactor);
            double crossProdMag = crossProd.magnitude;

            // Store the triangle normal and area if they exist (i.e. if the triangle has a normal and therefore an area).
            if(double.IsNaN(crossProdMag)) {
                this.undeformedTriangleNormals[triangleId] = Vec3D.zero;
                this.undeformedTriangleAreas[triangleId] = 0;
            } else {
                this.undeformedTriangleNormals[triangleId] = crossProd / crossProdMag;
                this.undeformedTriangleAreas[triangleId] = crossProdMag / 2d; // Triangle area is half of the cross product of any two of its edges.
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

        // Compute triangle normals and areas.
        this.recalcTriangleNormalsAndAreas(triangles, this.vertexPositions);

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
            case TimeSteppingMethod.OPTIMIZATION_INTEGRATOR: {
                this.doOptimizationIntegratorStep(deltaTime);
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
        VecD vertexEnergyGradient = this.getSystemEnergyGradient(triangles, this.vertexPositions);
        VecD vertexWindForce = this.getVertexWindForce(triangles, this.vertexPositions);
        VecD vertexCoordMasses = this.getVertexCoordinateMasses();
        VecD gravityForce = this.getVertexGravityForce(vertexCoordMasses);
        VecD step = this.kGradientDescent * (vertexWindForce + gravityForce - vertexEnergyGradient);
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
        VecD vertexEnergyGradient = this.getSystemEnergyGradient(triangles, this.vertexPositions);
        VecD vertexWindForce = this.getVertexWindForce(triangles, this.vertexPositions);

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
        VecD vertexEnergyGradient = this.getSystemEnergyGradient(triangles, this.vertexPositions);
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
     * Perform a simulation step based on this Optimization Integrator paper: https://www.math.ucla.edu/~jteran/papers/GSSJT15.pdf
     */
    private void doOptimizationIntegratorStep(double deltaTime) {
        // TODO - Implement.
        /*
         * Forces:
         *     f = -d_fi_d_x where x are the vertex coordinates.
         *  
         * Forces equilibrium:
         *     h(xNew) = massMatrix * (xNew - xOld - deltaTime * vOld) / deltaTime^2 + d_fi_d_x = 0
         *     Goal: Obtain xNew, being the new vertex coordinates.
         * 
         * Energy E:
         *     E(xNew) = 1 / (2 * deltaTime^2) * transpose(xNew - xOld - deltaTime * vOld) * massMatrix * (xNew - xOld - deltaTime * vOld) + fi // double
         *         This value will never be needed, only its gradient h(x) and its Hessian h'(x).
         * 
         * Energy E gradient:
         *     Approximation: Use constant mass within one step.
         *     d_E(x)_dx = h(x) = massMatrix * (x - xOld - deltaTime * vOld) / deltaTime^2 + d_fi_dx // VecD
         * 
         * Energy E Hessian:
         *     Approximation: Use constant mass within one step.
         *     dd_E(x)_dx_dx = d_h(x)_dx = d_(massMatrix * (x - xOld - deltaTime * vOld) / deltaTime^2 + d_fi_d_x)_dx
         *         = massMatrix / deltaTime^2 + dd_fi_dx_dx // MatD
         * 
         * 
         */
        
        // Declare constants.
        double terminationThreshold = 0.5d; // TODO - Set sensible value.
        double kappa = 0.01d; // Value as proposed by Optimization Integrator paper.
        double maxStepMagnitude = this.vertexPositions.Length * 0.001; // TODO - Set sensible value. Optimization Integrator paper uses 1000 (mesh size dependent?).
        int numVertices = this.vertexPositions.Length;
        int[] triangles = this.getMesh().triangles;

        // Perform Newton's method, setting steps until the termination criterion has been met.
        Vec3D[] newVertexPositions = new Vec3D[this.vertexPositions.Length]; // TODO - Use pos + deltaTime * velocity as initial guess.
        for(int i = 0; i < this.vertexPositions.Length; i++) {
            newVertexPositions[i] = this.vertexPositions[i].clone();
        }
        VecD vertexCoordMasses = this.getVertexCoordinateMasses(); // Masses per vertex coordinate. Format: {m_v1x, m_v1y, m_v1z, ...}.
        int maxNumIterations = 1000;
        int iteration = 0;
        while(true) {

            // Limit amount of iterations to prevent endless loops.
            if(++iteration > maxNumIterations) {
                print("Maximum number of iterations reached in Optimization Integrator update: " + maxNumIterations + ". Returning without taking a step.");
                return;
            }

            // Get system-wide energy gradient.
            VecD energyGradient = this.getSystemEnergyGradient(triangles, newVertexPositions);

            // Get wind force.
            VecD windForce = this.getVertexWindForce(triangles, newVertexPositions);

            // Get gravity force.
            VecD gravityForce = this.getVertexGravityForce(vertexCoordMasses);

            // Set energy gradient and wind force to zero for constrained vertices.
            // This causes the E gradient to be zero for them as well, causing it not to get in the way of the minimization problem.
            for(int i = 0; i < this.verticesMovementConstraints.Length; i++) {
                if(this.verticesMovementConstraints[i]) {
                    for(int coord = 0; coord < 3; coord++) {
                        energyGradient[3 * i + coord] = 0;
                        windForce[3 * i + coord] = 0;
                        gravityForce[3 * i + coord] = 0;
                    }
                }
            }

            // Get E gradient.
            VecD eGradient = VecD.multiplyElementWise(vertexCoordMasses,
                    new VecD(newVertexPositions) - new VecD(this.vertexPositions) - deltaTime * this.vertexVelocities)
                    / (deltaTime * deltaTime) + energyGradient - windForce - gravityForce;

            // Terminate when the termination criterion has been met.
            double eGradientMagnitude = eGradient.magnitude;
            //print("E gradient magnitude: " + eGradientMagnitude + " (threshold: " + terminationThreshold + ")");
            if(eGradientMagnitude < terminationThreshold) {
                break;
            }

            // Assemble system-wide energy Hessian.
            // TODO - Represent this using triples instead of a full matrix.
            MatD energyHess = new MatD(numVertices * 3, numVertices * 3);
            foreach(Edge edge in this.edges) {

                // The length energy Hessian consists of 4 (3x3) parts that have to be inserted into the matrix.
                MatD lengthEnergyHess = this.getEdgeLengthEnergyHess(newVertexPositions, edge.ve1, edge.ve2);
                makeHessPositiveDefinite(lengthEnergyHess);
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {

                        // ddLengthEnergy_dv1_dv1.
                        energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] = lengthEnergyHess[i, j];

                        // ddLengthEnergy_dv2_dv2.
                        energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] = lengthEnergyHess[i + 3, j + 3];

                        // ddLengthEnergy_dv1_dv2.
                        energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] = lengthEnergyHess[i, j + 3];

                        // ddLengthEnergy_dv2_dv1.
                        energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] = lengthEnergyHess[i + 3, j];
                    }
                }

                // The bending energy Hessian consists of 16 (3x3) parts that have to be inserted into the matrix.
                if(edge.hasSideFlaps()) {
                    MatD bendEnergyHess = this.getEdgeBendEnergyHess(newVertexPositions, edge);
                    makeHessPositiveDefinite(bendEnergyHess);
                    for(int i = 0; i < 3; i++) {
                        for(int j = 0; j < 3; j++) {
                            energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] = bendEnergyHess[i, j];
                            energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] = bendEnergyHess[i, j + 3];
                            energyHess[edge.ve1 * 3 + i, edge.vf1 * 3 + j] = bendEnergyHess[i, j + 6];
                            energyHess[edge.ve1 * 3 + i, edge.vf2 * 3 + j] = bendEnergyHess[i, j + 9];
                        
                            energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] = bendEnergyHess[i + 3, j];
                            energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] = bendEnergyHess[i + 3, j + 3];
                            energyHess[edge.ve2 * 3 + i, edge.vf1 * 3 + j] = bendEnergyHess[i + 3, j + 6];
                            energyHess[edge.ve2 * 3 + i, edge.vf2 * 3 + j] = bendEnergyHess[i + 3, j + 9];
                        
                            energyHess[edge.vf1 * 3 + i, edge.ve1 * 3 + j] = bendEnergyHess[i + 6, j];
                            energyHess[edge.vf1 * 3 + i, edge.ve2 * 3 + j] = bendEnergyHess[i + 6, j + 3];
                            energyHess[edge.vf1 * 3 + i, edge.vf1 * 3 + j] = bendEnergyHess[i + 6, j + 6];
                            energyHess[edge.vf1 * 3 + i, edge.vf2 * 3 + j] = bendEnergyHess[i + 6, j + 9];
                        
                            energyHess[edge.vf2 * 3 + i, edge.ve1 * 3 + j] = bendEnergyHess[i + 9, j];
                            energyHess[edge.vf2 * 3 + i, edge.ve2 * 3 + j] = bendEnergyHess[i + 9, j + 3];
                            energyHess[edge.vf2 * 3 + i, edge.vf1 * 3 + j] = bendEnergyHess[i + 9, j + 6];
                            energyHess[edge.vf2 * 3 + i, edge.vf2 * 3 + j] = bendEnergyHess[i + 9, j + 9];
                        }
                    }
                }
            }
            for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId += 3) {
                int v1 = triangles[triangleId];
                int v2 = triangles[triangleId + 1];
                int v3 = triangles[triangleId + 2];
                MatD areaEnergyHess = this.getTriangleAreaEnergyHessian(newVertexPositions, triangleId, v1, v2, v3);
                makeHessPositiveDefinite(areaEnergyHess);
                for(int i = 0; i < 3; i++) {
                    for(int j = 0; j < 3; j++) {

                        // ddAreaEnergy_dvi_dvi.
                        energyHess[v1 * 3 + i, v1 * 3 + j] = areaEnergyHess[i, j];
                        energyHess[v2 * 3 + i, v2 * 3 + j] = areaEnergyHess[i + 3, j + 3];
                        energyHess[v3 * 3 + i, v3 * 3 + j] = areaEnergyHess[i + 6, j + 6];

                        // ddAreaEnergy_dv1_dv2 & ddAreaEnergy_dv2_dv1.
                        energyHess[v1 * 3 + i, v2 * 3 + j] = areaEnergyHess[i, j + 3];
                        energyHess[v2 * 3 + i, v1 * 3 + j] = areaEnergyHess[i + 3, j];

                        // ddAreaEnergy_dv1_dv3 & ddAreaEnergy_dv3_dv1.
                        energyHess[v1 * 3 + i, v3 * 3 + j] = areaEnergyHess[i, j + 6];
                        energyHess[v3 * 3 + i, v1 * 3 + j] = areaEnergyHess[i + 6, j];

                        // ddAreaEnergy_dv2_dv3 & ddAreaEnergy_dv3_dv2.
                        energyHess[v2 * 3 + i, v3 * 3 + j] = areaEnergyHess[i + 3, j + 6];
                        energyHess[v3 * 3 + i, v2 * 3 + j] = areaEnergyHess[i + 6, j + 3];
                    }
                }
            }

            // Get E Hessian.
            MatD eHess = energyHess;
            double deltaTimeSquare = deltaTime * deltaTime;
            for(int i = 0; i < eHess.numRows; i++) {
                eHess[i, i] += vertexCoordMasses[i] / deltaTimeSquare;
            }

            // Compute Newton step.
            // TODO - Rewrite to step = -inverse(eHess) * eGradient, being equivalent to eHess * step = eGradient (use sparse direct solver).
            VecD step = -eGradient;

            // Ensure that the step is in downhill direction.
            // If a < b, then the step is suitable. Otherwise try -a < b. As a last resort, fall back to gradient descent.
            double a = VecD.dot(step, eGradient);
            double b = -kappa * step.magnitude * eGradient.magnitude;
            if(a >= b) {
                if(-a < b) {
                    step = -step;
                } else {
                    step = -eGradient;
                }
            }

            // Clamp max step magnitude.
            double stepMagnitude = step.magnitude;
            if(stepMagnitude > maxStepMagnitude) {
                step = step / stepMagnitude * maxStepMagnitude;
            }

            // TODO - Choose step size alpha in direction delta_x using a line search.
            double alpha = 0.01d;

            // Take step: x += alpha * step.
            for(int i = 0; i < newVertexPositions.Length; i++) {
                newVertexPositions[i].x += alpha * step[3 * i];
                newVertexPositions[i].y += alpha * step[3 * i + 1];
                newVertexPositions[i].z += alpha * step[3 * i + 2];
            }
        }
        
        // Update vertex velocity.
        this.vertexVelocities = (new VecD(newVertexPositions) - new VecD(this.vertexPositions)) / deltaTime;

        // Update vertex positions.
        this.vertexPositions = newVertexPositions;

        // Update mesh.
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

    private void recalcTriangleNormalsAndAreas(int[] triangles, Vec3D[] vertices) {
        for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
            int triangleBaseIndex = triangleId * 3;
            int v1 = triangles[triangleBaseIndex];
            int v2 = triangles[triangleBaseIndex + 1];
            int v3 = triangles[triangleBaseIndex + 2];
            Vec3D e12 = vertices[v2] - vertices[v1];
            Vec3D e13 = vertices[v3] - vertices[v1];
            Vec3D cross_e12_e13 = Vec3D.cross(e12, e13);
            double crossprod_length = cross_e12_e13.magnitude; // Length is the same, regardless of which edges are used.

            // Store the triangle normal and area if they exist (i.e. if the triangle has a normal and therefore an area).
            if(double.IsNaN(crossprod_length)) {
                print("Encountered zero-area triangle.");

                // The triangle area is 0 and the triangle has infinitely many normals in a circle (two edges parallel) or sphere (all vertices on the same point).
                this.triangleNormals[triangleId] = Vec3D.zero;
                this.triangleAreas[triangleId] = 0;

                // Technically, there are infinitely many options for the normal to change, and the area will almust always grow.
                // This simplification might introduce a small error in one timestep, but only in this exact zero-area case.
                this.dTriangleNormals_dv1[triangleId] = new Vec3D[] {Vec3D.zero, Vec3D.zero, Vec3D.zero};
                this.dTriangleNormals_dv2[triangleId] = new Vec3D[] {Vec3D.zero, Vec3D.zero, Vec3D.zero};
                this.dTriangleNormals_dv3[triangleId] = new Vec3D[] {Vec3D.zero, Vec3D.zero, Vec3D.zero};
                this.dTriangleAreas_dv1[triangleId] = Vec3D.zero;
                this.dTriangleAreas_dv2[triangleId] = Vec3D.zero;
                this.dTriangleAreas_dv3[triangleId] = Vec3D.zero;
                continue;
            }
            this.triangleNormals[triangleId] = cross_e12_e13 / crossprod_length;
            this.triangleAreas[triangleId] = crossprod_length / 2d; // Triangle area is half of the cross product of any two of its edges.

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
            Vec3D d_cross_e12_e13_length_d_e12 = 1d / crossprod_length * new Vec3D(
                    (e12.x * e13.z - e13.x * e12.z) * e13.z + (e12.x * e13.y - e13.x * e12.y) * e13.y,
                    (e12.y * e13.z - e13.y * e12.z) * e13.z + (e12.x * e13.y - e13.x * e12.y) * -e13.x,
                    (e12.y * e13.z - e13.y * e12.z) * -e13.y + (e12.x * e13.z - e13.x * e12.z) * -e13.x);
            Vec3D d_cross_e12_e13_length_d_e13 = 1f / crossprod_length * new Vec3D(
                    (e12.y * e13.z - e13.y * e12.z) * 0      + (e12.x * e13.z - e13.x * e12.z) * -e12.z + (e12.x * e13.y - e13.x * e12.y) * -e12.y,
                    (e12.y * e13.z - e13.y * e12.z) * -e12.z + (e12.x * e13.z - e13.x * e12.z) * 0      + (e12.x * e13.y - e13.x * e12.y) * e12.x,
                    (e12.y * e13.z - e13.y * e12.z) * e12.y  + (e12.x * e13.z - e13.x * e12.z) * e12.x  + (e12.x * e13.y - e13.x * e12.y) * 0);
            Vec3D d_cross_e12_e13_d_e12x = new Vec3D(0f, e13.z, e13.y);
            Vec3D d_cross_e12_e13_d_e12y = new Vec3D(e13.z, 0f, -e13.x);
            Vec3D d_cross_e12_e13_d_e12z = new Vec3D(-e13.y, -e13.x, 0f);
            Vec3D d_cross_e12_e13_d_e13x = new Vec3D(0f, -e12.z, -e12.y);
            Vec3D d_cross_e12_e13_d_e13y = new Vec3D(-e12.z, 0f, e12.x);
            Vec3D d_cross_e12_e13_d_e13z = new Vec3D(e12.y, e12.x, 0f);
            Vec3D d_n_d_e12x = (crossprod_length * d_cross_e12_e13_d_e12x - cross_e12_e13 * d_cross_e12_e13_length_d_e12.x) / (crossprod_length * crossprod_length);
            Vec3D d_n_d_e12y = (crossprod_length * d_cross_e12_e13_d_e12y - cross_e12_e13 * d_cross_e12_e13_length_d_e12.y) / (crossprod_length * crossprod_length);
            Vec3D d_n_d_e12z = (crossprod_length * d_cross_e12_e13_d_e12z - cross_e12_e13 * d_cross_e12_e13_length_d_e12.z) / (crossprod_length * crossprod_length);
            Vec3D d_n_d_e13x = (crossprod_length * d_cross_e12_e13_d_e13x - cross_e12_e13 * d_cross_e12_e13_length_d_e13.x) / (crossprod_length * crossprod_length);
            Vec3D d_n_d_e13y = (crossprod_length * d_cross_e12_e13_d_e13y - cross_e12_e13 * d_cross_e12_e13_length_d_e13.y) / (crossprod_length * crossprod_length);
            Vec3D d_n_d_e13z = (crossprod_length * d_cross_e12_e13_d_e13z - cross_e12_e13 * d_cross_e12_e13_length_d_e13.z) / (crossprod_length * crossprod_length);
            Vec3D[] d_n_d_e12 = new Vec3D[] {d_n_d_e12x, d_n_d_e12y, d_n_d_e12z};
            Vec3D[] d_n_d_e13 = new Vec3D[] {d_n_d_e13x, d_n_d_e13y, d_n_d_e13z};
            Vec3D d_e12_d_v1 = new Vec3D(-1f, -1f, -1f);
            Vec3D d_e12_d_v2 = new Vec3D(1f, 1f, 1f);
            Vec3D d_e13_d_v3 = new Vec3D(1f, 1f, 1f);
            this.dTriangleNormals_dv1[triangleId] = new Vec3D[] {d_n_d_e12[0] * d_e12_d_v1.x, d_n_d_e12[1] * d_e12_d_v1.y, d_n_d_e12[2] * d_e12_d_v1.z};
            this.dTriangleNormals_dv2[triangleId] = new Vec3D[] {d_n_d_e12[0] * d_e12_d_v2.x, d_n_d_e12[1] * d_e12_d_v2.y, d_n_d_e12[2] * d_e12_d_v2.z};
            this.dTriangleNormals_dv3[triangleId] = new Vec3D[] {d_n_d_e13[0] * d_e13_d_v3.x, d_n_d_e13[1] * d_e13_d_v3.y, d_n_d_e13[2] * d_e13_d_v3.z};

            /*
             * Triangle area partial derivatives:
             * e1 === e12 and e2 === e13 and crossprod_length === cross_e1_e2_length.
             * triangleArea = crossprod_length / 2f
             * dTriangleArea_dv1 = {d_cross_e12_e13_length_d_e12x / 2f, d_cross_e12_e13_length_d_e12x / 2f, d_cross_e12_e13_length_d_e12z / 2f} .* d_e12_d_v1
             * dTriangleArea_dv2 = {d_cross_e12_e13_length_d_e12x / 2f, d_cross_e12_e13_length_d_e12y / 2f, d_cross_e12_e13_length_d_e12z / 2f} .* d_e12_d_v2
             * dTriangleArea_dv3 = {d_cross_e12_e13_length_d_e13x / 2f, d_cross_e12_e13_length_d_e13y / 2f, d_cross_e12_e13_length_d_e13z / 2f} .* d_e13_d_v3
             */
            this.dTriangleAreas_dv1[triangleId] = new Vec3D(
                    d_cross_e12_e13_length_d_e12.x * d_e12_d_v1.x,
                    d_cross_e12_e13_length_d_e12.y * d_e12_d_v1.y,
                    d_cross_e12_e13_length_d_e12.z * d_e12_d_v1.z) / 2f;
            this.dTriangleAreas_dv2[triangleId] = new Vec3D(
                    d_cross_e12_e13_length_d_e12.x * d_e12_d_v2.x,
                    d_cross_e12_e13_length_d_e12.y * d_e12_d_v2.y,
                    d_cross_e12_e13_length_d_e12.z * d_e12_d_v2.z) / 2f;
            this.dTriangleAreas_dv3[triangleId] = new Vec3D(
                    d_cross_e12_e13_length_d_e13.x * d_e13_d_v3.x,
                    d_cross_e12_e13_length_d_e13.y * d_e13_d_v3.y,
                    d_cross_e12_e13_length_d_e13.z * d_e13_d_v3.z) / 2f;
        }
    }

    private double getSystemEnergy(int[] triangles, Vec3D[] vertices) {

        // Initialize system energy.
        double systemEnergy = 0d;

        // Compute edge length and bending energy.
        foreach(Edge edge in this.edges) {

            // Compute edge length energy.
            if(this.kLength != 0f) {
                systemEnergy += this.kLength * this.getEdgeLengthEnergy(vertices, edge.ve1, edge.ve2);
            }

            // Compute edge bending energy.
            if(edge.hasSideFlaps() && this.kBend != 0f) {
                systemEnergy += this.kBend * this.getEdgeBendEnergy(vertices, edge);
            }
        }

        // Compute triangle area energy.
        if(this.kArea != 0f) {
            for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
                int v1 = triangles[3 * triangleId];
                int v2 = triangles[3 * triangleId + 1];
                int v3 = triangles[3 * triangleId + 2];
                systemEnergy += this.kArea * this.getTriangleAreaEnergy(vertices, triangleId, v1, v2, v3);
            }
        }

        // Return the result.
        return systemEnergy;
    }

    private VecD getSystemEnergyGradient(int[] triangles, Vec3D[] vertices) {

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
                VecD edgeBendEnergyGrad = this.kBend * this.getEdgeBendEnergyGradient(vertices, edge);
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

        // Return the result.
        return vertexEnergyGradient;
    }

    /**
     * Calculate lumped vertex mass per vertex coordinate. This will return a vector in format: {m_v1x, m_v1y, m_v1z, m_v2x, ...}.
     * Note that m_vix = m_viy = m_viz = m_vi for any vertex i.
     */
    private VecD getVertexCoordinateMasses() {
        int numVertices = this.vertexPositions.Length;
        VecD vertexCoordinateMasses = new VecD(numVertices * 3);
        for(int i = 0; i < numVertices; i++) {
            double vertexArea = 0d;
            foreach(int triangleId in this.sortedVertexTriangles[i]) {
                vertexArea += this.triangleAreas[triangleId];
            }
            vertexArea /= 3d;
            double mass = vertexArea * this.shellThickness * this.shellMaterialDensity;
            vertexCoordinateMasses[3 * i] = mass;
            vertexCoordinateMasses[3 * i + 1] = mass;
            vertexCoordinateMasses[3 * i + 2] = mass;
        }
        return vertexCoordinateMasses;
    }

    /*
     * Calculates the wind force acting on each vertex.
     * Returns the wind force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
     */
    private VecD getVertexWindForce(int[] triangles, VecD[] vertices) {

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
            double triangleArea = this.triangleAreas[triangleId];
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

    /*
     * Calculates the gravity force acting on each vertex.
     * Expects the vertex coordinate masses to be in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
     * Returns the gravity force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
     */
    private VecD getVertexGravityForce(VecD vertexCoordMasses) {
        VecD gravityForce = new VecD(vertexCoordMasses.length);
        for(int yCoord = 1; yCoord < vertexCoordMasses.length; yCoord += 3) {
            gravityForce[yCoord] = -vertexCoordMasses[yCoord] * this.gravityConstant;
        }
        return gravityForce;
    }

    private double getEdgeLengthEnergy(Vec3D[] vertexPositions, int v1, int v2) {
        Vec3D edge = vertexPositions[v2] - vertexPositions[v1]; // Vector from v1 to v2.
        Vec3D undeformedEdge = new Vec3D(this.originalVertices[v2] - this.originalVertices[v1]);
        double edgeLength = edge.magnitude;
        double undeformedEdgeLength = undeformedEdge.magnitude * this.undeformedEdgeLengthFactor;
        double a = 1d - edgeLength / undeformedEdgeLength;
        return a * a * undeformedEdgeLength;
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
        double undeformedEdgeLength = (this.originalVertices[v2] - this.originalVertices[v1]).magnitude * this.undeformedEdgeLengthFactor;
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
    private MatD getEdgeLengthEnergyHess(Vec3D[] vertices, int v1, int v2) {
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
        Vec3D edge = vertexPositions[v2] - vertexPositions[v1]; // Vector from v1 to v2.
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
        double undeformedEdgeLength = (this.originalVertices[v2] - this.originalVertices[v1]).magnitude * this.undeformedEdgeLengthFactor;
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

    private double getTriangleAreaEnergy(Vec3D[] vertexPositions, int triangleId, int v1Ind, int v2Ind, int v3Ind) {
        
        // Get triangle area and undeformed triangle area.
        double triangleArea = this.triangleAreas[triangleId];
        double undeformedTriangleArea = this.undeformedTriangleAreas[triangleId];

        // Return zero if a zero undeformed area was found.
        if(undeformedTriangleArea == 0d) {
            return 0d;
        }

        // Calculate triangle energy.
        double a = 1d - triangleArea / undeformedTriangleArea;
        return a * a * undeformedTriangleArea;
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

    /*
     * Computes the triangle area energy Hessian for the triangle defined by vertices v1, v2 and v3.
     * Returns the 9x9 Hessian towards all combinations of {v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z}.
     */
    private MatD getTriangleAreaEnergyHessian(Vec3D[] vertices, int triangleId, int v1Ind, int v2Ind, int v3Ind) {
        // TODO - Test this implementation and protect against possible infinite/NaN cases if necessary.
        
        // TODO - This is a copy from the energy gradient code. Combine this in a way to prevent double calculations.
        // Return if the area is zero.
        if(this.triangleAreas[triangleId] == 0d) {
            return new MatD(9, 9); // Area is 0 m^2.
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
        //VecD d_triangleEnergy_dv1v2v3 = (2d * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2d) * d_triangleArea_dv1v2v3;
        // TODO - Copied gradient code ends here (See TODO above).

        // Compute triangle energy Hessian.
        // hess[i, j] = (crossProdLength * (a[i] * a[j] + d_a[j]_di * d + b[i] * b[j] + d_b[j]_di * e + c[i] * c[j] + d_c[j]_di * f) - (d * a[j] + e * b[j] + f * c[j]) * d_crossProdLength_di) / crossProdLength^2 / 2
        MatD areaEnergyHessian = new MatD(9, 9);
        MatD d_a_dv1v2v3 = new MatD(new double[,] {
            {0,  0,  0, 0,  0,  0, 0,  0,  0}, // d_a_dv1x
			{0,  0,  0, 0,  0, -1, 0,  0,  1}, // d_a_dv1y
			{0,  0,  0, 0,  1,  0, 0, -1,  0}, // d_a_dv1z
			{0,  0,  0, 0,  0,  0, 0,  0,  0}, // d_a_dv2x
			{0,  0,  1, 0,  0,  0, 0,  0, -1}, // d_a_dv2y
			{0, -1,  0, 0,  0,  0, 0,  1,  0}, // d_a_dv2z
			{0,  0,  0, 0,  0,  0, 0,  0,  0}, // d_a_dv3x
			{0,  0, -1, 0,  0,  1, 0,  0,  0}, // d_a_dv3y
			{0,  1,  0, 0, -1,  0, 0,  0,  0}  // d_a_dv3z
        });
        MatD d_b_dv1v2v3 = new MatD(new double[,] {
            { 0, 0,  0,  0, 0,  1,  0, 0, -1}, // d_b_dv1x
			{ 0, 0,  0,  0, 0,  0,  0, 0,  0}, // d_b_dv1y
			{ 0, 0,  0, -1, 0,  0,  1, 0,  0}, // d_b_dv1z
			{ 0, 0, -1,  0, 0,  0,  0, 0,  1}, // d_b_dv2x
			{ 0, 0,  0,  0, 0,  0,  0, 0,  0}, // d_b_dv2y
			{ 1, 0,  0,  0, 0,  0, -1, 0,  0}, // d_b_dv2z
			{ 0, 0,  1,  0, 0, -1,  0, 0,  0}, // d_b_dv3x
			{ 0, 0,  0,  0, 0,  0,  0, 0,  0}, // d_b_dv3y
			{-1, 0,  0,  1, 0,  0,  0, 0,  0}  // d_b_dv3z
        });
        MatD d_c_dv1v2v3 = new MatD(new double[,] {
            { 0,  0, 0,  0, -1, 0,  0,  1, 0}, // d_c_dv1x
			{ 0,  0, 0,  1,  0, 0, -1,  0, 0}, // d_c_dv1y
			{ 0,  0, 0,  0,  0, 0,  0,  0, 0}, // d_c_dv1z
			{ 0,  1, 0,  0,  0, 0,  0, -1, 0}, // d_c_dv2x
			{-1,  0, 0,  0,  0, 0,  1,  0, 0}, // d_c_dv2y
			{ 0,  0, 0,  0,  0, 0,  0,  0, 0}, // d_c_dv2z
			{ 0, -1, 0,  0,  1, 0,  0,  0, 0}, // d_c_dv3x
			{ 1,  0, 0, -1,  0, 0,  0,  0, 0}, // d_c_dv3y
			{ 0,  0, 0,  0,  0, 0,  0,  0, 0}  // d_c_dv3z
        });
        for(int i = 0; i < areaEnergyHessian.numRows; i++) {
            for(int j = 0; j < areaEnergyHessian.numColumns; j++) {
                areaEnergyHessian[i, j] = (
                        (a[i] * a[j] + d_a_dv1v2v3[j, i] * d + b[i] * b[j] + d_b_dv1v2v3[j, i] * e + c[i] * c[j] + d_b_dv1v2v3[j, i] * f)
                        - (d * a[j] + e * b[j] + f * c[j]) * d_crossProdLength_dv1v2v3[i] / crossProdLength
                        ) / crossProdLength / 2d;
            }
        }
        return areaEnergyHessian;
    }

    private double getEdgeBendEnergy(Vec3D[] vertexPositions, Edge edge) {

        // Define required constants.
        double undeformedEdgeLength = (this.originalVertices[edge.ve2] - this.originalVertices[edge.ve1]).magnitude * this.undeformedEdgeLengthFactor;
        // h_e_undeformed is a third of the average triangle height, where the height is twice the triangle area divided by the triangle width.
        double h_e_undeformed = (this.undeformedTriangleAreas[edge.triangleId1] + this.undeformedTriangleAreas[edge.triangleId2]) / undeformedEdgeLength / 3d;
        double teta_e_undeformed = 0; // TODO - Add option to use a non-flat rest state depending on the useFlatUndeformedBendState field.

        // Get triangle normals.
        Vec3D n1 = this.triangleNormals[edge.triangleId1];
        Vec3D n2 = this.triangleNormals[edge.triangleId2];

        // Return if at least one of the triangles has a zero area and no normal.
        if(n1 == null || n2 == null) {
            return 0d;
        }
        
        // Calculate hinge angle.
        double teta_e_sign = Mathf.Sign((float) VecD.dot(n1, vertexPositions[edge.vf2] - vertexPositions[edge.ve1])); // 1 if teta_e positive, -1 if negative.
        double dot_n1_n2 = VecD.dot(n1, n2);
        double teta_e = -Mathf.Acos((float) (dot_n1_n2 > 1d ? 1d : (dot_n1_n2 < -1d ? -1d : dot_n1_n2))) * teta_e_sign; // Limit the dot product of the normals at 1 (fix precision errors).

        double a = teta_e - teta_e_undeformed;
        return a * a * undeformedEdgeLength / h_e_undeformed;
    }

    /**
     * Computes the bending energy gradient of the given edge. This gradient is taken towards the edge-defining vertices ve1 and ve2, as well
     * as the vertices vf1 and vf2, completing triangles 1 and 2 connected to the edge respectively.
     * The return values are the bending energy gradient towards: {ve1x. ve1y, ve1z, ve2x. ve2y, ve2z, vf1x. vf1y, vf1z, vf2x. vf2y, vf2z}.
     */
    private VecD getEdgeBendEnergyGradient(Vec3D[] vertices, Edge edge) {

        // Define required constants.
        double undeformedEdgeLength = (this.originalVertices[edge.ve2] - this.originalVertices[edge.ve1]).magnitude * this.undeformedEdgeLengthFactor;
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
     * Computes the bending energy Hessian of the given edge. This Hessian is taken towards the edge-defining vertices ve1 and ve2, as well
     * as the vertices vf1 and vf2, completing triangles 1 and 2 connected to the edge respectively.
     * The returned Hessian contains partial derivatives in rows and columns in order: {ve1x. ve1y, ve1z, ve2x. ve2y, ve2z, vf1x. vf1y, vf1z, vf2x. vf2y, vf2z}.
     */
    private MatD getEdgeBendEnergyHess(Vec3D[] vertices, Edge edge) {

        // Define required constants.
        double undeformedEdgeLength = (this.originalVertices[edge.ve2] - this.originalVertices[edge.ve1]).magnitude * this.undeformedEdgeLengthFactor;
        // h_e_undeformed is a third of the average triangle height, where the height is twice the triangle area divided by the triangle width.
        double h_e_undeformed = (this.undeformedTriangleAreas[edge.triangleId1] + this.undeformedTriangleAreas[edge.triangleId2]) / undeformedEdgeLength / 3d;
        double teta_e_undeformed = 0; // TODO - Add option to use a non-flat rest state depending on the useFlatUndeformedBendState field.

        // Get triangle normals.
        Vec3D n1 = this.triangleNormals[edge.triangleId1];
        Vec3D n2 = this.triangleNormals[edge.triangleId2];

        // Return if at least one of the triangles has a zero area and no normal.
        if(n1 == null || n2 == null) {
            return new MatD(12, 12);
        }
        
        // Calculate hinge angle.
        double teta_e_sign = Mathf.Sign((float) VecD.dot(n1, vertices[edge.vf2] - vertices[edge.ve1])); // 1 if teta_e positive, -1 if negative.
        double dot_n1_n2 = VecD.dot(n1, n2);
        double teta_e = -Mathf.Acos((float) (dot_n1_n2 > 1d ? 1d : (dot_n1_n2 < -1d ? -1d : dot_n1_n2))) * teta_e_sign; // Limit the dot product of the normals at 1 (fix precision errors).

        // TODO - Remove this check if it can't happen. Check what happens for non-existing normals though.
        if(double.IsNaN(teta_e)) {
            print("NaN teta_e. sign: " + teta_e_sign + ", dot(n1, n2): " + VecD.dot(n1, n2));
            // teta_e = Mathf.PI; // The triangles are on top of each other, which is both 180 and -180 degrees.
            return new MatD(12, 12);
        }
        
        // Calculate energy derivative towards hinge angle teta.
        double d_fi_d_teta = 2d * (teta_e - teta_e_undeformed) * undeformedEdgeLength / h_e_undeformed;

        // Return if the energy gradient is zero.
        if(d_fi_d_teta == 0d) {
            return new MatD(12, 12);
        }

        // Define variables following paper: https://studios.disneyresearch.com/wp-content/uploads/2019/03/Discrete-Bending-Forces-and-Their-Jacobians-Paper.pdf.
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
        
        double t1_omega_00 = 1 / (t1_h0 * t1_h0);
        double t1_omega_01 = 1 / (t1_h0 * t1_h1);
        double t1_omega_02 = 1 / (t1_h0 * t1_h2);
        double t1_omega_10 = t1_omega_01;
        double t1_omega_11 = 1 / (t1_h1 * t1_h1);
        double t1_omega_12 = 1 / (t1_h1 * t1_h2);
        double t1_omega_20 = t1_omega_02;
        double t1_omega_21 = t1_omega_12;
        double t1_omega_22 = 1 / (t1_h2 * t1_h2);
        double t2_omega_00 = 1 / (t2_h0 * t2_h0);
        double t2_omega_01 = 1 / (t2_h0 * t2_h1);
        double t2_omega_02 = 1 / (t2_h0 * t2_h2);
        double t2_omega_10 = t2_omega_01;
        double t2_omega_11 = 1 / (t2_h1 * t2_h1);
        double t2_omega_12 = 1 / (t2_h1 * t2_h2);
        double t2_omega_20 = t2_omega_02;
        double t2_omega_21 = t2_omega_12;
        double t2_omega_22 = 1 / (t2_h2 * t2_h2);
        
        Vec3D t1_e0_normal = Vec3D.cross(e0.unit, n1); // "m_0" in paper.
        Vec3D t1_e1_normal = Vec3D.cross(t1_e1.unit, n1); // "m_1" in paper.
        Vec3D t1_e2_normal = Vec3D.cross(t1_e2.unit, n1); // "m_2" in paper.
        Vec3D t2_e0_normal = Vec3D.cross(e0.unit, n2); // "~m_0" in paper.
        Vec3D t2_e1_normal = Vec3D.cross(t2_e1.unit, n2); // "~m_1" in paper.
        Vec3D t2_e2_normal = Vec3D.cross(t2_e2.unit, n2); // "~m_2" in paper.

        MatD t1_M0 = MatD.fromVecMultiplication(n1, t1_e0_normal);
        MatD t1_M1 = MatD.fromVecMultiplication(n1, t1_e1_normal);
        MatD t1_M2 = MatD.fromVecMultiplication(n1, t1_e2_normal);
        MatD t2_M0 = MatD.fromVecMultiplication(n2, t2_e0_normal);
        MatD t2_M1 = MatD.fromVecMultiplication(n2, t2_e1_normal);
        MatD t2_M2 = MatD.fromVecMultiplication(n2, t2_e2_normal);
        
        MatD t1_Q0 = t1_omega_00 * t1_M0;
        MatD t1_Q1 = t1_omega_01 * t1_M1;
        MatD t1_Q2 = t1_omega_02 * t1_M2;
        MatD t2_Q0 = t2_omega_00 * t2_M0;
        MatD t2_Q1 = t2_omega_01 * t2_M1;
        MatD t2_Q2 = t2_omega_02 * t2_M2;

        double e0_magnitude = e0.magnitude;
        MatD t1_N0 = t1_M0 / (e0_magnitude * e0_magnitude);
        MatD t2_N0 = t2_M0 / (e0_magnitude * e0_magnitude);
        
        MatD t1_P10 = t1_omega_10 * Mathf.Cos(t1_alpha1) * t1_M0.transpose;
        MatD t1_P20 = t1_omega_20 * Mathf.Cos(t1_alpha2) * t1_M0.transpose;
        MatD t1_P12 = t1_omega_12 * Mathf.Cos(t1_alpha1) * t1_M2.transpose;
        MatD t1_P21 = t1_omega_21 * Mathf.Cos(t1_alpha2) * t1_M1.transpose;
        MatD t1_P11 = t1_omega_11 * Mathf.Cos(t1_alpha1) * t1_M1.transpose;
        MatD t1_P22 = t1_omega_22 * Mathf.Cos(t1_alpha2) * t1_M2.transpose;
        MatD t2_P10 = t2_omega_10 * Mathf.Cos(t2_alpha1) * t2_M0.transpose;
        MatD t2_P20 = t2_omega_20 * Mathf.Cos(t2_alpha2) * t2_M0.transpose;
        MatD t2_P12 = t2_omega_12 * Mathf.Cos(t2_alpha1) * t2_M2.transpose;
        MatD t2_P21 = t2_omega_21 * Mathf.Cos(t2_alpha2) * t2_M1.transpose;
        MatD t2_P11 = t2_omega_11 * Mathf.Cos(t2_alpha1) * t2_M1.transpose;
        MatD t2_P22 = t2_omega_22 * Mathf.Cos(t2_alpha2) * t2_M2.transpose;

        // Construct 3x3 building blocks for the teta Hessian.
        MatD teta_hess_00 = -getMatPlusTransposedMat(t1_Q0);
        MatD teta_hess_33 = -getMatPlusTransposedMat(t2_Q0);
        MatD teta_hess_11 = getMatPlusTransposedMat(t1_P11) - t1_N0 + getMatPlusTransposedMat(t2_P11) - t2_N0;
        MatD teta_hess_22 = getMatPlusTransposedMat(t1_P22) - t1_N0 + getMatPlusTransposedMat(t2_P22) - t2_N0;
        MatD teta_hess_10 = getMatPlusTransposedMat(t1_P10) - t1_Q1;
        MatD teta_hess_20 = getMatPlusTransposedMat(t1_P20) - t1_Q2;
        MatD teta_hess_13 = getMatPlusTransposedMat(t2_P10) - t2_Q1;
        MatD teta_hess_23 = getMatPlusTransposedMat(t2_P20) - t2_Q2;
        MatD teta_hess_12 = t1_P12 + t1_P21.transpose + t1_N0 + t2_P12 + t2_P21.transpose + t2_N0;
        MatD teta_hess_03 = new MatD(3, 3);

        // Construct teta Hessian.
        MatD teta_hess = new MatD(12, 12);
        for(int row = 0; row < 3; row++) {
            for(int col = 0; col < 3; col++) {

                // Set Hessian diagonal building blocks.
                teta_hess[row, col] = teta_hess_00[row, col];
                teta_hess[row + 3, col + 3] = teta_hess_11[row, col];
                teta_hess[row + 6, col + 6] = teta_hess_22[row, col];
                teta_hess[row + 9, col + 9] = teta_hess_33[row, col];

                // Set remaining symmetric non-diagonal building blocks.
                teta_hess[row + 3, col] = teta_hess_10[row, col];
                teta_hess[row, col + 3] = teta_hess_10[col, row];
                
                teta_hess[row + 6, col] = teta_hess_20[row, col];
                teta_hess[row, col + 6] = teta_hess_20[col, row];
                
                teta_hess[row + 3, col + 9] = teta_hess_13[row, col];
                teta_hess[row + 9, col + 3] = teta_hess_13[col, row];
                
                teta_hess[row + 6, col + 9] = teta_hess_23[row, col];
                teta_hess[row + 9, col + 6] = teta_hess_23[col, row];
                
                teta_hess[row + 3, col + 6] = teta_hess_12[row, col];
                teta_hess[row + 6, col + 3] = teta_hess_12[col, row];
                
                teta_hess[row, col + 9] = teta_hess_03[row, col];
                teta_hess[row + 9, col] = teta_hess_03[col, row];
            }
        }

        // Calculate derivatives of teta towards all 4 vertices.
        VecD d_teta_d_ve1 = Mathf.Cos(t1_alpha2) / t1_h1 * n1 + Mathf.Cos(t2_alpha2) / t2_h1 * n2;
        VecD d_teta_d_ve2 = Mathf.Cos(t1_alpha1) / t1_h2 * n1 + Mathf.Cos(t2_alpha1) / t2_h2 * n2;
        VecD d_teta_d_vf1 = -n1 / t1_h0;
        VecD d_teta_d_vf2 = -n2 / t2_h0;
        
        // Assemble hinge angle gradient.
        VecD d_teta_d_ve1ve2vf1vf2 = new VecD(d_teta_d_ve1, d_teta_d_ve2, d_teta_d_vf1, d_teta_d_vf2);

        //// Calculate the bending energy gradient.
        //VecD bendEnergyGrad = d_fi_d_teta * d_teta_d_ve1ve2vf1vf2;

        // Calculate bending energy Hessian.
        /*
         * bendEnergyHess(fi) = fi' * Hess(teta) + fi'' * grad_teta_trans * grad_teta
         * Discrete shells bending energy: fi_i(teta_i) = ||e|| / h_e_undeformed * (teta_i - teta_e_undeformed)^2
         * fi' = d_fi_d_teta = ||e|| / h_e_undeformed * 2 * (teta - teta_e_undeformed)
         * fi'' = dd_fi_d_teta_d_teta = ||e|| / h_e_undeformed * 2
         */
        return d_fi_d_teta * teta_hess + 2d * undeformedEdgeLength / h_e_undeformed * MatD.fromVecMultiplication(d_teta_d_ve1ve2vf1vf2, d_teta_d_ve1ve2vf1vf2);
    }

    private static MatD getMatPlusTransposedMat(MatD mat) {
        if(mat.numRows != mat.numColumns) {
            throw new Exception("Given matrix is not a square matrix.");
        }
        MatD ret = new MatD(mat.numRows, mat.numColumns);
        for(int row = 0; row < mat.numRows; row++) {
            for(int col = 0; col < mat.numColumns; col++) {
                ret[row, col] = mat[col, row] + mat[row, col];
            }
        }
        return ret;
    }

    /**
     * Makes the given Hessian positive definite by adding a multiple of the identity matrix to it.
     */
    private static void makeHessPositiveDefinite(MatD hess) {
        if(hess.numRows != hess.numColumns) {
            throw new Exception("Given matrix is not a square matrix.");
        }
        MatD mat = hess.Clone();
        double diagIncrement = 0;
        for(int i = 0; i < mat.numRows; i++) {

            // Make pivot positive.
            if(mat[i, i] <= 0) {
                double inc = -mat[i, i] * 1.0001d; // Add small delta to ensure a positive non-zero value. Multiply to prevent rounding to 0.
                if(inc == 0d) {
                    inc = 0.0001d; // Set small delta to ensure a positive non-zero value.
                }
                diagIncrement += inc;
                mat.addDiag(inc);
            }

            // Sweep lower rows.
            for(int row = i + 1; row < mat.numRows; row++) {
                double factor = mat[row, i] / mat[i, i];
                for(int col = 0; col < mat.numColumns; col++) {
                    mat[row, col] -= factor * mat[i, col];
                }
            }
        }
        if(diagIncrement != 0d) {
            hess.addDiag(diagIncrement);
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
                if(mat[row, col] != 0d) {
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
