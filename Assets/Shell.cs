using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading.Tasks;
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
	public double timeScale = 1f;
	public Vector3 windPressure; // [N/m^2]. TODO - Could also apply scalar pressure in triangle normal directions.
	public double gravityConstant = 9.81d;

	public float measurementsGenerateFactor = 0.1f; // Defines the number of generated measurements by multiplying this with the amount of vertices.
	private Vec3D[] measurements = null; // Measurements. One element for each vertex, null meaning that there is no measurement for that vertex.

	// Gradient descent specific settings.
	public double kGradientDescent;
	public double maxGradientDescentStep;

	// Explicit integration specific settings.
	public float dampingFactor = 0.99f; // [F * s / m] = [kg / s] ? Factor applied to vertex velocity per time step. TODO - Replace with proper energy dissipation.

	// Optimization integrator specific settings.
	public double dampingConstant = 1d;
	public double kMeasurementsPenalty = 1d;
	public double eGradientMagnitudeTerminationThreshold = 0.5d;
	private double lastLineSearchAlpha = 1d;
	public double MaxLineSearchTimeMS = 10000;
	public double MinLineSearchAlpha = 0.000001d;
	public double MaxNewtonsMethodLoopTimeMS = 10000;

	// Cached vertex/triangle properties.
	private Vec3D[] triangleNormals;
	private Vec3D[] undeformedTriangleNormals;
	private double[] triangleAreas;
	private double[] undeformedTriangleAreas;

	// Simulation recording.
	private MeshRecorder meshRecorder = null;
	private Boolean isRecording = false;

	void Awake() {
		QualitySettings.vSyncCount = 0; // Disable V-sync.
		Application.targetFrameRate = 100; // Set max framerate.
		System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US"); // To print decimal points instead of commas.
	}

	// Start is called before the first frame update.
	void Start() {

		// Initialize fields. Doing this overwrites the values set in Unity's inspector.
		this.kLength = 10f;
		this.kArea = 10f;
		this.kBend = 0.1f;
		this.kGradientDescent = 1f;
		this.maxGradientDescentStep = 0.01f;
		this.windPressure = new Vector3(0f, 0f, 10f);

		// Create the new object in the scene.
		Mesh mesh;
		double undeformedInnerEdgeLengthFactor;
		if(this.shellObj == null) {
			this.shellObj = new GameObject();
			this.shellObj.AddComponent<MeshRenderer>();
			this.shellObj.AddComponent<MeshFilter>();
			//mesh = MeshHelper.createTriangleMesh(5, 5, 0);
			//mesh = MeshHelper.createSquareMesh(5, 5, 1);
			mesh = MeshHelper.createTriangleMesh(5, 5, 5); // 5 subdivisions leads to 561 vertices and 3072 triangles.
			undeformedInnerEdgeLengthFactor = 1.1d;
		} else {
			mesh = this.shellObj.GetComponent<MeshFilter>().mesh;
			undeformedInnerEdgeLengthFactor = 1d;
		}

		// Set the mesh renderer.
		MeshRenderer meshRenderer = this.shellObj.GetComponent<MeshRenderer>();
		meshRenderer.sharedMaterial = new Material(Shader.Find("Custom/StandardTwoSides"));

		// Load the mesh.
		this.loadMesh(mesh, undeformedInnerEdgeLengthFactor);

		// Set the shell position.
		this.shellObj.transform.position = new Vector3(0, 1, 0);
	}

	private void loadMesh(Mesh mesh, double undeformedInnerEdgeLengthFactor) {

		// Set the mesh in the mesh filter.
		MeshFilter meshFilter = this.shellObj.GetComponent<MeshFilter>();
		meshFilter.mesh = mesh;

		// Store undeformed vertices.
		this.originalVertices = (Vector3[]) mesh.vertices.Clone();

		// Get mesh variables for convenience.
		int[] triangles = mesh.triangles;
		Vector3[] vertices = mesh.vertices;
		int numVertices = mesh.vertexCount;

		// Create clockwise sorted triangle list per vertex.
		this.sortedVertexTriangles = MeshUtils.getSortedVertexTriangles(triangles, numVertices);

		// Cache mesh edges.
		this.edges = MeshUtils.getEdges(this.sortedVertexTriangles, triangles);

		// Set undeformed edge lengths.
		foreach(Edge edge in this.edges) {
			edge.undeformedLength = (new VecD(this.originalVertices[edge.ve2]) - new VecD(this.originalVertices[edge.ve1])).magnitude;

			// Apply undeformed edge length factor on inner edges.
			if(edge.hasSideFlaps()) {
				edge.undeformedLength *= undeformedInnerEdgeLengthFactor;
			}
		}

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
		this.triangleAreas = new double[numTriangles];
		this.undeformedTriangleNormals = new Vec3D[numTriangles];
		this.undeformedTriangleAreas = new double[numTriangles];
		for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
			int triangleBaseIndex = triangleId * 3;
			int v1 = triangles[triangleBaseIndex];
			int v2 = triangles[triangleBaseIndex + 1];
			int v3 = triangles[triangleBaseIndex + 2];
			Vec3D crossProd = Vec3D.cross(
					(new Vec3D(this.originalVertices[v2]) - new Vec3D(this.originalVertices[v1])) * undeformedInnerEdgeLengthFactor,
					(new Vec3D(this.originalVertices[v3]) - new Vec3D(this.originalVertices[v1])) * undeformedInnerEdgeLengthFactor); // TODO - Don't apply factor to outer edges.
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
		//	Vector3 vertex = vertices[i];
		//	if(vertex.x > vertices[vertMaxX].x) {
		//		vertMaxX = i;
		//	} else if(vertex.x < vertices[vertMinX].x) {
		//		vertMinX = i;
		//	}
		//	if(vertex.y > vertices[vertMaxY].y) {
		//		vertMaxY = i;
		//	} else if(vertex.y < vertices[vertMinY].y) {
		//		vertMinY = i;
		//	}
		//	if(vertex.z > vertices[vertMaxZ].z) {
		//		vertMaxZ = i;
		//	} else if(vertex.z < vertices[vertMinZ].z) {
		//		vertMinZ = i;
		//	}
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
		//	if(this.verticesMovementConstraints[i]) {
		//		vertices[i].z += 2;
		//	}
		//}
		//mesh.vertices = vertices;
		//mesh.RecalculateNormals();

		print("Loaded mesh with " + mesh.vertices.Length + " vertices, " + this.edges.Count + " edges and " + mesh.triangles.Length + " triangles.");
	}

	private void reset() {
		Mesh mesh = this.getMesh();
		Vector3[] vertices = mesh.vertices;
		for(int i = 0; i < mesh.vertices.Length; i++) {
			vertices[i] = new Vector3(this.originalVertices[i].x, this.originalVertices[i].y, this.originalVertices[i].z);

			this.vertexPositions[i] = new Vec3D(this.originalVertices[i].x, this.originalVertices[i].y, this.originalVertices[i].z);
			for(int coord = 0; coord < 3; coord++) {
				this.vertexVelocities[3 * i + coord] = 0;
				this.vertexAccelerations[3 * i + coord] = 0;
			}
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
		if(Input.GetKeyDown(KeyCode.Q)) {
			this.doUpdate = false;
			print("Simulation paused.");
		}

		// Handle single step hotkey.
		if(Input.GetKeyDown(KeyCode.G)) {
			double deltaTime = Time.fixedDeltaTime * this.timeScale;
			this.simulationStep(deltaTime);
			print("Performed single step (" + (deltaTime < 1d ? ((deltaTime * 1000d) + "ms") : deltaTime + "s") + ").");
		}

		// Handle record start/stop/replay.
		if(Input.GetKeyDown(KeyCode.Y)) {
			if(this.isRecording) {
				print("Recording already active.");
			} else {
				this.meshRecorder = new MeshRecorder(this.getMesh().triangles, this.vertexPositions);
				this.isRecording = true;
				print("New recording started.");
			}
		}
		if(Input.GetKeyDown(KeyCode.U)) {
			if(!this.isRecording) {
				print("No recording running.");
			} else {
				this.isRecording = false;
				print("Recording stopped.");
			}
		}
		if(Input.GetKeyDown(KeyCode.I)) {
			if(this.meshRecorder == null) {
				print("No recording available for replay.");
			} else if(this.meshRecorder.isPlaying()) {
				print("Recording is already playing.");
			} else {
				this.meshRecorder.replay();
			}
		}

		// Update recording.
		if(this.meshRecorder != null) {
			this.meshRecorder.update();
		}
	}

	// FixedUpdate is called every fixed interval (Edit -> Project Settings -> Time -> Fixed Timestep).
	void FixedUpdate() {

		// Don't update if no noticable time has passed or when the simulation has been paused.
		if(Time.deltaTime == 0 || Time.timeScale == 0 || (!this.doUpdate && !Input.GetKey(KeyCode.F))) {
			return;
		}

		// Perform a simulation step.
		this.simulationStep(Time.deltaTime * this.timeScale);
	}

	private void simulationStep(double deltaTime) {

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

		// Update mesh recorder.
		if(this.isRecording) {
			this.meshRecorder.record(deltaTime, this.vertexPositions);
		}
	}

	private void doGradientDescentStep() {
		int[] triangles = this.getMesh().triangles;
		VecD vertexEnergyGradient = this.getSystemEnergyGradient(triangles, this.vertexPositions);
		VecD vertexWindForce = this.getVertexWindForce(triangles, this.vertexPositions, new VecD(this.windPressure));
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
		VecD vertexWindForce = this.getVertexWindForce(triangles, this.vertexPositions, new VecD(this.windPressure));

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
			MatD lengthEnergyHess = this.getEdgeLengthEnergyHess(this.vertexPositions, edge);
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
		 */
		
		// Start timer.
		Stopwatch stopWatch = Stopwatch.StartNew();
		
		// Declare constants.
		double terminationThreshold = this.eGradientMagnitudeTerminationThreshold;
		double kappa = 0.01d; // Value as proposed by Optimization Integrator paper.
		int numVertices = this.vertexPositions.Length;
		double maxStepMagnitude = 1000d; //numVertices * 0.001; // TODO - Set sensible value. Optimization Integrator paper uses 1000 (mesh size dependent?).
		int[] triangles = this.getMesh().triangles;

		// Get vertex masses. These are assumed to be constant, which is an okay approximation for small steps.
		VecD vertexCoordMasses = this.getVertexCoordinateMasses(); // Masses per vertex coordinate. Format: {m_v1x, m_v1y, m_v1z, ...}.

		// Assemble system-wide energy Hessian.
		MathNet.Numerics.LinearAlgebra.Double.SparseMatrix energyHess =
			this.getSystemEnergyHessianSparseRepresentationMultiThreadedTriplets(this.vertexPositions, triangles);

		// Get E Hessian. This is a simplified version without measurements penalty, wind force and gravity force.
		/**
		  * Get E Hessian.
		  * This is a simplified version without measurements penalty and wind force.
		  * The gravity force is not dependent on the positions, so the gravity Hessian is a zero matrix.
		  * The wind force is dependent on the rotation and area of the triangles, but for small steps, this is constant enough to ignore.
		  * TODO - Check whether the measurements penalty Hessian should be included here. The length energy Hessian code can be used for this.
		  */
		MathNet.Numerics.LinearAlgebra.Double.SparseMatrix eHess = MathNet.Numerics.LinearAlgebra.Double.SparseMatrix.OfMatrix(energyHess);
		double deltaTimeSquare = deltaTime * deltaTime;
		for(int i = 0; i < eHess.RowCount; i++) {
			eHess[i, i] += vertexCoordMasses[i] / deltaTimeSquare;
		}
		
		// Perform Newton's method, setting steps until the termination criterion has been met.
		VecD vertexPositionsFlat = new VecD(this.vertexPositions);
		Vec3D[] newVertexPositions = new Vec3D[numVertices];
		for(int i = 0; i < numVertices; i++) {

			// Set initial guess: pos + deltaTime * velocity.
			newVertexPositions[i] = new Vec3D();
			for(int coord = 0; coord < 3; coord++) {
				newVertexPositions[i][coord] = this.vertexPositions[i][coord] + deltaTime * this.vertexVelocities[i * 3 + coord];
			}
		}
		int iteration = 0;
		while(true) {

			// Limit amount of time spent to prevent endless loops.
			iteration++;
			if(stopWatch.ElapsedMilliseconds > this.MaxNewtonsMethodLoopTimeMS) {
				print(stopWatch.ElapsedMilliseconds + "ms: Maximum time reached in Optimization Integrator update after "
						+ iteration + " iterations. Returning without taking a step.");
				return;
			}

			// Get system-wide energy gradient.
			this.recalcTriangleNormalsAndAreas(triangles, newVertexPositions);
			VecD energyGradient = this.getSystemEnergyGradient(triangles, newVertexPositions);

			// Get wind force.
			VecD windForce = this.getVertexWindForce(triangles, newVertexPositions, new VecD(this.windPressure));

			// Get gravity force.
			VecD gravityForce = this.getVertexGravityForce(vertexCoordMasses);

			// Get damping force.
			VecD nextVertexVelocities = new VecD(newVertexPositions).sub(vertexPositionsFlat);
			VecD dampingForce = -this.dampingConstant * (energyHess * nextVertexVelocities);

			// Get measurement penalty gradient.
			// Energy fi_penalty = kPenalty * sum(measurements k) {(x - k)^2}
			// Energy gradient d_fi_penalty_dx = kPenalty * sum(measurements k) {2 * (x - k)}
			VecD penaltyEnergyGradient = new VecD(3 * numVertices);
			if(this.measurements != null) {
				for(int vertexInd = 0; vertexInd < numVertices; vertexInd++) {
					if(this.measurements[vertexInd] != null) {
						for(int coord = 0; coord < 3; coord++) {
							penaltyEnergyGradient[3 * vertexInd + coord] =
								this.kMeasurementsPenalty * 2d * (this.vertexPositions[vertexInd][coord] - this.measurements[vertexInd][coord]);
						}
					}
				}
			}

			// Set forces to zero for constrained vertices.
			// This causes the E gradient to be zero for them as well, causing it not to get in the way of the minimization problem.
			for(int i = 0; i < numVertices; i++) {
				if(this.verticesMovementConstraints[i]) {
					for(int coord = 0; coord < 3; coord++) {
						energyGradient[3 * i + coord] = 0;
						windForce[3 * i + coord] = 0;
						gravityForce[3 * i + coord] = 0;
						dampingForce[3 * i + coord] = 0;
						penaltyEnergyGradient[3 * i + coord] = 0;
					}
				}
			}

			// Get E gradient.
			VecD eGradient = new VecD(newVertexPositions).sub(vertexPositionsFlat).sub(deltaTime * this.vertexVelocities)
					.multiplyElementWise(vertexCoordMasses).div(deltaTimeSquare)
					.add(energyGradient).add(penaltyEnergyGradient).sub(windForce).sub(gravityForce).sub(dampingForce);

			// Terminate when the termination criterion has been met.
			double eGradientMagnitude = eGradient.magnitude;
			print(stopWatch.ElapsedMilliseconds + "ms: E gradient magnitude: " + eGradientMagnitude + " (threshold: " + terminationThreshold
					+ ", iteration: " + iteration + ")");
			if(eGradientMagnitude < terminationThreshold) {
				print(stopWatch.ElapsedMilliseconds + "ms: Finished on iteration: " + iteration + ".");
				if(iteration == 1) {
					this.doUpdate = false;
					print("Finished on the first iteration. Pausing simulation.");
				}
				break;
			}

			// Compute Newton step.
			// step = -inverse(eHess) * eGradient === eHess * step = -eGradient
			//VecD step = -eGradient;
			VecD step = this.sparseDirectSolve(eHess, -eGradient);

			// Prevent constrained vertices from moving.
			for(int i = 0; i < numVertices; i++) {
				if(this.verticesMovementConstraints[i]) {
					for(int coord = 0; coord < 3; coord++) {
						step[3 * i + coord] = 0;
					}
				}
			}

			// Ensure that the step is in downhill direction.
			// If a < b, then the step is suitable. Otherwise try -a < b. As a last resort, fall back to gradient descent.
			double stepMagnitude = step.magnitude;
			double a = VecD.dot(step, eGradient);
			double b = -kappa * stepMagnitude * eGradientMagnitude;
			if(a >= b) {
				if(-a < b) {
					step = -step;
				} else {
					step = -eGradient;
					stepMagnitude = eGradientMagnitude;
				}
			}

			// Clamp max step magnitude.
			if(stepMagnitude > maxStepMagnitude) {
				step.mul(maxStepMagnitude / stepMagnitude); // Equivalent of: step = step / stepMagnitude * maxStepMagnitude;
			}

			// Choose step size alpha in direction delta_x using a line search.

			// Define new vertex positions for the line search, starting at the current new vertex positions.
			Vec3D[] newNewVertexPositions = new Vec3D[numVertices];
			for(int i = 0; i < numVertices; i++) {
				newNewVertexPositions[i] = new Vec3D();
			}

			// Take steps with adjusting alpha until sufficient decrease has been achieved.
			double alpha = this.lastLineSearchAlpha;
			double bestAlpha = double.NaN;
			double e = this.getE(triangles, vertexPositionsFlat, newVertexPositions,
					vertexCoordMasses, new VecD(this.windPressure), energyHess, deltaTime);
			double bestE = e;
			double c = 1d;
			long lineSearchLoopStartTime = stopWatch.ElapsedMilliseconds;
			while(true) {
				if(stopWatch.ElapsedMilliseconds - lineSearchLoopStartTime > this.MaxLineSearchTimeMS) {
					print(stopWatch.ElapsedMilliseconds + "ms: Spent too long in line search. Breaking with alpha = " + alpha);
					break;
				}

				// Take step: x += alpha * step.
				for(int i = 0; i < newNewVertexPositions.Length; i++) {
					newNewVertexPositions[i].x = newVertexPositions[i].x + alpha * step[3 * i];
					newNewVertexPositions[i].y = newVertexPositions[i].y + alpha * step[3 * i + 1];
					newNewVertexPositions[i].z = newVertexPositions[i].z + alpha * step[3 * i + 2];
				}

				// Terminate when there is sufficient E decrease. Adjust alpha otherwise.
				this.recalcTriangleNormalsAndAreas(triangles, newNewVertexPositions);
				double newE = this.getE(triangles, vertexPositionsFlat, newNewVertexPositions,
						vertexCoordMasses, new VecD(this.windPressure), energyHess, deltaTime);
				if(newE <= bestE) {
				//if(newE <= e + c * alpha * VecD.dot(eGradient, step)) {
					
					// Alpha is suitable, but there might be a higher value of alpha that is still suitable. Increase alpha and store the current best alpha.
					bestAlpha = alpha;
					bestE = newE;
					alpha *= 2;
				} else if(!double.IsNaN(bestAlpha)) {

					// A best alpha was set, but a higher value of alpha didn't make it. Return the best value.
					alpha = bestAlpha;
					print(stopWatch.ElapsedMilliseconds + "ms: Terminating with alpha: " + alpha
							+ " after " + (stopWatch.ElapsedMilliseconds - lineSearchLoopStartTime) + "ms (bestE: " + bestE + ").");
					break;
				} else {

					// Alpha is too high to be suitable. Decrease alpha.
					alpha /= 2d;

					// Just take the step if alpha gets too small.
					if(alpha <= this.MinLineSearchAlpha) {
						alpha = this.MinLineSearchAlpha;
						print(stopWatch.ElapsedMilliseconds + "ms: Alpha is getting too small. Setting alpha: " + alpha);
						break;
					}
				}
			}
			this.lastLineSearchAlpha = alpha;

			// Take step: x += alpha * step.
			for(int i = 0; i < newVertexPositions.Length; i++) {
				newVertexPositions[i].x += alpha * step[3 * i];
				newVertexPositions[i].y += alpha * step[3 * i + 1];
				newVertexPositions[i].z += alpha * step[3 * i + 2];
			}

		}
		
		// Update vertex velocity.
		this.vertexVelocities = new VecD(newVertexPositions).sub(vertexPositionsFlat).div(deltaTime);

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
				continue;
			}
			this.triangleNormals[triangleId] = cross_e12_e13 / crossprod_length;
			this.triangleAreas[triangleId] = crossprod_length / 2d; // Triangle area is half of the cross product of any two of its edges.
		}
	}

	private double getSystemEnergy(int[] triangles, Vec3D[] vertices) {

		// Initialize system energy.
		double systemEnergy = 0d;

		// Compute edge length and bending energy.
		foreach(Edge edge in this.edges) {

			// Compute edge length energy.
			if(this.kLength != 0f) {
				systemEnergy += this.kLength * this.getEdgeLengthEnergy(vertices, edge);
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
				VecD edgeLengthEnergyGrad = this.kLength * this.getEdgeLengthEnergyGradient(vertices, edge);
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

	private MatD getSystemEnergyHessian(Vec3D[] vertexPositions, int[] triangles) {
		
		// Assemble system-wide energy Hessian.
		int numVertices = vertexPositions.Length;
		MatD energyHess = new MatD(numVertices * 3, numVertices * 3);
		foreach(Edge edge in this.edges) {

			// The length energy Hessian consists of 4 (3x3) parts that have to be inserted into the matrix.
			if(this.kLength != 0d) {
				MatD lengthEnergyHess = this.getEdgeLengthEnergyHess(vertexPositions, edge).mul(this.kLength);
				makeHessPositiveDefinite(lengthEnergyHess);
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddLengthEnergy_dv1_dv1.
						energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] += lengthEnergyHess[i, j];

						// ddLengthEnergy_dv2_dv2.
						energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] += lengthEnergyHess[i + 3, j + 3];

						// ddLengthEnergy_dv1_dv2.
						energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] += lengthEnergyHess[i, j + 3];

						// ddLengthEnergy_dv2_dv1.
						energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] += lengthEnergyHess[i + 3, j];
					}
				}
			}

			// The bending energy Hessian consists of 16 (3x3) parts that have to be inserted into the matrix.
			if(this.kBend != 0d && edge.hasSideFlaps()) {
				MatD bendEnergyHess = this.getEdgeBendEnergyHess(vertexPositions, edge).mul(this.kBend);
				makeHessPositiveDefinite(bendEnergyHess);
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {
						energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i, j];
						energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i, j + 3];
						energyHess[edge.ve1 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i, j + 6];
						energyHess[edge.ve1 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i, j + 9];
						
						energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 3, j];
						energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 3, j + 3];
						energyHess[edge.ve2 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 3, j + 6];
						energyHess[edge.ve2 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 3, j + 9];
						
						energyHess[edge.vf1 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 6, j];
						energyHess[edge.vf1 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 6, j + 3];
						energyHess[edge.vf1 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 6, j + 6];
						energyHess[edge.vf1 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 6, j + 9];
						
						energyHess[edge.vf2 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 9, j];
						energyHess[edge.vf2 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 9, j + 3];
						energyHess[edge.vf2 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 9, j + 6];
						energyHess[edge.vf2 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 9, j + 9];
					}
				}
			}
		}
		if(this.kArea != 0d) {
			for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId += 3) {
				int v1 = triangles[triangleId];
				int v2 = triangles[triangleId + 1];
				int v3 = triangles[triangleId + 2];
				MatD areaEnergyHess = this.getTriangleAreaEnergyHessian(vertexPositions, triangleId, v1, v2, v3).mul(this.kArea);
				makeHessPositiveDefinite(areaEnergyHess);
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddAreaEnergy_dvi_dvi.
						energyHess[v1 * 3 + i, v1 * 3 + j] += areaEnergyHess[i, j];
						energyHess[v2 * 3 + i, v2 * 3 + j] += areaEnergyHess[i + 3, j + 3];
						energyHess[v3 * 3 + i, v3 * 3 + j] += areaEnergyHess[i + 6, j + 6];

						// ddAreaEnergy_dv1_dv2 & ddAreaEnergy_dv2_dv1.
						energyHess[v1 * 3 + i, v2 * 3 + j] += areaEnergyHess[i, j + 3];
						energyHess[v2 * 3 + i, v1 * 3 + j] += areaEnergyHess[i + 3, j];

						// ddAreaEnergy_dv1_dv3 & ddAreaEnergy_dv3_dv1.
						energyHess[v1 * 3 + i, v3 * 3 + j] += areaEnergyHess[i, j + 6];
						energyHess[v3 * 3 + i, v1 * 3 + j] += areaEnergyHess[i + 6, j];

						// ddAreaEnergy_dv2_dv3 & ddAreaEnergy_dv3_dv2.
						energyHess[v2 * 3 + i, v3 * 3 + j] += areaEnergyHess[i + 3, j + 6];
						energyHess[v3 * 3 + i, v2 * 3 + j] += areaEnergyHess[i + 6, j + 3];
					}
				}
			}
		}
		return energyHess;
	}

	private MathNet.Numerics.LinearAlgebra.Double.SparseMatrix getSystemEnergyHessianSparseRepresentation(Vec3D[] vertexPositions, int[] triangles) {
		
		// Assemble system-wide energy Hessian.
		int numVertices = vertexPositions.Length;
		MathNet.Numerics.LinearAlgebra.Double.SparseMatrix energyHess = new MathNet.Numerics.LinearAlgebra.Double.SparseMatrix(numVertices * 3, numVertices * 3);
		foreach(Edge edge in this.edges) {

			// The length energy Hessian consists of 4 (3x3) parts that have to be inserted into the matrix.
			if(this.kLength != 0d) {
				MatD lengthEnergyHess = this.getEdgeLengthEnergyHess(vertexPositions, edge).mul(this.kLength);
				makeHessPositiveDefinite(lengthEnergyHess);
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddLengthEnergy_dv1_dv1.
						energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] += lengthEnergyHess[i, j];

						// ddLengthEnergy_dv2_dv2.
						energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] += lengthEnergyHess[i + 3, j + 3];

						// ddLengthEnergy_dv1_dv2.
						energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] += lengthEnergyHess[i, j + 3];

						// ddLengthEnergy_dv2_dv1.
						energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] += lengthEnergyHess[i + 3, j];
					}
				}
			}

			// The bending energy Hessian consists of 16 (3x3) parts that have to be inserted into the matrix.
			if(this.kBend != 0d && edge.hasSideFlaps()) {
				MatD bendEnergyHess = this.getEdgeBendEnergyHess(vertexPositions, edge).mul(this.kBend);
				makeHessPositiveDefinite(bendEnergyHess);
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {
						energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i, j];
						energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i, j + 3];
						energyHess[edge.ve1 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i, j + 6];
						energyHess[edge.ve1 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i, j + 9];
						
						energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 3, j];
						energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 3, j + 3];
						energyHess[edge.ve2 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 3, j + 6];
						energyHess[edge.ve2 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 3, j + 9];
						
						energyHess[edge.vf1 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 6, j];
						energyHess[edge.vf1 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 6, j + 3];
						energyHess[edge.vf1 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 6, j + 6];
						energyHess[edge.vf1 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 6, j + 9];
						
						energyHess[edge.vf2 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 9, j];
						energyHess[edge.vf2 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 9, j + 3];
						energyHess[edge.vf2 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 9, j + 6];
						energyHess[edge.vf2 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 9, j + 9];
					}
				}
			}
		}
		if(this.kArea != 0d) {
			for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
				int triangleBaseIndex = triangleId * 3;
				int v1 = triangles[triangleBaseIndex];
				int v2 = triangles[triangleBaseIndex + 1];
				int v3 = triangles[triangleBaseIndex + 2];
				MatD areaEnergyHess = this.getTriangleAreaEnergyHessian(vertexPositions, triangleId, v1, v2, v3).mul(this.kArea);
				makeHessPositiveDefinite(areaEnergyHess);
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddAreaEnergy_dvi_dvi.
						energyHess[v1 * 3 + i, v1 * 3 + j] += areaEnergyHess[i, j];
						energyHess[v2 * 3 + i, v2 * 3 + j] += areaEnergyHess[i + 3, j + 3];
						energyHess[v3 * 3 + i, v3 * 3 + j] += areaEnergyHess[i + 6, j + 6];

						// ddAreaEnergy_dv1_dv2 & ddAreaEnergy_dv2_dv1.
						energyHess[v1 * 3 + i, v2 * 3 + j] += areaEnergyHess[i, j + 3];
						energyHess[v2 * 3 + i, v1 * 3 + j] += areaEnergyHess[i + 3, j];

						// ddAreaEnergy_dv1_dv3 & ddAreaEnergy_dv3_dv1.
						energyHess[v1 * 3 + i, v3 * 3 + j] += areaEnergyHess[i, j + 6];
						energyHess[v3 * 3 + i, v1 * 3 + j] += areaEnergyHess[i + 6, j];

						// ddAreaEnergy_dv2_dv3 & ddAreaEnergy_dv3_dv2.
						energyHess[v2 * 3 + i, v3 * 3 + j] += areaEnergyHess[i + 3, j + 6];
						energyHess[v3 * 3 + i, v2 * 3 + j] += areaEnergyHess[i + 6, j + 3];
					}
				}
			}
		}
		return energyHess;
	}

	private MathNet.Numerics.LinearAlgebra.Double.SparseMatrix getSystemEnergyHessianSparseRepresentationMultiThreaded(Vec3D[] vertexPositions, int[] triangles) {
		
		// Assemble system-wide energy Hessian.
		int numVertices = vertexPositions.Length;
		MathNet.Numerics.LinearAlgebra.Double.SparseMatrix energyHess = new MathNet.Numerics.LinearAlgebra.Double.SparseMatrix(numVertices * 3, numVertices * 3);
		
		// The length energy Hessian consists of 4 (3x3) parts that have to be inserted into the matrix.
		if(this.kLength != 0d) {
			MatD[] lengthEnergyHessians = new MatD[this.edges.Count];
			Parallel.For(0, this.edges.Count, (i) => {
				lengthEnergyHessians[i] = this.getEdgeLengthEnergyHess(vertexPositions, this.edges[i]).mul(this.kLength);
				makeHessPositiveDefinite(lengthEnergyHessians[i]);
			});

			for(int edgeInd = 0; edgeInd < this.edges.Count; edgeInd++) {
				Edge edge = this.edges[edgeInd];
				MatD lengthEnergyHess = lengthEnergyHessians[edgeInd];
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddLengthEnergy_dv1_dv1.
						energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] += lengthEnergyHess[i, j];

						// ddLengthEnergy_dv2_dv2.
						energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] += lengthEnergyHess[i + 3, j + 3];

						// ddLengthEnergy_dv1_dv2.
						energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] += lengthEnergyHess[i, j + 3];

						// ddLengthEnergy_dv2_dv1.
						energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] += lengthEnergyHess[i + 3, j];
					}
				}
			}
		}
		
		// The bending energy Hessian consists of 16 (3x3) parts that have to be inserted into the matrix.
		if(this.kBend != 0d) {
			MatD[] bendEnergyHessians = new MatD[this.edges.Count];
			Parallel.For(0, this.edges.Count, (i) => {
				if(this.edges[i].hasSideFlaps()) {
					bendEnergyHessians[i] = this.getEdgeBendEnergyHess(vertexPositions, this.edges[i]).mul(this.kBend);
					makeHessPositiveDefinite(bendEnergyHessians[i]);
				} else {
					bendEnergyHessians[i] = null;
				}
			});

			for(int edgeInd = 0; edgeInd < this.edges.Count; edgeInd++) {
				MatD bendEnergyHess = bendEnergyHessians[edgeInd];
				if(bendEnergyHess == null) {
					continue;
				}
				Edge edge = this.edges[edgeInd];
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {
						energyHess[edge.ve1 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i, j];
						energyHess[edge.ve1 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i, j + 3];
						energyHess[edge.ve1 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i, j + 6];
						energyHess[edge.ve1 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i, j + 9];
						
						energyHess[edge.ve2 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 3, j];
						energyHess[edge.ve2 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 3, j + 3];
						energyHess[edge.ve2 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 3, j + 6];
						energyHess[edge.ve2 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 3, j + 9];
						
						energyHess[edge.vf1 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 6, j];
						energyHess[edge.vf1 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 6, j + 3];
						energyHess[edge.vf1 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 6, j + 6];
						energyHess[edge.vf1 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 6, j + 9];
						
						energyHess[edge.vf2 * 3 + i, edge.ve1 * 3 + j] += bendEnergyHess[i + 9, j];
						energyHess[edge.vf2 * 3 + i, edge.ve2 * 3 + j] += bendEnergyHess[i + 9, j + 3];
						energyHess[edge.vf2 * 3 + i, edge.vf1 * 3 + j] += bendEnergyHess[i + 9, j + 6];
						energyHess[edge.vf2 * 3 + i, edge.vf2 * 3 + j] += bendEnergyHess[i + 9, j + 9];
					}
				}
			}
		}
		
		// The area energy Hessian consists of 9 (3x3) parts that have to be inserted into the matrix.
		if(this.kArea != 0d) {
			int numTriangles = triangles.Length / 3;
			MatD[] areaEnergyHessians = new MatD[numTriangles];
			Parallel.For(0, numTriangles, (triangleId) => {
				int triangleBaseIndex = 3 * triangleId;
				int v1 = triangles[triangleBaseIndex];
				int v2 = triangles[triangleBaseIndex + 1];
				int v3 = triangles[triangleBaseIndex + 2];
				areaEnergyHessians[triangleId] = this.getTriangleAreaEnergyHessian(vertexPositions, triangleId, v1, v2, v3).mul(this.kArea);
				makeHessPositiveDefinite(areaEnergyHessians[triangleId]);
			});

			for(int triangleId = 0; triangleId < numTriangles; triangleId++) {
				int trianglebaseIndex = 3 * triangleId;
				int v1 = triangles[trianglebaseIndex];
				int v2 = triangles[trianglebaseIndex + 1];
				int v3 = triangles[trianglebaseIndex + 2];
				MatD areaEnergyHess = areaEnergyHessians[triangleId];
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddAreaEnergy_dvi_dvi.
						energyHess[v1 * 3 + i, v1 * 3 + j] += areaEnergyHess[i, j];
						energyHess[v2 * 3 + i, v2 * 3 + j] += areaEnergyHess[i + 3, j + 3];
						energyHess[v3 * 3 + i, v3 * 3 + j] += areaEnergyHess[i + 6, j + 6];

						// ddAreaEnergy_dv1_dv2 & ddAreaEnergy_dv2_dv1.
						energyHess[v1 * 3 + i, v2 * 3 + j] += areaEnergyHess[i, j + 3];
						energyHess[v2 * 3 + i, v1 * 3 + j] += areaEnergyHess[i + 3, j];

						// ddAreaEnergy_dv1_dv3 & ddAreaEnergy_dv3_dv1.
						energyHess[v1 * 3 + i, v3 * 3 + j] += areaEnergyHess[i, j + 6];
						energyHess[v3 * 3 + i, v1 * 3 + j] += areaEnergyHess[i + 6, j];

						// ddAreaEnergy_dv2_dv3 & ddAreaEnergy_dv3_dv2.
						energyHess[v2 * 3 + i, v3 * 3 + j] += areaEnergyHess[i + 3, j + 6];
						energyHess[v3 * 3 + i, v2 * 3 + j] += areaEnergyHess[i + 6, j + 3];
					}
				}
			}
		}
		return energyHess;
	}

	private MathNet.Numerics.LinearAlgebra.Double.SparseMatrix getSystemEnergyHessianSparseRepresentationMultiThreadedTriplets(Vec3D[] vertexPositions, int[] triangles) {
		Stopwatch stopwatch = Stopwatch.StartNew();

		// Define constants and allocate memory for all Hessian part triplets.
		int numVertices = vertexPositions.Length;
		int numEdges = this.edges.Count;
		int numTriangles = triangles.Length / 3;

		int numLengthHessParts = (this.kLength != 0d ? numEdges : 0);
		int numAreaHessParts = (this.kArea != 0d ? numTriangles : 0);
		int numBendHessParts = (this.kBend != 0d ? numEdges : 0); // TODO - Get number of inner edges instead, and cache an array of inner edges for iteration.
		
		int numTriplets = numLengthHessParts * 6 * 6 + numAreaHessParts * 9 * 9 + numBendHessParts * 12 * 12;
		int lengthHessPartsOffset = 0;
		int areaHessPartsOffset = lengthHessPartsOffset + numLengthHessParts * 6 * 6;
		int bendHessPartsOffset = areaHessPartsOffset + numAreaHessParts * 9 * 9;
		Tuple<int, int, double>[] energyHessTriplets = new Tuple<int, int, double>[numTriplets];

		// The length energy Hessian consists of 4 (3x3) parts that have to be inserted into the matrix.
		long startTime = stopwatch.ElapsedMilliseconds;
		if(this.kLength != 0d) {
			Parallel.For(0, numEdges, (edgeInd) => {
				Edge edge = this.edges[edgeInd];
				MatD lengthEnergyHess = this.getEdgeLengthEnergyHess(vertexPositions, edge).mul(this.kLength);
				makeHessPositiveDefinite(lengthEnergyHess);
				
				int tripletsOffset = lengthHessPartsOffset + 6 * 6 * edgeInd;
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddLengthEnergy_dv1_dv1.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve1 * 3 + i, edge.ve1 * 3 + j, lengthEnergyHess[i, j]);

						// ddLengthEnergy_dv2_dv2.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve2 * 3 + i, edge.ve2 * 3 + j, lengthEnergyHess[i + 3, j + 3]);

						// ddLengthEnergy_dv1_dv2.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve1 * 3 + i, edge.ve2 * 3 + j, lengthEnergyHess[i, j + 3]);

						// ddLengthEnergy_dv2_dv1.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve2 * 3 + i, edge.ve1 * 3 + j, lengthEnergyHess[i + 3, j]);
					}
				}
			});
		}
		print(stopwatch.ElapsedMilliseconds + ": Length hess done: " + (stopwatch.ElapsedMilliseconds - startTime) + "ms.");

		// The bending energy Hessian consists of 16 (3x3) parts that have to be inserted into the matrix.
		startTime = stopwatch.ElapsedMilliseconds;
		if(this.kBend != 0d) {
			Parallel.For(0, numEdges, (edgeInd) => {
				Edge edge = this.edges[edgeInd];
				if(!edge.hasSideFlaps()) {
					return; // TODO - Instead, cache all inner edges so that we don't have to allocate memory for outer edges without bending Hessian.
				}
				MatD bendEnergyHess = this.getEdgeBendEnergyHess(vertexPositions, edge).mul(this.kBend);
				makeHessPositiveDefinite(bendEnergyHess);

				int tripletsOffset = bendHessPartsOffset + 12 * 12 * edgeInd;
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve1 * 3 + i, edge.ve1 * 3 + j, bendEnergyHess[i, j]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve1 * 3 + i, edge.ve2 * 3 + j, bendEnergyHess[i, j + 3]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve1 * 3 + i, edge.vf1 * 3 + j, bendEnergyHess[i, j + 6]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve1 * 3 + i, edge.vf2 * 3 + j, bendEnergyHess[i, j + 9]);
						
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve2 * 3 + i, edge.ve1 * 3 + j, bendEnergyHess[i + 3, j]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve2 * 3 + i, edge.ve2 * 3 + j, bendEnergyHess[i + 3, j + 3]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve2 * 3 + i, edge.vf1 * 3 + j, bendEnergyHess[i + 3, j + 6]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.ve2 * 3 + i, edge.vf2 * 3 + j, bendEnergyHess[i + 3, j + 9]);
						
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf1 * 3 + i, edge.ve1 * 3 + j, bendEnergyHess[i + 6, j]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf1 * 3 + i, edge.ve2 * 3 + j, bendEnergyHess[i + 6, j + 3]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf1 * 3 + i, edge.vf1 * 3 + j, bendEnergyHess[i + 6, j + 6]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf1 * 3 + i, edge.vf2 * 3 + j, bendEnergyHess[i + 6, j + 9]);
						
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf2 * 3 + i, edge.ve1 * 3 + j, bendEnergyHess[i + 9, j]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf2 * 3 + i, edge.ve2 * 3 + j, bendEnergyHess[i + 9, j + 3]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf2 * 3 + i, edge.vf1 * 3 + j, bendEnergyHess[i + 9, j + 6]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(edge.vf2 * 3 + i, edge.vf2 * 3 + j, bendEnergyHess[i + 9, j + 9]);
					}
				}
			});
		}
		print(stopwatch.ElapsedMilliseconds + ": Bend hess done: " + (stopwatch.ElapsedMilliseconds - startTime) + "ms.");
		
		// The area energy Hessian consists of 9 (3x3) parts that have to be inserted into the matrix.
		startTime = stopwatch.ElapsedMilliseconds;
		if(this.kArea != 0d) {
			Parallel.For(0, numTriangles, (triangleId) => {
				int trianglebaseIndex = 3 * triangleId;
				int v1 = triangles[trianglebaseIndex];
				int v2 = triangles[trianglebaseIndex + 1];
				int v3 = triangles[trianglebaseIndex + 2];
				MatD areaEnergyHess = this.getTriangleAreaEnergyHessian(vertexPositions, triangleId, v1, v2, v3).mul(this.kArea);
				makeHessPositiveDefinite(areaEnergyHess);
				
				int tripletsOffset = areaHessPartsOffset + 9 * 9 * triangleId;
				for(int i = 0; i < 3; i++) {
					for(int j = 0; j < 3; j++) {

						// ddAreaEnergy_dvi_dvi.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v1 * 3 + i, v1 * 3 + j, areaEnergyHess[i, j]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v2 * 3 + i, v2 * 3 + j, areaEnergyHess[i + 3, j + 3]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v3 * 3 + i, v3 * 3 + j, areaEnergyHess[i + 6, j + 6]);

						// ddAreaEnergy_dv1_dv2 & ddAreaEnergy_dv2_dv1.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v1 * 3 + i, v2 * 3 + j, areaEnergyHess[i, j + 3]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v2 * 3 + i, v1 * 3 + j, areaEnergyHess[i + 3, j]);

						// ddAreaEnergy_dv1_dv3 & ddAreaEnergy_dv3_dv1.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v1 * 3 + i, v3 * 3 + j, areaEnergyHess[i, j + 6]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v3 * 3 + i, v1 * 3 + j, areaEnergyHess[i + 6, j]);

						// ddAreaEnergy_dv2_dv3 & ddAreaEnergy_dv3_dv2.
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v2 * 3 + i, v3 * 3 + j, areaEnergyHess[i + 3, j + 6]);
						energyHessTriplets[tripletsOffset++] = new Tuple<int, int, double>(v3 * 3 + i, v2 * 3 + j, areaEnergyHess[i + 6, j + 3]);
					}
				}
			});
		}
		print(stopwatch.ElapsedMilliseconds + ": Area hess done: " + (stopwatch.ElapsedMilliseconds - startTime) + "ms.");

		// Sort triplets first on row and then on column.
		startTime = stopwatch.ElapsedMilliseconds;
		Array.Sort(energyHessTriplets, (o1, o2) =>
				(o1 == null ? (o2 == null ? 0 : -1) : o2 == null ? 1
						: (o1.Item1 > o2.Item1 ? 1 : (o1.Item1 < o2.Item1 ? -1
								: (o1.Item2 > o2.Item2 ? 1 : (o1.Item2 < o2.Item2 ? -1 : 0))))));
		print(stopwatch.ElapsedMilliseconds + ": Sorting triplets done: " + (stopwatch.ElapsedMilliseconds - startTime) + "ms.");

		// Sum up triplets at the same index and move them to the left, resulting in trailing nulls.
		startTime = stopwatch.ElapsedMilliseconds;
		Tuple<int, int, double> lastTriplet = null;
		int lastTripletInd = -1;
		for(int i = 0; i < energyHessTriplets.Length; i++) {
			Tuple<int, int, double> triplet = energyHessTriplets[i];
			if(triplet == null) {
				continue;
			} else if(triplet.Item3 == 0d) {
				energyHessTriplets[i] = null;
			} else if(lastTriplet == null) {

				// Found first triplet. Put it in the first index of the triplets array.
				energyHessTriplets[i] = null;
				energyHessTriplets[0] = triplet;
				lastTriplet = triplet;
				lastTripletInd = 0;
			} else if(triplet.Item1 == lastTriplet.Item1 && triplet.Item2 == lastTriplet.Item2) {

				// Merge triplet into lastTriplet.
				energyHessTriplets[lastTripletInd] = new Tuple<int, int, double>(lastTriplet.Item1, lastTriplet.Item2, lastTriplet.Item3 + triplet.Item3);
				lastTriplet = energyHessTriplets[lastTripletInd];
				energyHessTriplets[i] = null;
			} else {

				// Move the non-null triplet to the leftmost free location in the array.
				energyHessTriplets[i] = null;
				energyHessTriplets[++lastTripletInd] = triplet;
				lastTriplet = triplet;
			}
		}
		print(stopwatch.ElapsedMilliseconds + ": Merging triplets done: " + (stopwatch.ElapsedMilliseconds - startTime) + "ms.");

		// Assemble and return system-wide energy Hessian.
		startTime = stopwatch.ElapsedMilliseconds;
		MathNet.Numerics.LinearAlgebra.Double.SparseMatrix energyHess = new MathNet.Numerics.LinearAlgebra.Double.SparseMatrix(numVertices * 3, numVertices * 3);
		foreach(Tuple<int, int, double> triplet in energyHessTriplets) {
			if(triplet == null) {
				break; // All nulls are trailing, so there won't be more data.
			}
			energyHess[triplet.Item1, triplet.Item2] = triplet.Item3;
		}
		print(stopwatch.ElapsedMilliseconds + ": Creating SparseMatrix done: " + (stopwatch.ElapsedMilliseconds - startTime) + "ms.");
		return energyHess;
	}

	private double getE(int[] triangles, VecD vertexPositionsFlat, Vec3D[] newVertexPositions, VecD vertexCoordMasses,
			VecD windVelocity, MathNet.Numerics.LinearAlgebra.Double.SparseMatrix energyHess, double deltaTime) {
		int numVertices = newVertexPositions.Length;

		VecD windForce = this.getVertexWindForce(triangles, newVertexPositions, windVelocity);
		for(int i = 0; i < numVertices; i++) {
			if(this.verticesMovementConstraints[i]) {
				for(int coord = 0; coord < 3; coord++) {
					windForce[3 * i + coord] = 0;
				}
			}
		}

		VecD deltaVertexPositions = new VecD(newVertexPositions).sub(vertexPositionsFlat); // These are also the new vertex velocities.
		double systemEnergy = this.getSystemEnergy(triangles, newVertexPositions);
		double gravityWork = 0;
		double windWork = -VecD.dot(windForce, deltaVertexPositions); // Approximation: Consider triangle normals and areas constant.
		VecD dampingForce = -this.dampingConstant * (energyHess * deltaVertexPositions);
		double dampingWork = -VecD.dot(dampingForce, deltaVertexPositions);
		for(int i = 0; i < numVertices; i++) {
			if(!this.verticesMovementConstraints[i]) { // Not necessary when comparing energy, but lets consider them part of the outside world.
				int yCoordIndex = 3 * i + 1;
				gravityWork += vertexCoordMasses[yCoordIndex] * this.gravityConstant * deltaVertexPositions[yCoordIndex];
			}
		}
		VecD squareTerm = new VecD(newVertexPositions).sub(vertexPositionsFlat).sub(deltaTime * this.vertexVelocities);
		return VecD.dot(VecD.multiplyElementWise(squareTerm, squareTerm), vertexCoordMasses) / (2d * deltaTime * deltaTime)
				+ systemEnergy + gravityWork + windWork + dampingWork;
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
	private VecD getVertexWindForce(int[] triangles, VecD[] vertices, VecD windPressure) {

		// Initialize vertex wind force array.
		VecD vertexWindForce = new VecD(vertices.Length * 3);

		// Compute vertex wind force.
		for(int triangleId = 0; triangleId < triangles.Length / 3; triangleId++) {
			int triangleBaseIndex = triangleId * 3;
			int v1 = triangles[triangleBaseIndex];
			int v2 = triangles[triangleBaseIndex + 1];
			int v3 = triangles[triangleBaseIndex + 2];

			// Compute a third of the projection of the wind pressure vector on the triangle normal.
			// The factor of a third is caused by the triangle wind force being divided over the 3 vertices of that triangle.
			VecD triangleNormal = this.triangleNormals[triangleId];
			if(triangleNormal == null) {
				continue; // Triangle has a zero-area and no normal, so the projected wind force is zero as well.
			}
			double triangleArea = this.triangleAreas[triangleId];
			VecD triangleVertexWindForce = (VecD.dot(windPressure, triangleNormal) * triangleArea / 3d) * triangleNormal;

			// Add a third of the total triangle wind force to each of its vertices.
			for(int coord = 0; coord < 3; coord++) {
				vertexWindForce[3 * v1 + coord] += triangleVertexWindForce[coord];
				vertexWindForce[3 * v2 + coord] += triangleVertexWindForce[coord];
				vertexWindForce[3 * v3 + coord] += triangleVertexWindForce[coord];
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

	private double getEdgeLengthEnergy(Vec3D[] vertexPositions, Edge edge) {
		Vec3D e = vertexPositions[edge.ve2] - vertexPositions[edge.ve1]; // Vector from v1 to v2.
		double edgeLength = e.magnitude;
		double a = 1d - edgeLength / edge.undeformedLength;
		return a * a * edge.undeformedLength;
	}

	/**
	 * Computes the edge length energy gradient of the edge between vertices v1 and v2, towards {v1x, v1y, v1z, v2x, v2y, v2z}.
	 * The result is a vector of length 6.
	 */
	private VecD getEdgeLengthEnergyGradient(Vec3D[] vertexPositions, Edge edge) {
		int v1 = edge.ve1;
		int v2 = edge.ve2;
		Vec3D e = vertexPositions[v2] - vertexPositions[v1]; // Vector from v1 to v2.
		double edgeLength = e.magnitude;
		if(double.IsNaN(edgeLength)) {
			return new VecD(0, 0, 0, 0, 0, 0); // Edge is zero-length, so the gradient is 0.
		}
		Vec3D dEdgeLength_dv1 = (vertexPositions[v1] - vertexPositions[v2]).div(edgeLength);

		/*
		 * Perform:
		 * Vec3D dEdgeLength_dv2 = -dEdgeLength_dv1;
		 * VecD dEdgeLength_dv1v2 = new VecD(dEdgeLength_dv1, dEdgeLength_dv2); // Partial derivative towards {v1x, v1y, v1z, v2x, v2y, v2z}.
		 */
		VecD dEdgeLength_dv1v2 = new VecD(6);
		for(int i = 0; i < 3; i++) {
			dEdgeLength_dv1v2[i] = dEdgeLength_dv1[i];
			dEdgeLength_dv1v2[i + 3] = -dEdgeLength_dv1[i];
		}

		return dEdgeLength_dv1v2.mul(2 * edgeLength / edge.undeformedLength - 2);
	}

	/**
	 * Computes the edge length energy Hessian of the edge between vertices v1 and v2, towards {v1x, v1y, v1z, v2x, v2y, v2z}.
	 * The result is a 6x6 matrix containing all combinations of double partial derivatives towards {v1x, v1y, v1z, v2x, v2y, v2z}.
	 */
	private MatD getEdgeLengthEnergyHess(Vec3D[] vertices, Edge edge) {
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
		int v1 = edge.ve1;
		int v2 = edge.ve2;
		Vec3D e = vertexPositions[v2] - vertexPositions[v1]; // Vector from v1 to v2.
		double edgeLength = e.magnitude;
		if(double.IsNaN(edgeLength)) {
			return new MatD(6, 6);
		}
		Vec3D dEdgeLength_dv1 = (vertexPositions[v1] - vertexPositions[v2]).div(edgeLength);

		/*
		 * Perform:
		 * Vec3D dEdgeLength_dv2 = -dEdgeLength_dv1;
		 * VecD dEdgeLength_dv1v2 = new VecD(dEdgeLength_dv1, dEdgeLength_dv2); // Partial derivative towards {v1x, v1y, v1z, v2x, v2y, v2z}.
		 */
		VecD dEdgeLength_dv1v2 = new VecD(6);
		for(int i = 0; i < 3; i++) {
			dEdgeLength_dv1v2[i] = dEdgeLength_dv1[i];
			dEdgeLength_dv1v2[i + 3] = -dEdgeLength_dv1[i];
		}

		//VecD dEdgeEnergy_dv1v2 = (2 * edgeLength / edge.undeformedLength - 2) * dEdgeLength_dv1v2;
		// TODO - Copied gradient code ends here (See TODO above).

		// Calculate edge length Hessian.
		double edgeLengthSquare = edgeLength * edgeLength;
		double edgeLengthCube = edgeLengthSquare * edgeLength;
		MatD ddEdgeLength_dv1_dv1 = new MatD(new double[,] {
			{edgeLengthSquare - e[0] * e[0],                  - e[1] * e[0],                  - e[2] * e[0]},
			{                 - e[0] * e[1], edgeLengthSquare - e[1] * e[1],                  - e[2] * e[1]},
			{                 - e[0] * e[2],                  - e[1] * e[2], edgeLengthSquare - e[2] * e[2]}
		}).div(edgeLengthCube);
		
		MatD ddEdgeLength_dv1_dv2 = ddEdgeLength_dv1_dv1.Clone();
		ddEdgeLength_dv1_dv2[0, 0] = -ddEdgeLength_dv1_dv2[0, 0];
		ddEdgeLength_dv1_dv2[1, 1] = -ddEdgeLength_dv1_dv2[1, 1];
		ddEdgeLength_dv1_dv2[2, 2] = -ddEdgeLength_dv1_dv2[2, 2];
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
		// This calculation uses the fact that the length energy Hessian is symmetric.
		MatD lengthEnergyHess = new MatD(6, 6);
		double a = (2 * edgeLength / edge.undeformedLength - 2);
		for(int i = 0; i < 6; i++) {
			lengthEnergyHess[i, i] = 2d / edge.undeformedLength * dEdgeLength_dv1v2[i] * dEdgeLength_dv1v2[i] + a * edgeLengthHess[i, i];
			for(int j = i + 1; j < 6; j++) {
				lengthEnergyHess[i, j] = lengthEnergyHess[j, i] = 2d / edge.undeformedLength * dEdgeLength_dv1v2[i] * dEdgeLength_dv1v2[j] + a * edgeLengthHess[i, j];
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

		/*
		 * Perform:
		 * VecD d_crossProdLength_dv1v2v3 = (d * a + e * b + f * c) / crossProdLength;
		 * VecD d_triangleArea_dv1v2v3 = d_crossProdLength_dv1v2v3 / 2d;
		 * VecD d_triangleEnergy_dv1v2v3 = (2d * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2d) * d_triangleArea_dv1v2v3;
		 */
		return new VecD(a.mul(d)).add(b.mul(e)).add(c.mul(f)).div(crossProdLength * 2d)
				.mul(2d * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2d);
	}

	// Area Hessian constants.
	private static readonly MatD d_a_dv1v2v3 = new MatD(new double[,] {
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
	private static readonly MatD d_b_dv1v2v3 = new MatD(new double[,] {
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
	private static readonly MatD d_c_dv1v2v3 = new MatD(new double[,] {
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
	/*
	 * Computes the triangle area energy Hessian for the triangle defined by vertices v1, v2 and v3.
	 * Returns the 9x9 Hessian towards all combinations of {v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z}.
	 */
	private MatD getTriangleAreaEnergyHessian(Vec3D[] vertices, int triangleId, int v1Ind, int v2Ind, int v3Ind) {
		
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
		VecD a = new VecD(0		  , v3.z - v2.z, v2.y - v3.y, 0                        , v2.z - v3.z + v1.z - v2.z, v2.y - v1.y + v3.y - v2.y, 0		  , v2.z - v1.z, v1.y - v2.y);
		VecD b = new VecD(v2.z - v3.z, 0		  , v3.x - v2.x, v2.z - v1.z + v3.z - v2.z, 0                        , v2.x - v3.x + v1.x - v2.x, v1.z - v2.z, 0		  , v2.x - v1.x);
		VecD c = new VecD(v3.y - v2.y, v2.x - v3.x, 0		  , v2.y - v3.y + v1.y - v2.y, v2.x - v1.x + v3.x - v2.x, 0                        , v2.y - v1.y, v1.x - v2.x, 0		  );
		double d = ((v1.y - v2.y) * (v3.z - v2.z) - (v3.y - v2.y) * (v1.z - v2.z));
		double e = ((v2.x - v1.x) * (v3.z - v2.z) + (v3.x - v2.x) * (v1.z - v2.z));
		double f = ((v1.x - v2.x) * (v3.y - v2.y) - (v3.x - v2.x) * (v1.y - v2.y));

		// Perform: VecD d_crossProdLength_dv1v2v3 = (d * a + e * b + f * c).div(crossProdLength);
		VecD d_crossProdLength_dv1v2v3 = (d * a).add(e * b).add(f * c).div(crossProdLength);
		// TODO - Copied gradient code ends here (See TODO above).

		// Compute triangle energy Hessian.
		// hess[i, j] = (crossProdLength * (a[i] * a[j] + d_a[j]_di * d + b[i] * b[j] + d_b[j]_di * e + c[i] * c[j] + d_c[j]_di * f) - (d * a[j] + e * b[j] + f * c[j]) * d_crossProdLength_di) / crossProdLength^2 / 2
		MatD areaEnergyHessian = new MatD(9, 9);
		for(int i = 0; i < areaEnergyHessian.numRows; i++) {
			for(int j = i; j < areaEnergyHessian.numColumns; j++) { // Use knowledge that the Hessian is symmetric.
				areaEnergyHessian[i, j] = areaEnergyHessian[j, i] = (
						(a[i] * a[j] + d_a_dv1v2v3[j, i] * d + b[i] * b[j] + d_b_dv1v2v3[j, i] * e + c[i] * c[j] + d_c_dv1v2v3[j, i] * f)
						- (d * a[j] + e * b[j] + f * c[j]) * d_crossProdLength_dv1v2v3[i] / crossProdLength
						) / crossProdLength / 2d;
			}
		}
		return areaEnergyHessian;
	}

	private double getEdgeBendEnergy(Vec3D[] vertexPositions, Edge edge) {

		// Define required constants.
		// h_e_undeformed is a third of the average triangle height, where the height is twice the triangle area divided by the triangle width.
		double h_e_undeformed = (this.undeformedTriangleAreas[edge.triangleId1] + this.undeformedTriangleAreas[edge.triangleId2]) / edge.undeformedLength / 3d;
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
		return a * a * edge.undeformedLength / h_e_undeformed;
	}

	/**
	 * Computes the bending energy gradient of the given edge. This gradient is taken towards the edge-defining vertices ve1 and ve2, as well
	 * as the vertices vf1 and vf2, completing triangles 1 and 2 connected to the edge respectively.
	 * The return values are the bending energy gradient towards: {ve1x. ve1y, ve1z, ve2x. ve2y, ve2z, vf1x. vf1y, vf1z, vf2x. vf2y, vf2z}.
	 */
	private VecD getEdgeBendEnergyGradient(Vec3D[] vertices, Edge edge) {

		// Define required constants.
		// h_e_undeformed is a third of the average triangle height, where the height is twice the triangle area divided by the triangle width.
		double h_e_undeformed = (this.undeformedTriangleAreas[edge.triangleId1] + this.undeformedTriangleAreas[edge.triangleId2]) / edge.undeformedLength / 3d;
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
		double d_fi_d_teta = 2d * (teta_e - teta_e_undeformed) * edge.undeformedLength / h_e_undeformed;

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

		double e0_magnitude = e0.magnitude;
		double t1_e1_magnitude = t1_e1.magnitude;
		double t1_e2_magnitude = t1_e2.magnitude;
		double t2_e1_magnitude = t2_e1.magnitude;
		double t2_e2_magnitude = t2_e2.magnitude;

		float t1_alpha1 = Mathf.Acos((float) (VecD.dot(e0, t1_e2) / (e0.magnitude * t1_e2_magnitude)));
		float t1_alpha2 = Mathf.Acos((float) (VecD.dot(-e0, t1_e1) / (e0.magnitude * t1_e1_magnitude)));
		float t2_alpha1 = Mathf.Acos((float) (VecD.dot(e0, t2_e2) / (e0.magnitude * t2_e2_magnitude)));
		float t2_alpha2 = Mathf.Acos((float) (VecD.dot(-e0, t2_e1) / (e0.magnitude * t2_e1_magnitude)));

		double t1_h0 = 2d * this.triangleAreas[edge.triangleId1] / e0_magnitude;
		double t1_h1 = 2d * this.triangleAreas[edge.triangleId1] / t1_e1_magnitude;
		double t1_h2 = 2d * this.triangleAreas[edge.triangleId1] / t1_e2_magnitude;
		double t2_h0 = 2d * this.triangleAreas[edge.triangleId2] / e0_magnitude;
		double t2_h1 = 2d * this.triangleAreas[edge.triangleId2] / t2_e1_magnitude;
		double t2_h2 = 2d * this.triangleAreas[edge.triangleId2] / t2_e2_magnitude;

		// Calculate derivatives of teta towards all 4 vertices.
		VecD d_teta_d_ve1 = (Mathf.Cos(t1_alpha2) / t1_h1) * n1 + (Mathf.Cos(t2_alpha2) / t2_h1) * n2;
		VecD d_teta_d_ve2 = (Mathf.Cos(t1_alpha1) / t1_h2) * n1 + (Mathf.Cos(t2_alpha1) / t2_h2) * n2;
		VecD d_teta_d_vf1 = (-n1).div(t1_h0); // Perform: -n1 / t1_h0;
		VecD d_teta_d_vf2 = (-n2).div(t2_h0); // Perform: -n2 / t2_h0;
		
		// Assemble hinge angle gradient.
		VecD d_teta_d_ve1ve2vf1vf2 = new VecD(d_teta_d_ve1, d_teta_d_ve2, d_teta_d_vf1, d_teta_d_vf2);

		// Return the bending energy gradient.
		return d_teta_d_ve1ve2vf1vf2.mul(d_fi_d_teta);
	}

	/**
	 * Computes the bending energy Hessian of the given edge. This Hessian is taken towards the edge-defining vertices ve1 and ve2, as well
	 * as the vertices vf1 and vf2, completing triangles 1 and 2 connected to the edge respectively.
	 * The returned Hessian contains partial derivatives in rows and columns in order: {ve1x. ve1y, ve1z, ve2x. ve2y, ve2z, vf1x. vf1y, vf1z, vf2x. vf2y, vf2z}.
	 */
	private MatD getEdgeBendEnergyHess(Vec3D[] vertices, Edge edge) {

		// Define required constants.
		// h_e_undeformed is a third of the average triangle height, where the height is twice the triangle area divided by the triangle width.
		double h_e_undeformed = (this.undeformedTriangleAreas[edge.triangleId1] + this.undeformedTriangleAreas[edge.triangleId2]) / edge.undeformedLength / 3d;
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
		double d_fi_d_teta = 2d * (teta_e - teta_e_undeformed) * edge.undeformedLength / h_e_undeformed;

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

		double e0_magnitude = e0.magnitude;
		double t1_e1_magnitude = t1_e1.magnitude;
		double t1_e2_magnitude = t1_e2.magnitude;
		double t2_e1_magnitude = t2_e1.magnitude;
		double t2_e2_magnitude = t2_e2.magnitude;

		double cos_t1_alpha1 = VecD.dot(e0, t1_e2) / (e0_magnitude * t1_e2_magnitude);
		double cos_t1_alpha2 = VecD.dot(-e0, t1_e1) / (e0_magnitude * t1_e1_magnitude);
		double cos_t2_alpha1 = VecD.dot(e0, t2_e2) / (e0_magnitude * t2_e2_magnitude);
		double cos_t2_alpha2 = VecD.dot(-e0, t2_e1) / (e0_magnitude * t2_e1_magnitude);

		double t1_h0 = 2d * this.triangleAreas[edge.triangleId1] / e0_magnitude;
		double t1_h1 = 2d * this.triangleAreas[edge.triangleId1] / t1_e1_magnitude;
		double t1_h2 = 2d * this.triangleAreas[edge.triangleId1] / t1_e2_magnitude;
		double t2_h0 = 2d * this.triangleAreas[edge.triangleId2] / e0_magnitude;
		double t2_h1 = 2d * this.triangleAreas[edge.triangleId2] / t2_e1_magnitude;
		double t2_h2 = 2d * this.triangleAreas[edge.triangleId2] / t2_e2_magnitude;
		
		double t1_omega_00 = 1d / (t1_h0 * t1_h0);
		double t1_omega_01 = 1d / (t1_h0 * t1_h1);
		double t1_omega_02 = 1d / (t1_h0 * t1_h2);
		double t1_omega_10 = t1_omega_01;
		double t1_omega_11 = 1d / (t1_h1 * t1_h1);
		double t1_omega_12 = 1d / (t1_h1 * t1_h2);
		double t1_omega_20 = t1_omega_02;
		double t1_omega_21 = t1_omega_12;
		double t1_omega_22 = 1d / (t1_h2 * t1_h2);
		double t2_omega_00 = 1d / (t2_h0 * t2_h0);
		double t2_omega_01 = 1d / (t2_h0 * t2_h1);
		double t2_omega_02 = 1d / (t2_h0 * t2_h2);
		double t2_omega_10 = t2_omega_01;
		double t2_omega_11 = 1d / (t2_h1 * t2_h1);
		double t2_omega_12 = 1d / (t2_h1 * t2_h2);
		double t2_omega_20 = t2_omega_02;
		double t2_omega_21 = t2_omega_12;
		double t2_omega_22 = 1d / (t2_h2 * t2_h2);
		
		Vec3D e0_unit = e0 / e0_magnitude;
		Vec3D t1_e0_normal = Vec3D.cross(e0_unit, n1); // "m_0" in paper.
		Vec3D t1_e1_normal = Vec3D.cross(t1_e1 / t1_e1_magnitude, n1); // "m_1" in paper.
		Vec3D t1_e2_normal = Vec3D.cross(n1, t1_e2 / t1_e2_magnitude); // "m_2" in paper.
		Vec3D t2_e0_normal = Vec3D.cross(n2, e0_unit); // "~m_0" in paper.
		Vec3D t2_e1_normal = Vec3D.cross(n2, t2_e1 / t2_e1_magnitude); // "~m_1" in paper.
		Vec3D t2_e2_normal = Vec3D.cross(t2_e2 / t2_e2_magnitude, n2); // "~m_2" in paper.

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

		MatD t1_N0 = t1_M0 / (e0_magnitude * e0_magnitude);
		MatD t2_N0 = t2_M0 / (e0_magnitude * e0_magnitude);
		
		MatD t1_P10 = t1_M0.transpose.mul(t1_omega_10 * cos_t1_alpha1); // Perform: MatD t1_P10 = t1_omega_10 * Mathf.Cos(t1_alpha1) * t1_M0.transpose;
		MatD t1_P20 = t1_M0.transpose.mul(t1_omega_20 * cos_t1_alpha2);
		MatD t1_P12 = t1_M2.transpose.mul(t1_omega_12 * cos_t1_alpha1);
		MatD t1_P21 = t1_M1.transpose.mul(t1_omega_21 * cos_t1_alpha2);
		MatD t1_P11 = t1_M1.transpose.mul(t1_omega_11 * cos_t1_alpha1);
		MatD t1_P22 = t1_M2.transpose.mul(t1_omega_22 * cos_t1_alpha2);
		MatD t2_P10 = t2_M0.transpose.mul(t2_omega_10 * cos_t2_alpha1);
		MatD t2_P20 = t2_M0.transpose.mul(t2_omega_20 * cos_t2_alpha2);
		MatD t2_P12 = t2_M2.transpose.mul(t2_omega_12 * cos_t2_alpha1);
		MatD t2_P21 = t2_M1.transpose.mul(t2_omega_21 * cos_t2_alpha2);
		MatD t2_P11 = t2_M1.transpose.mul(t2_omega_11 * cos_t2_alpha1);
		MatD t2_P22 = t2_M2.transpose.mul(t2_omega_22 * cos_t2_alpha2);

		// Construct 3x3 building blocks for the teta Hessian.
		MatD teta_hess_00 = -getMatPlusTransposedMat(t1_Q0);
		MatD teta_hess_33 = -getMatPlusTransposedMat(t2_Q0);
		MatD teta_hess_11 = getMatPlusTransposedMat(t1_P11).sub(t1_N0).add(getMatPlusTransposedMat(t2_P11)).sub(t2_N0);
		MatD teta_hess_22 = getMatPlusTransposedMat(t1_P22).sub(t1_N0).add(getMatPlusTransposedMat(t2_P22)).sub(t2_N0);
		MatD teta_hess_10 = getMatPlusTransposedMat(t1_P10).sub(t1_Q1);
		MatD teta_hess_20 = getMatPlusTransposedMat(t1_P20).sub(t1_Q2);
		MatD teta_hess_13 = getMatPlusTransposedMat(t2_P10).sub(t2_Q1);
		MatD teta_hess_23 = getMatPlusTransposedMat(t2_P20).sub(t2_Q2);
		MatD teta_hess_12 = t1_P21.transpose.add(t1_P12).add(t1_N0).add(t2_P12).add(t2_P21.transpose).add(t2_N0);
		MatD teta_hess_03 = new MatD(3, 3);

		// Construct teta Hessian.
		MatD teta_hess = new MatD(12, 12);
		for(int row = 0; row < 3; row++) {
			for(int col = row; col < 3; col++) {

				// Set Hessian diagonal building blocks.
				// Use only the top triangle of the blocks to ensure that there are no rounding errors in the result that will make it non-symmetrical.
				teta_hess[row, col] = teta_hess[col, row] = teta_hess_00[row, col];
				teta_hess[row + 3, col + 3] = teta_hess[col + 3, row + 3] = teta_hess_11[row, col];
				teta_hess[row + 6, col + 6] = teta_hess[col + 6, row + 6] = teta_hess_22[row, col];
				teta_hess[row + 9, col + 9] = teta_hess[col + 9, row + 9] = teta_hess_33[row, col];
			}
			for(int col = 0; col < 3; col++) {

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
		VecD d_teta_d_ve1 = (cos_t1_alpha2 / t1_h1) * n1 + (cos_t2_alpha2 / t2_h1) * n2;
		VecD d_teta_d_ve2 = (cos_t1_alpha1 / t1_h2) * n1 + (cos_t2_alpha1 / t2_h2) * n2;
		VecD d_teta_d_vf1 = (-n1).div(t1_h0); // Perform: -n1 / t1_h0;
		VecD d_teta_d_vf2 = (-n2).div(t2_h0); // Perform: -n2 / t2_h0;
		
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
		 * 
		 * Perform: return d_fi_d_teta * teta_hess + 2d * edge.undeformedLength / h_e_undeformed * MatD.fromVecMultiplication(d_teta_d_ve1ve2vf1vf2, d_teta_d_ve1ve2vf1vf2);
		 */
		return teta_hess.mul(d_fi_d_teta).add(MatD.fromVecMultiplication(d_teta_d_ve1ve2vf1vf2, d_teta_d_ve1ve2vf1vf2).mul(2d * edge.undeformedLength / h_e_undeformed));
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

		// TODO - Remove old implementation if the new implementation works properly.
		//MatD mat = hess.Clone();
		//double diagIncrement = 0;
		//for(int i = 0; i < mat.numRows; i++) {

		//	// Make pivot positive.
		//	if(mat[i, i] <= 0) {
		//		double inc = -mat[i, i] * 1.0001d; // Add small delta to ensure a positive non-zero value. Multiply to prevent rounding to 0.
		//		if(inc == 0d) {
		//			inc = 0.0001d; // Set small delta to ensure a positive non-zero value.
		//		}
		//		diagIncrement += inc;
		//		mat.addDiag(inc);
		//	}

		//	// Sweep lower rows.
		//	for(int row = i + 1; row < mat.numRows; row++) {
		//		double factor = mat[row, i] / mat[i, i];
		//		for(int col = 0; col < mat.numColumns; col++) {
		//			mat[row, col] -= factor * mat[i, col];
		//		}
		//	}
		//}
		//if(diagIncrement != 0d) {
		//	hess.addDiag(diagIncrement);
		//}

		double minEigenValue = 0.1d;

		MathNet.Numerics.LinearAlgebra.Matrix<double> mat = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.DenseOfArray(hess.asDoubleArray());
		MathNet.Numerics.LinearAlgebra.Factorization.Evd<double> eigValDecomp = mat.Evd();
		MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex> eigenValues = eigValDecomp.EigenValues;
		MathNet.Numerics.LinearAlgebra.Matrix<double> eigenVectors = eigValDecomp.EigenVectors;
		//print("Hess: " + hess);
		for(int i = 0; i < eigenValues.Count; i++) {
			//print("Eigenvalue: " + eigenValues[i]);
			if(eigenValues[i].Imaginary != 0d) {
				print("Imaginary eigenvalue found: " + eigenValues[i]);

				// Return the identity matrix as fallback.
				for(int row = 0; row < hess.numRows; row++) {
					for(int col = 0; col < hess.numColumns; col++) {
						hess[row, col] = (row == col ? 1d : 0d);
					}
				}
				return;
				//throw new Exception("Imaginary eigenvalue found: " + eigenValues[i]);
			}
			if(eigenValues[i].Real < minEigenValue) {

				// Get eigen vector.
				VecD eigenVector = new VecD(eigenVectors.RowCount);
				for(int j = 0; j < eigenVectors.RowCount; j++) {
					eigenVector[j] = eigenVectors[j, i];
				}
				//print("EigenValues[" + i + "]: " + eigenValues[i] + ", Eigenvec: " + eigenVector);

				// Create vec * vec' matrix.
				MatD eigVecMat = MatD.fromVecMultiplication(eigenVector, eigenVector);
				//print("EigVecMat: " + eigVecMat);

				// Add multiple of the vec * vec' matrix to the Hessian to make the corresponding eigenvalue positive.
				double factor = -eigenValues[i].Real + minEigenValue; // Add small delta to ensure a positive non-zero value.
				//print("Factor: " + factor);
				//if(factor < minEigenValue) {
				//	print("Applying rounding fix. Factor: " + factor);
				//	factor = minEigenValue;
				//}
				if(factor + eigenValues[i].Real <= 0d) {
					print("Factor is not high enough. Rounding error? Resulting eigenvalue: " + (factor + eigenValues[i].Real));
				}
				hess.add(eigVecMat.mul(factor));

				//mat = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.DenseOfArray(hess.asDoubleArray());
				//eigValDecomp = mat.Evd();
				//eigenValues = eigValDecomp.EigenValues;
				//eigenVectors = eigValDecomp.EigenVectors;
			}
		}
		//print("Hess (after): " + hess);

		//// TODO - Remove test below after validating that it always succeeds.
		//mat = MathNet.Numerics.LinearAlgebra.Matrix<double>.Build.DenseOfArray(hess.asDoubleArray());
		//eigValDecomp = mat.Evd();
		//eigenValues = eigValDecomp.EigenValues;
		//eigenVectors = eigValDecomp.EigenVectors;
		//for(int i = 0; i < eigenValues.Count; i++) {
		//	if(eigenValues[i].Real <= 0d) {
		//		print("Negative eigenvalue! eigenValues[" + i + "] = " + eigenValues[i] + " Matrix: " + hess);

		//		// Return the identity matrix as fallback.
		//		for(int row = 0; row < hess.numRows; row++) {
		//			for(int col = 0; col < hess.numColumns; col++) {
		//				mat[row, col] = (row == col ? 1d : 0d);
		//			}
		//		}
		//		return;
		//		//throw new Exception("Negative eigenvalue! eigenValues[" + i + "] = " + eigenValues[i] + " Matrix: " + hess);
		//	}
		//}
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

	/*
	 * Solves the given linear equation using a sparse direct solver.
	 * Equation: mat * x = vec
	 * Returns vector x.
	 */
	private VecD sparseDirectSolve(MatD mat, VecD vec) {
		return ParDiSoLib.sparseDirectSymmetricSolve(mat, vec);
	}

	/*
	 * Solves the given linear equation using a sparse direct solver.
	 * Equation: mat * x = vec
	 * Returns vector x.
	 */
	private VecD sparseDirectSolve(MathNet.Numerics.LinearAlgebra.Double.SparseMatrix mat, VecD vec) {
		return ParDiSoLib.sparseDirectSymmetricSolve(mat, vec);
	}

	public void onSaveSailShapeButtonPress() {
		Mesh mesh = this.getMesh();
		new SailConfiguration(this.vertexPositions, mesh.triangles).storeToFile("sail");
	}

	public void onLoadSailShapeButtonPress() {
		SailConfiguration sailConfiguration = SailConfiguration.loadFromFile("sail");
		Mesh mesh = this.getMesh();
		mesh.vertices = vecToVec(sailConfiguration.vertexPositions);
		mesh.triangles = sailConfiguration.triangles;
		mesh.normals = new Vector3[sailConfiguration.vertexPositions.Length];
		mesh.RecalculateNormals();
		this.loadMesh(mesh, 1d);
	}

	public void onSaveSailMeasurementsButtonPress() {
		Mesh mesh = this.getMesh();
		if(this.measurementsGenerateFactor < 0f) {
			this.measurementsGenerateFactor = 0f;
		} else if(this.measurementsGenerateFactor > 1f) {
			this.measurementsGenerateFactor = 1f;
		}
		int numMeasurements = (int) (this.vertexPositions.Length * this.measurementsGenerateFactor);
		new SailMeasurements(this.vertexPositions, MeshUtils.generateMeasurements(this.vertexPositions, numMeasurements)).storeToFile("measurements");
	}

	public void onLoadSailMeasurementsButtonPress() {
		SailMeasurements sailMeasurements = SailMeasurements.loadFromFile("measurements");
		this.measurements = sailMeasurements.measurements;
	}

	private static Vector3[] vecToVec(Vec3D[] vec) {
		Vector3[] ret = new Vector3[vec.Length];
		for(int i = 0; i < vec.Length; i++) {
			ret[i] = vecToVec(vec[i]);
		}
		return ret;
	}

	private static Vector3 vecToVec(Vec3D vec) {
		return new Vector3((float) vec.x, (float) vec.y, (float) vec.z);
	}
}
