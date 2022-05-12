using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Threading.Tasks;
using UnityEditor;
using UnityEngine;
using UnityEngine.UI;

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
	private Vector3 shellObjPosition = new Vector3(0f, 1f, 0f);

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
	public Vector3 windPressureVec; // [N/m^2] Wind pressure vector, applied to triangles through triangle normal direction.
	public double windPressure; // [N/m^2] Static wind pressure, applied to triangles in triangle normal direction.
	public double gravityConstant = 9.81d;

	public float measurementsGenerateFactor = 0.1f; // Defines the number of generated measurements by multiplying this with the amount of vertices.
	private SailMeasurements measurements = null;
	private Vec3D[] realDataTruthMeasurements = null; // Contains a subset of measurements to estimate the reconstruction distance with.

	// Gradient descent specific settings.
	public double kGradientDescent;
	public double maxGradientDescentStep;

	// Explicit integration specific settings.
	public float dampingFactor = 0.99f; // [F * s / m] = [kg / s] ? Factor applied to vertex velocity per time step.

	// Optimization integrator specific settings.
	public bool DoStaticMinimization = false; // When set to true, velocity preserving and damping terms will be ignored.
	public double VelocityTerminationThreshold = 0.01d; // [m/s]. When the maximum vertex velocity after a step is below this threshold, simulation is stopped.
	public double dampingConstant = 0.001d;
	public double directVelocityDampingFactor = 1d;
	public double maxWindSpeed = 50d;
	public double maxDeltaWindSpeed = 5d;
	public double kMeasurementSprings = 0d;
	public double eGradientMagnitudeTerminationThreshold = 0.5d;
	private double lastLineSearchAlpha = 1d;
	public double MaxLineSearchTimeMS = 10000;
	public double MinLineSearchAlpha = 0.000001d;
	public int MinNumNewtonIterations = 0; // Minimum number of Newton steps to do before allowing termination.
	public double MaxNewtonsMethodLoopTimeMS = 10000;
	public double windNoiseFuncMagnitudeSlope = 0d;

	// Cached vertex/triangle properties.
	private Vec3D[] triangleNormals;
	private Vec3D[] undeformedTriangleNormals;
	private double[] triangleAreas;
	private double[] undeformedTriangleAreas;

	// Simulation recording.
	public static string storageBaseDirPath;
	private MeshRecorder meshRecorder = null;
	private bool isRecording = false;

	// Reconstruction.
	private List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>(); // A list of reconstruction setups for automated reconstruction.
	private int reconstructionSetupsIndex = -1;
	public ReconstructionStage ReconstructionStage = ReconstructionStage.DISABLED;
	private int stepCount = 0;
	public double MaxVertexMoveToMeasurementsStep = 0.05; // Max distance per step to move vertices towards their corresponding measurements with to snap to measurements.
	public int NumWindReconstructionSteps1 = 100; // The number of steps to take for the first wind reconstruction phase.
	public int NumWindReconstructionSteps2 = 100; // The number of steps to take for the second wind reconstruction phase.

	// Debugging.
	private Visualizer visualizer;
	public VisualizationType visualizationType = VisualizationType.NONE;
	private Text avgReconDistTextObj;
	private Text maxReconDistTextObj;

	void Awake() {
		storageBaseDirPath = Application.dataPath + "/StoredData";
		QualitySettings.vSyncCount = 0; // Disable V-sync.
		Application.targetFrameRate = 100; // Set max framerate.
		System.Threading.Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US"); // To print decimal points instead of commas.
	}

	// Start is called before the first frame update.
	void Start() {
		//MeasurementGeneration.generateFarthestFirstMeasurements();
		//MeasurementGeneration.generateTest3FarthestFirstMeasurements();
		//MeasurementGeneration.generateTest4FarthestFirstMeasurements();

		// Initialize game object references.
		this.avgReconDistTextObj = GameObject.Find("AvgReconDistText").GetComponent<Text>();
		this.maxReconDistTextObj = GameObject.Find("MaxReconDistText").GetComponent<Text>();

		// Initialize visualizer.
		this.visualizer = new Visualizer(this.shellObjPosition);

		// Create the new object in the scene.
		Mesh mesh;
		double undeformedInnerEdgeLengthFactor;
		if(this.shellObj == null) {
			this.shellObj = new GameObject();
			this.shellObj.AddComponent<MeshRenderer>();
			this.shellObj.AddComponent<MeshFilter>();
			//mesh = MeshHelper.createSquareMesh(5, 5, 1);
			//mesh = MeshHelper.createTriangleMesh(5, 5, 1);
			//undeformedInnerEdgeLengthFactor = 1d;
			mesh = MeshHelper.createTriangleMesh(3.5f, 7f, 5); // 5 subdivisions leads to 561 vertices and 1024 triangles.
			undeformedInnerEdgeLengthFactor = 1.01d;
			//mesh = MeshHelper.createFlatRectangleMesh(0.1f, 0.05f, 4);
			//mesh = MeshHelper.createFlatRectangleMesh(5f, 5f, 1);
			//undeformedInnerEdgeLengthFactor = 1d;
		} else {
			mesh = this.shellObj.GetComponent<MeshFilter>().mesh;
			undeformedInnerEdgeLengthFactor = 1d;
		}

		// Set the mesh material.
		MeshRenderer meshRenderer = this.shellObj.GetComponent<MeshRenderer>();
		meshRenderer.sharedMaterial = new Material(Shader.Find("Custom/StandardTwoSides"));

		// Load the mesh.
		this.loadMesh(mesh, undeformedInnerEdgeLengthFactor);

		// Constrain the mesh.
		this.verticesMovementConstraints = this.createVertexContraints();

		// Set the shell position.
		this.shellObj.transform.position = this.shellObjPosition;

		//// Load the default sail and open a measurements load dialog. TODO - Remove temporary direct load.
		//this.loadSailConfiguration(SailConfiguration.loadFromFile(storageBaseDirPath
		//		+ "/SailData/numVerts=561, kLength=500, kArea=500, kBend=0.01, thickness=0.002, density=40, steady state with no external forces.sailshapedata"));
		//this.onLoadSailMeasurementsButtonPress();

		// Initialize automated reconstruction setups.
		if(this.ReconstructionStage != ReconstructionStage.DISABLED) {
			this.reconstructionSetups.AddRange(this.getReconstructionSetupsTest1FarthestFirst());
			this.reconstructionSetups.AddRange(this.getReconstructionSetupsTest2FarthestFirst());
			this.reconstructionSetups.AddRange(this.getReconstructionSetupsTest3FarthestFirst());
			this.reconstructionSetups.AddRange(this.getReconstructionSetupsTest4FarthestFirst());
			this.reconstructionSetups.AddRange(this.getReconstructionSetupsTest5FarthestFirst());
			this.reconstructionSetups.AddRange(this.getReconstructionSetupsTest6FarthestFirst());
			this.reconstructionSetups.AddRange(this.getReconstructionSetupsRealDataTest());
			this.ReconstructionStage = ReconstructionStage.DONE;
		}

		// TODO - Remove temp test.
		//this.onLoadSailShapeButtonPress();
		//this.recalcTriangleNormalsAndAreas(mesh.triangles, this.vertexPositions);
		//int numRuns = 10;
		//Stopwatch stopWatch = Stopwatch.StartNew();
		//for(int i = 0; i < numRuns; i++) {
		//	double startTime = stopWatch.ElapsedMilliseconds;
		//	MathNet.Numerics.LinearAlgebra.Double.SparseMatrix mat = this.getSystemEnergyHessianSparseRepresentation(this.vertexPositions, this.getMesh().triangles);
		//	print("getSystemEnergyHessianSparseRepresentation returned after " + (stopWatch.ElapsedMilliseconds - startTime) + "ms. " + mat.RowCount + "x" + mat.ColumnCount);
		//}
		//print("getSystemEnergyHessianSparseRepresentation returned after average " + (stopWatch.ElapsedMilliseconds / numRuns) + "ms.");

		//// TODO - Remove temp test.
		//this.kLength = 1;
		//this.kArea = 1;
		//this.kBend = 1;
		//this.onLoadSailShapeButtonPress();
		//this.recalcTriangleNormalsAndAreas(mesh.triangles, this.vertexPositions);
		//MathNet.Numerics.LinearAlgebra.Double.SparseMatrix mat1 = this.getSystemEnergyHessianSparseRepresentation(this.vertexPositions, this.getMesh().triangles);
		//MathNet.Numerics.LinearAlgebra.Double.SparseMatrix mat2 = this.getSystemEnergyHessianSparseRepresentationMultiThreadedTriplets(this.vertexPositions, this.getMesh().triangles);
		//MathNet.Numerics.LinearAlgebra.Double.SparseMatrix mat3 = this.getSystemEnergyHessianSparseRepresentationMultiThreaded(this.vertexPositions, this.getMesh().triangles);
		//for(int row = 0; row < mat1.RowCount; row++) {
		//	for(int col = 0; col < mat1.ColumnCount; col++) {
		//		double v1 = mat1[row, col];
		//		double v2 = mat2[row, col];
		//		double v3 = mat3[row, col];
		//		//if(v1 != v2 || v2 != v3) {
		//		if(v3 != v2) {
		//			print("Mismatching value in Hessian implementations: row = " + row + ", col = " + col + ", v1 = " + v1 + ", v2 = " + v2 + ", v3 = " + v3
		//					+ ", v3 - v2 = " + (v3 - v2));
		//		}
		//	}
		//}
		//print("mat1: " + new MatD(mat1.ToArray()).toFancyString());
		//print("mat2: " + new MatD(mat2.ToArray()).toFancyString());
		//print("Test finished.");
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest1() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double windMag in new double[] {64.6963, 1035.1}) {
			foreach(double windDeg in new double[] {30, 45, 60}) {
				foreach(int n in new int[] {1, 5, 10, 15}) {
					foreach(double m in new double[] {1, 2.5, 5}) {
						string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
								+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
						reconstructionSetups.Add(new ReconstructionSetup {
							sailStartConfigurationRelPath = restConfigurationSailRelPath,
							sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ",m=" + m + ".measurements",
							resultsStorageRelPath = "test1/" + fileNameNoEx + "/n=" + n + ",m=" + m + ".results",
							kLength = 15621f,
							kArea = 8125f,
							kBend = 10f,
							shellThickness = 0.00025f,
							shellMaterialDensity = 920f,
							useFlatUndeformedBendState = true,
							initialWindPressureVec = new Vec3D(0, 0, 0),
							initialWindPressure = 0d,
							gravityConstant = 9.81d,
							doStaticMinimization = true,
							maxWindSpeed = this.maxWindSpeed,
							maxDeltaWindSpeed = this.maxDeltaWindSpeed,
							minNumNewtonIterations = this.MinNumNewtonIterations,
							numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
							numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
						});
					}
				}
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest1FarthestFirst() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double windMag in new double[] {64.6963, 1035.1}) {
			foreach(double windDeg in new double[] {30, 45, 60}) {
				foreach(int n in new int[] {5, 10, 20, 30, 40, 50, 70, 90, 110, 130, 140, 150, 160}) {
					string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
							+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
					reconstructionSetups.Add(new ReconstructionSetup {
						sailStartConfigurationRelPath = restConfigurationSailRelPath,
						sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ".measurements",
						resultsStorageRelPath = "test1/" + fileNameNoEx + "/n=" + n + ".results",
						kLength = 15621f,
						kArea = 8125f,
						kBend = 10f,
						shellThickness = 0.00025f,
						shellMaterialDensity = 920f,
						useFlatUndeformedBendState = true,
						initialWindPressureVec = new Vec3D(0, 0, 0),
						initialWindPressure = 0d,
						gravityConstant = 9.81d,
						doStaticMinimization = true,
						maxWindSpeed = this.maxWindSpeed,
						maxDeltaWindSpeed = this.maxDeltaWindSpeed,
						minNumNewtonIterations = this.MinNumNewtonIterations,
						numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
						numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
					});
				}
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest2() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double windMag in new double[] {1035.1}) {
			foreach(double windDeg in new double[] {45}) {
				foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
					int n = (int) nm[0];
					double m = nm[1];
					foreach(double stiffnessFactor in new double[] {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
							1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0}) {
						string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
								+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
						reconstructionSetups.Add(new ReconstructionSetup {
							sailStartConfigurationRelPath = restConfigurationSailRelPath,
							sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ",m=" + m + ".measurements",
							resultsStorageRelPath = "test2/" + fileNameNoEx + "/n=" + n + ",m=" + m + ",stiffFac=" + stiffnessFactor + ".results",
							kLength = (float) (15621d * stiffnessFactor),
							kArea = (float) (8125d * stiffnessFactor),
							kBend = 10f,
							shellThickness = 0.00025f,
							shellMaterialDensity = 920f,
							useFlatUndeformedBendState = true,
							initialWindPressureVec = new Vec3D(0, 0, 0),
							initialWindPressure = 0d,
							gravityConstant = 9.81d,
							doStaticMinimization = true,
							maxWindSpeed = this.maxWindSpeed,
							maxDeltaWindSpeed = this.maxDeltaWindSpeed,
							minNumNewtonIterations = this.MinNumNewtonIterations,
							numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
							numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
						});
					}
				}
			}
		}
		return reconstructionSetups;
	}
	private List<ReconstructionSetup> getReconstructionSetupsTest2FarthestFirst() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double windMag in new double[] {1035.1}) {
			foreach(double windDeg in new double[] {45}) {
				foreach(int n in new int[] {30, 50, 110}) {
					foreach(double stiffnessFactor in new double[] {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
							1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 3.0, 4.0, 5.0}) {
						string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
								+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
						reconstructionSetups.Add(new ReconstructionSetup {
							sailStartConfigurationRelPath = restConfigurationSailRelPath,
							sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ".measurements",
							resultsStorageRelPath = "test2/" + fileNameNoEx + "/n=" + n + ",stiffFac=" + stiffnessFactor + ".results",
							kLength = (float) (15621d * stiffnessFactor),
							kArea = (float) (8125d * stiffnessFactor),
							kBend = 10f,
							shellThickness = 0.00025f,
							shellMaterialDensity = 920f,
							useFlatUndeformedBendState = true,
							initialWindPressureVec = new Vec3D(0, 0, 0),
							initialWindPressure = 0d,
							gravityConstant = 9.81d,
							doStaticMinimization = true,
							maxWindSpeed = this.maxWindSpeed,
							maxDeltaWindSpeed = this.maxDeltaWindSpeed,
							minNumNewtonIterations = this.MinNumNewtonIterations,
							numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
							numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
						});
					}
				}
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest3() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
			int n = (int) nm[0];
			double m = nm[1];
			foreach(double windPressFac in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
				string measurementsDirName = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=45,windPressFac=" + windPressFac;
				reconstructionSetups.Add(new ReconstructionSetup {
					sailStartConfigurationRelPath = restConfigurationSailRelPath,
					sailMeasurementsRelPath = "test3/" + measurementsDirName + "/n=" + n + ",m=" + m + ".measurements",
					resultsStorageRelPath = "test3/" + measurementsDirName + "/n=" + n + ",m=" + m + ".results",
					kLength = 15621f,
					kArea = 8125f,
					kBend = 10f,
					shellThickness = 0.00025f,
					shellMaterialDensity = 920f,
					useFlatUndeformedBendState = true,
					initialWindPressureVec = new Vec3D(0, 0, 0),
					initialWindPressure = 0d,
					gravityConstant = 9.81d,
					doStaticMinimization = true,
					maxWindSpeed = this.maxWindSpeed,
					maxDeltaWindSpeed = this.maxDeltaWindSpeed,
					minNumNewtonIterations = this.MinNumNewtonIterations,
					numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
					numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
				});
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest3FarthestFirst() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(int n in new int[] {30, 50, 110}) {
			foreach(double windPressFac in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
				string measurementsDirName = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=45,windPressFac=" + windPressFac;
				reconstructionSetups.Add(new ReconstructionSetup {
					sailStartConfigurationRelPath = restConfigurationSailRelPath,
					sailMeasurementsRelPath = "test3/" + measurementsDirName + "/n=" + n + ".measurements",
					resultsStorageRelPath = "test3/" + measurementsDirName + "/n=" + n + ".results",
					kLength = 15621f,
					kArea = 8125f,
					kBend = 10f,
					shellThickness = 0.00025f,
					shellMaterialDensity = 920f,
					useFlatUndeformedBendState = true,
					initialWindPressureVec = new Vec3D(0, 0, 0),
					initialWindPressure = 0d,
					gravityConstant = 9.81d,
					doStaticMinimization = true,
					maxWindSpeed = this.maxWindSpeed,
					maxDeltaWindSpeed = this.maxDeltaWindSpeed,
					minNumNewtonIterations = this.MinNumNewtonIterations,
					numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
					numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
				});
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest4() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
			int n = (int) nm[0];
			double m = nm[1];
			foreach(double noiseMagSlope in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
				string measurementsDirName = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=60,noiseMagSlope=" + noiseMagSlope;
				reconstructionSetups.Add(new ReconstructionSetup {
					sailStartConfigurationRelPath = restConfigurationSailRelPath,
					sailMeasurementsRelPath = "test4/" + measurementsDirName + "/n=" + n + ",m=" + m + ".measurements",
					resultsStorageRelPath = "test4/" + measurementsDirName + "/n=" + n + ",m=" + m + ".results",
					kLength = 15621f,
					kArea = 8125f,
					kBend = 10f,
					shellThickness = 0.00025f,
					shellMaterialDensity = 920f,
					useFlatUndeformedBendState = true,
					initialWindPressureVec = new Vec3D(0, 0, 0),
					initialWindPressure = 0d,
					gravityConstant = 9.81d,
					doStaticMinimization = true,
					maxWindSpeed = this.maxWindSpeed,
					maxDeltaWindSpeed = this.maxDeltaWindSpeed,
					minNumNewtonIterations = this.MinNumNewtonIterations,
					numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
					numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
				});
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest4FarthestFirst() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(int n in new int[] {30, 50, 110}) {
			foreach(double noiseMagSlope in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {
				string measurementsDirName = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=60,noiseMagSlope=" + noiseMagSlope;
				reconstructionSetups.Add(new ReconstructionSetup {
					sailStartConfigurationRelPath = restConfigurationSailRelPath,
					sailMeasurementsRelPath = "test4/" + measurementsDirName + "/n=" + n + ".measurements",
					resultsStorageRelPath = "test4/" + measurementsDirName + "/n=" + n + ".results",
					kLength = 15621f,
					kArea = 8125f,
					kBend = 10f,
					shellThickness = 0.00025f,
					shellMaterialDensity = 920f,
					useFlatUndeformedBendState = true,
					initialWindPressureVec = new Vec3D(0, 0, 0),
					initialWindPressure = 0d,
					gravityConstant = 9.81d,
					doStaticMinimization = true,
					maxWindSpeed = this.maxWindSpeed,
					maxDeltaWindSpeed = this.maxDeltaWindSpeed,
					minNumNewtonIterations = this.MinNumNewtonIterations,
					numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
					numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
				});
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest5() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double windMag in new double[] {1035.1}) {
			foreach(double windDeg in new double[] {45}) {
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
				foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
					int n = (int) nm[0];
					double m = nm[1];
					foreach(Tuple<Vec3D, double>[] measurementsIgnoreSpheres in new Tuple<Vec3D, double>[][] {
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(0d, 0d, 0d), 1.5d)},
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(0d, 0d, 0d), 2.5d)},
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(3.5d, 0d, 0d), 1.5d)},
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(3.5d, 0d, 0d), 2.5d)},
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(0d, 7d, 0d), 2.5d)},
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(3.5d / 3d, 7d / 3d, 0d), 1d)},
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(3.5d / 3d, 7d / 3d, 0d), 1.5d)}}) {
						string ignoreSpheresStr = "{";
						foreach(Tuple<Vec3D, double> measurementsIgnoreSphere in measurementsIgnoreSpheres) {
							Vec3D origin = measurementsIgnoreSphere.Item1;
							double radius = measurementsIgnoreSphere.Item2;
							if(ignoreSpheresStr.Length > 1) {
								ignoreSpheresStr += ", ";
							}
							ignoreSpheresStr += "{" + origin.ToString() + ", " + radius + "}";
						}
						ignoreSpheresStr += "}";
						reconstructionSetups.Add(new ReconstructionSetup {
							sailStartConfigurationRelPath = restConfigurationSailRelPath,
							sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ",m=" + m + ".measurements",
							resultsStorageRelPath = "test5/" + fileNameNoEx + "/n=" + n + ",m=" + m
									+ ",ignoreSpheres=" + ignoreSpheresStr.Replace(", ", ",") + ".results",
							kLength = 15621f,
							kArea = 8125f,
							kBend = 10f,
							shellThickness = 0.00025f,
							shellMaterialDensity = 920f,
							useFlatUndeformedBendState = true,
							initialWindPressureVec = new Vec3D(0, 0, 0),
							initialWindPressure = 0d,
							gravityConstant = 9.81d,
							doStaticMinimization = true,
							maxWindSpeed = this.maxWindSpeed,
							maxDeltaWindSpeed = this.maxDeltaWindSpeed,
							minNumNewtonIterations = this.MinNumNewtonIterations,
							numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
							numWindReconstructionSteps2 = this.NumWindReconstructionSteps2,
							measurementsIgnoreSpheres = measurementsIgnoreSpheres
						});
					}
				}
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest5FarthestFirst() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double windMag in new double[] {1035.1}) {
			foreach(double windDeg in new double[] {45}) {
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
				foreach(int n in new int[] {30, 50, 110}) {
					foreach(Tuple<Vec3D, double>[] measurementsIgnoreSpheres in new Tuple<Vec3D, double>[][] {
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(0d, 0d, 0d), 1.75d)}, // Mast-boom corner, half sail width radius.
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(3.5d, 0d, 0d), 1.75d)}, // End of boom corner, half sail width radius.
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(0d, 7d, 0d), 3.5d)}, // Top of mast corner, half sail height radius.

							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(1d, 2.5d, 0d), 1.75d)}, // Sail middle, very large radius
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(1d, 2.5d, 0d), 1.25d)}, // Sail middle, likely the most important vertices.

							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(1d, 3d, 0d), 1d)}, // Sail left bottom area small.
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(2.5d, 3d, 0d), 1d)}, // Sail right bottom area small.
							new Tuple<Vec3D, double>[] {new Tuple<Vec3D, double>(new Vec3D(1.75d, 1.5d, 0d), 1d)} // Sail middle small.
						}) {
						string ignoreSpheresStr = "{";
						foreach(Tuple<Vec3D, double> measurementsIgnoreSphere in measurementsIgnoreSpheres) {
							Vec3D origin = measurementsIgnoreSphere.Item1;
							double radius = measurementsIgnoreSphere.Item2;
							if(ignoreSpheresStr.Length > 1) {
								ignoreSpheresStr += ", ";
							}
							ignoreSpheresStr += "{" + origin.ToString() + ", " + radius + "}";
						}
						ignoreSpheresStr += "}";
						reconstructionSetups.Add(new ReconstructionSetup {
							sailStartConfigurationRelPath = restConfigurationSailRelPath,
							sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ".measurements",
							resultsStorageRelPath = "test5/" + fileNameNoEx + "/n=" + n
									+ ",ignoreSpheres=" + ignoreSpheresStr.Replace(", ", ",") + ".results",
							kLength = 15621f,
							kArea = 8125f,
							kBend = 10f,
							shellThickness = 0.00025f,
							shellMaterialDensity = 920f,
							useFlatUndeformedBendState = true,
							initialWindPressureVec = new Vec3D(0, 0, 0),
							initialWindPressure = 0d,
							gravityConstant = 9.81d,
							doStaticMinimization = true,
							maxWindSpeed = this.maxWindSpeed,
							maxDeltaWindSpeed = this.maxDeltaWindSpeed,
							minNumNewtonIterations = this.MinNumNewtonIterations,
							numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
							numWindReconstructionSteps2 = this.NumWindReconstructionSteps2,
							measurementsIgnoreSpheres = measurementsIgnoreSpheres
						});
					}
				}
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest6() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
			int n = (int) nm[0];
			double m = nm[1];
			for(int randomSeed = 0; randomSeed < 10; randomSeed++) {
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920 ss with g=9.81,wm=1035.1,wd=60";
				reconstructionSetups.Add(new ReconstructionSetup {
					sailStartConfigurationRelPath = restConfigurationSailRelPath,
					sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ",m=" + m + ".measurements",
					resultsStorageRelPath = "test6/" + fileNameNoEx + ",randomSeed=" + randomSeed + "/n=" + n + ",m=" + m + ".results",
					kLength = 15621f,
					kArea = 8125f,
					kBend = 10f,
					shellThickness = 0.00025f,
					shellMaterialDensity = 920f,
					useFlatUndeformedBendState = true,
					initialWindPressureVec = new Vec3D(0, 0, 0),
					initialWindPressure = 0d,
					gravityConstant = 9.81d,
					doStaticMinimization = true,
					maxWindSpeed = this.maxWindSpeed,
					maxDeltaWindSpeed = this.maxDeltaWindSpeed,
					minNumNewtonIterations = this.MinNumNewtonIterations,
					numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
					numWindReconstructionSteps2 = this.NumWindReconstructionSteps2,
					initialPositionNoiseMagnitude = 0.01,
					initialPositionNoiseRandomSeed = randomSeed
				});
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsTest6FarthestFirst() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		string restConfigurationSailRelPath = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920,"
				+ " steady state without external forces.sailshapedata";
		foreach(int n in new int[] {30, 50, 110}) {
			for(int randomSeed = 0; randomSeed < 20; randomSeed++) {
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920 ss with g=9.81,wm=1035.1,wd=60";
				reconstructionSetups.Add(new ReconstructionSetup {
					sailStartConfigurationRelPath = restConfigurationSailRelPath,
					sailMeasurementsRelPath = fileNameNoEx + "/n=" + n + ".measurements",
					resultsStorageRelPath = "test6/" + fileNameNoEx + ",randomSeed=" + randomSeed + "/n=" + n + ".results",
					kLength = 15621f,
					kArea = 8125f,
					kBend = 10f,
					shellThickness = 0.00025f,
					shellMaterialDensity = 920f,
					useFlatUndeformedBendState = true,
					initialWindPressureVec = new Vec3D(0, 0, 0),
					initialWindPressure = 0d,
					gravityConstant = 9.81d,
					doStaticMinimization = true,
					maxWindSpeed = this.maxWindSpeed,
					maxDeltaWindSpeed = this.maxDeltaWindSpeed,
					minNumNewtonIterations = this.MinNumNewtonIterations,
					numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
					numWindReconstructionSteps2 = this.NumWindReconstructionSteps2,
					initialPositionNoiseMagnitude = 0.01,
					initialPositionNoiseRandomSeed = randomSeed
				});
			}
		}
		return reconstructionSetups;
	}

	private List<ReconstructionSetup> getReconstructionSetupsRealDataTest() {
		List<ReconstructionSetup> reconstructionSetups = new List<ReconstructionSetup>();
		//foreach(int n in new int[] {10, 20, 30, 40, 50, 70, 90, 110, 130}) { // TODO - Restore (or pick useful measurements anyways).
		foreach(int n in new int[] {30, 50, 80, 110}) {
			//foreach(int n in new int[] {110}) {
			foreach(int measurementsFileId in new int[] {1, 2, 3, 4}) { // Measurements 1, 2 and 4 best match the rest shape mesh (from measurements 2).
			//for(int measurementsFileId = 1; measurementsFileId <= 6; measurementsFileId++) {
				reconstructionSetups.Add(new ReconstructionSetup {
					sailStartConfigurationRelPath = "realDataTest/restshape.sailshapedata",
					sailMeasurementsRelPath = "realDataTest/markerConfig_" + measurementsFileId + ".measurements",
					resultsStorageRelPath = "realDataTest/markerConfig_" + measurementsFileId + "/n=" + n + ".results",
					kLength = 15621f,
					kArea = 8125f,
					kBend = 10f,
					shellThickness = 0.00025f,
					shellMaterialDensity = 920f,
					useFlatUndeformedBendState = true,
					numRestShapeSubdivisions = 1, // Slightly enlarge the rest shape mesh.
					isRealWorldDataTest = true,
					numRealWorldDataFarthestFirstMarkers = n,
					initialWindPressureVec = new Vec3D(0, 0, 0),
					initialWindPressure = 0d,
					gravityConstant = 9.81d,
					doStaticMinimization = true,
					maxWindSpeed = this.maxWindSpeed,
					maxDeltaWindSpeed = this.maxDeltaWindSpeed,
					minNumNewtonIterations = this.MinNumNewtonIterations,
					numWindReconstructionSteps1 = this.NumWindReconstructionSteps1,
					numWindReconstructionSteps2 = this.NumWindReconstructionSteps2
				});
			}
		}
		return reconstructionSetups;
	}

	private bool[] createVertexContraints() {
		//return MeshHelper.createOuterEdgeVertexContraints(this.vertexPositions, this.edges);

		// Constrain the mast and the boom vertices.
		bool[] verticesMovementConstraints = new bool[this.vertexPositions.Length];
		for(int i = 0; i < this.vertexPositions.Length; i++) {
			verticesMovementConstraints[i] = (this.vertexPositions[i].x <= 0d || this.vertexPositions[i].y <= 0d);
		}
		return verticesMovementConstraints;
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
		Dictionary<int, Dictionary<int, Edge>> vertexEdgeMap = MeshUtils.getVertexEdgeMap(this.sortedVertexTriangles, triangles);
		this.edges = MeshUtils.getEdges(vertexEdgeMap);

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
					new Vec3D(this.originalVertices[v2]) - new Vec3D(this.originalVertices[v1]),
					new Vec3D(this.originalVertices[v3]) - new Vec3D(this.originalVertices[v1]));
			double crossProdMag = crossProd.magnitude;

			// Store the triangle normal and area if they exist (i.e. if the triangle has a normal and therefore an area).
			if(double.IsNaN(crossProdMag)) {
				this.undeformedTriangleNormals[triangleId] = Vec3D.zero;
				this.undeformedTriangleAreas[triangleId] = 0;
			} else {
				this.undeformedTriangleNormals[triangleId] = crossProd / crossProdMag;

				// Get the triangle area using the edge lengths. These might have been adjusted, making the vertex positions useless to determine the undeformed edge length.
				// This uses Heron's formula, which should be equivalent to an area of crossProdMag / 2d.
				double a = vertexEdgeMap[v1][v2].undeformedLength;
				double b = vertexEdgeMap[v2][v3].undeformedLength;
				double c = vertexEdgeMap[v3][v1].undeformedLength;
				double p = (a + b + c) / 2d;
				this.undeformedTriangleAreas[triangleId] = Math.Sqrt(p * (p - a) * (p - b) * (p - c));
			}
		}

		print("Loaded mesh with " + mesh.vertices.Length + " vertices, " + this.edges.Count + " edges and " + (mesh.triangles.Length / 3) + " triangles.");
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

		// Recreate vertex constraints.
		this.verticesMovementConstraints = this.createVertexContraints();

		// Clear visualizer.
		this.visualizer.clear();

		// Reset reconstruction error parameters.
		this.stepCount = 0;
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
				this.meshRecorder.stop();
				this.meshRecorder.replay();
				print("Recording replay interrupted and restarted.");
			} else {
				this.meshRecorder.replay();
				print("Recording replay started.");
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
		try {
			this.simulationStep(Time.deltaTime * this.timeScale);
		} catch (Exception e) {
			print(e.ToString());
		}
	}

	private void simulationStep(double deltaTime) {

		// Load a new reconstruction setup if we are running automated reconstruction.
		while(this.ReconstructionStage == ReconstructionStage.DONE && this.reconstructionSetupsIndex + 1 < this.reconstructionSetups.Count) {
			this.reconstructionSetupsIndex++;
			print("Loading reconstruction setup " + (this.reconstructionSetupsIndex + 1) + "/" + this.reconstructionSetups.Count);
			ReconstructionSetup reconSetup = this.reconstructionSetups[this.reconstructionSetupsIndex];
			if(reconSetup.resultsStorageRelPath != null) {
				string resultsFilePath = storageBaseDirPath + "/Results/" + reconSetup.resultsStorageRelPath;
				if(File.Exists(resultsFilePath)) {
					print("Skipping reconstruction setup " + (this.reconstructionSetupsIndex + 1) + "/" + this.reconstructionSetups.Count
							+ " because a result file for this reconstruction already exists.");
					continue;
				}
			}
			this.loadSailConfiguration(SailConfiguration.loadFromFile(
					storageBaseDirPath + "/SailData/" + reconSetup.sailStartConfigurationRelPath), reconSetup.numRestShapeSubdivisions);
			this.measurements = SailMeasurements.loadFromFile(storageBaseDirPath + "/SailData/" + reconSetup.sailMeasurementsRelPath);
			this.kLength = reconSetup.kLength;
			this.kArea = reconSetup.kArea;
			this.kBend = reconSetup.kBend;
			this.shellThickness = reconSetup.shellThickness;
			this.shellMaterialDensity = reconSetup.shellMaterialDensity;
			this.useFlatUndeformedBendState = reconSetup.useFlatUndeformedBendState;
			this.windPressureVec = vecToVec(reconSetup.initialWindPressureVec);
			this.windPressure = reconSetup.initialWindPressure;
			this.gravityConstant = reconSetup.gravityConstant;
			this.DoStaticMinimization = reconSetup.doStaticMinimization;
			this.maxWindSpeed = reconSetup.maxWindSpeed;
			this.maxDeltaWindSpeed = reconSetup.maxDeltaWindSpeed;
			this.MinNumNewtonIterations = reconSetup.minNumNewtonIterations;
			this.NumWindReconstructionSteps1 = reconSetup.numWindReconstructionSteps1;
			this.NumWindReconstructionSteps2 = reconSetup.numWindReconstructionSteps2;

			// Enlarge measurements array if necessary if the mesh was subdivided.
			if(reconSetup.numRestShapeSubdivisions > 0 && this.measurements.measurements.Length < this.vertexPositions.Length) {
				Vec3D[] measurements = new Vec3D[this.vertexPositions.Length];
				for(int i = 0; i < this.measurements.measurements.Length; i++) {
					measurements[i] = this.measurements.measurements[i];
				}
				for(int i = this.measurements.measurements.Length; i < measurements.Length; i++) {
					measurements[i] = null;
				}
				this.measurements = new SailMeasurements(null, measurements);
			}

			// Perform real-world data specific operations. This code is written specifically for one mesh.
			if(reconSetup.isRealWorldDataTest) {

				// Generate constraints along the mast.
				for(int i = 0; i < this.vertexPositions.Length; i++) {
					this.verticesMovementConstraints[i] = false;
				}
				foreach(Edge edge in this.edges) {
					if(!edge.hasSideFlaps()) {
						if(this.vertexPositions[edge.ve1].x <= 0 || this.vertexPositions[edge.ve1].y <= 0.5) {
							this.verticesMovementConstraints[edge.ve1] = true;
						}
						if(this.vertexPositions[edge.ve2].x <= 0 || this.vertexPositions[edge.ve2].y <= 0.5) {
							this.verticesMovementConstraints[edge.ve2] = true;
						}
					}
				}

				// Visualize vertex constraints.
				for(int i = 0; i < this.vertexPositions.Length; i++) {
					if(this.verticesMovementConstraints[i]) {
						this.visualizer.visualizePoints(this.vertexPositions[i]);
					}
				}

				// Sample n measurements from the measurements array.
				Vec3D[] originalMeasurements = this.measurements.measurements;
				int n = reconSetup.numRealWorldDataFarthestFirstMarkers;
				Vec3D[] measurements = MeshUtils.generateFarthestPointSamplingSailMeasurements(originalMeasurements, n, this.verticesMovementConstraints);
				this.measurements = new SailMeasurements(null, measurements);

				// Use the measurements to estimate the reconstruction distance.
				this.realDataTruthMeasurements = new Vec3D[this.vertexPositions.Length];
				for(int i = 0; i < this.vertexPositions.Length; i++) {
					if(originalMeasurements[i] != null) {
						this.realDataTruthMeasurements[i] = originalMeasurements[i].clone();
					}
				}
			}

			int numMeasurements = 0;
			for(int i = 0; i < this.measurements.measurements.Length; i++) {
				if(this.measurements.measurements[i] != null) {
					numMeasurements++;
				}
			}
			print("Loaded reconstruction setup " + (this.reconstructionSetupsIndex + 1) + "/" + this.reconstructionSetups.Count
					+ "\n\tSail start configuration: " + reconSetup.sailStartConfigurationRelPath
					+ "\n\tMeasurements and truth configuration: " + reconSetup.sailMeasurementsRelPath
					+ "\n\tNumber of measurements: " + numMeasurements);

			// Remove measurements if they should be ignored.
			if(reconSetup.measurementsIgnoreSpheres != null) {
				Vec3D[] measurements = this.measurements.measurements;
				for(int i = 0; i < measurements.Length; i++) {
					if(measurements[i] != null) {
						foreach(Tuple<Vec3D, double> measurementsIgnoreSphere in reconSetup.measurementsIgnoreSpheres) {
							if((measurementsIgnoreSphere.Item1 - measurements[i]).squareMagnitude <= measurementsIgnoreSphere.Item2 * measurementsIgnoreSphere.Item2) {
								measurements[i] = null;
							}
						}
					}
				}
			}

			// Constrain vertices that have corresponding measurements. They will be smoothly snapped to these measurements later.
			int numVerticesConstrained = 0;
			for(int i = 0; i < this.vertexPositions.Length; i++) {
				if(this.measurements.measurements[i] != null) {
					this.verticesMovementConstraints[i] = true;
					numVerticesConstrained++;
				}
			}
			print("Measurement vertices locked: " + numVerticesConstrained);

			// Apply noise to all individual vertices of the loaded mesh to allow for giving slightly different initial sail configurations.
			if(reconSetup.initialPositionNoiseMagnitude != 0d) {
				System.Random random = new System.Random(reconSetup.initialPositionNoiseRandomSeed);
				for(int i = 0; i < this.vertexPositions.Length; i++) {
					if(!this.verticesMovementConstraints[i]) {
						Vec3D randomVec = new Vec3D(random.Next(-10000, 10000), random.Next(-10000, 10000), random.Next(-10000, 10000));
						randomVec *= reconSetup.initialPositionNoiseMagnitude / randomVec.magnitude;
						this.vertexPositions[i].add(randomVec);
					}
				}
			}

			this.ReconstructionStage = ReconstructionStage.MOVING_VERTICES_TO_MEASUREMENTS;
			this.stepCount = 0;
			return;
		}
		if(this.ReconstructionStage == ReconstructionStage.DONE) {
			print("Reconstruction(s) done. Pausing simulation.");
			this.doUpdate = false;
			return;
		}

		// Get the mesh.
		Mesh mesh = this.getMesh();
		int[] triangles = mesh.triangles;

		// Compute triangle normals and areas.
		this.recalcTriangleNormalsAndAreas(triangles, this.vertexPositions);

		// Clear previous visualizations.
		this.visualizer.clear();

		// Update the vertices using the chosen time stepping method.
		switch(this.timeSteppingMethod) {
			case TimeSteppingMethod.GRADIENT_DESCENT: {
				this.doGradientDescentStep();
				break;
			}
			case TimeSteppingMethod.EXPLICIT: {
				this.doExplititIntegrationStep(deltaTime);
				break;
			}
			case TimeSteppingMethod.OPTIMIZATION_INTEGRATOR: {
				this.doOptimizationIntegratorStep(deltaTime);
				break;
			}
		}

		// Move the vertices towards their corresponding measurements. This is the more smooth version of snapping vertices to measurements.
		if(this.ReconstructionStage == ReconstructionStage.MOVING_VERTICES_TO_MEASUREMENTS) {

			// Perform steps.
			bool allTargetsReached = true;
			for(int i = 0; i < this.vertexPositions.Length; i++) {
				if(this.measurements.measurements[i] != null && this.verticesMovementConstraints[i]) {

					// Get limited step from vertex towards its corresponding measurement.
					Vec3D step = this.measurements.measurements[i] - this.vertexPositions[i];
					if(step.x == 0 && step.y == 0 && step.z == 0) {
						continue; // Already at target.
					}
					double stepMag = step.magnitude;
					if(stepMag > MaxVertexMoveToMeasurementsStep) {
						step.div(stepMag).mul(MaxVertexMoveToMeasurementsStep);
					}

					// Apply step, or snap to the measurement position if the step might be non-zero due to double precision.
					if(stepMag < 0.00001) {
						this.vertexPositions[i] = this.measurements.measurements[i].clone();
					} else {
						this.vertexPositions[i].add(step);
						allTargetsReached = false;
					}
				}
			}

			// Continue with the next reconstruction stage when all vertices have reached their corresponding measurements.
			if(allTargetsReached) {
				this.ReconstructionStage = ReconstructionStage.RECONSTRUCT_WIND;
				this.stepCount = 0;
			}
		}

		// Update and print best reconstruction error. The step method is responsible for the wind reconstruction itself.
		// Also update the max delta wind speed to at first move quickly and later more precise towards the target.
		if(this.ReconstructionStage == ReconstructionStage.RECONSTRUCT_WIND && this.measurements != null) {
			this.stepCount++;

			double measurementsSquaredError = this.getSquaredMeasurementsError(this.vertexPositions);
			double reconstructionDistance = this.getReconstructionDistance(this.vertexPositions);
			double averageReconstructionDistance = this.getAverageReconstructionDistance(this.vertexPositions);
			double maxReconstructionDistance = this.getMaxReconstructionDistance(this.vertexPositions);

			print("Step: " + this.stepCount + ", Squared measurements error: " + measurementsSquaredError
					+ ", average reconstruction distance: " + averageReconstructionDistance
					+ ", max reconstruction distance: " + maxReconstructionDistance);
			this.avgReconDistTextObj.text = this.avgReconDistTextObj.text.Split(':')[0] + ": " + averageReconstructionDistance + "m";
			this.maxReconDistTextObj.text = this.maxReconDistTextObj.text.Split(':')[0] + ": " + maxReconstructionDistance + "m";
			
			// Set max delta wind speed.
			if(this.stepCount >= this.NumWindReconstructionSteps1) {
				this.maxDeltaWindSpeed = 1d;
			} else {
				this.maxDeltaWindSpeed = 1d + 49d * (1d - this.stepCount / (double) this.NumWindReconstructionSteps1); // Linear from 50 to 1.
			}
			
			// Stop simulation when the desired amount of reconstruction steps for this phase has been reached.
			if(this.stepCount >= this.NumWindReconstructionSteps1 + this.NumWindReconstructionSteps2) {
				print("Step threshold has been reached. Wind reconstruction step complete.");
				if(this.reconstructionSetupsIndex >= this.reconstructionSetups.Count || this.reconstructionSetupsIndex == -1) {
					this.doUpdate = false;
				} else {

					// We are running automated reconstruction. Store the results and continue.
					print("Reconstruction complete.");
					ReconstructionSetup reconSetup = this.reconstructionSetups[this.reconstructionSetupsIndex];
					if(reconSetup.resultsStorageRelPath != null) {
						string filePath = storageBaseDirPath + "/Results/" + reconSetup.resultsStorageRelPath;
						string dirPath = Path.GetDirectoryName(filePath);
						print("Storing reconstruction results to file: " + filePath);
						if(!Directory.Exists(dirPath)) {
							Directory.CreateDirectory(dirPath);
						}
						File.WriteAllText(filePath, "reconstructionDistance = " + reconstructionDistance
								+ "\naverageReconstructionDistance = " + averageReconstructionDistance
								+ "\nmaxReconstructionDistance = " + maxReconstructionDistance
								+ "\nnumMeasurements = " + this.measurements.getNumMeasurements());
					} else {
						print("Not storing reconstruction results because no result file was provided. Results:"
								+ "\n\treconstructionDistance = " + reconstructionDistance
								+ "\n\taverageReconstructionDistance = " + averageReconstructionDistance
								+ "\n\tmaxReconstructionDistance = " + maxReconstructionDistance
								+ "\n\tnumMeasurements = " + this.measurements.getNumMeasurements());
					}
				}
				this.ReconstructionStage = ReconstructionStage.DONE;
				this.stepCount = 0;
			}
		} else {
			this.avgReconDistTextObj.text = this.avgReconDistTextObj.text.Split(':')[0] + ": -";
			this.maxReconDistTextObj.text = this.maxReconDistTextObj.text.Split(':')[0] + ": -";
		}

		// Update visualizations.
		switch(this.visualizationType) {
			case VisualizationType.MEASUREMENTS: {
				if(this.measurements != null) {
					Vec3D[] measurements = this.measurements.measurements;
					List<Vec3D> measurementPositions = new List<Vec3D>();
					for(int i = 0; i < measurements.Length; i++) {
						if(measurements[i] != null) {
							measurementPositions.Add(measurements[i]);
						}
					}
					this.visualizer.visualizePoints(measurementPositions.ToArray(), 0.1f);
				}
				break;
			}
			case VisualizationType.VERTEX_POSITIONS: {
				this.visualizer.visualizePoints(this.vertexPositions, 0.1f);
				break;
			}
			case VisualizationType.VERTEX_VELOCITIES: {
				this.visualizer.visualizeVectors(this.vertexPositions, this.vertexVelocities, 0.2f);
				break;
			}
			case VisualizationType.VERTEX_MEASUREMENT_DIFF: {
				if(this.measurements != null) {
					Vec3D[] measurements = this.measurements.measurements;
					VecD measurementsError = new VecD(this.vertexPositions.Length * 3);
					for(int vertexInd = 0; vertexInd < this.vertexPositions.Length; vertexInd++) {
						if(measurements[vertexInd] != null) {
							for(int coord = 0; coord < 3; coord++) {
								measurementsError[3 * vertexInd + coord] = measurements[vertexInd][coord] - this.vertexPositions[vertexInd][coord];
							}
						}
					}
					this.visualizer.visualizeVectors(this.vertexPositions, measurementsError, 1f);
				}
				break;
			}
			case VisualizationType.WIND_FORCE: {
				this.recalcTriangleNormalsAndAreas(triangles, this.vertexPositions);
				VecD windForce = this.getVertexWindForce(triangles, this.vertexPositions, new Vec3D(this.windPressureVec), this.windPressure);
				this.visualizer.visualizeVectors(this.vertexPositions, windForce);
				break;
			}
			case VisualizationType.ENERGY_GRADIENT: {
				VecD energyGradient = getSystemEnergyGradient(triangles, this.vertexPositions);
				this.visualizer.visualizeVectors(this.vertexPositions, energyGradient);
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
		VecD vertexWindForce = this.getVertexWindForce(triangles, this.vertexPositions, new Vec3D(this.windPressureVec), this.windPressure);
		VecD vertexCoordMasses = this.getVertexCoordinateMasses();
		VecD gravityForce = this.getVertexGravityForce(vertexCoordMasses);
		VecD step = this.kGradientDescent * (new VecD(vertexWindForce).add(gravityForce).sub(vertexEnergyGradient));
		for(int vertexInd = 0; vertexInd < this.vertexPositions.Length; vertexInd++) {
			if(!this.verticesMovementConstraints[vertexInd]) {
				Vec3D stepVec = new Vec3D(step[3 * vertexInd], step[3 * vertexInd + 1], step[3 * vertexInd + 2]);

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
				double stepVecMag = stepVec.magnitude;
				if(stepVecMag > this.maxGradientDescentStep) {
					stepVec = stepVec.mul(this.maxGradientDescentStep / stepVecMag);
				}

				// Apply step.
				this.vertexPositions[vertexInd].add(stepVec);
			}
		}
		this.updateMesh();
	}

	private void doExplititIntegrationStep(double deltaTime) {
		
		// Calculate vertex forces.
		int[] triangles = this.getMesh().triangles;
		VecD vertexEnergyGradient = this.getSystemEnergyGradient(triangles, this.vertexPositions);
		VecD vertexWindForce = this.getVertexWindForce(triangles, this.vertexPositions, new Vec3D(this.windPressureVec), this.windPressure);
		VecD vertexCoordMasses = this.getVertexCoordinateMasses();
		VecD gravityForce = this.getVertexGravityForce(vertexCoordMasses);

		// Calculate new vertex accelerations.
		this.vertexAccelerations = new VecD(vertexWindForce).add(gravityForce).sub(vertexEnergyGradient).divideElementWise(vertexCoordMasses);

		// Update vertex velocities.
		this.vertexVelocities.add(new VecD(this.vertexAccelerations).mul(deltaTime));

		// Apply velocity-based damping (simple heuristic where we assume that a constant percentage of the velocity is lost due to friction).
		this.vertexVelocities.mul(this.dampingFactor);

		// Update vertex positions and apply movement constraints.
		for(int i = 0; i < this.vertexPositions.Length; i++) {
			if(this.verticesMovementConstraints[i]) {
				this.vertexVelocities[3 * i] = 0;
				this.vertexVelocities[3 * i + 1] = 0;
				this.vertexVelocities[3 * i + 2] = 0;
				this.vertexAccelerations[3 * i] = 0;
				this.vertexAccelerations[3 * i + 1] = 0;
				this.vertexAccelerations[3 * i + 2] = 0;
				continue;
			}
			for(int coord = 0; coord < 3; coord++) {
				this.vertexPositions[i][coord] += deltaTime * this.vertexVelocities[3 * i + coord]
						+ deltaTime * deltaTime * this.vertexAccelerations[3 * i + coord] / 2d;
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

		// Get vertex masses.
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

		// Update wind force if measurements are set.
		Vec3D windPressureVec = new Vec3D(this.windPressureVec);
		if(this.measurements != null && this.ReconstructionStage == ReconstructionStage.RECONSTRUCT_WIND) {

			// Get system-wide energy gradient.
			this.recalcTriangleNormalsAndAreas(triangles, this.vertexPositions);
			VecD energyGradient = this.getSystemEnergyGradient(triangles, this.vertexPositions);

			// Get gravity force.
			VecD gravityForce = this.getVertexGravityForce(vertexCoordMasses);

			// Compute the force that should be put onto the sail by the wind in order to result in static force equilibrium.
			VecD virtMeasurementsErrorForce = new VecD(energyGradient).sub(gravityForce);

			// Compute the wind pressure vector that best represents the desired force.
			Vec3D newWindPressureVec = this.getWindPressureVec(triangles, this.vertexPositions, virtMeasurementsErrorForce);

			// Limit new wind pressure vector magnitude.
			double newWindPressureVecMag = newWindPressureVec.magnitude;
			if(newWindPressureVecMag > this.maxWindSpeed) {
				newWindPressureVec.div(newWindPressureVecMag).mul(this.maxWindSpeed);
			}

			// Compute and limit delta wind pressure vector magnitude.
			Vec3D deltaWindPressureVec = newWindPressureVec - windPressureVec;
			double deltaWindPressureVecMag = deltaWindPressureVec.magnitude;
			if(deltaWindPressureVecMag > this.maxDeltaWindSpeed) {
				deltaWindPressureVec.div(deltaWindPressureVecMag).mul(this.maxDeltaWindSpeed);
			}

			// Apply delta wind velocity.
			windPressureVec += deltaWindPressureVec;
		}

		// Perform Newton's method, setting steps until the termination criterion has been met.
		VecD vertexPositionsFlat = new VecD(this.vertexPositions);
		Vec3D[] newVertexPositions = new Vec3D[numVertices];
		for(int i = 0; i < numVertices; i++) {

			// Reset vertex velocity for constrained vertices. This is necessary when constraints are dynamically added while the mesh had a velocity.
			if(this.verticesMovementConstraints[i]) {
				for(int coord = 0; coord < 3; coord++) {
					this.vertexVelocities[i * 3 + coord] = 0d;
				}
			}

			// Set initial guess: pos + deltaTime * velocity.
			newVertexPositions[i] = new Vec3D();
			for(int coord = 0; coord < 3; coord++) {
				newVertexPositions[i][coord] = this.vertexPositions[i][coord] + deltaTime * this.vertexVelocities[i * 3 + coord];
			}
		}
		int iteration = 0;
		while(true) {

			// Get system-wide energy gradient.
			this.recalcTriangleNormalsAndAreas(triangles, newVertexPositions);
			VecD energyGradient = this.getSystemEnergyGradient(triangles, newVertexPositions);

			// Get wind force.
			VecD windForce = this.getVertexWindForce(triangles, newVertexPositions, windPressureVec, this.windPressure);

			// Get gravity force.
			VecD gravityForce = this.getVertexGravityForce(vertexCoordMasses);

			// Get measurements force.
			VecD measurementsForce = this.getVertexMeasurementSpringForce(this.vertexPositions, this.measurements);

			// Get damping force.
			VecD nextVertexVelocities = new VecD(newVertexPositions).sub(vertexPositionsFlat).div(deltaTime);
			//VecD dampingForce = -this.dampingConstant * (energyHess * nextVertexVelocities);
			VecD dampingForce = (this.DoStaticMinimization
				? new VecD(numVertices * 3)
				: -this.dampingConstant * nextVertexVelocities); // Damping force using an identity matrix as energy Hessian.

			// Set forces to zero for constrained vertices.
			// This causes the E gradient to be zero for them as well, causing it not to get in the way of the minimization problem.
			for(int i = 0; i < numVertices; i++) {
				if(this.verticesMovementConstraints[i]) {
					for(int coord = 0; coord < 3; coord++) {
						energyGradient[3 * i + coord] = 0;
						windForce[3 * i + coord] = 0;
						gravityForce[3 * i + coord] = 0;
						measurementsForce[3 * i + coord] = 0;
						dampingForce[3 * i + coord] = 0;
					}
				}
			}

			// Get E gradient.
			VecD velocityPreservingTermPart = (this.DoStaticMinimization
					? new VecD(numVertices * 3)
					: new VecD(this.vertexVelocities).multiplyElementWise(vertexCoordMasses).div(deltaTime));
			if(!this.DoStaticMinimization) {
				for(int i = 0; i < numVertices * 3; i++) {

					// Limit the damping force magnitude to the velocity preserving force magnitude.
					double velocityPreservingTermPartMag = Math.Abs(velocityPreservingTermPart[i]);
					if(Math.Abs(dampingForce[i]) > velocityPreservingTermPartMag) {
						dampingForce[i] = Math.Sign(dampingForce[i]) * velocityPreservingTermPartMag;
					}

					// Set the damping force to 0 when it does not go against the velocity preserving force.
					if(Math.Sign(dampingForce[i]) == Math.Sign(velocityPreservingTermPart[i])) {
						dampingForce[i] = 0d;
					}
				}
			}
			VecD eGradient = new VecD(newVertexPositions).sub(vertexPositionsFlat)
					.multiplyElementWise(vertexCoordMasses).div(deltaTimeSquare).sub(velocityPreservingTermPart)
					.add(energyGradient).sub(windForce).sub(gravityForce).sub(measurementsForce).sub(dampingForce);

			// Terminate when the termination criterion has been met.
			double eGradientMagnitude = eGradient.magnitude;
			print(stopWatch.ElapsedMilliseconds + "ms: E gradient magnitude: " + eGradientMagnitude + " (threshold: " + terminationThreshold
					+ ", iteration: " + iteration + ")");
			if(eGradientMagnitude < terminationThreshold && iteration > this.MinNumNewtonIterations) {
				print(stopWatch.ElapsedMilliseconds + "ms: Finished on iteration: " + iteration + ".");
				break;
			}

			// Limit amount of time spent to prevent endless loops.
			iteration++;
			if(stopWatch.ElapsedMilliseconds > this.MaxNewtonsMethodLoopTimeMS) {
				print(stopWatch.ElapsedMilliseconds + "ms: Maximum time reached in Optimization Integrator update after "
						+ iteration + " iterations. Taking non-optimal step with E gradient magnitude: "
						+ eGradientMagnitude + " (threshold: " + terminationThreshold + ")");
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
					vertexCoordMasses, windPressureVec, this.windPressure, energyHess, deltaTime);
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
						vertexCoordMasses, windPressureVec, this.windPressure, energyHess, deltaTime);
				if(newE <= bestE) {
				//if(newE <= e + c * alpha * VecD.dot(eGradient, step)) {
					
					// Alpha is suitable, but there might be a higher value of alpha that is still suitable. Increase alpha and store the current best alpha.
					bestAlpha = alpha;
					bestE = newE;
					alpha *= 2;

					// Limit alpha. This prevents stable zero-step configurations from resulting in an infinite alpha.
					if(alpha >= 10) {
						alpha = 10;
						print(stopWatch.ElapsedMilliseconds + "ms: Alpha is getting too large. Setting alpha: " + alpha);
						break;
					}
				} else if(!double.IsNaN(bestAlpha)) {

					// A best alpha was set, but a higher value of alpha didn't make it. Return the best value.
					alpha = bestAlpha;
					print(stopWatch.ElapsedMilliseconds + "ms: Terminating with alpha: " + alpha
							+ " after " + (stopWatch.ElapsedMilliseconds - lineSearchLoopStartTime) + "ms (bestE: " + bestE + ").");
					break;
				} else {

					// Alpha is too high to be suitable. Decrease alpha.
					alpha /= 2d;

					// Just take the step if alpha gets too small. We're either at perfect reconstruction or close to it.
					// In the case that this does occasionally introduce an error, next steps will still go towards a valid solution.
					if(alpha <= this.MinLineSearchAlpha) {
						alpha = 0d;
						print(stopWatch.ElapsedMilliseconds + "ms: Alpha is getting too small. Setting alpha: " + alpha);
						goto newtLoopEnd;
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
		newtLoopEnd:

		// Update vertex velocity.
		this.vertexVelocities = new VecD(newVertexPositions).sub(vertexPositionsFlat).div(deltaTime);

		// Pause simulation when the maximum vertex velocity is below the set threshold.
		double maxVertexVelocity = 0d;
		double averageVertexVelocity = 0d;
		for(int i = 0; i < numVertices; i++) {
			int baseInd = 3 * i;
			double vertexVelocity = new Vec3D(this.vertexVelocities[baseInd], this.vertexVelocities[baseInd + 1], this.vertexVelocities[baseInd + 2]).magnitude;
			averageVertexVelocity += vertexVelocity;
			if(vertexVelocity > maxVertexVelocity) {
				maxVertexVelocity = vertexVelocity;
			}
		}
		averageVertexVelocity /= numVertices;
		print("Vertex velocity after step, before direct velocity damping: Max = " + maxVertexVelocity + ", average = " + averageVertexVelocity);
		if(maxVertexVelocity < this.VelocityTerminationThreshold
				&& (this.ReconstructionStage == ReconstructionStage.DISABLED || this.ReconstructionStage == ReconstructionStage.RECONSTRUCT_WIND)) {
			this.doUpdate = false;
			print("Simulation finished. Maximum vertex velocity is below the velocity termination threshold: "
					+ maxVertexVelocity + " m/s < ." + this.VelocityTerminationThreshold + " m/s.");
		}

		// Apply direct velocity damping.
		if(this.directVelocityDampingFactor != 1d) {
			this.vertexVelocities.mul(this.directVelocityDampingFactor);
		}

		// Update vertex positions.
		this.vertexPositions = newVertexPositions;

		// Update wind velocity.
		this.windPressureVec = new Vector3((float) windPressureVec[0], (float) windPressureVec[1], (float) windPressureVec[2]);

		// Update mesh.
		this.updateMesh();

		// Visualize step.
		if(this.visualizationType == VisualizationType.STEP) {
			this.visualizer.visualizeVectors(this.vertexPositions, new VecD(newVertexPositions).sub(vertexPositionsFlat), 10f);
		}
	}

	private Vec3D getWindPressureVec(int[] triangles, Vec3D[] vertexPositions, VecD virtMeasurementsErrorForce) {
		int numVertices = vertexPositions.Length;

		// Set the virtual measurements error force to 0 for points without measurements.
		// This is necessary to not cause the wind force to attempt to cancel acceleration and therefore keep non-measurements vertices in place.
		int numMeasurements = 0;
		for(int i = 0; i < numVertices; i++) {
			if(this.measurements == null || this.measurements.measurements[i] == null) {
				virtMeasurementsErrorForce[i * 3] = 0d;
				virtMeasurementsErrorForce[i * 3 + 1] = 0d;
				virtMeasurementsErrorForce[i * 3 + 2] = 0d;
			} else {
				numMeasurements++;
			}
		}
		if(numMeasurements == 0) {
			return new Vec3D();
		}

		// Determine the wind pressure direction.
		// For the direction, the direction of the sum of virtual measurements error forces is used.
		// For the magnitude, we first determine the wind force produced by the unit wind direction and then scale it to on average match the virtual error forces.
		Vec3D windPressureVec = new Vec3D();
		double averageVirtMeasurementsErrorForceMagnitude = 0d;
		for(int i = 0; i < numVertices; i++) {
			Vec3D virtMeasurementErrorForce = new Vec3D(
					virtMeasurementsErrorForce[3 * i], virtMeasurementsErrorForce[3 * i + 1], virtMeasurementsErrorForce[3 * i + 2]);
			windPressureVec += virtMeasurementErrorForce;
			averageVirtMeasurementsErrorForceMagnitude += virtMeasurementErrorForce.magnitude;
		}
		double magnitude = windPressureVec.magnitude;
		if(magnitude == 0d) {
			return new Vec3D();
		}
		windPressureVec /= magnitude;
		averageVirtMeasurementsErrorForceMagnitude /= numMeasurements;

		VecD deltaVertexUnitWindForce = this.getVertexWindForce(triangles, vertexPositions, windPressureVec, 0d);
		double averageUnitDeltaWindForceMagnitude = 0d;
		for(int i = 0; i < numVertices; i++) {
			Vec3D deltaWindForce = new Vec3D(
					deltaVertexUnitWindForce[3 * i], deltaVertexUnitWindForce[3 * i + 1], deltaVertexUnitWindForce[3 * i + 2]);
			averageUnitDeltaWindForceMagnitude += deltaWindForce.magnitude;
		}
		averageUnitDeltaWindForceMagnitude /= numVertices;
		double deltaWindForceMagnitude = averageVirtMeasurementsErrorForceMagnitude / averageUnitDeltaWindForceMagnitude;
		windPressureVec.mul(deltaWindForceMagnitude);

		return windPressureVec;
	}

	/*
	 * Gets the sum of per-coordinate squared errors between the current vertex positions and the vertex positions at the time the measurements were taken.
	 * Returns 0 if no measurement is set.
	 */
	private double getSquaredMeasurementsError(Vec3D[] vertexPositions) {
		double error = 0;
		if(this.measurements != null) {
			Vec3D[] measurements = this.measurements.measurements;
			for(int i = 0; i < vertexPositions.Length; i++) {
				if(measurements[i] != null) {
					Vec3D diff = measurements[i] - vertexPositions[i];
					error += VecD.dot(diff, diff);
				}
			}
		}
		return error;
	}

	/*
	 * Gets the sum of distances between the current vertex positions and the vertex positions at the time the measurements were taken.
	 * Returns 0 if no measurement is set.
	 */
	private double getReconstructionDistance(Vec3D[] vertexPositions) {
		double error = 0;
		if(this.realDataTruthMeasurements == null) {
			if(this.measurements != null && this.measurements.vertexPositions != null) {
				Vec3D[] correctPositions = this.measurements.vertexPositions;
				for(int i = 0; i < vertexPositions.Length; i++) {
					Vec3D diff = correctPositions[i] - vertexPositions[i];
					error += diff.magnitude;
				}
			}
		} else {
			Vec3D[] correctPositions = this.realDataTruthMeasurements;
			for(int i = 0; i < vertexPositions.Length; i++) {
				if(correctPositions[i] != null) {
					Vec3D diff = correctPositions[i] - vertexPositions[i];
					error += diff.magnitude;
				}
			}
		}
		return error;
	}

	/*
	 * Gets the average distance between the current vertex positions and the vertex positions at the time the measurements were taken.
	 * Extrapolates the average distance to unknown vertex positions if vertex position truths are unknown.
	 * Returns 0 if no measurement is set.
	 */
	private double getAverageReconstructionDistance(Vec3D[] vertexPositions) {
		if(this.realDataTruthMeasurements == null) {
			double error = 0;
			if(this.measurements != null && this.measurements.vertexPositions != null) {
				Vec3D[] correctPositions = this.measurements.vertexPositions;
				for(int i = 0; i < vertexPositions.Length; i++) {
					Vec3D diff = correctPositions[i] - vertexPositions[i];
					error += diff.magnitude;
				}
			}
			return error / vertexPositions.Length;
		} else {
			int numConstraints = 0;
			int numConstrainedMeasurements = 0;
			int numUnconstrainedMeasurements = 0;
			double unconstrainedError = 0;
			double constrainedError = 0;
			Vec3D[] correctPositions = this.realDataTruthMeasurements;
			for(int i = 0; i < vertexPositions.Length; i++) {
				if(this.verticesMovementConstraints[i]) {
					numConstraints++;
				}
				if(correctPositions[i] != null) {
					double diffMag = (correctPositions[i] - vertexPositions[i]).magnitude;
					if(this.verticesMovementConstraints[i]) {
						constrainedError += diffMag;
						numConstrainedMeasurements++;
					} else {
						unconstrainedError += diffMag;
						numUnconstrainedMeasurements++;
					}
				}
			}
			double constrainedMeasurementFactor = (double) numConstrainedMeasurements / (double) numConstraints;
			return (1d - constrainedMeasurementFactor) * unconstrainedError / numUnconstrainedMeasurements
					+ constrainedMeasurementFactor * constrainedError / numConstrainedMeasurements;
		}
	}

	/*
	 * Gets the maximum of the distances between the current vertex positions and the vertex positions at the time the measurements were taken.
	 * Returns 0 if no measurement is set.
	 */
	private double getMaxReconstructionDistance(Vec3D[] vertexPositions) {
		if(this.realDataTruthMeasurements == null) {
			if(this.measurements == null || this.measurements.vertexPositions == null) {
				return 0d;
			}
			double maxDiffMag = 0d;
			Vec3D[] correctPositions = this.measurements.vertexPositions;
			for(int i = 0; i < vertexPositions.Length; i++) {
				Vec3D diff = correctPositions[i] - vertexPositions[i];
				double diffMag = diff.magnitude;
				if(diffMag > maxDiffMag) {
					maxDiffMag = diffMag;
				}
			}
			return maxDiffMag;
		} else {
			double maxDiffMag = 0d;
			Vec3D[] correctPositions = this.realDataTruthMeasurements;
			for(int i = 0; i < vertexPositions.Length; i++) {
				if(correctPositions[i] != null) {
					Vec3D diff = correctPositions[i] - vertexPositions[i];
					double diffMag = diff.magnitude;
					if(diffMag > maxDiffMag) {
						maxDiffMag = diffMag;
					}
				}
			}
			return maxDiffMag;
		}
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

		// The bending energy Hessian consists of 16 (3x3) parts that have to be inserted into the matrix.
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
		
		// The area energy Hessian consists of 9 (3x3) parts that have to be inserted into the matrix.
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

		// Sort triplets first on row and then on column.
		Array.Sort(energyHessTriplets, (o1, o2) =>
				(o1 == null ? (o2 == null ? 0 : -1) : o2 == null ? 1
						: (o1.Item1 > o2.Item1 ? 1 : (o1.Item1 < o2.Item1 ? -1
								: (o1.Item2 > o2.Item2 ? 1 : (o1.Item2 < o2.Item2 ? -1 : 0))))));

		// Sum up triplets at the same index and move them to the left, resulting in trailing nulls.
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

		// Assemble and return system-wide energy Hessian.
		MathNet.Numerics.LinearAlgebra.Double.SparseMatrix energyHess = new MathNet.Numerics.LinearAlgebra.Double.SparseMatrix(numVertices * 3, numVertices * 3);
		foreach(Tuple<int, int, double> triplet in energyHessTriplets) {
			if(triplet == null) {
				break; // All nulls are trailing, so there won't be more data.
			}
			energyHess[triplet.Item1, triplet.Item2] = triplet.Item3;
		}
		return energyHess;
	}

	private double getE(int[] triangles, VecD vertexPositionsFlat, Vec3D[] newVertexPositions, VecD vertexCoordMasses,
			Vec3D windPressureVec, double windPressure, MathNet.Numerics.LinearAlgebra.Double.SparseMatrix energyHess, double deltaTime) {
		int numVertices = newVertexPositions.Length;

		VecD windForce = this.getVertexWindForce(triangles, newVertexPositions, windPressureVec, windPressure);
		for(int i = 0; i < numVertices; i++) {
			if(this.verticesMovementConstraints[i]) {
				for(int coord = 0; coord < 3; coord++) {
					windForce[3 * i + coord] = 0;
				}
			}
		}

		double measurementSpringEnergy = this.getVertexMeasurementSpringEnergy(newVertexPositions, this.measurements);

		VecD deltaVertexPositions = new VecD(newVertexPositions).sub(vertexPositionsFlat);
		double systemEnergy = this.getSystemEnergy(triangles, newVertexPositions);
		double gravityWork = 0;
		double windWork = -VecD.dot(windForce, deltaVertexPositions); // Approximation: Consider triangle normals and areas constant.
		//VecD dampingForce = -(this.dampingConstant / deltaTime) * (energyHess * deltaVertexPositions);
		VecD dampingForce = (this.DoStaticMinimization
				? new VecD(numVertices * 3)
				: -(this.dampingConstant / deltaTime) * deltaVertexPositions); // Damping force using an identity matrix as energy Hessian.
		double dampingWork = -VecD.dot(dampingForce, deltaVertexPositions) / 2d;
		for(int i = 0; i < numVertices; i++) {
			if(!this.verticesMovementConstraints[i]) { // Not necessary when comparing energy, but lets consider them part of the outside world.
				int yCoordIndex = 3 * i + 1;
				gravityWork += vertexCoordMasses[yCoordIndex] * this.gravityConstant * deltaVertexPositions[yCoordIndex];
			}
		}
		VecD squareTerm = new VecD(newVertexPositions).sub(vertexPositionsFlat);
		if(!this.DoStaticMinimization) {
			squareTerm.sub(deltaTime * this.vertexVelocities);
		}
		return VecD.dot(VecD.multiplyElementWise(squareTerm, squareTerm), vertexCoordMasses) / (2d * deltaTime * deltaTime)
				+ systemEnergy + gravityWork + measurementSpringEnergy + windWork + dampingWork;
	}

	/**
	 * Calculate lumped vertex mass per vertex coordinate. This will return a vector in format: {m_v1x, m_v1y, m_v1z, m_v2x, ...}.
	 * Note that m_vix = m_viy = m_viz = m_vi for any vertex i, and that the masses are defined by the undeformed triangle areas.
	 */
	private VecD getVertexCoordinateMasses() {
		int numVertices = this.vertexPositions.Length;
		VecD vertexCoordinateMasses = new VecD(numVertices * 3);
		for(int i = 0; i < numVertices; i++) {
			double vertexArea = 0d;
			foreach(int triangleId in this.sortedVertexTriangles[i]) {
				vertexArea += this.undeformedTriangleAreas[triangleId];
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
	 * Calculates the wind forces acting on the vertices.
	 * These forces are the sum of:
	 * - The projection of the wind pressure vector onto each triangle normal, multiplied by the triangle area.
	 * - The wind vector scalar multiplied by each triangle area.
	 * Each triangle force is then distributed evenly over the vertices defining that triangle.
	 * Returns the wind force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
	 */
	private VecD getVertexWindForce(int[] triangles, Vec3D[] vertices, Vec3D windPressureVec, double windPressure) {

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
			double windNoiseMagnitudeFactor = 1d;
			if(this.windNoiseFuncMagnitudeSlope != 0d) {
				double xTriangleMiddle = (vertices[v1].x + vertices[v2].x + vertices[v3].x) / 3d;
				double sailWidth = 3.5d;
				double c = sailWidth / 2;
				windNoiseMagnitudeFactor = 1d + this.windNoiseFuncMagnitudeSlope * (xTriangleMiddle / c - 1d);
			}
			VecD triangleVertexWindForce = ((VecD.dot(windPressureVec, triangleNormal) + windPressure) * windNoiseMagnitudeFactor * triangleArea / 3d) * triangleNormal;

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

	/*
	 * Calculates the measurements spring force acting on each vertex. This is an artificial force that should only be used during reconstruction.
	 * Measurement spring forces are equivalent to placing springs between measurements and their corresponding vertices.
	 * springForce[i] = kMeasurementSprings * (vertexPositions[i] - measurements[i]) for all coordinates of measurement i if it exists.
	 * Returns the measurement spring force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
	 */
	private VecD getVertexMeasurementSpringForce(Vec3D[] vertexPositions, SailMeasurements measurements) {
		if(measurements == null) {
			return new VecD(3 * vertexPositions.Length);
		}
		return this.getVertexMeasurementSpringForce(vertexPositions, measurements.measurements);
	}

	/*
	 * Calculates the measurements spring force acting on each vertex. This is an artificial force that should only be used during reconstruction.
	 * Measurement spring forces are equivalent to placing springs between measurements and their corresponding vertices.
	 * springForce[i] = kMeasurementSprings * (measurements[i] - vertexPositions[i]) for all coordinates of measurement i if it exists.
	 * Returns the measurement spring force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
	 */
	private VecD getVertexMeasurementSpringForce(Vec3D[] vertexPositions, Vec3D[] measurements) {
		VecD measurementSpringForces = new VecD(3 * vertexPositions.Length);
		if(this.kMeasurementSprings != 0d) {
			for(int vertexInd = 0; vertexInd < vertexPositions.Length; vertexInd++) {
				if(measurements[vertexInd] != null) {
					for(int coord = 0; coord < 3; coord++) {
						measurementSpringForces[3 * vertexInd + coord] =
							this.kMeasurementSprings * (measurements[vertexInd][coord] - vertexPositions[vertexInd][coord]);
					}
				}
			}
		}
		return measurementSpringForces;
	}

	/*
	 * Calculates the total measurements spring energy (sum of measurement spring energies). This is an artificial energy that should only be used during reconstruction.
	 * springEnergy[i] = 1/2 * kMeasurementSprings * (vertexPositions[i] - measurements[i])^2 for all coordinates of measurement i if it exists.
	 * Returns the measurement spring force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
	 */
	private double getVertexMeasurementSpringEnergy(Vec3D[] vertexPositions, SailMeasurements measurements) {
		if(measurements == null) {
			return 0d;
		}
		return this.getVertexMeasurementSpringEnergy(vertexPositions, measurements.measurements);
	}

	/*
	 * Calculates the total measurements spring energy (sum of measurement spring energies). This is an artificial energy that should only be used during reconstruction.
	 * springEnergy[i] = 1/2 * kMeasurementSprings * (vertexPositions[i] - measurements[i])^2 for all coordinates of measurement i if it exists.
	 * Returns the measurement spring force in format: {v1x, v1y, v1z, v2x, v2y, v2z, ...}.
	 */
	private double getVertexMeasurementSpringEnergy(Vec3D[] vertexPositions, Vec3D[] measurements) {
		double measurementSpringEnergy = 0d;
		for(int vertexInd = 0; vertexInd < vertexPositions.Length; vertexInd++) {
			if(measurements[vertexInd] != null) {
				for(int coord = 0; coord < 3; coord++) {
					double diff = vertexPositions[vertexInd][coord] - measurements[vertexInd][coord];
					measurementSpringEnergy += this.kMeasurementSprings / 2d * diff * diff;
				}
			}
		}
		return measurementSpringEnergy;
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
		string filePath = EditorUtility.SaveFilePanel("Save sail shape", storageBaseDirPath + "/SailData", "", "sailshapedata");
		if(filePath.Length != 0) {
			Mesh mesh = this.getMesh();
			new SailConfiguration(this.vertexPositions, mesh.triangles, this.verticesMovementConstraints).storeToFile(filePath);
		}
	}

	public void onLoadSailShapeButtonPress() {
		string filePath = EditorUtility.OpenFilePanel("Load sail shape", storageBaseDirPath + "/SailData", "sailshapedata");
		if(filePath.Length != 0) {
			this.loadSailConfiguration(SailConfiguration.loadFromFile(filePath));
		}
	}

	private void loadSailConfiguration(SailConfiguration sailConfiguration) {
		this.loadSailConfiguration(sailConfiguration, 0);
	}

	private void loadSailConfiguration(SailConfiguration sailConfiguration, int numSubDivisions) {
		Mesh mesh = new Mesh();
		mesh.vertices = vecToVec(sailConfiguration.vertexPositions);
		mesh.triangles = sailConfiguration.triangles;

		// Perform additional mesh subdivisions.
		if(numSubDivisions > 0) {
			mesh.RecalculateNormals();
			for(int i = 0; i < numSubDivisions; i++) {
				MeshHelper.Subdivide(mesh);
			}

			// Mark the newly added vertices as not constrained.
			this.verticesMovementConstraints = new bool[mesh.vertexCount];
			for(int i = 0; i < sailConfiguration.vertexConstraints.Length; i++) {
				this.verticesMovementConstraints[i] = sailConfiguration.vertexConstraints[i];
			}
			for(int i = sailConfiguration.vertexConstraints.Length; i < this.verticesMovementConstraints.Length; i++) {
				this.verticesMovementConstraints[i] = false;
			}
		} else {
			this.verticesMovementConstraints = sailConfiguration.vertexConstraints;
		}

		// Load the mesh.
		mesh.RecalculateNormals();
		this.loadMesh(mesh, 1d);
	}

	public void onSaveSailMeasurementsButtonPress() {
		string filePath = EditorUtility.SaveFilePanel("Save sail measurements", storageBaseDirPath + "/SailData", "", "measurements");
		if(filePath.Length != 0) {
			if(this.measurementsGenerateFactor < 0f) {
				this.measurementsGenerateFactor = 0f;
			} else if(this.measurementsGenerateFactor > 1f) {
				this.measurementsGenerateFactor = 1f;
			}
			int numMeasurements = (int) (this.vertexPositions.Length * this.measurementsGenerateFactor);
			new SailMeasurements(this.vertexPositions, MeshUtils.generateFarthestPointSamplingSailMeasurements(
					this.vertexPositions, numMeasurements, this.verticesMovementConstraints)).storeToFile(filePath);
		}
	}

	public void onLoadSailMeasurementsButtonPress() {
		string filePath = EditorUtility.OpenFilePanel("Load sail measurements", storageBaseDirPath + "/SailData", "measurements");
		if(filePath.Length != 0) {
			this.measurements = SailMeasurements.loadFromFile(filePath);
			int count = 0;
			for(int i = 0; i < this.measurements.measurements.Length; i++) {
				if(this.measurements.measurements[i] != null) {
					count++;
				}
			}
			print("Loaded " + count + " measurement positions.");
		}
	}

	public void onSaveRecordingButtonPress() {
		if(this.meshRecorder == null) {
			print("No recording available to save.");
			return;
		}
		string filePath = EditorUtility.SaveFilePanel("Save recording", storageBaseDirPath + "/Recordings", "", "rec");
		if(filePath.Length != 0) {
			this.meshRecorder.storeToFile(filePath);
			print("Recording saved.");
		}
	}

	public void onLoadRecordingButtonPress() {
		string filePath = EditorUtility.OpenFilePanel("Load recording", storageBaseDirPath + "/Recordings", "rec");
		if(filePath.Length != 0) {
			if(this.meshRecorder != null) {
				this.meshRecorder.stop();
			}
			this.meshRecorder = MeshRecorder.loadFromFile(filePath);
			print("Recording loaded.");
		}
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
