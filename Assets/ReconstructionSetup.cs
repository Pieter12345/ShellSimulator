﻿using System;

public class ReconstructionSetup {

	public string sailStartConfigurationRelPath;
	public string sailMeasurementsRelPath;
	public string resultsStorageRelPath;
	public float kLength;
	public float kArea;
	public float kBend;
	public float shellThickness;
	public float shellMaterialDensity;
	public bool useFlatUndeformedBendState;
	public int numRestShapeSubdivisions = 0; // Additional subdivisions applied to the sail rest shape.
	public bool isRealWorldDataTest = false; // Convenient property to perform additional real-world-data specific actions.
	public int numRealWorldDataFarthestFirstMarkers = 0; // Amount of farthest-first selected markers to use from the measurements data.
	public Vec3D initialWindPressureVec;
	public double initialWindPressure;
	public int initialPositionNoiseRandomSeed = 0;
	public double initialPositionNoiseMagnitude = 0d; // Applied in random directions to all non-constrained vertex positions on mesh load.
	public double gravityConstant;
	public bool doStaticMinimization;
	public double maxWindSpeed;
	public double maxDeltaWindSpeed;
	public int minNumNewtonIterations;
	public int numWindReconstructionSteps1;
	public int numWindReconstructionSteps2;
	public Tuple<Vec3D, double>[] measurementsIgnoreSpheres = null; // Used to ignore measurements within a radius from a position.
}
