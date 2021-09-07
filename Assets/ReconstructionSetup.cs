using System;

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
	public Vec3D initialWindPressure;
	public double gravityConstant;
	public bool doStaticMinimization;
	public double maxWindSpeed;
	public double maxDeltaWindSpeed;
	public int minNumNewtonIterations;
	public int numWindReconstructionSteps;
	public int numSnapReconstructionSteps;
}
