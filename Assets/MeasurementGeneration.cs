using System.Collections;
using UnityEngine;

public class MeasurementGeneration {

	public static void generateLineDensityMeasurements() {
		double sailWidth = 3.5d;
		double sailHeight = 7d;
		foreach(string windMag in new string[] {"64.6963", "1035.1"}) {
			foreach(string windDeg in new string[] {"30", "45", "60"}) {

				// Load the sail configuration.
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx;
				foreach(int n in new int[] {1, 5, 10, 15}) {
					foreach(double m in new double[] {1, 2.5, 5}) {
						SailMeasurements measurements = new SailMeasurements(sailConf.vertexPositions,
								MeshUtils.generateTriangularSailMeasurements(
										sailConf.vertexPositions, sailWidth, sailHeight, n, m));
						string measurementsFilePath = measurementsStorageDir + "/" + "n=" + n + ",m=" + m + ".measurements";
						measurements.storeToFile(measurementsFilePath);
						MonoBehaviour.print("Generated measurements: " + measurementsFilePath);
					}
				}
			}
		}
	}

	public static void generateFarthestFirstMeasurements() {
		foreach(string windMag in new string[] {"64.6963", "1035.1"}) {
			foreach(string windDeg in new string[] {"30", "45", "60"}) {

				// Load the sail configuration.
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=" + windMag + ",wd=" + windDeg;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx;
				foreach(int n in new int[] {5, 10, 20, 30, 40, 50, 70, 90, 110, 130, 140, 150, 160}) {
					SailMeasurements measurements = new SailMeasurements(sailConf.vertexPositions,
							MeshUtils.generateFarthestPointSamplingSailMeasurements(
									sailConf.vertexPositions, n, sailConf.vertexConstraints));
					string measurementsFilePath = measurementsStorageDir + "/" + "n=" + n + ".measurements";
					measurements.storeToFile(measurementsFilePath);
					MonoBehaviour.print("Generated measurements: " + measurementsFilePath);
				}
			}
		}
	}

	public static void generateTest3LineDensityMeasurements() {
		double sailWidth = 3.5d;
		double sailHeight = 7d;
		foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
			int n = (int) nm[0];
			double m = nm[1];
			foreach(double windPressFac in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {

				// Load the sail configuration.
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=45,windPressFac=" + windPressFac;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/test3/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/test3/" + fileNameNoEx;
				SailMeasurements measurements = new SailMeasurements(sailConf.vertexPositions,
						MeshUtils.generateTriangularSailMeasurements(
								sailConf.vertexPositions, sailWidth, sailHeight, n, m));
				string measurementsFilePath = measurementsStorageDir + "/" + "n=" + n + ",m=" + m + ".measurements";
				measurements.storeToFile(measurementsFilePath);
				MonoBehaviour.print("Generated measurements: " + measurementsFilePath);
			}
		}
	}

	public static void generateTest3FarthestFirstMeasurements() {
		foreach(int n in new int[] {30, 50, 110}) {
			foreach(double windPressFac in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {

				// Load the sail configuration.
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=45,windPressFac=" + windPressFac;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/test3/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/test3/" + fileNameNoEx;
				SailMeasurements measurements = new SailMeasurements(sailConf.vertexPositions,
						MeshUtils.generateFarthestPointSamplingSailMeasurements(
								sailConf.vertexPositions, n, sailConf.vertexConstraints));
				string measurementsFilePath = measurementsStorageDir + "/" + "n=" + n + ".measurements";
				measurements.storeToFile(measurementsFilePath);
				MonoBehaviour.print("Generated measurements: " + measurementsFilePath);
			}
		}
	}

	public static void generateTest4LineDensityMeasurements() {
		double sailWidth = 3.5d;
		double sailHeight = 7d;
		foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
			int n = (int) nm[0];
			double m = nm[1];
			foreach(double noiseMagSlope in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {

				// Load the sail configuration.
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=60,noiseMagSlope=" + noiseMagSlope;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/test4/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/test4/" + fileNameNoEx;
				SailMeasurements measurements = new SailMeasurements(sailConf.vertexPositions,
						MeshUtils.generateTriangularSailMeasurements(
								sailConf.vertexPositions, sailWidth, sailHeight, n, m));
				string measurementsFilePath = measurementsStorageDir + "/" + "n=" + n + ",m=" + m + ".measurements";
				measurements.storeToFile(measurementsFilePath);
				MonoBehaviour.print("Generated measurements: " + measurementsFilePath);
			}
		}
	}

	public static void generateTest4FarthestFirstMeasurements() {
		foreach(int n in new int[] {30, 50, 110}) {
			foreach(double noiseMagSlope in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {

				// Load the sail configuration.
				string fileNameNoEx = "nv=561,kL=15621,kA=8125,kB=10,t=0.00025,d=920"
						+ " ss with g=9.81,wm=1035.1,wd=60,noiseMagSlope=" + noiseMagSlope;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/test4/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/test4/" + fileNameNoEx;
				SailMeasurements measurements = new SailMeasurements(sailConf.vertexPositions,
						MeshUtils.generateFarthestPointSamplingSailMeasurements(
								sailConf.vertexPositions, n, sailConf.vertexConstraints));
				string measurementsFilePath = measurementsStorageDir + "/" + "n=" + n + ".measurements";
				measurements.storeToFile(measurementsFilePath);
				MonoBehaviour.print("Generated measurements: " + measurementsFilePath);
			}
		}
	}
}
