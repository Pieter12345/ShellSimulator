using System.Collections;
using UnityEngine;

public class MeasurementGeneration {

	public static void generateMeasurements() {
		double sailWidth = 3.5d;
		double sailHeight = 7d;
		foreach(string windMag in new string[] {"9.6", "38.4"}) {
			foreach(string windDeg in new string[] {"30", "45", "60"}) {

				// Load the sail configuration.
				string fileNameNoEx = "numVerts=561, kLength=500, kArea=500, kBend=0.01, thickness=0.002, density=40,"
						+ " steady state with g=9.81, windMag=" + windMag + ", windDeg=" + windDeg;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx;
				foreach(int n in new int[] {1, 5, 10, 15}) {
					foreach(double m in new double[] {1, 2.5, 5}) {
						SailMeasurements measurements = new SailMeasurements(sailConf.vertexPositions,
								MeshUtils.generateTriangularSailMeasurements(
										sailConf.vertexPositions, sailWidth, sailHeight, n, m));
						string measurementsFilePath = measurementsStorageDir + "/" + "n=" + n + ", m=" + m + ".measurements";
						measurements.storeToFile(measurementsFilePath);
						MonoBehaviour.print("Generated measurements: " + measurementsFilePath);
					}
				}
			}
		}
	}

	public static void generateTest3Measurements() {
		double sailWidth = 3.5d;
		double sailHeight = 7d;
		foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
			int n = (int) nm[0];
			double m = nm[1];
			foreach(double windPressFac in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {

				// Load the sail configuration.
				string fileNameNoEx = "numVerts=561,kLength=500,kArea=500,kBend=0.01,thickness=0.002,density=40,"
						+ "steady state with g=9.81,windMag=38.4,windDeg=45,windPressFac=" + windPressFac;
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

	public static void generateTest4Measurements() {
		double sailWidth = 3.5d;
		double sailHeight = 7d;
		foreach(double[] nm in new double[][] {new double[] {5, 2.5}, new double[] {15, 5}}) {
			int n = (int) nm[0];
			double m = nm[1];
			foreach(double noiseMagSlope in new double[] {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}) {

				// Load the sail configuration.
				string fileNameNoEx = "numVerts=561,kLength=500,kArea=500,kBend=0.01,thickness=0.002,density=40,"
						+ "steady state with g=9.81,windMag=38.4,windDeg=60,noiseMagSlope=" + noiseMagSlope;
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
}
