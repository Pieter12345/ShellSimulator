using System.Collections;
using UnityEngine;

public class MeasurementGeneration {

	public static void generateMeasurements() {
		double sailWidth = 3.5d;
		double sailHeight = 7d;
		foreach(string windMag in new string[] {"9.6", "38.4"}) {
			foreach(string windDeg in new string[] {"30", "45", "60"}) {

				// Load the sail configuration.
				string fileNameNoEx = "numVerts=561, kLength=500, kArea=500, kBend=0.01, thickness=0.002, density=40," +
						" steady state with g=9.81, windMag=" + windMag + ", windDeg=" + windDeg;
				SailConfiguration sailConf = SailConfiguration.loadFromFile(
						Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx + ".sailshapedata");

				// Generate and store the measurements.
				string measurementsStorageDir = Shell.storageBaseDirPath + "/SailData/" + fileNameNoEx;
				foreach(int n in new int[] {1, 5, 10}) {
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
}