using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class SailMeasurements {
	
	/*
	 * Vertex positions at the time the measurements were taken. This can be used as ground truth.
	 * This is null when the ground truth configuration is unknown.
	 */
	public Vec3D[] vertexPositions { get; }

	/*
	 * The measurements. One element for each vertex, null meaning that there is no measurement for that vertex.
	 */
	public Vec3D[] measurements { get; }

	public SailMeasurements(Vec3D[] vertexPositions, Vec3D[] measurements) {
		this.vertexPositions = (vertexPositions == null ? null : (Vec3D[]) vertexPositions.Clone());
		this.measurements = (Vec3D[]) measurements.Clone();
	}

	public void storeToFile(string filePath) {
		string dirPath = Path.GetDirectoryName(filePath);
		MonoBehaviour.print("Storing sail measurements to file: " + filePath);
		if(!Directory.Exists(dirPath)) {
			Directory.CreateDirectory(dirPath);
		}
		FileStream fs = (File.Exists(filePath) ? File.OpenWrite(filePath) : File.Create(filePath));
		BinaryWriter writer = new BinaryWriter(fs);

		writer.Write((int) this.vertexPositions.Length);

		foreach(Vec3D pos in this.vertexPositions) {
			writer.Write((double) pos.x);
			writer.Write((double) pos.y);
			writer.Write((double) pos.z);
		}

		foreach(Vec3D pos in this.measurements) {
			if(pos != null) {
				writer.Write(true);
				writer.Write((double) pos.x);
				writer.Write((double) pos.y);
				writer.Write((double) pos.z);
			} else {
				writer.Write(false);
			}
		}
		writer.Close();
	}

	public static SailMeasurements loadFromFile(string filePath) {
		MonoBehaviour.print("Loading sail measurements from file: " + filePath);
		FileStream fs = File.OpenRead(filePath);
		BinaryReader reader = new BinaryReader(fs);

		int numVertices = reader.ReadInt32();

		Vec3D[] vertexPositions = new Vec3D[numVertices];
		for(int i = 0; i < numVertices; i++) {
			vertexPositions[i] = new Vec3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
		}

		Vec3D[] measurements = new Vec3D[numVertices];
		for(int i = 0; i < numVertices; i++) {
			if(reader.ReadBoolean()) {
				measurements[i] = new Vec3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
			} else {
				measurements[i] = null;
			}
		}
		reader.Close();
		return new SailMeasurements(vertexPositions, measurements);
	}

	public int getNumMeasurements() {
		int numMeasurements = 0;
		foreach(VecD measurement in this.measurements) {
			if(measurement != null) {
				numMeasurements++;
			}
		}
		return numMeasurements;
	}
}
