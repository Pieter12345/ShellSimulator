using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class MeshRecorder {
	
	// Record variables.
	private int[] triangles;
	private List<Tuple<double, Vec3D[]>> vertexPositionsTimeTuples;

	// Replay variables.
	private GameObject gameObject;
	private double replayTime;
	private int replayIndex = -1;
	private double lastUpdateTime;

	public MeshRecorder(int[] triangles, Vec3D[] initialVertexPositions) {
		this.triangles = triangles;
		this.vertexPositionsTimeTuples = new List<Tuple<double, Vec3D[]>> {
			new Tuple<double, Vec3D[]>(0, Vec3D.clone(initialVertexPositions))
		};
	}

	private MeshRecorder(int[] triangles, List<Tuple<double, Vec3D[]>> vertexPositionsTimeTuples) {
		this.triangles = triangles;
		this.vertexPositionsTimeTuples = vertexPositionsTimeTuples;
	}

	public void record(double deltaTime, Vec3D[] vertexPositions) {
		if(this.isPlaying()) {
			throw new Exception("Cannot record during playback.");
		}
		this.vertexPositionsTimeTuples.Add(new Tuple<double, Vec3D[]>(deltaTime, Vec3D.clone(vertexPositions)));
	}

	public void replay() {

		// Create the GameObject used to render the mesh.
		this.gameObject = new GameObject();
		this.gameObject.AddComponent<MeshRenderer>();
		this.gameObject.AddComponent<MeshFilter>();

		// Set the mesh renderer.
		MeshRenderer meshRenderer = this.gameObject.GetComponent<MeshRenderer>();
		meshRenderer.sharedMaterial = new Material(Shader.Find("Custom/StandardTwoSides"));

		// Set GameObject position.
		this.gameObject.transform.Translate(10, 1, 0);

		// Create mesh.
		Mesh mesh = new Mesh();
		mesh.vertices = vecToVec(this.vertexPositionsTimeTuples[0].Item2);
		mesh.triangles = this.triangles;
		mesh.RecalculateNormals();
		this.gameObject.GetComponent<MeshFilter>().mesh = mesh;

		// Reset replay variables.
		this.replayTime = 0;
		this.replayIndex = 0;
		this.lastUpdateTime = 0;
	}

	public bool isPlaying() {
		return this.replayIndex >= 0;
	}

	public void stop() {
		this.replayIndex = -1;
		UnityEngine.Object.Destroy(this.gameObject);
	}

	public void update() {
		if(this.replayIndex < 0) {
			return;
		}
		if(this.replayIndex >= this.vertexPositionsTimeTuples.Count) {
			this.stop();
			return;
		}
		Tuple<double, Vec3D[]> tuple = this.vertexPositionsTimeTuples[this.replayIndex];
		if(this.replayTime >= this.lastUpdateTime + tuple.Item1) {
			Mesh mesh = this.gameObject.GetComponent<MeshFilter>().mesh;
			mesh.vertices = vecToVec(tuple.Item2);
			mesh.RecalculateNormals();
			this.lastUpdateTime += tuple.Item1;
			this.replayIndex++;
		}
		this.replayTime += Time.deltaTime;
	}

	public void storeToFile(string filePath) {
		MonoBehaviour.print("Storing sail measurements to file: " + filePath);
		string dirPath = Path.GetDirectoryName(filePath);
		if(!Directory.Exists(dirPath)) {
			Directory.CreateDirectory(dirPath);
		}
		FileStream fs = (File.Exists(filePath) ? File.OpenWrite(filePath) : File.Create(filePath));
		BinaryWriter writer = new BinaryWriter(fs);

		// Write triangles.
		writer.Write(this.triangles.Length);
		foreach(int val in this.triangles) {
			writer.Write(val);
		}

		// Write time position tuples.
		writer.Write(this.vertexPositionsTimeTuples.Count);
		foreach(Tuple<double, Vec3D[]> tuple in this.vertexPositionsTimeTuples) {
			double deltaTime = tuple.Item1;
			Vec3D[] vertexPositions = tuple.Item2;
			writer.Write(deltaTime);
			writer.Write(vertexPositions.Length);
			foreach(Vec3D vertexPosition in vertexPositions) {
				writer.Write(vertexPosition.x);
				writer.Write(vertexPosition.y);
				writer.Write(vertexPosition.z);
			}
		}
		writer.Close();
	}

	public static MeshRecorder loadFromFile(string filePath) {
		MonoBehaviour.print("Loading sail measurements from file: " + filePath);
		FileStream fs = File.OpenRead(filePath);
		BinaryReader reader = new BinaryReader(fs);

		// Read triangles.
		int triangleLength = reader.ReadInt32();
		int[] triangles = new int[triangleLength];
		for(int i = 0; i < triangleLength; i++) {
			triangles[i] = reader.ReadInt32();
		}

		// Read time position tuples.
		int numVertexPositionTuples = reader.ReadInt32();
		List<Tuple<double, Vec3D[]>> vertexPositionsTimeTuples = new List<Tuple<double, Vec3D[]>>();
		for(int i = 0; i < numVertexPositionTuples; i++) {
			double deltaTime = reader.ReadDouble();
			int numVertices = reader.ReadInt32();
			Vec3D[] vertexPositions = new Vec3D[numVertices];
			for(int j = 0; j < numVertices; j++) {
				vertexPositions[j] = new Vec3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
			}
			vertexPositionsTimeTuples.Add(new Tuple<double, Vec3D[]>(deltaTime, vertexPositions));
		}
		reader.Close();
		return new MeshRecorder(triangles, vertexPositionsTimeTuples);
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
