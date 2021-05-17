using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MeshRecorder {
	
	// Record variables.
	private int[] triangles;
	private List<Edge> edges;
	private List<Tuple<double, Vec3D[]>> vertexPositionsTimeTuples = new List<Tuple<double, Vec3D[]>>();

	// Replay variables.
	private GameObject gameObject;
	private double replayTime;
	private int replayIndex = -1;
	private double lastUpdateTime;

	public MeshRecorder(int[] triangles, List<Edge> edges, Vec3D[] initialVertexPositions) {
		this.triangles = triangles;
		this.edges = edges;
		this.vertexPositionsTimeTuples.Add(new Tuple<double, Vec3D[]>(0, Vec3D.clone(initialVertexPositions)));
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
