﻿using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Visualizer {

	private List<GameObject> visualizationGameObjects = new List<GameObject>();
	private Vector3 origin;

	public Visualizer(Vector3 origin) {
		this.origin = origin;
	}

	public void visualizeVectors(Vec3D[] relPositions, VecD vectors, Color color, float scale = 1f, float vecThickness = 0.02f) {
		this.visualizeVectors(new VecD(relPositions), vectors, color, scale, vecThickness);
	}

	public void visualizeVectors(VecD relPositions, VecD vectors, Color color, float scale = 1f, float vecThickness = 0.02f) {
		for(int i = 0; i < relPositions.length / 3; i++) {
			Vector3 relStartPos = new Vector3((float) relPositions[3 * i], (float) relPositions[3 * i + 1], (float) relPositions[3 * i + 2]);
			Vector3 vec = new Vector3((float) vectors[3 * i], (float) vectors[3 * i + 1], (float) vectors[3 * i + 2]);
			this.visualizeVectors(relStartPos, vec, color, scale, vecThickness);
		}
	}

	public void visualizeVectors(Vector3 relPosition, Vector3 vector, Color color, float scale = 1f, float vecDiameter = 0.01f) {
		Vector3 position = relPosition + this.origin;
		Vector3 vec = vector * scale;
		float vecMag = vec.magnitude;
		if(vecMag != 0f) {
			GameObject gameObj = GameObject.CreatePrimitive(PrimitiveType.Cylinder); // Defaults to 2m height, 1m diameter, origin at center.
			Renderer renderer = gameObj.GetComponent<Renderer>();
			renderer.material.color = color;
			Vector3 objMiddle = position + vec / 2f;
			gameObj.transform.position = objMiddle;
			Vector3 rotAxis = Vector3.Cross(new Vector3(0f, 1f, 0f), vec);
			float angle = Vector3.SignedAngle(new Vector3(0f, 1f, 0f), vec, rotAxis);
			gameObj.transform.rotation = Quaternion.AngleAxis(angle, rotAxis);
			gameObj.transform.localScale = new Vector3(vecDiameter, 0.5f * vecMag, vecDiameter);
			this.visualizationGameObjects.Add(gameObj);
		}
	}

	public void visualizePoints(Vec3D[] relPositions, Color color, float radius = 0.1f) {
		this.visualizePoints(new VecD(relPositions), color, radius);
	}

	public void visualizePoints(VecD relPositions, Color color, float radius = 0.1f) {
		for(int i = 0; i < relPositions.length / 3; i++) {
			Vector3 relPosition = new Vector3((float) relPositions[3 * i], (float) relPositions[3 * i + 1], (float) relPositions[3 * i + 2]);
			this.visualizePoints(relPosition, color, radius);
		}
	}

	public void visualizePoints(Vector3 relPosition, Color color, float radius = 0.1f) {
		Vector3 position = relPosition + this.origin;
		GameObject gameObj = GameObject.CreatePrimitive(PrimitiveType.Sphere); // Defaults to 1m diameter, origin at center.
		Renderer renderer = gameObj.GetComponent<Renderer>();
		renderer.material.color = color;
		gameObj.transform.position = position;
		gameObj.transform.localScale = new Vector3(radius * 2f, radius * 2f, radius * 2f); // Factor 2 to scale the radius, rather than the diameter.
		this.visualizationGameObjects.Add(gameObj);
	}

	public void clear() {
		foreach(GameObject visualizationGameObj in this.visualizationGameObjects) {
			MonoBehaviour.Destroy(visualizationGameObj);
		}
		this.visualizationGameObjects.Clear();
	}
}
