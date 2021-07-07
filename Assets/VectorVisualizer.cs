using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class VectorVisualizer {

	private List<GameObject> vectorGameObjects = new List<GameObject>();
	private Vector3 offset;

	public VectorVisualizer(Vector3 offset) {
		this.offset = offset;
	}

	public void visualize(Vec3D[] positions, VecD vectors, float scale = 1f, float vecThickness = 0.02f) {
		this.visualize(new VecD(positions), vectors, scale, vecThickness);
	}

	public void visualize(VecD positions, VecD vectors, float scale = 1f, float vecThickness = 0.02f) {
		for(int i = 0; i < positions.length / 3; i++) {
			Vector3 origin = new Vector3((float) positions[3 * i], (float) positions[3 * i + 1], (float) positions[3 * i + 2]);
			Vector3 vec = new Vector3((float) vectors[3 * i], (float) vectors[3 * i + 1], (float) vectors[3 * i + 2]);
			this.visualize(origin, vec, scale, vecThickness);
		}
	}

	public void visualize(Vector3 position, Vector3 vector, float scale = 1f, float vecDiameter = 0.01f) {
		Vector3 origin = position + this.offset;
		Vector3 vec = vector * scale;
		float vecMag = vec.magnitude;
		if(vecMag != 0f) {
			GameObject gameObj = GameObject.CreatePrimitive(PrimitiveType.Cylinder); // Defaults to 2m height, 1m diameter, origin at center.
			Vector3 objMiddle = origin + vec / 2f;
			gameObj.transform.position = objMiddle;
			Vector3 rotAxis = Vector3.Cross(new Vector3(0f, 1f, 0f), vec);
			float angle = Vector3.SignedAngle(new Vector3(0f, 1f, 0f), vec, rotAxis);
			gameObj.transform.rotation = Quaternion.AngleAxis(angle, rotAxis);
			gameObj.transform.localScale = new Vector3(vecDiameter, 0.5f * vecMag, vecDiameter);
			this.vectorGameObjects.Add(gameObj);
		}
	}

	public void clear() {
		foreach(GameObject vectorGameObj in this.vectorGameObjects) {
			MonoBehaviour.Destroy(vectorGameObj);
		}
		this.vectorGameObjects.Clear();
	}
}
