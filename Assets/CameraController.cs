using UnityEngine;
using System.Collections;

public class CameraController : MonoBehaviour {

	// Camera movement settings.
	public float movementSpeed = 10f;

	// Camera rotation settings.
	public float mouseSensitivity = 10f;
	public bool doRotateWhenMouseNotCaptured = false;
	private float minRotationX = 10f; // 0 = down.
	private float maxRotationX = 170f; // 180 = up.

	private bool isLeftMouseBtnDown = false;

	void Update() {

		// Get mouse input.
		float horizontal = Input.GetAxis("Mouse X");
		float vertical = Input.GetAxis("Mouse Y");
		if(vertical > 10f) { vertical = 10f; } // Limit y so it cant glitch by looking up or down too fast.
		if(vertical < -10f) { vertical = -10f; } // Limit y so it cant glitch by looking up or down too fast.

		// Set camera rotation.
		if(this.doRotateWhenMouseNotCaptured || Cursor.lockState == CursorLockMode.Locked) {
			Vector3 cameraRotation = this.transform.rotation.eulerAngles;
			cameraRotation += new Vector3(-vertical, horizontal, 0f) * this.mouseSensitivity;
			cameraRotation.x = (cameraRotation.x + 360f) % 360f; // Maps [-360, inf] to [0, 360].
			if(cameraRotation.x < 180f && cameraRotation.x > 90f - this.minRotationX) {
				cameraRotation.x = 90f - this.minRotationX;
			}
			if(cameraRotation.x > 180f && cameraRotation.x < 270f + (180f - this.maxRotationX)) {
				cameraRotation.x = 270f + (180f - this.maxRotationX);
			}
			this.transform.rotation = Quaternion.Euler(cameraRotation);
		}

		// Update camera position.
		Vector3 normalizedInputVelocity = this.getNormalizedInputVelocity();
		this.transform.Translate(normalizedInputVelocity * Time.deltaTime * this.movementSpeed);

		// Handle cursor capture. Capture on left mouse button release to allow UI buttons to function.
		if(Input.GetKeyDown(KeyCode.Escape)) {
			Cursor.lockState = CursorLockMode.None;
		}
		if(this.isLeftMouseBtnDown != Input.GetMouseButton(0)) {
			this.isLeftMouseBtnDown = !this.isLeftMouseBtnDown;
			if(!this.isLeftMouseBtnDown) { // On left mouse button release.
				Cursor.lockState = CursorLockMode.Locked;
			}
		}
	}

	private Vector3 getNormalizedInputVelocity() {
		Vector3 velocity = new Vector3();
		if(Input.GetKey(KeyCode.W)) {
			velocity += new Vector3(0, 0, 1);
		}
		if(Input.GetKey(KeyCode.S)) {
			velocity += new Vector3(0, 0, -1);
		}
		if(Input.GetKey(KeyCode.A)) {
			velocity += new Vector3(-1, 0, 0);
		}
		if(Input.GetKey(KeyCode.D)) {
			velocity += new Vector3(1, 0, 0);
		}
		if(Input.GetKey(KeyCode.LeftShift)) {
			velocity += new Vector3(0, -1, 0);
		}
		if(Input.GetKey(KeyCode.Space)) {
			velocity += new Vector3(0, 1, 0);
		}
		if(velocity.magnitude != 0f) {
			velocity = velocity.normalized;
		}
		return velocity;
	}
}
