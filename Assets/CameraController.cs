using UnityEngine;
using System.Collections;

public class CameraController : MonoBehaviour {

    // Camera movement settings.
    float movementSpeed = 10f;

    float mainSpeed = 10.0f; //regular speed
    float shiftAdd = 25.0f; //multiplied by how long shift is held.  Basically running
    float maxShift = 100.0f; //Maximum speed when holdin gshift
    float camSens = 0.15f; //How sensitive it with mouse
    private Vector3 lastMouse = new Vector3(255, 255, 255); //kind of in the middle of the screen, rather than at the top (play)
    private float totalRun = 1.0f;

    // Camera rotation settings.
    public float mouseSensitivity = 10f;
    private float minRotationX = 10f; // 0 = down.
    private float maxRotationX = 170f; // 180 = up.

    void Update() {

        // Get mouse input.
        float horizontal = Input.GetAxis("Mouse X");
        float vertical = Input.GetAxis("Mouse Y");
        if(vertical > 10f) { vertical = 10f; } // Limit y so it cant glitch by looking up or down too fast.
        if(vertical < -10f) { vertical = -10f; } // Limit y so it cant glitch by looking up or down too fast.

        // Set camera rotation.
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

        // Update camera position.
        Vector3 normalizedInputVelocity = this.getNormalizedInputVelocity();
        this.transform.Translate(normalizedInputVelocity * Time.deltaTime * this.movementSpeed);

        // Handle cursor capture.
        if(Input.GetKeyDown(KeyCode.Escape)) {
            Cursor.lockState = CursorLockMode.None;
        }
        if(Input.GetMouseButtonDown(0)) {
            Cursor.lockState = CursorLockMode.Locked;
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
