using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ShellHelper {















    
    //public static float getEdgeLengthEnergy(Vector3[] vertices, Vector3[] originalVertices, int v1, int v2) {
    //    Vector3 edge = vertices[v2] - vertices[v1]; // Vector from v1 to v2.
    //    Vector3 undeformedEdge = originalVertices[v2] - originalVertices[v1];
    //    float edgeLength = edge.magnitude;
    //    float undeformedEdgeLength = undeformedEdge.magnitude;
    //    return Mathf.Pow(1f - edgeLength / undeformedEdgeLength, 2) * undeformedEdgeLength;
    //}

    ///**
    // * Computes the edge length energy gradient of vertex v1, for the edge between vertices v1 and v2.
    // */
    //private static Vector3 getEdgeLengthEnergyGradient(Vector3[] vertices, Vector3[] originalVertices, int v1, int v2) {
    //    Vector3 edge = vertices[v2] - vertices[v1]; // Vector from v1 to v2.
    //    Vector3 undeformedEdge = originalVertices[v2] - originalVertices[v1];
    //    float edgeLength = edge.magnitude;
    //    if(float.IsNaN(edgeLength)) {
    //        return Vector3.zero; // Edge is zero-length, so the gradient is 0.
    //    }
    //    float undeformedEdgeLength = undeformedEdge.magnitude;
    //    Vector3 dEdgeLength = (vertices[v1] - vertices[v2]) / edgeLength;
    //    Vector3 result = (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength;
    //    if(float.IsNaN(result.x) || float.IsNaN(result.y) || float.IsNaN(result.z)) {
    //        MonoBehaviour.print("NaN length gradient: " + result + " undeformedEdgeLength: " + undeformedEdgeLength
    //                + " edgeLength: " + edgeLength);
    //        return Vector3.zero;
    //    }
    //    if(float.IsInfinity(result.x) || float.IsInfinity(result.y) || float.IsInfinity(result.z)) {
    //        MonoBehaviour.print("Infinite length gradient: " + result + " undeformedEdgeLength: " + undeformedEdgeLength
    //                + " edgeLength: " + edgeLength);
    //        return Vector3.zero;
    //    }
    //    return (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength;
    //}

    ///**
    // * Computes the edge length energy Hessian of vertex v1, for the edge between vertices v1 and v2.
    // */
    //private Vector3[] getEdgeLengthEnergyHess(Vector3[] vertices, Vector3[] originalVertices, int v1, int v2) {
    //    /*
    //     * EdgeLength (float):
    //     *     sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
    //     * 
    //     * dEdgeLength_dv1x (float):
    //     *     (v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2) = (v1x - v2x) / edgeLength
    //     * 
    //     * dEdgeLength_dv1 (Vector3):
    //     *     (v1 - v2) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2) = (v1 - v2) / edgeLength
    //     * 
    //     * ddEdgeLength_dv1x_dv1y (float):
    //     *    ((v1x - v2x) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
    //     *    (edgeLength * (v1x - v2x)' - (v1x - v2x) * edgeLength') / edgeLength^2
    //     *    (edgeLength * 0 - (v1x - v2x) * dEdgeLength.y) / edgeLength^2
    //     *    (v2x - v1x) * dEdgeLength.y / edgeLength^2
    //     *
    //     * ddEdgeLength_dv1y_dv1x (float):
    //     *    ((v1y - v2y) / sqrt((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2))' = ((v1x - v2x) / edgeLength)'
    //     *    (edgeLength * (v1y - v2y)' - (v1y - v2y) * edgeLength') / edgeLength^2
    //     *    (edgeLength * 0 - (v1y - v2y) * dEdgeLength.x) / edgeLength^2
    //     *    (v2y - v1y) * dEdgeLength.x / edgeLength^2 = (v2y - v1y) * (v1x - v2x) / edgeLength / edgeLength^2 = (v2y - v1y) * (v1x - v2x) / edgeLength^3
    //     * 
    //     * ddEdgeLength_dv1x_dv1x (float):
    //     *     ((v1x - v2x) / edgeLength)'
    //     *     (edgeLength * (v1x - v2x)' - (v1x - v2x) * dEdgeLength.x) / edgeLength^2
    //     *     (edgeLength - (v1x - v2x) * dEdgeLength.x) / ((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
    //     *         = (edgeLength - (v1x - v2x) * dEdgeLength.x) / edgeLength^2
    //     *         = (edgeLength - (v1x - v2x) * (v1x - v2x) / edgeLength) / edgeLength^2
    //     *         = (edgeLength^2 - (v1x - v2x)^2) / edgeLength^3
    //     * 
    //     * ddEdgeLength_dv1y_dv1y (float):
    //     *     ((v1y - v2y) / edgeLength)'
    //     *     (edgeLength * (v1y - v2y)' - (v1y - v2y) * dEdgeLength.y) / edgeLength^2
    //     *     (edgeLength - (v1y - v2y) * dEdgeLength.y) / ((v2x - v1x)^2 + (v2y - v1y)^2 + (v2z - v1z)^2)
    //     *     (edgeLength - (v1y - v2y) * dEdgeLength.y) / edgeLength^2
    //     *     (edgeLength - (v1y - v2y) * (v1 - v2) / edgeLength) / edgeLength^2
    //     *     (edgeLength^2 - (v1y - v2y)^2) / edgeLength^3
    //     * 
    //     * ddEdgeLength_dv1_dv1 (symmetric 3x3 matrix, Hessian):
    //     *     | ddEdgeLength_dv1x_dv1x, ddEdgeLength_dv1y_dv1x, ddEdgeLength_dv1z_dv1x |
    //     *     | ddEdgeLength_dv1x_dv1y, ddEdgeLength_dv1y_dv1y, ddEdgeLength_dv1z_dv1y |
    //     *     | ddEdgeLength_dv1x_dv1z, ddEdgeLength_dv1y_dv1z, ddEdgeLength_dv1z_dv1z |
    //     *     =
    //     *     | edgeLength^2 - (v1x - v2x)^2, (v2y - v1y) * (v1x - v2x)   , (v2z - v1z) * (v1x - v2x) |
    //     *     | (v2x - v1x) * (v1y - v2y)   , edgeLength^2 - (v2y - v1y)^2, (v2z - v1z) * (v1x - v2x) | / edgeLength^3
    //     *     | (v2x - v1x) * (v1z - v2z)   , (v2y - v1y) * (v1z - v2z), edgeLength^2 - (v2z - v1z)^2 |
    //     * 
    //     * Length energy gradient (Vector3):
    //     *     (2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_dv1
    //     * 
    //     * Length energy Hessian (3x3 matrix):
    //     *     ((2 * edgeLength / undeformedEdgeLength - 2) * dEdgeLength_dv1)'
    //     *     (2 * edgeLength / undeformedEdgeLength - 2)' * dEdgeLength_dv1 + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1
    //     *     2 / undeformedEdgeLength * dEdgeLength_dv1 + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1
    //     *     =
    //     *     | 2 / undeformedEdgeLength * dEdgeLength_dv1x + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1x |
    //     *     | 2 / undeformedEdgeLength * dEdgeLength_dv1y + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1y |
    //     *     | 2 / undeformedEdgeLength * dEdgeLength_dv1z + (2 * edgeLength / undeformedEdgeLength - 2) * ddEdgeLength_dv1_dv1z |
    //     */

    //    Vector3 e_v1_v2 = vertices[v2] - vertices[v1];
    //    Vector3 undeformedEdge = originalVertices[v2] - originalVertices[v1];
    //    float edgeLength = e_v1_v2.magnitude;
    //    if(float.IsNaN(edgeLength)) {
    //        return new Vector3[] { Vector3.zero, Vector3.zero, Vector3.zero };
    //    }
    //    float undeformedEdgeLength = undeformedEdge.magnitude;
    //    Vector3 dEdgeLength_dv1 = (vertices[v1] - vertices[v2]) / edgeLength;

    //    float edgeLengthSquare = edgeLength * edgeLength;
    //    float edgeLengthCube = edgeLengthSquare * edgeLength;
    //    Vector3[] ddEdgeLength_dv1_dv1 = new Vector3[] {
    //        new Vector3(edgeLengthSquare - e_v1_v2.x * e_v1_v2.x,                  - e_v1_v2.y * e_v1_v2.x,                  - e_v1_v2.z * e_v1_v2.x) / edgeLengthCube,
    //        new Vector3(                 - e_v1_v2.x * e_v1_v2.y, edgeLengthSquare - e_v1_v2.y * e_v1_v2.y,                  - e_v1_v2.z * e_v1_v2.y) / edgeLengthCube,
    //        new Vector3(                 - e_v1_v2.x * e_v1_v2.z,                  - e_v1_v2.y * e_v1_v2.z, edgeLengthSquare - e_v1_v2.z * e_v1_v2.z) / edgeLengthCube
    //    };

    //    // Calculate length energy Hessian. This matrix is symmetrical, so transposing it has no effect.
    //    Vector3 a = 2 / undeformedEdgeLength * dEdgeLength_dv1;
    //    float b = 2 * edgeLength / undeformedEdgeLength - 2;
    //    Vector3[] lengthEnergyHess = new Vector3[] {
    //        a + b * ddEdgeLength_dv1_dv1[0],
    //        a + b * ddEdgeLength_dv1_dv1[1],
    //        a + b * ddEdgeLength_dv1_dv1[2]
    //    };
    //    return lengthEnergyHess;
    //}

    ///**
    // * Computes the triangle area energy gradient of vertex v1, for the triangle defined by vertices v1, v2 and v3.
    // * Note that 1/3 of the area energy is used for v1 (the other two parts will be used for v2 and v3).
    // */
    //private Vector3 getTriangleAreaEnergyGradient(Vector3[] vertices, Vector3[] originalVertices, int triangleId, int v1, int v2, int v3) {

    //    // Get two edges, where only one is dependent on vertex v1.
    //    Vector3 edge21 = vertices[v1] - vertices[v2]; // dEdge21 / dv1 = {1, 1, 1}.
    //    Vector3 edge23 = vertices[v3] - vertices[v2]; // dEdge23 / dv1 = {0, 0, 0}.

    //    // Calculate the triangle area gradient.
    //    if(this.triangleAreas[triangleId] == 0f) {
    //        return Vector3.zero; // Area is 0 m^2, so the gradient is 0.
    //    }
    //    float crossProdLength = this.triangleAreas[triangleId] * 2f;
    //    Vector3 dCrossProdLength = 1f / crossProdLength * new Vector3(
    //            (edge21.x * edge23.z - edge23.x * edge21.z) * edge23.z + (edge21.x * edge23.y - edge23.x * edge21.y) * edge23.y,
    //            (edge21.y * edge23.z - edge23.y * edge21.z) * edge23.z + (edge21.x * edge23.y - edge23.x * edge21.y) * -edge23.x,
    //            (edge21.y * edge23.z - edge23.y * edge21.z) * -edge23.y + (edge21.x * edge23.z - edge23.x * edge21.z) * -edge23.x);
    //    Vector3 dTriangleArea = dCrossProdLength / 6f; // Area of triangle is half the cross product length, and we only look at a third.

    //    // Calculate the area energy gradient.
    //    Vector3 result = (2 * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2) * dTriangleArea;
    //    if(float.IsNaN(result.x) || float.IsNaN(result.y) || float.IsNaN(result.z)) {
    //        print("Nan area gradient: " + result + " triangleArea: " + this.triangleAreas[triangleId]
    //                + " dTriangleArea: " + dTriangleArea + " crossProdLength: " + crossProdLength);
    //        return Vector3.zero;
    //    }
    //    if(float.IsInfinity(result.x) || float.IsInfinity(result.y) || float.IsInfinity(result.z)) {
    //        print("Infinite area gradient: " + result + " triangleArea: " + this.triangleAreas[triangleId]
    //                + " dTriangleArea: " + dTriangleArea + " crossProdLength: " + crossProdLength);
    //        return Vector3.zero;
    //    }
    //    return (2 * this.triangleAreas[triangleId] / this.undeformedTriangleAreas[triangleId] - 2) * dTriangleArea;
    //}
}
