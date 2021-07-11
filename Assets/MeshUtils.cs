using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MeshUtils {

	public static List<int>[] getSortedVertexTriangles(int[] triangles, int numVertices) {
		
		// Create triangle list per vertex.
		List<int>[] sortedVertexTriangles = new List<int>[numVertices];
		for(int i = 0; i < sortedVertexTriangles.Length; i++) {
			sortedVertexTriangles[i] = new List<int>();
		}
		for(int i = 0; i < triangles.Length; i += 3) {
			int triangleId = i / 3;
			sortedVertexTriangles[triangles[i]].Add(triangleId);
			sortedVertexTriangles[triangles[i + 1]].Add(triangleId);
			sortedVertexTriangles[triangles[i + 2]].Add(triangleId);
		}

		// Sort triangle list per vertex in clockwise order.
		for(int vertexInd = 0; vertexInd < sortedVertexTriangles.Length; vertexInd++) {
			List<int> triangleList = sortedVertexTriangles[vertexInd];
			if(triangleList.Count > 0) {
				List<int> sortedTriangleList = new List<int>(triangleList.Count);
				sortedTriangleList.Add(triangleList[0]);
				for(int i = 1; i < triangleList.Count; i++) {

					// Get triangle vertices.
					int triangleBaseInd1 = triangleList[i] * 3;
					int v11 = triangles[triangleBaseInd1];
					int v12 = triangles[triangleBaseInd1 + 1];
					int v13 = triangles[triangleBaseInd1 + 2];

					// Get the vertex indices of the other vertices that are connected to this vertex's edges.
					int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
					int otherVertexAntiClockwiseInd1 = (vertexInd == v11 ? v12 : (vertexInd == v12 ? v13 : v11));

					// Insert the triangle into the sorted triangles list.
					bool inserted = false;
					for(int j = 0; j < sortedTriangleList.Count; j++) {

						// Get triangle vertices.
						int triangleBaseInd2 = sortedTriangleList[j] * 3;
						int v21 = triangles[triangleBaseInd2];
						int v22 = triangles[triangleBaseInd2 + 1];
						int v23 = triangles[triangleBaseInd2 + 2];

						// Get the vertex indices of the other vertices that are connected to this vertex's edges.
						int otherVertexClockwiseInd2 = (vertexInd == v21 ? v23 : (vertexInd == v23 ? v22 : v21));
						int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));

						// Insert the triangle before the current triangle if it is placed next to it in anti-clockwise rotation.
						if(otherVertexClockwiseInd1 == otherVertexAntiClockwiseInd2) {
							sortedTriangleList.Insert(j, triangleList[i]);
							inserted = true;
							break;
						}

						// Insert the triangle after the current triangle if it is placed next to it in clockwise rotation.
						if(otherVertexAntiClockwiseInd1 == otherVertexClockwiseInd2) {
							sortedTriangleList.Insert(j + 1, triangleList[i]);
							inserted = true;
							break;
						}
					}
					if(!inserted) {
						sortedTriangleList.Add(triangleList[i]);
					}
				}

				// Replace vertex triangles list with sorted vertex triangles list.
				sortedVertexTriangles[vertexInd] = sortedTriangleList;
			}
		}
		return sortedVertexTriangles;
	}

	public static List<Edge> getEdges(List<int>[] sortedVertexTriangles, int[] triangles) {
		List<Edge> edges = new List<Edge>();
		for(int vertexInd = 0; vertexInd < sortedVertexTriangles.Length; vertexInd++) {
			List<int> triangleList = sortedVertexTriangles[vertexInd];
			for(int i = 0; i < triangleList.Count; i++) {

				// Get two possibly adjacent triangles in clockwise direction.
				int triangleId = triangleList[i];
				int nextTriangleId = triangleList[(i + 1) % triangleList.Count];

				// Get triangle vertices.
				int triangleBaseInd1 = triangleId * 3;
				int triangleBaseInd2 = nextTriangleId * 3;
				int v11 = triangles[triangleBaseInd1];
				int v12 = triangles[triangleBaseInd1 + 1];
				int v13 = triangles[triangleBaseInd1 + 2];
				int v21 = triangles[triangleBaseInd2];
				int v22 = triangles[triangleBaseInd2 + 1];
				int v23 = triangles[triangleBaseInd2 + 2];

				// Get the vertex indices of the other vertices that are connected to the possibly shared triangle edge.
				int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
				int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));
				
				// Define edge vertices.
				int ve1 = vertexInd;
				int ve2 = otherVertexClockwiseInd1;
				
				// Add vertices defining the side flaps if the edge is shared between two triangles.
				bool edgeSharedByTriangles = (otherVertexClockwiseInd1 == otherVertexAntiClockwiseInd2);
				if(edgeSharedByTriangles) {

					// Skip shared edges where v1 > v2. This prevents shared edges from being added twice.
					if(ve1 > ve2) {
						continue;
					}

					// Initialize edge and flap vertices.
					int vf1 = (vertexInd == v11 ? v12 : (vertexInd == v12 ? v13 : v11));
					int vf2 = (vertexInd == v21 ? v23 : (vertexInd == v22 ? v21 : v22));

					// Add edge with side flaps.
					edges.Add(new Edge(ve1, ve2, vf1, vf2, triangleId, nextTriangleId));
				} else {

					// Add edge without side flaps.
					edges.Add(new Edge(ve1, ve2));
				}
			}
		}
		return edges;
	}

	/**
	 * Gets a list of all edges within the given nested vertex edges map
	 * (as returned by MeshUtils.getVertexEdgesMap(List<int>[], int[])).
	 * Returns a list containing all unique edges in the nested vertex edge map.
	 * Note that uniqueness is only guaranteed if the input mapping is in the correct format.
	 */
	public static List<Edge> getEdges(Dictionary<int, Dictionary<int, Edge>> vertexEdgeMap) {
		List<Edge> edges = new List<Edge>();
		foreach(KeyValuePair<int, Dictionary<int, Edge>> entry in vertexEdgeMap) {
			int v1 = entry.Key;
			Dictionary<int, Edge> v1EdgeMap = entry.Value;
			foreach(Edge edge in v1EdgeMap.Values) {
				if(edge.ve1 == v1) { // Only add edges in one direction to prevent duplicates.
					edges.Add(edge);
				}
			}
		}
		return edges;
	}

	/**
	 * Creates edges based on the sorted vertex triangles (list of triangle ids per vertex) and triangles
	 * (mapping from triangle id to triangle vertices).
	 * Returns a nested map containing the edge based on their defining vertices.
	 * Example: Both map[v1][v2] and map[v2][v1] return the edge between vertices v1 and v2.
	 */
	public static Dictionary<int, Dictionary<int, Edge>> getVertexEdgeMap(List<int>[] sortedVertexTriangles, int[] triangles) {
		Dictionary<int, Dictionary<int, Edge>> edgeMap = new Dictionary<int, Dictionary<int, Edge>>();
		for(int vertexInd = 0; vertexInd < sortedVertexTriangles.Length; vertexInd++) {
			List<int> triangleList = sortedVertexTriangles[vertexInd];
			for(int i = 0; i < triangleList.Count; i++) {

				// Get two possibly adjacent triangles in clockwise direction.
				int triangleId = triangleList[i];
				int nextTriangleId = triangleList[(i + 1) % triangleList.Count];

				// Get triangle vertices.
				int triangleBaseInd1 = triangleId * 3;
				int triangleBaseInd2 = nextTriangleId * 3;
				int v11 = triangles[triangleBaseInd1];
				int v12 = triangles[triangleBaseInd1 + 1];
				int v13 = triangles[triangleBaseInd1 + 2];
				int v21 = triangles[triangleBaseInd2];
				int v22 = triangles[triangleBaseInd2 + 1];
				int v23 = triangles[triangleBaseInd2 + 2];

				// Get the vertex indices of the other vertices that are connected to the possibly shared triangle edge.
				int otherVertexClockwiseInd1 = (vertexInd == v11 ? v13 : (vertexInd == v13 ? v12 : v11));
				int otherVertexAntiClockwiseInd2 = (vertexInd == v21 ? v22 : (vertexInd == v22 ? v23 : v21));
				
				// Define edge vertices.
				int ve1 = vertexInd;
				int ve2 = otherVertexClockwiseInd1;
				
				// Add vertices defining the side flaps if the edge is shared between two triangles.
				Edge edge;
				bool edgeSharedByTriangles = (otherVertexClockwiseInd1 == otherVertexAntiClockwiseInd2);
				if(edgeSharedByTriangles) {

					// Skip shared edges where v1 > v2. This prevents shared edges from being added twice.
					if(ve1 > ve2) {
						continue;
					}

					// Initialize edge and flap vertices.
					int vf1 = (vertexInd == v11 ? v12 : (vertexInd == v12 ? v13 : v11));
					int vf2 = (vertexInd == v21 ? v23 : (vertexInd == v22 ? v21 : v22));

					// Add edge with side flaps.
					edge = new Edge(ve1, ve2, vf1, vf2, triangleId, nextTriangleId);
				} else {

					// Add edge without side flaps.
					edge = new Edge(ve1, ve2);
				}

				// Add edge to both ve1 and ve2 edge maps.
				if(!edgeMap.TryGetValue(edge.ve1, out Dictionary<int, Edge> ve1EdgeMap)) {
					ve1EdgeMap = new Dictionary<int, Edge>();
					edgeMap.Add(edge.ve1, ve1EdgeMap);
				}
				ve1EdgeMap.Add(edge.ve2, edge);
				if(!edgeMap.TryGetValue(edge.ve2, out Dictionary<int, Edge> ve2EdgeMap)) {
					ve2EdgeMap = new Dictionary<int, Edge>();
					edgeMap.Add(edge.ve2, ve2EdgeMap);
				}
				ve2EdgeMap.Add(edge.ve1, edge);
			}
		}
		return edgeMap;
	}

	public static Vec3D[] generateMeasurements(Vec3D[] vertexPositions, int numMeasurements) {
		if(numMeasurements > vertexPositions.Length) {
			throw new Exception("Cannot take more measurements than the amount of vertex positions.");
		}

		// Sample measurements following a heuristic to get a fairly good spread over the mesh while still being deterministic.
		Vec3D[] measurements = new Vec3D[vertexPositions.Length];
		int inc = vertexPositions.Length / numMeasurements;
		for(int i = 0; i < vertexPositions.Length; i += inc) {
			measurements[i] = vertexPositions[i].clone();
		}
		return measurements;
	}

	/**
	 * Generates measurements for a triangular sail with the given width and height.
	 * The sail is expected to have its 90 degree corner at vertex position (0, 0, 0),
	 * going into the positive x and y directions.
	 * The measurements will be generated by creating the given number of horizontal lines on the sail,
	 * filling these lines in with measurementDensity measurements per meter.
	 * The measurements will be spaced as equally as possibly, attempting to leaving a gap between the first/last
	 * measurements and their corresponding sides. This mostly just prevents measurements from being placed on
	 * constrained vertices, becoming useless for reconstruction.
	 */
	public static Vec3D[] generateTriangularSailMeasurements(Vec3D[] vertexPositions,
			double sailWidth, double sailHeight, int numMeasurementLines, double measurementDensity) {
		Vec3D[] measurements = new Vec3D[vertexPositions.Length];

		// Populate a triangle with the given number of measurement lines and measurement density.
		double deltaMeasurementHeight = sailHeight / (double) (numMeasurementLines + 1);
		double y = deltaMeasurementHeight;
		for(int i = 0; i < numMeasurementLines; i++) {
			double deltaMeasurementWidth = 1d / measurementDensity;
			double lineWidth = sailWidth - (y / sailHeight) * sailWidth; // Only valid for a triangle.
			int numLineMeasurements = (int) (lineWidth / deltaMeasurementWidth); // Floor to int.
			double x = deltaMeasurementWidth;
			for(int j = 0; j < numLineMeasurements; j++) {

				// Create a measurement for the closest match to the given x-y coordinate in the x-y plane (ignoring z coordinates).
				double bestDistSquare = double.MaxValue;
				int bestDistSquareInd = -1;
				for(int k = 0; k < vertexPositions.Length; k++) {
					Vec3D vertexPos = vertexPositions[k];
					double distSquare = (vertexPos.x - x) * (vertexPos.x - x) + (vertexPos.y - y) * (vertexPos.y - y);
					if(distSquare < bestDistSquare) {
						bestDistSquare = distSquare;
						bestDistSquareInd = k;
					}
				}
				measurements[bestDistSquareInd] = vertexPositions[bestDistSquareInd].clone();
				
				x += deltaMeasurementWidth;
			}
			y += deltaMeasurementHeight;
		}
		return measurements;
	}
}
