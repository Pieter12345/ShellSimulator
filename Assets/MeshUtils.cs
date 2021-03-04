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
}
