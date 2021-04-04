using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class SailConfiguration {
    
    public Vec3D[] vertexPositions { get; }
    public int[] triangles { get; }

    public SailConfiguration(Vec3D[] vertexPositions, int[] triangles) {
        this.vertexPositions = (Vec3D[]) vertexPositions.Clone();
        this.triangles = (int[]) triangles.Clone();
    }

    public void storeToFile(string fileName) {
        string dirPath = Application.dataPath + "/SailData";
        string filePath = dirPath + "/" + fileName + ".sailshapedata";
        MonoBehaviour.print("Storing sail shape to file: " + filePath);
        if(!Directory.Exists(dirPath)) {
            Directory.CreateDirectory(dirPath);
        }
        FileStream fs = (File.Exists(filePath) ? File.OpenWrite(filePath) : File.Create(filePath));
        BinaryWriter writer = new BinaryWriter(fs);

        writer.Write((int) this.triangles.Length);
        foreach(int triangle in this.triangles) {
            writer.Write((int) triangle);
        }

        writer.Write((int) this.vertexPositions.Length);
        foreach(Vec3D pos in this.vertexPositions) {
            writer.Write((double) pos.x);
            writer.Write((double) pos.y);
            writer.Write((double) pos.z);
        }
    }

    public static SailConfiguration loadFromFile(string fileName) {
        string filePath = Application.dataPath + "/SailData/" + fileName + ".sailshapedata";
        MonoBehaviour.print("Loading sail shape from file: " + filePath);
        FileStream fs = File.OpenRead(filePath);
        BinaryReader reader = new BinaryReader(fs);

        int numTriangles = reader.ReadInt32();
        int[] triangles = new int[numTriangles];
        for(int i = 0; i < numTriangles; i++) {
            triangles[i] = reader.ReadInt32();
        }

        int numVertices = reader.ReadInt32();
        Vec3D[] vertexPositions = new Vec3D[numVertices];
        for(int i = 0; i < numVertices; i++) {
            vertexPositions[i] = new Vec3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
        }
        return new SailConfiguration(vertexPositions, triangles);
    }
}
