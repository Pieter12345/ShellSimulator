using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;

public class SailMeasurements {
    
    public Vec3D[] measurements { get; }

    public SailMeasurements(Vec3D[] measurements) {
        this.measurements = (Vec3D[]) measurements.Clone();
    }

    public void storeToFile(string fileName) {
        string dirPath = Application.dataPath + "/SailData";
        string filePath = dirPath + "/" + fileName + ".sailmeasurements";
        MonoBehaviour.print("Storing sail measurements to file: " + filePath);
        if(!Directory.Exists(dirPath)) {
            Directory.CreateDirectory(dirPath);
        }
        FileStream fs = (File.Exists(filePath) ? File.OpenWrite(filePath) : File.Create(filePath));
        BinaryWriter writer = new BinaryWriter(fs);

        writer.Write((int) this.measurements.Length);
        foreach(Vec3D pos in this.measurements) {
            if(pos != null) {
                writer.Write(true);
                writer.Write((double) pos.x);
                writer.Write((double) pos.y);
                writer.Write((double) pos.z);
            } else {
                writer.Write(false);
            }
        }
    }

    public static SailMeasurements loadFromFile(string fileName) {
        string filePath = Application.dataPath + "/SailData/" + fileName + ".sailmeasurements";
        MonoBehaviour.print("Loading sail measurements from file: " + filePath);
        FileStream fs = File.OpenRead(filePath);
        BinaryReader reader = new BinaryReader(fs);

        int numVertices = reader.ReadInt32();
        Vec3D[] measurements = new Vec3D[numVertices];
        for(int i = 0; i < numVertices; i++) {
            if(reader.ReadBoolean()) {
                measurements[i] = new Vec3D(reader.ReadDouble(), reader.ReadDouble(), reader.ReadDouble());
            } else {
                measurements[i] = null;
            }
        }
        return new SailMeasurements(measurements);
    }
}
