using System;

public class Edge {

    public int ve1 { get; } // Edge vertex 1.
    public int ve2 { get; } // Edge vertex 2.
    public int vf1 { get; } // Flap vertex from triangle 1.
    public int vf2 { get; } // Flap vertex from triangle 2.
    
    public int triangleId1 { get; }
    public int triangleId2 { get; }

    public double undeformedLength { set; get; }

    public Edge(int v1, int v2) : this(v1, v2, -1) {
    }

    public Edge(int v1, int v2, double undeformedLength) : this(v1, v2, -1, -1, -1, -1, undeformedLength) {
    }

    public Edge(int ve1, int ve2, int vf1, int vf2, int triangleId1, int triangleId2) : this(ve1, ve2, vf1, vf2, triangleId1, triangleId2, -1) {
    }

    public Edge(int ve1, int ve2, int vf1, int vf2, int triangleId1, int triangleId2, double undeformedLength) {
        this.ve1 = ve1;
        this.ve2 = ve2;
        this.vf1 = vf1;
        this.vf2 = vf2;
        this.triangleId1 = triangleId1;
        this.triangleId2 = triangleId2;
        this.undeformedLength = undeformedLength;

        if(ve1 == ve2 || ve1 == vf1 || ve1 == vf2 || ve2 == vf1 || ve2 == vf2
                || (vf1 == vf2 && vf1 >= 0) || (triangleId1 == triangleId2 && triangleId1 >= 0)) {
            throw new Exception("An edge cannot be defined using the same vertex twice. Edge: " + this);
        }
    }

    public bool hasSideFlaps() {
        return this.vf1 >= 0 && this.vf2 >= 0;
    }

    public override string ToString() {
        return base.ToString() + "{ve1 = " + this.ve1 + ", ve2 = " + this.ve2 + ", vf1 = " + this.vf1 + ", vf2 = " + this.vf2
            + ", triangleId1 = " + this.triangleId1 + ", triangleId2 = " + this.triangleId2 + "}";
    }
}
