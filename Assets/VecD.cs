using System;
using UnityEngine;

public class VecD {

    private double[] vector;
    
    public int length { get { return this.vector.Length; } }
    public double magnitude {
        get {
            double sum = 0;
            for(int i = 0; i < this.vector.Length; i++) {
                sum += this.vector[i] * this.vector[i];
            }
            return Math.Sqrt(sum);
        }
    }

    public VecD(int size) {
        this.vector = new double[size];
    }

    public VecD(params double[] elements) {
        this.vector = elements;
    }

    public VecD(Vector3 vec) {
        this.vector = new double[] {vec.x, vec.y, vec.z};
    }

    public VecD(params VecD[] vecs) {
        int size = 0;
        for(int i = 0; i < vecs.Length; i++) {
            size += vecs[i].length;
        }
        this.vector = new double[size];
        int ind = 0;
        for(int i = 0; i < vecs.Length; i++) {
            for(int j = 0; j < vecs[i].length; j++) {
                this.vector[ind++] = vecs[i][j];
            }
        }
    }

    public double this[int i] {
        get { return this.vector[i]; }
        set { this.vector[i] = value; }
    }
    
    public static VecD operator +(VecD v) => v;
    public static VecD operator -(VecD v) => new VecD(v.length).sub(v);
    public static VecD operator +(VecD v1, VecD v2) => new VecD((double[]) v1.vector.Clone()).add(v2);
    public static VecD operator +(VecD v, double val) => new VecD((double[]) v.vector.Clone()).add(val);
    public static VecD operator +(double val, VecD v) => new VecD((double[]) v.vector.Clone()).add(val);
    public static VecD operator -(VecD v1, VecD v2) => new VecD((double[]) v1.vector.Clone()).sub(v2);
    public static VecD operator -(VecD v, double val) => new VecD((double[]) v.vector.Clone()).sub(val);
    public static VecD operator -(double val, VecD v) => new VecD((double[]) v.vector.Clone()).sub(val);
    public static VecD operator *(VecD v, double val) => new VecD((double[]) v.vector.Clone()).mul(val);
    public static VecD operator *(double val, VecD v) => new VecD((double[]) v.vector.Clone()).mul(val);
    public static VecD operator /(VecD v, double val) => new VecD((double[]) v.vector.Clone()).div(val);

    public VecD add(VecD v) {
        ensureVectorsSameLength(this, v);
        for(int i = 0; i < this.length; i++) {
            this.vector[i] += v.vector[i];
        }
        return this;
    }

    public VecD add(double val) {
        for(int i = 0; i < this.length; i++) {
            this.vector[i] += val;
        }
        return this;
    }

    public VecD sub(VecD v) {
        ensureVectorsSameLength(this, v);
        for(int i = 0; i < this.length; i++) {
            this.vector[i] -= v.vector[i];
        }
        return this;
    }

    public VecD sub(double val) {
        for(int i = 0; i < this.length; i++) {
            this.vector[i] -= val;
        }
        return this;
    }

    public VecD mul(double val) {
        for(int i = 0; i < this.vector.Length; i++) {
            this.vector[i] *= val;
        }
        return this;
    }

    public VecD div(double val) {
        for(int i = 0; i < this.vector.Length; i++) {
            this.vector[i] /= val;
        }
        return this;
    }

    public double[] asDoubleArray() {
        return (double[]) this.vector.Clone();
    }

    public override string ToString() {
        string str = (this.vector.Length > 0 ? this.vector[0].ToString() : "");
        for(int i = 1; i < this.vector.Length; i++) {
            str += ", " + this.vector[i];
        }
        return base.ToString() + "{" + str + "}";
    }

    public static double dot(VecD v1, VecD v2) {
        ensureVectorsSameLength(v1, v2);
        double ret = 0;
        for(int i = 0; i < v1.length; i++) {
            ret += v1[i] * v2[i];
        }
        return ret;
    }

    private static void ensureVectorsSameLength(VecD v1, VecD v2) {
        if(v1.length != v2.length) {
            throw new ArithmeticException("Cannot perform binary operation on different size vectors. Sizes: " + v1.length + " and " + v2.length + ".");
        }
    }
}
