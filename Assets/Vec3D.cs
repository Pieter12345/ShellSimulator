using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Vec3D : VecD {
    
    public double x { get { return this[0]; } set { this[0] = value; } }
    public double y { get { return this[1]; } set { this[1] = value; } }
    public double z { get { return this[2]; } set { this[2] = value; } }

    public static Vec3D zero { get { return new Vec3D(0, 0, 0); } }
    
    public Vec3D() : base(3) {
    }
    
    public Vec3D(double x, double y, double z) : base(x, y, z) {
    }

    public Vec3D(Vec3D vec) : base(vec.x, vec.y, vec.z) {
    }

    public Vec3D(Vector3 vec) : base(vec.x, vec.y, vec.z) {
    }
    
    public static Vec3D operator +(Vec3D v) => v;
    public static Vec3D operator -(Vec3D v) => new Vec3D().sub(v);
    public static Vec3D operator +(Vec3D v1, Vec3D v2) => new Vec3D(v1).add(v2);
    public static VecD operator +(Vec3D v, double val) => new Vec3D(v).add(val);
    public static VecD operator +(double val, Vec3D v) => new Vec3D(v).add(val);
    public static Vec3D operator -(Vec3D v1, Vec3D v2) => new Vec3D(v1).sub(v2);
    public static VecD operator -(Vec3D v, double val) => new Vec3D(v).sub(val);
    public static VecD operator -(double val, Vec3D v) => new Vec3D(v).sub(val);
    public static Vec3D operator *(Vec3D v, double val) => new Vec3D(v).mul(val);
    public static Vec3D operator *(double val, Vec3D v) => new Vec3D(v).mul(val);
    public static Vec3D operator /(Vec3D v, double val) => new Vec3D(v).div(val);

    public Vec3D add(Vec3D v) {
        this.x += v.x;
        this.y += v.y;
        this.z += v.z;
        return this;
    }

    public new Vec3D add(double val) {
        this.x += val;
        this.y += val;
        this.z += val;
        return this;
    }

    public Vec3D sub(Vec3D v) {
        this.x -= v.x;
        this.y -= v.y;
        this.z -= v.z;
        return this;
    }

    public new Vec3D sub(double val) {
        this.x -= val;
        this.y -= val;
        this.z -= val;
        return this;
    }

    public new Vec3D mul(double val) {
        this.x *= val;
        this.y *= val;
        this.z *= val;
        return this;
    }

    public new Vec3D div(double val) {
        this.x /= val;
        this.y /= val;
        this.z /= val;
        return this;
    }

    public static Vec3D cross(Vec3D v1, Vec3D v2) {
        return new Vec3D(v1.y * v2.z - v2.y * v1.z, -v1.x * v2.z + v2.x * v1.z, v1.x * v2.y - v2.x * v1.y);
    }
}
