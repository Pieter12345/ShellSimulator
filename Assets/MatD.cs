using System;

public class MatD {

    private double[,] matrix;
    
    public int numRows { get { return matrix.GetLength(0); } }
    public int numColumns { get { return matrix.GetLength(1); } }

    public MatD(int rows, int columns) {
        this.matrix = new double[rows, columns];
        // TODO - Need to initialize to 0?
    }

    public MatD(double[,] data) {
        this.matrix = data;
    }

    public MatD(params VecD[] rows) {
        int numColumns = (rows.Length == 0 ? 0 : rows[0].length);
        this.matrix = new double[rows.Length, numColumns];
        for(int row = 0; row < rows.Length; row++) {
            if(rows[row].length != numColumns) {
                throw new ArgumentException("All rows that compose a matrix must have the same length.");
            }
            for(int col = 0; col < numColumns; col++) {
                this.matrix[row, col] = rows[row][col];
            }
        }
    }

    public double this[int i, int j] {
        get { return this.matrix[i, j]; }
        set { this.matrix[i, j] = value; }
    }
    
    public static MatD operator +(MatD m) => m;
    public static MatD operator -(MatD m) => new MatD(m.numRows, m.numColumns).sub(m);
    public static MatD operator +(MatD m1, MatD m2) => new MatD(m1.numRows, m1.numColumns).add(m1).add(m2);
    public static MatD operator -(MatD m1, MatD m2) => new MatD(m1.numRows, m1.numColumns).add(m1).sub(m2);
    public static MatD operator *(MatD m, double v) => new MatD(m.numRows, m.numColumns).add(m).mul(v);
    public static MatD operator *(double v, MatD m) => new MatD(m.numRows, m.numColumns).add(m).mul(v);
    public static MatD operator *(MatD m1, MatD m2) => mul(m1, m2);
    public static VecD operator *(VecD v, MatD m) => mul(v, m);
    public static VecD operator *(MatD m, VecD v) => mul(m, v);
    public static MatD operator /(MatD m, double v) => new MatD(m.numRows, m.numColumns).add(m).div(v);

    public MatD add(MatD m) {
        ensureMatricesSameLength(this, m);
        for(int i = 0; i < this.numRows; i++) {
            for(int j = 0; j < this.numColumns; j++) {
                this.matrix[i, j] += m.matrix[i, j];
            }
        }
        return this;
    }

    public MatD sub(MatD m) {
        ensureMatricesSameLength(this, m);
        for(int i = 0; i < this.numRows; i++) {
            for(int j = 0; j < this.numColumns; j++) {
                this.matrix[i, j] -= m.matrix[i, j];
            }
        }
        return this;
    }

    public MatD mul(double val) {
        for(int i = 0; i < this.numRows; i++) {
            for(int j = 0; j < this.numColumns; j++) {
                this.matrix[i, j] *= val;
            }
        }
        return this;
    }

    public static MatD mul(MatD m1, MatD m2) {
        if(m1.numColumns != m2.numRows) {
            throw new ArithmeticException("Cannot perform matrix multiplication. Dimensions not compatible:"
                + " (" + m1.numRows + " x " + m1.numColumns + ") * (" + m2.numRows + " x " + m2.numColumns + ").");
        }
        MatD res = new MatD(m1.numRows, m2.numColumns);
        for(int i = 0; i < res.numRows; i++) {
            for(int j = 0; j < res.numColumns; j++) {
                res[i, j] = 0;
                for(int col1 = 0; col1 < m1.numColumns; col1++) {
                    for(int row2 = 0; row2 < m2.numRows; row2++) {
                        res[i, j] += m1[i, col1] * m2[row2, j];
                    }
                }
            }
        }
        return res;
    }

    public static VecD mul(MatD m, VecD v) {
        if(m.numColumns != v.length) {
            throw new ArithmeticException("Cannot perform matrix multiplication. Dimensions not compatible:"
                    + " (" + m.numRows + " x " + m.numColumns + ") * " + v.length + ".");
        }
        VecD res = new VecD(m.numColumns);
        for(int i = 0; i < res.length; i++) {
            res[i] = 0;
            for(int matCol = 0; matCol < m.numColumns; matCol++) {
                for(int vecInd = 0; vecInd < v.length; vecInd++) {
                    res[i] += m[i, matCol] * v[vecInd];
                }
            }
        }
        return res;
    }

    public static VecD mul(VecD v, MatD m) {
        if(m.numRows != v.length) {
            throw new ArithmeticException("Cannot perform matrix multiplication. Dimensions not compatible: "
                    + v.length + " * (" + m.numRows + " x " + m.numColumns + ").");
        }
        VecD res = new VecD(m.numRows);
        for(int i = 0; i < res.length; i++) {
            res[i] = 0;
            for(int matRow = 0; matRow < m.numRows; matRow++) {
                for(int vecInd = 0; vecInd < v.length; vecInd++) {
                    res[i] += m[matRow, i] * v[vecInd];
                }
            }
        }
        return res;
    }

    public MatD div(double val) {
        for(int i = 0; i < this.numRows; i++) {
            for(int j = 0; j < this.numColumns; j++) {
                this.matrix[i, j] /= val;
            }
        }
        return this;
    }

    public MatD addDiag(double val) {
        if(this.numRows != this.numColumns) {
            throw new ArithmeticException("Cannot perform diagonal addition on a non-square matrix. Dimensions:"
                + " (" + this.numRows + " x " + this.numColumns + ").");
        }
        for(int i = 0; i < this.numRows; i++) {
            this.matrix[i, i] += val;
        }
        return this;
    }

    private static void ensureMatricesSameLength(MatD m1, MatD m2) {
        if(m1.numRows != m2.numRows || m1.numColumns != m2.numColumns) {
            throw new ArithmeticException("Cannot perform binary operation on different size matrices. Dimensions:"
                + " (" + m1.numRows + " x " + m1.numColumns + ") and (" + m2.numRows + " x " + m2.numColumns + ").");
        }
    }

    public bool containsNaN() {
        for(int row = 0; row < this.numRows; row++) {
            for(int col = 0; col < this.numColumns; col++) {
                if(double.IsNaN(this[row, col])) {
                    return true;
                }
            }
        }
        return false;
    }

    public bool containsInf() {
        for(int row = 0; row < this.numRows; row++) {
            for(int col = 0; col < this.numColumns; col++) {
                if(double.IsInfinity(this[row, col])) {
                    return true;
                }
            }
        }
        return false;
    }

    public MatD Clone() {
        return new MatD((double[,]) this.matrix.Clone());
    }

    public override string ToString() {
        string str = (this.numRows > 0 ? this.getRowStr(0) : "");
        for(int i = 1; i < this.numRows; i++) {
            str += "; " + this.getRowStr(i);
        }
        return base.ToString() + "{" + str + "}";
    }

    private string getRowStr(int row) {
        string str = (this.numColumns > 0 ? this.matrix[row, 0].ToString() : "");
        for(int i = 1; i < this.numColumns; i++) {
            str += ", " + this.matrix[row, i];
        }
        return str;
    }
}
