using pardiso;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

class ParDiSoLib {
    
    public static VecD sparseDirectSymmetricSolve(MathNet.Numerics.LinearAlgebra.Double.SparseMatrix mat, VecD vec) {

        // Documentation: https://pardiso-project.org/manual/manual.pdf

        //// Example from documentation.
        //MatD matd = new MatD(new double[,] {
        //    {7, 0, 1, 0, 0, 2, 7, 0},
        //    {0, -4, 8, 0, 2, 0, 0, 0},
        //    {0, 0, 1, 0, 0, 0, 0, 5},
        //    {0, 0, 0, 7, 0, 0, 9, 0},
        //    {0, 0, 0, 0, 5, -1, 5, 0},
        //    {0, 0, 0, 0, 0, 0, 0, 5},
        //    {0, 0, 0, 0, 0, 0, 11, 0},
        //    {0, 0, 0, 0, 0, 0, 0, 5},
        //});
        //vec = new VecD(1, 1, 1, 1, 1, 1, 1, 1);
        //mat = new MathNet.Numerics.LinearAlgebra.Double.SparseMatrix(matd.numRows, matd.numColumns);
        //for(int row = 0; row < matd.numRows; row++) {
        //    for(int col = 0; col < matd.numColumns; col++) {
        //        if(matd[row, col] != 0d) {
        //            mat[row, col] = matd[row, col];
        //        }
        //    }
        //}

        // Input validation.
        if(mat.RowCount != mat.ColumnCount) {
            throw new Exception("Matrix must be square.");
        }
        if(vec.length != mat.RowCount) {
            throw new Exception("Vector size must match matrix dimensions.");
        }
        if(vec.containsInf() || vec.containsNaN()) {
            throw new Exception("Inf and NaN values are not supported in the vector.");
        }

        // Create compressed sparse row format.
        IEnumerable<Tuple<int, int, double>> enumerable = mat.Storage.EnumerateNonZeroIndexed();
        IEnumerator<Tuple<int, int, double>> enumerator = enumerable.GetEnumerator();
        List<int> iaList = new List<int>();
        List<int> jaList = new List<int>();
        List<double> aList = new List<double>();
        int mtype = -2; // Real symmetric indefinite, diagonal or Buch-Kaufman pivoting.
        while(enumerator.MoveNext()) {
            Tuple<int, int, double> tuple = enumerator.Current;
            int row = tuple.Item1;
            int col = tuple.Item2;
            double val = tuple.Item3;

            // Skip left bottom triangle matrix. We're telling the algorithm that it's symmetric, so it won't need this.
            if(col < row) {
                continue;
            }

            // Insert diagonal zeros if needed.
            while(row >= iaList.Count && col != iaList.Count) {
                aList.Add(0d);
                iaList.Add(aList.Count); // Store 1-indexed diagonal value index.
                jaList.Add(iaList.Count); // Store 1-indexed column index.
            }

            // Insert value.
            aList.Add(val);
            jaList.Add(col + 1); // 1-indexed.
            if(row == col) {
                iaList.Add(aList.Count); // Store 1-indexed diagonal value index.
            }
        }
        iaList.Add(aList.Count + 1); // This would be the starting index of the next unexisting row. Without it, PARDISO crashes.
        int n = mat.RowCount;
        int[] ia = iaList.ToArray();
        int[] ja = jaList.ToArray();
        double[] a = aList.ToArray();

        // Define b and x in linear system: mat * x = b.
        double[] b = vec.asDoubleArray();
        double[] x = new double[n];
        int nrhs = 1; // Number of right hand sides.
        // Internal solver memory pointer pt, 32-bit: int pt[64]; 64-bit: long int pt[64] or void *pt[64] should be OK on both architectures void *pt[64];
        IntPtr[] pt = new IntPtr[64];

        // Pardiso control parameters.
        int[] iparm = new int[64];
        int maxfct, mnum, phase, error, msglvl;

        // Auxiliary variables.
        int i;
        double[] ddum = new double[1]; /* Double dummy */
        int[] idum = new int[1]; /* Integer dummy. */
                                    /* ----------------------------------------------------------------- */
                                    /* .. Setup Pardiso control parameters. */
                                    /* ----------------------------------------------------------------- */
        for(i = 0; i < 64; i++) {
            iparm[i] = 0;
        }
        iparm[0] = 1; /* No solver default */
        iparm[1] = 2; /* Fill-in reordering from METIS */
                        /* Numbers of processors, value of OMP_NUM_THREADS */
        iparm[2] = 1;
        iparm[3] = 0; /* No iterative-direct algorithm */
        iparm[4] = 0; /* No user fill-in reducing permutation */
        iparm[5] = 0; /* Write solution into x */
        iparm[6] = 0; /* Not in use */
        iparm[7] = 2; /* Max numbers of iterative refinement steps */
        iparm[8] = 0; /* Not in use */
        iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
        iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
        iparm[11] = 0; /* Not in use */
        iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off
                    * (default for symmetric). Try iparm[12] = 1 in case of
                    *  inappropriate accuracy */
        iparm[13] = 0; /* Output: Number of perturbed pivots */
        iparm[14] = 0; /* Not in use */
        iparm[15] = 0; /* Not in use */
        iparm[16] = 0; /* Not in use */
        iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1; /* Output: Mflops for LU factorization */
        iparm[19] = 0; /* Output: Numbers of CG Iterations */
        maxfct = 1; /* Maximum number of numerical factorizations. */
        mnum = 1; /* Which factorization to use. */
        msglvl = 1; /* Print statistical information in file */
        error = 0; /* Initialize error flag */
                    /* ----------------------------------------------------------------- */
                    /* .. Initialize the internal solver memory pointer. This is only */
                    /* necessary for the FIRST call of the PARDISO solver. */
                    /* ----------------------------------------------------------------- */
        for(i = 0; i < 64; i++) {
            pt[i] = IntPtr.Zero;
        }
        /* ----------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        /* ----------------------------------------------------------------- */
        phase = 11;
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        if(error != 0) {
            throw new Exception("ERROR during symbolic factorization: " + error);
        }
        //MonoBehaviour.print("\nReordering completed ... ");
        //MonoBehaviour.print("Number of nonzeros in factors = " + iparm[17]);
        //MonoBehaviour.print("Number of factorization MFLOPS = " + iparm[18]);

        /* ----------------------------------------------------------------- */
        /* .. Numerical factorization. */
        /* ----------------------------------------------------------------- */
        phase = 22;
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        if(error != 0) {
            throw new Exception("ERROR during numerical factorization: " + error);
        }
        //MonoBehaviour.print("\nFactorization completed ... ");

        /* ----------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* ----------------------------------------------------------------- */
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, b, x, ref error);
        if(error != 0) {
            throw new Exception("ERROR during solution: " + error);
        }

        //MonoBehaviour.print("\nSolve completed ... ");
        //MonoBehaviour.print("\nThe solution of the system is: ");
        //for(i = 0; i < n; i++) {
        //    MonoBehaviour.print(" x [" + i + "] = " + x[i]);
        //}

        /* ----------------------------------------------------------------- */
        /* .. Termination and release of memory. */
        /* ----------------------------------------------------------------- */
        phase = -1; /* Release internal memory. */
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, ddum, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        
        // Return the result.
        return new VecD(x);
    }
    
    public static VecD sparseDirectSymmetricSolve(MatD mat, VecD vec) {

        // Documentation: https://pardiso-project.org/manual/manual.pdf

        // Example from documentation.
        //mat = new MatD(new double[,] {
        //    {7, 0, 1, 0, 0, 2, 7, 0},
        //    {0, -4, 8, 0, 2, 0, 0, 0},
        //    {0, 0, 1, 0, 0, 0, 0, 5},
        //    {0, 0, 0, 7, 0, 0, 9, 0},
        //    {0, 0, 0, 0, 5, -1, 5, 0},
        //    {0, 0, 0, 0, 0, 0, 0, 5},
        //    {0, 0, 0, 0, 0, 0, 11, 0},
        //    {0, 0, 0, 0, 0, 0, 0, 5},
        //});
        //vec = new VecD(1, 1, 1, 1, 1, 1, 1, 1);

        // Input validation.
        if(mat.numRows != mat.numColumns) {
            throw new Exception("Matrix must be square.");
        }
        if(vec.length != mat.numRows) {
            throw new Exception("Vector size must match matrix dimensions.");
        }
        if(mat.containsInf() || mat.containsNaN()) {
            throw new Exception("Inf and NaN values are not supported in the matrix.");
        }
        if(vec.containsInf() || vec.containsNaN()) {
            throw new Exception("Inf and NaN values are not supported in the vector.");
        }

        // Convert symmatric matrix to sparse matrix format.
        List<int> iaList = new List<int>();
        List<int> jaList = new List<int>();
        List<double> aList = new List<double>();
        int mtype = -2; // Real symmetric indefinite, diagonal or Buch-Kaufman pivoting.
        for(int row = 0; row < mat.numRows; row++) {
            for(int col = row; col < mat.numRows; col++) {
                double val = mat[row, col];
                if(val != 0d || row == col) { // Create entry for non-zero and matrix diagonal elements.
                    jaList.Add(col + 1); // 1-indexed.
                    aList.Add(val);
                    if(iaList.Count == row) {
                        iaList.Add(aList.Count); // 1-indexed.
                    }
                }
            }
        }
        iaList.Add(aList.Count + 1); // This would be the starting index of the next unexisting row. Without it, PARDISO crashes.
        int n = mat.numRows;
        int[] ia = iaList.ToArray();
        int[] ja = jaList.ToArray();
        double[] a = aList.ToArray();

        // Define b and x in linear system: mat * x = b.
        double[] b = vec.asDoubleArray();
        double[] x = new double[n];
        int nrhs = 1; // Number of right hand sides.
        // Internal solver memory pointer pt, 32-bit: int pt[64]; 64-bit: long int pt[64] or void *pt[64] should be OK on both architectures void *pt[64];
        IntPtr[] pt = new IntPtr[64];

        // Pardiso control parameters.
        int[] iparm = new int[64];
        int maxfct, mnum, phase, error, msglvl;

        // Auxiliary variables.
        int i;
        double[] ddum = new double[1]; /* Double dummy */
        int[] idum = new int[1]; /* Integer dummy. */
                                    /* ----------------------------------------------------------------- */
                                    /* .. Setup Pardiso control parameters. */
                                    /* ----------------------------------------------------------------- */
        for(i = 0; i < 64; i++) {
            iparm[i] = 0;
        }
        iparm[0] = 1; /* No solver default */
        iparm[1] = 2; /* Fill-in reordering from METIS */
                        /* Numbers of processors, value of OMP_NUM_THREADS */
        iparm[2] = 1;
        iparm[3] = 0; /* No iterative-direct algorithm */
        iparm[4] = 0; /* No user fill-in reducing permutation */
        iparm[5] = 0; /* Write solution into x */
        iparm[6] = 0; /* Not in use */
        iparm[7] = 2; /* Max numbers of iterative refinement steps */
        iparm[8] = 0; /* Not in use */
        iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
        iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
        iparm[11] = 0; /* Not in use */
        iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off
                    * (default for symmetric). Try iparm[12] = 1 in case of
                    *  inappropriate accuracy */
        iparm[13] = 0; /* Output: Number of perturbed pivots */
        iparm[14] = 0; /* Not in use */
        iparm[15] = 0; /* Not in use */
        iparm[16] = 0; /* Not in use */
        iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1; /* Output: Mflops for LU factorization */
        iparm[19] = 0; /* Output: Numbers of CG Iterations */
        maxfct = 1; /* Maximum number of numerical factorizations. */
        mnum = 1; /* Which factorization to use. */
        msglvl = 1; /* Print statistical information in file */
        error = 0; /* Initialize error flag */
                    /* ----------------------------------------------------------------- */
                    /* .. Initialize the internal solver memory pointer. This is only */
                    /* necessary for the FIRST call of the PARDISO solver. */
                    /* ----------------------------------------------------------------- */
        for(i = 0; i < 64; i++) {
            pt[i] = IntPtr.Zero;
        }
        /* ----------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        /* ----------------------------------------------------------------- */
        phase = 11;
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        if(error != 0) {
            throw new Exception("ERROR during symbolic factorization: " + error);
        }
        //MonoBehaviour.print("\nReordering completed ... ");
        //MonoBehaviour.print("Number of nonzeros in factors = " + iparm[17]);
        //MonoBehaviour.print("Number of factorization MFLOPS = " + iparm[18]);

        /* ----------------------------------------------------------------- */
        /* .. Numerical factorization. */
        /* ----------------------------------------------------------------- */
        phase = 22;
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        if(error != 0) {
            throw new Exception("ERROR during numerical factorization: " + error);
        }
        //MonoBehaviour.print("\nFactorization completed ... ");

        /* ----------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* ----------------------------------------------------------------- */
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, b, x, ref error);
        if(error != 0) {
            throw new Exception("ERROR during solution: " + error);
        }

        //MonoBehaviour.print("\nSolve completed ... ");
        //MonoBehaviour.print("\nThe solution of the system is: ");
        //for(i = 0; i < n; i++) {
        //    MonoBehaviour.print(" x [" + i + "] = " + x[i]);
        //}

        /* ----------------------------------------------------------------- */
        /* .. Termination and release of memory. */
        /* ----------------------------------------------------------------- */
        phase = -1; /* Release internal memory. */
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, ddum, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        
        // Return the result.
        return new VecD(x);
    }

    public static int runTest() {
        /* Matrix data. */
        int n = 8;
        int[] ia/*[9]*/ = new int[] { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
        int[] ja/*[18]*/ = new int[] { 1, 3, 6, 7,
                                        2, 3, 5,
                                        3, 8,
                                        4, 7,
                                        5, 6, 7,
                                        6, 8,
                                        7,
                                        8 };
        double[] a/*[18]*/ = new double[] { 7.0, 1.0, 2.0, 7.0,
                                            -4.0, 8.0, 2.0,
                                            1.0, 5.0,
                                            7.0, 9.0,
                                            5.0, 1.0, 5.0,
                                            -1.0, 5.0,
                                            11.0,
                                            5.0 };
        int mtype = -2; /* Real symmetric matrix */
                        /* RHS and solution vectors. */
        double[] b = new double[8];
        double[] x = new double[8];
        int nrhs = 1; /* Number of right hand sides. */
                        /* Internal solver memory pointer pt, */
                        /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
                        /* or void *pt[64] should be OK on both architectures */
                        /* void *pt[64]; */
        IntPtr[] pt = new IntPtr[64];
        /* Pardiso control parameters. */
        int[] iparm = new int[64];
        int maxfct, mnum, phase, error, msglvl;
        /* Auxiliary variables. */
        int i;
        double[] ddum = new double[1]; /* Double dummy */
        int[] idum = new int[1]; /* Integer dummy. */
                                    /* ----------------------------------------------------------------- */
                                    /* .. Setup Pardiso control parameters. */
                                    /* ----------------------------------------------------------------- */
        for(i = 0; i < 64; i++) {
            iparm[i] = 0;
        }
        iparm[0] = 1; /* No solver default */
        iparm[1] = 2; /* Fill-in reordering from METIS */
                        /* Numbers of processors, value of OMP_NUM_THREADS */
        iparm[2] = 1;
        iparm[3] = 0; /* No iterative-direct algorithm */
        iparm[4] = 0; /* No user fill-in reducing permutation */
        iparm[5] = 0; /* Write solution into x */
        iparm[6] = 0; /* Not in use */
        iparm[7] = 2; /* Max numbers of iterative refinement steps */
        iparm[8] = 0; /* Not in use */
        iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
        iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
        iparm[11] = 0; /* Not in use */
        iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off
                    * (default for symmetric). Try iparm[12] = 1 in case of
                    *  inappropriate accuracy */
        iparm[13] = 0; /* Output: Number of perturbed pivots */
        iparm[14] = 0; /* Not in use */
        iparm[15] = 0; /* Not in use */
        iparm[16] = 0; /* Not in use */
        iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1; /* Output: Mflops for LU factorization */
        iparm[19] = 0; /* Output: Numbers of CG Iterations */
        maxfct = 1; /* Maximum number of numerical factorizations. */
        mnum = 1; /* Which factorization to use. */
        msglvl = 1; /* Print statistical information in file */
        error = 0; /* Initialize error flag */
                    /* ----------------------------------------------------------------- */
                    /* .. Initialize the internal solver memory pointer. This is only */
                    /* necessary for the FIRST call of the PARDISO solver. */
                    /* ----------------------------------------------------------------- */
        for(i = 0; i < 64; i++) {
            pt[i] = IntPtr.Zero;
        }
        /* ----------------------------------------------------------------- */
        /* .. Reordering and Symbolic Factorization. This step also allocates */
        /* all memory that is necessary for the factorization. */
        /* ----------------------------------------------------------------- */
        phase = 11;
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        if(error != 0) {
            MonoBehaviour.print("\nERROR during symbolic factorization: " + error);
            return (1);
        }
        MonoBehaviour.print("\nReordering completed ... ");
        MonoBehaviour.print("Number of nonzeros in factors = " + iparm[17]);
        MonoBehaviour.print("Number of factorization MFLOPS = " + iparm[18]);
        /* ----------------------------------------------------------------- */
        /* .. Numerical factorization. */
        /* ----------------------------------------------------------------- */
        phase = 22;
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        if(error != 0) {
            MonoBehaviour.print("\nERROR during numerical factorization: " + error);
            return (2);
        }
        MonoBehaviour.print("\nFactorization completed ... ");
        /* ----------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* ----------------------------------------------------------------- */
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
                        /* Set right hand side to one. */
        for(i = 0; i < n; i++) {
            b[i] = 1;
        }
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, a, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, b, x, ref error);
        if(error != 0) {
            MonoBehaviour.print("\nERROR during solution: " + error);
            return (3);
        }
        MonoBehaviour.print("\nSolve completed ... ");
        MonoBehaviour.print("\nThe solution of the system is: ");
        for(i = 0; i < n; i++) {
            MonoBehaviour.print(" x [" + i + "] = " + x[i]);
        }
        /* ----------------------------------------------------------------- */
        /* .. Termination and release of memory. */
        /* ----------------------------------------------------------------- */
        phase = -1; /* Release internal memory. */
        Pardiso.pardiso(pt, ref maxfct, ref mnum, ref mtype, ref phase,
            ref n, ddum, ia, ja, idum, ref nrhs,
            iparm, ref msglvl, ddum, ddum, ref error);
        MonoBehaviour.print("TEST PASSED");
        return 0;
    }
}
