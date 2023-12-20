using Accord;
using Accord.Extensions;
using Accord.Extensions.BinaryTree;
using Accord.Math;
using Accord.Math.Decompositions;
using Accord.Math.Integration;
using FEM.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ScottPlot;
using System.Diagnostics;
using System.Drawing;
using Generate = MathNet.Numerics.Generate;

namespace FEM
{
    public static class FEMSolver1DLinear
    {
        //Local to global matrix
        public static int[,] LTOG(int N)
        {
            var T = new int[N - 1, 2];

            for (int i = 0; i < N - 1; ++i)
            {
                for (int j = 0; j < 2; ++j)
                    T[i, j] = i + j;
            }

            return T;
        }

        public static Dictionary<int, (int, int)> GenerateMesh(int N)
        {
            var T = LTOG(N);
            var mesh = new Dictionary<int, (int, int)>();

            for (int e = 0; e < T.GetLength(0); ++e)
            {
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        var k = (T[e, i], T[e, j]);

                        if (!mesh.ContainsKey(e))
                            mesh.Add(e, k);
                    }
                }
            }

            return mesh;
        }

        //Global to local coordinates transform
        public static double Local(double x0, int e, double[] x)
        {
            var h = x[1] - x[0];
            return 2 / h * (x0 - x[e]) - 1;
        }

        //Local to global coordinates transform
        public static double Global(double t, int e, double[] x)
        {
            var h = x[1] - x[0];
            return x[e] + 0.5 * h * (t + 1);
        }

        public static (double, double) Jacobian(double[] x)
        {
            var h = x[1] - x[0];
            return (2 / h, h / 2);
        }

        //Shape functions
        public static double N(double t, int i)
        {
            if (t < -1 || t > 1)
                return 0d;

            switch (i)
            {
                case 0:
                    return 0.5 - 0.5 * t;

                case 1:
                    return 0.5 * t + 0.5;
            }

            throw new ArgumentException();
        }

        //Shape functions' derivatives
        public static double D(double t, int i)
        {
            if (t < -1 || t > 1)
                return 0d;

            switch (i)
            {
                case 0:
                    return -0.5;

                case 1:
                    return 0.5;
            }

            throw new ArgumentException();
        }

        //Element stiffness matrix component
        public static double K(double detJ, double t, int i, int j)
        {
            return detJ * N(t, i) * N(t, j);
        }

        //First order element stiffness matrix component
        public static double K1(double invJ, double detJ, double t, int i, int j)
        {
            return detJ * invJ * N(t, i) * D(t, j);
        }

        //Second order element stiffness matrix component
        public static double K2(double invJ, double detJ, double t, int i, int j)
        {
            return -detJ * invJ * invJ * D(t, i) * D(t, j);
        }

        //Element load vector component
        public static double F(double detJ, double t, int i)
        {
            return detJ * N(t, i);
        }

        public static double Gauss(Func<double, double> f)
        {
            return Integrate.GaussLegendre(f, -1, 1, 1);
        }

        public static List<(double, Func<double, double>)> Solve(string[] equation, double x0, double L, int n)
        {
            var x = Generate.LinearSpaced(n, x0, L);
            var h = x[1] - x[0];

            var mesh = GenerateMesh(n);
            var T = LTOG(n);

            var A = new List<Matrix<double>>();
            var B = new List<Matrix<double>>();

            var H = CreateMatrix.Sparse<double>(n, n); 
            var D = CreateMatrix.Sparse<double>(n, n); 

            var a = new Expression[equation.Length];

            for (int i = 0; i < equation.Length; ++i)
                a[i] = new Expression(equation[i]);

            var jacobian = Jacobian(x);
            var detJ = jacobian.Item2;
            var invJ = jacobian.Item1;

            for (int e = 0; e < T.GetLength(0); ++e)
            {
                var Ae = CreateMatrix.Sparse<double>(2, 2);
                var Be = CreateMatrix.Sparse<double>(2, 2);

                for (int i = 0; i < 2; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        Ae[i, j] = Gauss(t =>
                        {
                            var x0 = Global(t, e, x);

                            return a[0].Calculate(x0) * K2(invJ, detJ, t, i, j) + a[1].Calculate(x0) * K1(invJ, detJ, t, i, j) + a[2].Calculate(x0) * K(detJ, t, i, j);
                        });

                        Be[i, j] = Gauss(t => K(detJ, t, i, j));
                    }
                }

                A.Add(Ae);
                B.Add(Be);
            }

            for (int i = 0; i < mesh.Count; ++i)
            {
                var element = mesh.ElementAt(i);
                var e = element.Key;
                var k = element.Value;

                H.SetSubMatrix(k.Item1, k.Item2, H.SubMatrix(k.Item1, 2, k.Item2, 2).Add(A[e]));
                D.SetSubMatrix(k.Item1, k.Item2, D.SubMatrix(k.Item1, 2, k.Item2, 2).Add(B[e]));
            }

            H = H.RemoveRow(n - 1).RemoveColumn(n - 1).RemoveRow(0).RemoveColumn(0);
            D = D.RemoveRow(n - 1).RemoveColumn(n - 1).RemoveRow(0).RemoveColumn(0);
            var evd = new GeneralizedEigenvalueDecomposition(H.ToArray(), D.ToArray());
            var E = evd.RealEigenvalues;
            var U = evd.Eigenvectors;
            var solutions = new List<(double, double[])>();

            for (int i = 0; i < E.Length; ++i)
                solutions.Add((E[i], U.GetColumn(i)));

            solutions = solutions.OrderBy(x => x.Item1).ToList();
            var result = new List<(double, Func<double, double>)>();

            for (int j = 0; j < solutions.Count; ++j)
            {
                var q_reduced = solutions[j].Item2;
                var q = new double[n];

                for (int i = 1; i < n - 1; ++i)
                    q[i] = q_reduced[i - 1];

                Func<double, double> u = x0 =>
                {
                    var y = 0d;

                    for (int e = 0; e < T.GetLength(0); ++e)
                    {
                        var t = Local(x0, e, x);

                        for (int i = 0; i < 2; ++i)
                            y += q[e + i] * N(t, i);
                    }

                    return y * y;
                };

                result.Add((solutions[j].Item1, u));
            }

            return result;
        }
    }
}
