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
    public static class FEMSolver1DQuadratic
    {
        public static List<double[]> PrepareMesh(ref double[] x, ref int N)
        {
            var new_x = new List<double>();
            new_x.AddRange(x);

            var mesh = new List<double[]>();

            for (int i = 0; i < N - 1; ++i)
            {
                var element = new List<double>();

                for (int j = 0; j < 2; ++j)
                    element.Add(x[i + j]);

                var center = (element[0] + element[1]) / 2d;

                if (!new_x.Contains(center))
                    new_x.Add(center);

                element.Insert(1, center);
                mesh.Add(element.ToArray());
            }

            new_x.Sort();
            x = new_x.ToArray();
            N = x.Length;

            return mesh;
        }

        //Local to global matrix
        public static int[,] LTOG(ref double[] x, ref int N)
        {
            var mesh = PrepareMesh(ref x, ref N);
            var T = new int[mesh.Count, 3];

            for (int i = 0; i < mesh.Count; ++i)
            {
                var element = mesh[i];

                for (int j = 0; j < 3; ++j)
                    T[i, j] = x.ToList().FindIndex(p => p == element[j]);
            }

            return T;
        }

        public static Dictionary<int, (int, int)> GenerateMesh(ref int[,] T)
        {
            var mesh = new Dictionary<int, (int, int)>();

            for (int e = 0; e < T.GetLength(0); ++e)
            {
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
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
        public static double Local(double x0, double xe, ref double[] x)
        {
            var h = x[1] - x[0];
            return 1 / h * (x0 - xe) - 1;
        }

        //Local to global coordinates transform
        public static double Global(double t, double xe, ref double[] x)
        {
            var h = x[1] - x[0];
            return xe + h * (t + 1);
        }

        public static (double, double) Jacobian(ref double[] x)
        {
            var h = x[1] - x[0];
            return (1 / h, h);
        }

        //Shape functions
        public static double N(double t, int i)
        {
            if (t < -1 || t > 1)
                return 0d;

            switch (i)
            {
                case 0:
                    return 1 / 3d * t * t - 0.5 * t + 1 / 6d;

                case 1:
                    return -4 / 3d * t * t + 4 / 3d;

                case 2:
                    return t * t + 0.5 * t - 0.5;
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
                    return 2 / 3d * t - 0.5;

                case 1:
                    return -8 / 3d * t;

                case 2:
                    return 2 * t + 0.5;
            }

            throw new ArgumentException();
        }

        //Element stiffness matrix component
        public static double K(double detJ, double t, int i, int j)
        {
            return detJ * N(t, i) * N(t, j);
        }

        //First order element stiffness matrix component
        public static double K1(double t, int i, int j)
        {
            return N(t, i) * D(t, j);
        }

        //Second order element stiffness matrix component
        public static double K2(double invJ, double t, int i, int j)
        {
            return -invJ * D(t, i) * D(t, j);
        }

        //Element load vector component
        public static double F(double detJ, double t, int i)
        {
            return detJ * N(t, i);
        }

        public static double Gauss(Func<double, double> f)
        {
            return Integrate.GaussLegendre(f, -1, 1, 2);
        }

        public static List<(double, Func<double, double>)> Solve(string[] equation, double x0, double L, int n)
        {
            var x = Generate.LinearSpaced(n, x0, L);
            var h = x[1] - x[0];

            var T = LTOG(ref x, ref n);
            var mesh = GenerateMesh(ref T);

            var A = new List<Matrix<double>>();
            var B = new List<Matrix<double>>();

            var H = CreateMatrix.Sparse<double>(n, n); 
            var D = CreateMatrix.Sparse<double>(n, n); 

            var a = new Expression[equation.Length];

            for (int i = 0; i < equation.Length; ++i)
                a[i] = new Expression(equation[i]);

            var jacobian = Jacobian(ref x);
            var detJ = jacobian.Item2;
            var invJ = jacobian.Item1;

            for (int e = 0; e < T.GetLength(0); ++e)
            {
                var Ae = CreateMatrix.Sparse<double>(3, 3);
                var Be = CreateMatrix.Sparse<double>(3, 3);

                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        Ae[i, j] = Gauss(t =>
                        {
                            var x0 = Global(t, x[T[e, 0]], ref x);

                            return a[0].Calculate(x0) * K2(invJ, t, i, j) + a[1].Calculate(x0) * K1(t, i, j) + a[2].Calculate(x0) * K(detJ, t, i, j);
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

                H.SetSubMatrix(k.Item1, k.Item2, H.SubMatrix(k.Item1, 3, k.Item2, 3).Add(A[e]));
                D.SetSubMatrix(k.Item1, k.Item2, D.SubMatrix(k.Item1, 3, k.Item2, 3).Add(B[e]));
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

                    for (int i = 0; i < mesh.Count; ++i)
                    {
                        var t = Local(x0, x[T[i, 0]], ref x);

                        for (int j = 0; j < 3; ++j)
                            y += q[T[i, j]] * N(t, j);
                    }

                    return y * y;
                };

                result.Add((solutions[j].Item1, u));
            }

            return result;
        }
    }
}
