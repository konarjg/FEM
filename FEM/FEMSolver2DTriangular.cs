using Accord;
using Accord.Extensions;
using Accord.Extensions.BinaryTree;
using Accord.Math;
using Accord.Math.Decompositions;
using Accord.Math.Integration;
using DelaunatorSharp;
using FEM.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using Meta.Numerics.Functions;
using Meta.Numerics.Matrices;
using ScottPlot;
using ScottPlot.MinMaxSearchStrategies;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Net.Http.Json;
using System.Numerics;
using System.Runtime.Intrinsics.Arm;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using Generate = MathNet.Numerics.Generate;
using Point = DelaunatorSharp.Point;

namespace FEM
{
    public static class FEMSolver2DTriangular
    {
        public static int[,] LTOG(double[] x, double[] y, int nx, int ny)
        {
            var domain = new List<IPoint>();

            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < ny; ++j)
                    domain.Add(new Point(x[i], y[j]));
            }

            var mesh = new Delaunator(domain.ToArray());
            var trigs = mesh.GetTriangles().ToList();
            var LTOG = new int[trigs.Count, 3];

            for (int i = 0; i < trigs.Count; ++i)
            {
                var trig = trigs[i];
                var points = trig.Points.ToList();

                for (int j = 0; j < points.Count; ++j)
                {
                    var x0 = x.ToList().FindIndex(p => p == points[j].X);
                    var y0 = y.ToList().FindIndex(p => p == points[j].Y);

                    var k = ny * x0 + y0;
                    LTOG[trig.Index, j] = k;
                }
            }

            return LTOG;
        }

        public static Dictionary<(int, int, int), (int, int)> GenerateMesh(int[,] T)
        {
            var mesh = new Dictionary<(int, int, int), (int, int)>();

            for (int e = 0; e < T.GetLength(0); ++e)
            {
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        var k = (T[e, i], T[e, j]);

                        if (!mesh.ContainsKey((e, i, j)))
                            mesh.Add((e, i, j), k);
                    }
                }
            }

            return mesh;
        }

        public static void PlotMesh(double[,] domain, int Nx, int Ny)
        {
            var x = Generate.LinearSpaced(Nx, domain[0, 0], domain[0, 1]);
            var y = Generate.LinearSpaced(Ny, domain[1, 0], domain[1, 1]);
            var mesh = LTOG(x, y, Nx, Ny);

            var plot = new Plot();
            var usedEdges = new List<(int, int)>();

            for (int e = 0; e < mesh.GetLength(0); ++e)
            {
                var k = mesh.GetRow(e);
                var i = k.DivideInteger(Ny);
                var j = k.Modulo(Ny);

                var r = new (double, double)[3];

                for (int n = 0; n < 3; ++n)
                    r[n] = (x[i[n]], y[j[n]]);

                if (!usedEdges.Contains((k[0], k[1])) && !usedEdges.Contains((k[1], k[0])))
                {
                    plot.AddLine(r[0].Item1, r[0].Item2, r[1].Item1, r[1].Item2, Color.Blue, 1);
                    usedEdges.Add((k[0], k[1]));
                }

                if (!usedEdges.Contains((k[1], k[2])) && !usedEdges.Contains((k[2], k[1])))
                {
                    plot.AddLine(r[1].Item1, r[1].Item2, r[2].Item1, r[2].Item2, Color.Blue, 1);
                    usedEdges.Add((k[1], k[2]));
                }

                if (!usedEdges.Contains((k[2], k[0])) && !usedEdges.Contains((k[0], k[2])))
                {
                    plot.AddLine(r[2].Item1, r[2].Item2, r[0].Item1, r[0].Item2, Color.Blue, 1);
                    usedEdges.Add((k[2], k[0]));
                }
            }

            var node = 0;

            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    plot.AddPoint(x[i], y[j], Color.Black, 1);
                    ++node;
                }
            }

            plot.SetAxisLimits(x[0] - 0.5, x[Nx - 1] + 0.5, y[0] - 0.5, y[Ny - 1] + 0.5);
            plot.SaveFig("plot.png");
            Process.Start("explorer.exe", "plot.png");
        }

        //Global to local coordinates transform
        public static (double, double) Local(double[,] invJ, (double, double) R0)
        {
            var T = invJ.Dot(new double[] { R0.Item1, R0.Item2 });
            return (T[0], T[1]);
        }

        //Local to global coordinates transform
        public static (double, double) Global(int[,] LTOG, (double, double) T, int e, (double[], double[]) R)
        {
            var x = R.Item1;
            var y = R.Item2;
            var ζ1 = T.Item1;
            var ζ2 = T.Item2;
            var Ny = y.Length;

            var x0 = ζ1 * x[LTOG[e, 0] / Ny] + ζ2 * x[LTOG[e, 1] / Ny] + (1 - ζ1 - ζ2) * x[LTOG[e, 2] / Ny];
            var y0 = ζ1 * y[LTOG[e, 0] % Ny] + ζ2 * y[LTOG[e, 1] % Ny] + (1 - ζ1 - ζ2) * y[LTOG[e, 2] % Ny];

            return (x0, y0);
        }

        //Shape functions
        public static double N(double ζ1, double ζ2, int i)
        {
            if (ζ1 < -1 || ζ1 > 1 || ζ2 < -1 || ζ2 > 1)
                return 0d;

            switch (i)
            {
                case 0:
                    return ζ1;

                case 1:
                    return ζ2;

                case 2:
                    return 1 - ζ1 - ζ2;
            }

            throw new ArgumentException();
        }

        //Shape functions' t partial derivatives
        public static double Dt(double ζ1, double ζ2, int i)
        {
            switch (i)
            {
                case 0:
                    return 1;

                case 1:
                    return 0;

                case 2:
                    return -1;
            }

            throw new ArgumentException();
        }

        //Shape functions' r partial derivatives
        public static double Dr(double ζ1, double ζ2, int i)
        {
            switch (i)
            {
                case 0:
                    return 0;

                case 1:
                    return 1;

                case 2:
                    return -1;
            }

            throw new ArgumentException();
        }

        //Calculates elemental inverse Jacobi matrix and the Jacobi matrix determinant
        public static (double[,], double) Jacobian(int[,] LTOG, int e, (double[], double[]) R)
        {
            var J = new double[2, 2];
            var x = R.Item1;
            var y = R.Item2;
            var Ny = y.Length;

            var xe = new double[3];
            var ye = new double[3];

            for (int i = 0; i < 3; ++i)
            {
                xe[i] = x[LTOG[e, i] / Ny];
                ye[i] = y[LTOG[e, i] % Ny];
            }

            J[0, 0] = xe[0] - xe[2];
            J[0, 1] = ye[0] - ye[2];
            J[1, 0] = xe[1] - xe[2];
            J[1, 1] = ye[1] - ye[2];

            return (J.Inverse(), J.Determinant());
        }

        //Element stiffness matrix component
        public static double K(double J, double t, double r, int i, int j)
        {
            return J * N(t, r, i) * N(t, r, j);
        }

        //First order element stiffness matrix component
        public static double K1(int var, double J, double[,] invJ, double t, double r, int i, int j)
        {
            var Dx = Dt(t, r, j) * invJ[0, 0] + Dr(t, r, j) * invJ[1, 0];
            var Dy = Dt(t, r, j) * invJ[1, 0] + Dr(t, r, j) * invJ[1, 1];

            if (var == 0)
                return J * N(t, r, i) * Dx;

            return J * N(t, r, i) * Dy;
        }

        //Second order element stiffness matrix component
        public static double K2(int var, double J, double[,] invJ, double t, double r, int i, int j)
        {
            var Dux = Dt(t, r, j) * invJ[0, 0] + Dr(t, r, j) * invJ[0, 1];
            var Duy = Dt(t, r, j) * invJ[1, 0] + Dr(t, r, j) * invJ[1, 1];
            var Dvx = Dt(t, r, i) * invJ[0, 0] + Dr(t, r, i) * invJ[0, 1];
            var Dvy = Dt(t, r, i) * invJ[1, 0] + Dr(t, r, i) * invJ[1, 1];

            if (var == 0)
                return -J * (Dux * Dvx);

            return -J * (Duy * Dvy);
        }

        //Element load vector component
        public static double F(double J, double t, double r, int i)
        {
            return J * N(t, r, i);
        }

        public static double Gauss(Func<double, double, double> f)
        {
            return f(1 / 3d, 1 / 3d) * 0.5;
        }

        public static List<(double, Func<double, double, double>)> Solve(string[] equation, double[,] domain, int nx, int ny)
        {
            var x = Generate.LinearSpaced(nx, domain[0, 0], domain[0, 1]);
            var y = Generate.LinearSpaced(ny, domain[1, 0], domain[1, 1]);
            var dx = x[1] - x[0];
            var dy = y[1] - y[0];

            var T = LTOG(x, y, nx, ny);
            var mesh = GenerateMesh(T);
            PlotMesh(domain, nx, ny);

            var A = new List<Matrix<double>>();
            var B = new List<Matrix<double>>();

            var H = CreateMatrix.Sparse<double>(nx * ny, nx * ny);
            var D = CreateMatrix.Sparse<double>(nx * ny, nx * ny);

            var a = new Expression[equation.Length];

            Console.WriteLine("System Size: {0}", nx * ny);

            for (int i = 0; i < equation.Length; ++i)
                a[i] = new Expression(equation[i]);

            for(int e = 0; e < mesh.Count / 9; ++e)
            {
                var jacobi = Jacobian(T, e, (x, y));
                var invJ = jacobi.Item1;
                var detJ = jacobi.Item2;

                var Ae = CreateMatrix.Sparse<double>(3, 3);
                var Be = CreateMatrix.Sparse<double>(3, 3);

                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        Ae[i, j] = Gauss((t, r) =>
                        {
                            var R0 = Global(T, (t, r), e, (x, y));
                            var ke = a[0].Calculate(R0.Item1, R0.Item2) * K2(0, detJ, invJ, t, r, i, j) + a[1].Calculate(R0.Item1, R0.Item2) * K2(1, detJ, invJ, t, r, i, j) + a[2].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[3].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[4].Calculate(R0.Item1, R0.Item2) * K(detJ, t, r, i, j);
                            return ke;
                        });

                        Be[i, j] = Gauss((t, r) =>
                        {
                            var R0 = Global(T, (t, r), e, (x, y));
                            var me = K(detJ, t, r, i, j);
                            return me;
                        });
                    }
                }

                A.Add(Ae);
                B.Add(Be);
            }

            var boundary = new List<int>();

            for (int i = 0; i < nx * ny; ++i)
            {
                if(i / ny == 0 || i / ny == nx - 1 || i % ny == 0 || i % ny == ny - 1)
                    boundary.Add(i);
            }

            var adjust = 0;
            boundary.Sort();

            for (int i = 0; i < mesh.Count; ++i)
            {
                var element = mesh.ElementAt(i);
                var e = element.Key.Item1;
                var I = element.Key.Item2;
                var J = element.Key.Item3;
                var k = element.Value;

                H[k.Item1, k.Item2] += A[e][I, J];
                D[k.Item1, k.Item2] += B[e][I, J];
            }

            for (int i = 0; i < boundary.Count; ++i)
            {
                H = H.RemoveRow(boundary[i] - adjust).RemoveColumn(boundary[i] - adjust);
                D = D.RemoveRow(boundary[i] - adjust).RemoveColumn(boundary[i] - adjust);
                ++adjust;
            }

            var evd = new GeneralizedEigenvalueDecomposition(H.ToArray(), D.ToArray());
            var E = evd.RealEigenvalues;
            var U = evd.Eigenvectors;

            var solutions = new List<(double, double[])>();

            for (int i = 0; i < E.Length; ++i)
                solutions.Add((E[i], U.GetColumn(i)));

            solutions = solutions.OrderBy(x => x.Item1).ToList();
            var result = new List<(double, Func<double, double, double>)>();

            for (int i = 0; i < solutions.Count; ++i)
            {
                var q_reduced = solutions[i].Item2;
                var q = new double[nx * ny];
                var m = 0;

                for (int k = 0; k < nx * ny; ++k)
                {
                    if (boundary.Contains(k))
                    {
                        q[k] = 0;
                        continue;
                    }

                    q[k] = q_reduced[m];
                    ++m;
                }

                Func<double, double, double> u = (x0, y0) =>
                {
                    var u = 0d;

                    for (int e = 0; e < T.GetLength(0); ++e)
                    {
                        var jacobi = Jacobian(T, e, (x, y));
                        var invJ = jacobi.Item1;

                        var R = Local(invJ, (x0, y0));

                        for (int j = 0; j < 3; ++j)
                            u += q[T[e, j]] * N(R.Item1, R.Item2, j);
                    }

                    return u * u;
                };

                var C = 1d / GaussLegendreRule.Integrate(u, x[0], x[x.Length - 1], y[0], y[y.Length - 1], 3);

                var f = new Func<double, double, double>((x, y) => C * u(x, y));
                result.Add((solutions[i].Item1, f));
            }

            return result;
        }
    }
}
