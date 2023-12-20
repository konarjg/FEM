using Accord;
using Accord.Extensions;
using Accord.Extensions.BinaryTree;
using Accord.Math;
using Accord.Math.Decompositions;
using Accord.Math.Integration;
using Accord.Statistics;
using Accord.Statistics.Distributions.DensityKernels;
using DelaunatorSharp;
using FEM.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
using Meta.Numerics.Functions;
using ScottPlot;
using ScottPlot.MinMaxSearchStrategies;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Net.Http.Headers;
using System.Net.Http.Json;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using Generate = MathNet.Numerics.Generate;

namespace FEM
{
    public static class FEMSolver2DPolar
    {
        //Local to global matrix
        public static int[,] LTOG(int Nx, int Ny)
        {
            var M = (Nx - 1) * (Ny - 1);
            var LTOG = new int[M, 4];

            var e = 0;

            for (int ey = 0; ey < Ny - 1; ++ey)
            {
                for (int ex = 0; ex < Nx - 1; ++ex)
                {
                    var i0 = Ny * ex + ey;

                    LTOG[e, 0] = i0;
                    LTOG[e, 1] = i0 + Ny;
                    LTOG[e, 2] = i0 + Ny + 1;
                    LTOG[e, 3] = i0 + 1;
                    ++e;
                }
            }

            return LTOG;
        }

        public static Dictionary<(int, int, int), (int, int)> GenerateMesh(int[,] T)
        {
            var mesh = new Dictionary<(int, int, int), (int, int)>();

            for (int e = 0; e < T.GetLength(0); ++e)
            {
                for (int i = 0; i < 4; ++i)
                {
                    for (int j = 0; j < 4; ++j)
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
            var mesh = LTOG(Nx, Ny);

            var plot = new Plot();
            var usedEdges = new List<(int, int)>();

            for (int e = 0; e < mesh.GetLength(0); ++e)
            {
                var k = mesh.GetRow(e);
                var i = k.DivideInteger(Ny);
                var j = k.Modulo(Ny);

                var r = new (double, double)[4];

                for (int n = 0; n < 4; ++n)
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

                if (!usedEdges.Contains((k[2], k[3])) && !usedEdges.Contains((k[3], k[2])))
                {
                    plot.AddLine(r[2].Item1, r[2].Item2, r[3].Item1, r[3].Item2, Color.Blue, 1);
                    usedEdges.Add((k[2], k[3]));
                }

                if (!usedEdges.Contains((k[3], k[0])) && !usedEdges.Contains((k[0], k[3])))
                {
                    plot.AddLine(r[3].Item1, r[3].Item2, r[0].Item1, r[0].Item2, Color.Blue, 1);
                    usedEdges.Add((k[3], k[0]));
                }
            }

            var node = 0;

            for (int j = 0; j < Ny; ++j)
            {
                for (int i = 0; i < Nx; ++i)
                {
                    plot.AddPoint(x[i], y[j], Color.Black, 3);
                    ++node;
                }
            }

            plot.SetAxisLimits(x[0] - 0.5, x[Nx - 1] + 0.5, y[0] - 0.5, y[Ny - 1] + 0.5);
            plot.SaveFig("plot.png");
            Process.Start("explorer.exe", "plot.png");
        }

        //Global to local coordinates transform
        public static (double, double) Local(int[,] LTOG,(double, double) R0, int e, (double[], double[]) R)
        {
            var x = R.Item1;
            var y = R.Item2;
            var x0 = R0.Item1;
            var y0 = R0.Item2;
            var dx = x[1] - x[0];
            var dy = y[1] - y[0];

            var xe = x[LTOG[e, 0] / y.Length];
            var ye = y[LTOG[e, 0] % y.Length];

            var t = 2 / dx * (x0 - xe) - 1;
            var r = 2 / dy * (y0 - ye) - 1;

            return (t, r);
        }

        //Local to global coordinates transform
        public static (double, double) Global(int[,] LTOG, (double, double) T, int e, (double[], double[]) R)
        {
            var x = R.Item1;
            var y = R.Item2;
            var t = T.Item1;
            var r = T.Item2;
            
            var Ny = y.Length;

            var dx = x[1] - x[0];
            var dy = y[1] - y[0];

            var x0 = x[LTOG[e, 0] / Ny] + dx / 2 * (t + 1);
            var y0 = y[LTOG[e, 0] % Ny] + dy / 2 * (r + 1);

            return (x0, y0);
        }

        //Shape functions
        public static double N(double t, double r, int i)
        {
            if (t < -1 || t > 1 || r < -1 || r > 1)
                return 0d;

            switch (i)
            {
                case 0:
                    return 0.25 * (1 - t) * (1 - r);

                case 1:
                    return 0.25 * (1 + t) * (1 - r);

                case 2:
                    return 0.25 * (1 + t) * (1 + r);

                case 3:
                    return 0.25 * (1 - t) * (1 + r);
            }

            throw new ArgumentException();
        }

        //Shape functions' t partial derivatives
        public static double Dt(double t, double r, int i)
        {
            switch (i)
            {
                case 0:
                    return -0.25 * (1 - r);

                case 1:
                    return 0.25 * (1 - r);

                case 2:
                    return 0.25 * (1 + r);

                case 3:
                    return -0.25 * (1 + r);
            }

            throw new ArgumentException();
        }

        //Shape functions' r partial derivatives
        public static double Dr(double t, double r, int i)
        {
            switch (i)
            {
                case 0:
                    return -0.25 * (1 - t);

                case 1:
                    return -0.25 * (1 + t);

                case 2:
                    return 0.25 * (1 + t);

                case 3:
                    return 0.25 * (1 - t);
            }

            throw new ArgumentException();
        }

        //Element stiffness matrix component
        public static double K(double J, double t, double r, int i, int j)
        {
            return J * N(t, r, i) * N(t, r, j);
        }

        //First order element stiffness matrix component
        public static double K1(int var, double J, double[] invJ, double t, double r, int i, int j)
        {
            var Dx = Dt(t, r, j) * invJ[0];
            var Dy = Dr(t, r, j) * invJ[1];

            if (var == 0)
                return J * N(t, r, i) * Dx;

            return J * N(t, r, i) * Dy;
        }

        //Second order element stiffness matrix component
        public static double K2(int var, double J, double[] invJ, double t, double r, int i, int j)
        {
            var Dux = Dt(t, r, j) * invJ[0];
            var Duy = Dr(t, r, j) * invJ[1];
            var Dvx = Dt(t, r, i) * invJ[0];
            var Dvy = Dr(t, r, i) * invJ[1];

            if (var == 0)
                return -J * (Dux * Dvx);

            return -J * (Duy * Dvy);
        }

        //Element load vector component
        public static double F(double J, double t, double r, int i, double dx, double dy)
        {
            return J * N(t, r, i);
        }

        public static double Gauss(Func<double, double, double> f)
        {
            return GaussLegendreRule.Integrate(f, -1, 1, -1, 1, 3);
        }

        public static List<(double, Func<double, double, double>)> Solve(string[] equation, double[,] domain, int nx, int ny)
        {
            var x = Generate.LinearSpaced(nx, domain[0, 0], domain[0, 1]);
            var y = Generate.LinearSpaced(ny, domain[1, 0], domain[1, 1]);
            var dx = x[1] - x[0];
            var dy = y[1] - y[0];

            var T = LTOG(nx, ny);
            var mesh = GenerateMesh(T);
            PlotMesh(domain, nx, ny);

            var H = CreateMatrix.Sparse<double>(nx * ny, nx * ny);
            var D = CreateMatrix.Sparse<double>(nx * ny, nx * ny);

            var A = new List<Matrix<double>>();
            var B = new List<Matrix<double>>();

            var a = new Expression[equation.Length];

            for (int i = 0; i < equation.Length; ++i)
                a[i] = new Expression(equation[i]);

            for (int e = 0; e < T.GetLength(0); ++e)
            {
                var invJ = new double[] { 2 / dx, 2 / dy };
                var detJ = dx * dy / 4;

                var Ae = CreateMatrix.Sparse<double>(4, 4);
                var Be = CreateMatrix.Sparse<double>(4, 4);

                for (int i = 0; i < 4; ++i)
                {
                    for (int j = 0; j < 4; ++j)
                    {
                        Ae[i, j] = Gauss((t, r) =>
                        {
                            var R0 = Global(T, (t, r), e, (x, y));

                            return a[0].Calculate(R0.Item1, R0.Item2) * K2(0, detJ, invJ, t, r, i, j) + a[1].Calculate(R0.Item1, R0.Item2) * K2(1, detJ, invJ, t, r, i, j) + a[2].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[3].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[4].Calculate(R0.Item1, R0.Item2) * K(detJ, t, r, i, j);
                        });

                        Be[i, j] = Gauss((t, r) =>
                        {
                            var R0 = Global(T, (t, r), e, (x, y));

                            return K(detJ, t, r, i, j);
                        });
                    }
                }

                A.Add(Ae);
                B.Add(Be);
            }

            var boundary = new List<int>();

            for (int i = 0; i < nx * ny; ++i)
            {
                if (i / ny == 0 || i / ny == nx - 1 || i % ny == 0 || i % ny == ny - 1)
                    boundary.Add(i);
            }

            boundary.Sort();
            var adjust = 0;

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
                var j = boundary[i];

                H = H.RemoveRow(j - adjust).RemoveColumn(j - adjust);
                D = D.RemoveRow(j - adjust).RemoveColumn(j - adjust);
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
                        var R = Local(T, (x0, y0), e, (x, y));

                        for (int j = 0; j < 4; ++j)
                            u += q[T[e, j]] * N(R.Item1, R.Item2, j);
                    }

                    return u * u;
                };

                var C = 1d / GaussLegendreRule.Integrate(u, x[0], x[x.Length - 1], y[0], y[y.Length - 1], 32);

                var f = new Func<double, double, double>((x, y) => C * u(x, y));
                result.Add((solutions[i].Item1, f));
            }

            return result;
        }
    }
}
