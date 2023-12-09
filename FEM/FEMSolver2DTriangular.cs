using Accord;
using Accord.Extensions;
using Accord.Extensions.BinaryTree;
using Accord.Math;
using Accord.Math.Decompositions;
using Accord.Math.Integration;
using FEM.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using ScottPlot;
using ScottPlot.MinMaxSearchStrategies;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Net.Http.Json;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using Generate = MathNet.Numerics.Generate;

namespace FEM
{
    public static class FEMSolver2DTriangular
    {
        //Local to global matrix
        public static int[,] LTOG(int Nx, int Ny)
        {
            var M = 2 * (Nx - 1) * (Ny - 1);
            var LTOG = new int[M, 3];

            var e = 1;

            for (int ey = 0; ey < Ny - 1; ++ey)
            {
                if (ey % 2 == 0)
                {
                    var i1 = ey;

                    for (int ex = 1; ex < 2 * (Nx - 1); ex += 2)
                    {
                        var i2 = i1 + Ny + 1;

                        LTOG[e - 1, 0] = i1;
                        LTOG[e - 1, 1] = i1 + Ny;
                        LTOG[e - 1, 2] = i1 + 1;

                        LTOG[e, 0] = i2;
                        LTOG[e, 1] = i2 - Ny;
                        LTOG[e, 2] = i2 - 1;

                        e += 2;
                        i1 += Ny;
                    }
                }
                else
                {
                    var i1 = Ny * (Nx - 1) + ey;

                    for (int ex = 1; ex < 2 * (Nx - 1); ex += 2)
                    {
                        var i2 = i1 - Ny + 1;

                        LTOG[e - 1, 0] = i1;
                        LTOG[e - 1, 1] = i1 + 1;
                        LTOG[e - 1, 2] = i1 - Ny;

                        LTOG[e, 0] = i2;
                        LTOG[e, 1] = i2 - 1;
                        LTOG[e, 2] = i2 + Ny;

                        e += 2;
                        i1 -= Ny;
                    }
                }
            }

            return LTOG;
        }

        public static void PlotMesh(double[,] domain, int Nx, int Ny)
        {
            var x = Generate.LinearSpaced(Nx, domain[0, 0], domain[0, 1]);
            var y = Generate.LinearSpaced(Ny, domain[1, 0], domain[1, 1]);
            var dx = x[1] - x[0];
            var dy = y[1] - y[0];
            var mesh = LTOG(Nx, Ny);

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
        public static (double, double) Local((double, double) R0, int e, (double[], double[]) R)
        {
            //todo
            return (0, 0);
        }

        //Local to global coordinates transform
        public static (double, double) Global(int[,] LTOG, (double, double) T, int e, (double[], double[]) R)
        {
            var x = R.Item1;
            var y = R.Item2;
            var t = T.Item1;
            var r = T.Item2;
            var Ny = y.Length;

            var x0 = 0d;
            var y0 = 0d;

            for (int j = 0; j < 3; ++j)
            {
                x0 += N(t, r, j) * x[LTOG[e, j] / Ny];
                y0 += N(t, r, j) * y[LTOG[e, j] % Ny];
            }

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
                    return 1 - t - r;

                case 1:
                    return t;

                case 2:
                    return r;
            }

            throw new ArgumentException();
        }

        //Shape functions' t partial derivatives
        public static double Dt(double t, double r, int i)
        {
            switch (i)
            {
                case 0:
                    return -1;

                case 1:
                    return 1;

                case 2:
                    return 0;
            }

            throw new ArgumentException();
        }

        //Shape functions' r partial derivatives
        public static double Dr(double t, double r, int i)
        {
            switch (i)
            {
                case 0:
                    return -1;

                case 1:
                    return 0;

                case 2:
                    return 1;
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

            var xe = new double[] { x[LTOG[e, 0] / Ny], x[LTOG[e, 1] / Ny], x[LTOG[e, 2] / Ny] };
            var ye = new double[] { y[LTOG[e, 0] % Ny], y[LTOG[e, 1] % Ny], y[LTOG[e, 2] % Ny] };

            J[0, 0] = xe[1] - xe[0];
            J[0, 1] = xe[2] - xe[0];
            J[1, 0] = ye[1] - ye[0];
            J[1, 1] = ye[2] - ye[0];

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
            var Dux = Dt(t, r, j) * invJ[0, 0] + Dr(t, r, j) * invJ[1, 0];
            var Duy = Dt(t, r, j) * invJ[1, 0] + Dr(t, r, j) * invJ[1, 1];
            var Dvx = Dt(t, r, i) * invJ[0, 0] + Dr(t, r, i) * invJ[1, 0];
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

        public static double/*List<(double, Func<double, double, double>)>*/ Solve(string[] equation, double[,] domain, int nx, int ny, bool polar = false)
        {
            var x = Generate.LinearSpaced(nx, domain[0, 0], domain[0, 1]);
            var y = Generate.LinearSpaced(ny, domain[1, 0], domain[1, 1]);
            var dx = x[1] - x[0];
            var dy = y[1] - y[0];

            var T = LTOG(nx, ny);
            PlotMesh(domain, nx, ny);

            var A = CreateMatrix.Sparse<double>(nx * ny, nx * ny);
            var B = CreateMatrix.Sparse<double>(nx * ny, nx * ny);

            var a = new Expression[equation.Length];

            for (int i = 0; i < equation.Length; ++i)
                a[i] = new Expression(equation[i]);

            Parallel.For(0, 2 * (nx - 1) * (ny - 1), e =>
            {
                var jacobi = Jacobian(T, e, (x, y));
                var invJ = jacobi.Item1;
                var detJ = jacobi.Item2;

                Parallel.For(0, 3, i =>
                {
                    Parallel.For(0, 3, j =>
                    {
                        lock (A)
                        {
                            lock (B)
                            {
                                var I = T[e, i];
                                var J = T[e, j];

                                A[I, J] += Gauss((t, r) =>
                                {
                                    var R0 = Global(T, (t, r), e, (x, y));
                                    var ke = a[0].Calculate(R0.Item1, R0.Item2) * K2(0, detJ, invJ, t, r, i, j) + a[1].Calculate(R0.Item1, R0.Item2) * K2(1, detJ, invJ, t, r, i, j) + a[2].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[3].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[4].Calculate(R0.Item1, R0.Item2) * K(detJ, t, r, i, j);
                                    return R0.Item1 * ke;
                                });

                                B[I, J] += Gauss((t, r) => 
                                {
                                    var R0 = Global(T, (t, r), e, (x, y));
                                    var me = K(detJ, t, r, i, j);
                                    return R0.Item1 * me;
                                });
                            }
                        }
                    });
                });
            });

            var boundary = new List<int>();

            Parallel.For(0, nx * ny, i =>
            {
                lock (boundary)
                { 
                    if (i / ny == 0 || i / ny == nx - 1 || i % ny == 0 || i % ny == ny - 1)
                        boundary.Add(i);
                }
            });

            var adjust = 0;
            boundary.Sort();

            for (int i = 0; i < boundary.Count; ++i)
            {
                A = A.RemoveRow(boundary[i] - adjust).RemoveColumn(boundary[i] - adjust);
                B = B.RemoveRow(boundary[i] - adjust).RemoveColumn(boundary[i] - adjust);
                ++adjust;
            }

            var evd = new GeneralizedEigenvalueDecomposition(A.ToArray(), B.ToArray());
            var E = evd.RealEigenvalues;
            var U = evd.Eigenvectors;

            var solutions = new List<(double, double[])>();

            for (int i = 0; i < E.Length; ++i)
                solutions.Add((E[i], U.GetColumn(i)));

            solutions = solutions.OrderBy(x => x.Item1).ToList();
            var result = new List<(double, Func<double, double>)>();

            var exact = new List<((int, int), double)>();
            var Lx = domain[0, 1];
            var Ly = domain[1, 1];

            for (int i = 1; i <= 4; ++i)
            {
                for (int j = 0; j <= 4; ++j)
                {
                    var jln = (i + 0.5 * j - 0.25) * Math.PI;
                    var Er = jln * jln;
                    exact.Add(((i, j), Er));
                }
            }

            exact = exact.OrderBy(x => x.Item2).ToList();
            var mse = 0d;
            var refer = 0d;

            for (int i = 0; i < exact.Count; ++i)
            {
                var measured = solutions[i].Item1;

                mse += 1 / 9d * Math.Pow(measured - exact[i].Item2, 2);
                refer += measured * measured;

                Console.WriteLine("|   Energy level {0} ~E={1:0.000} E={2:0.000}  |", exact[i].Item1, measured / 2, exact[i].Item2);
            }

            return (mse / refer) * 100;
        }
    }
}
