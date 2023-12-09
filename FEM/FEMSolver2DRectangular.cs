using Accord;
using Accord.Extensions;
using Accord.Extensions.BinaryTree;
using Accord.Math;
using Accord.Math.Decompositions;
using Accord.Math.Integration;
using Accord.Statistics;
using FEM.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
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
    public static class FEMSolver2DRectangular
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

            var xe = new double[] { x[LTOG[e, 0] / Ny], x[LTOG[e, 1] / Ny], x[LTOG[e, 2] / Ny], x[LTOG[e, 3] / Ny] };
            var ye = new double[] { y[LTOG[e, 0] % Ny], y[LTOG[e, 1] % Ny], y[LTOG[e, 2] % Ny], y[LTOG[e, 3] % Ny] };

            for (int i = 0; i < 4; ++i)
            {
                x0 += N(t, r, i) * xe[i];
                y0 += N(t, r, i) * ye[i];
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

        public static Func<double, double, (double[,], double)> Jacobian(int[,] LTOG, int e, (double[], double[]) R)
        {
            var J = new double[2, 2];
            var x = R.Item1;
            var y = R.Item2;
            var Ny = y.Length;

            var xye = new double[4, 2];
            var xe = new double[] { x[LTOG[e, 0] / Ny], x[LTOG[e, 1] / Ny], x[LTOG[e, 2] / Ny], x[LTOG[e, 3] / Ny] };
            var ye = new double[] { y[LTOG[e, 0] % Ny], y[LTOG[e, 1] % Ny], y[LTOG[e, 2] % Ny], x[LTOG[e, 3] % Ny] };
            xye = xye.SetColumn(0, xe).SetColumn(1, ye);

            return new Func<double, double, (double[,], double)>((t, r) =>
            {
                var Dxy = new double[2, 4];
                var Dx = new double[4] { Dt(t, r, 0), Dt(t, r, 1), Dt(t, r, 2), Dt(t, r, 3) };
                var Dy = new double[4] { Dr(t, r, 0), Dr(t, r, 1), Dr(t, r, 2), Dr(t, r, 3) };

                Dxy = Dxy.SetRow(0, Dx).SetRow(1, Dy);
                J = Dxy.Dot(xye);

                return (J.Inverse(), J.Determinant());
            });
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
        public static double F(double J, double t, double r, int i, double dx, double dy)
        {
            return J * N(t, r, i);
        }

        public static double Gauss(Func<double, double, double> f)
        {
            var t = new double[2] { -1 / Math.Sqrt(3), 1 / Math.Sqrt(3) };
            var r = new double[2] { -1 / Math.Sqrt(3), 1 / Math.Sqrt(3) };
            var w = new double[2] { 1, 1 };

            var y = 0d;

            for (int i = 0; i < 2; ++i)
            {
                for (int j = 0; j < 2; ++j)
                    y += f(t[i], r[j]) * w[i] * w[j];
            }

            return y;
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

            for (int e = 0; e < (nx - 1) * (ny - 1); ++e)
            {
                var jacobi = Jacobian(T, e, (x, y));

                for (int i = 0; i < 4; ++i)
                {
                    for (int j = 0; j < 4; ++j)
                    {
                        var I = T[e, i];
                        var J = T[e, j];

                        A[I, J] += Gauss((t, r) =>
                        {
                            var jacobi_e = jacobi(t, r);
                            var detJ = jacobi_e.Item2;
                            var invJ = jacobi_e.Item1;

                            var R0 = Global(T, (t, r), e, (x, y));
                            return R0.Item1 * (a[0].Calculate(R0.Item1, R0.Item2) * K2(0, detJ, invJ, t, r, i, j) + a[1].Calculate(R0.Item1, R0.Item2) * K2(1, detJ, invJ, t, r, i, j) + a[2].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[3].Calculate(R0.Item1, R0.Item2) * K1(0, detJ, invJ, t, r, i, j) + a[4].Calculate(R0.Item1, R0.Item2) * K(detJ, t, r, i, j));
                        });

                        B[I, J] += Gauss((t, r) =>
                        {
                            var jacobi_e = jacobi(t, r);
                            var detJ = jacobi_e.Item2;
                            var R0 = Global(T, (t, r), e, (x, y));

                            return R0.Item1 * K(detJ, t, r, i, j);
                        });
                    }
                }
            }

            var boundary = new List<int>();

            for (int i = 0; i < nx * ny; ++i)
            {
                if (i / ny == 0 || i / ny == nx - 1 || i % ny == 0 || i % ny == ny - 1)
                    boundary.Add(i);
            }

            boundary.Sort();
            var adjust = 0;

            for (int i = 0; i < boundary.Count; ++i)
            {
                var j = boundary[i];

                A = A.RemoveRow(j - adjust).RemoveColumn(j - adjust);
                B = B.RemoveRow(j - adjust).RemoveColumn(j - adjust);
                ++adjust;
            }

            var evd = new GeneralizedEigenvalueDecomposition(A.ToArray(), B.ToArray());
            var E = evd.RealEigenvalues;
            var U = evd.Eigenvectors;
            var solutions = new List<(double, double[])>();

            for (int i = 0; i < E.Length; ++i)
                solutions.Add((E[i], U.GetColumn(i)));

            solutions = solutions.OrderBy(x => x.Item1).Where(x => x.Item1 > 0).ToList();
            var result = new List<(double, Func<double, double>)>();

            var exact = new List<((int, int), double)>();
            var Lx = domain[0, 1];
            var Ly = domain[1, 1];

            for (int i = 1; i <= 3; ++i)
            {
                for (int j = 0; j <= 3; ++j)
                {
                    var Er = -Math.Pow(1 / Lx, 2) * (i * i + (j * j) / (double)(i * i));
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

                Console.WriteLine("|   Energy level {0} ~E={1:0.000} E={2:0.000}  |", exact[i].Item1, measured, exact[i].Item2);
            }

            return Math.Sqrt(mse / refer) * 100;
        }
    }
}
