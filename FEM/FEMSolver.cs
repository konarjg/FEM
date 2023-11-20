using Accord.Extensions;
using Accord.Extensions.BinaryTree;
using Accord.Math;
using Accord.Math.Decompositions;
using Accord.Math.Integration;
using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.LinearAlgebra;
using ScottPlot;
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
    public class FEMSolver
    {
        private static double PhiLinear(double x0, double[] x, int j)
        {
            var h = x[1] - x[0];

            if (j == 0)
            {
                if (x0 > x[1])
                    return 0d;

                return (x[1] - x0) / h;
            }
            else if (j == x.Length - 1)
            {
                if (x0 < x[x.Length - 2])
                    return 0d;

                return (x0 - x[x.Length - 2]) / h;
            }

            if (x0 < x[j - 1] || x0 > x[j + 1])
                return 0d;

            if (x0 <= x[j])
                return (x0 - x[j - 1]) / h;

            return (x[j + 1] - x0) / h;
        }

        private static double DPhiLinear(double x0, double[] x, int j)
        {
            var h = x[1] - x[0];

            if (j == 0)
            {
                if (x0 > x[1])
                    return 0d;

                return -1 / h;
            }
            else if (j == x.Length - 1)
            {
                if (x0 < x[x.Length - 2])
                    return 0d;

                return 1 / h;
            }

            if (x0 < x[j - 1] || x0 > x[j + 1])
                return 0d;

            if (x0 <= x[j])
                return 1 / h;

            return -1 / h;
        }

        private static double Integrate(Func<double, double> f, double a, double b)
        {
            return MathNet.Numerics.Integrate.GaussLegendre(f, a, b);
        }

        private static double Delta(int j, int N, double h)
        {
            if (j == 0 || j == N - 1)
                return 1 / h;

            return 0;
        }

        public static void Solve(double a, int N)
        {
            var x = Generate.LinearSpaced(N, -a, a);
            var h = x[1] - x[0];
            var A = CreateMatrix.Sparse<double>(N, N);
            var B = CreateMatrix.Sparse<double>(N, N);

            Parallel.For(0, N, j =>
            {
                lock (A)
                {
                    lock (B)
                    {
                        var a1 = 0d;
                        var a2 = 0d;
                        var b1 = 0d;
                        var b2 = 0d;

                        if (j - 1 >= 0)
                        {
                            A[j, j - 1] = -Integrate(x0 => DPhiLinear(x0, x, j) * DPhiLinear(x0, x, j - 1), x[j - 1], x[j]) - Integrate(x0 => x0 * x0 * PhiLinear(x0, x, j) * PhiLinear(x0, x, j - 1), x[j - 1], x[j]);
                            B[j, j - 1] = Integrate(x0 => PhiLinear(x0, x, j) * PhiLinear(x0, x, j - 1), x[j - 1], x[j]);
                            a1 = -Integrate(x0 => DPhiLinear(x0, x, j) * DPhiLinear(x0, x, j), x[j - 1], x[j]) - Integrate(x0 => x0 * x0 * PhiLinear(x0, x, j) * PhiLinear(x0, x, j), x[j - 1], x[j]); 
                            b1 = Integrate(x0 => PhiLinear(x0, x, j) * PhiLinear(x0, x, j), x[j - 1], x[j]);
                        }
                        if (j + 1 < N)
                        {
                            A[j, j + 1] = -Integrate(x0 => DPhiLinear(x0, x, j) * DPhiLinear(x0, x, j + 1), x[j], x[j + 1]) - Integrate(x0 => x0 * x0 * PhiLinear(x0, x, j) * PhiLinear(x0, x, j + 1), x[j], x[j + 1]);
                            B[j, j + 1] = Integrate(x0 => PhiLinear(x0, x, j) * PhiLinear(x0, x, j + 1), x[j], x[j + 1]);
                            a2 = -Integrate(x0 => DPhiLinear(x0, x, j) * DPhiLinear(x0, x, j), x[j], x[j + 1]) - Integrate(x0 => x0 * x0 * PhiLinear(x0, x, j) * PhiLinear(x0, x, j), x[j], x[j + 1]);
                            b2 = Integrate(x0 => PhiLinear(x0, x, j) * PhiLinear(x0, x, j), x[j], x[j + 1]);
                        }

                        A[j, j] = Delta(j, N, h) + a1 + a2;
                        B[j, j] = b1 + b2;
                    }
                }
            });

            A = A.Multiply(-1).RemoveRow(N - 1).RemoveColumn(N - 1).RemoveRow(0).RemoveColumn(0);
            B = B.RemoveRow(N - 1).RemoveColumn(N - 1).RemoveRow(0).RemoveColumn(0);

            var evd = new GeneralizedEigenvalueDecomposition(A.ToArray(), B.ToArray());
            var E = evd.RealEigenvalues;
            var U = evd.Eigenvectors;

            var solutions = new List<(double, double[])>();

            for (int i = 0; i < E.Length; ++i)
                solutions.Add((E[i], U.GetColumn(i)));

            solutions = solutions.Where(x => x.Item1 >= 0).OrderBy(x => x.Item1).ToList();
            var plot = new Plot();
            plot.SetAxisLimits(-a, a, -1, 1);
            plot.Title(string.Format("4 First Eigenstates Of Quantum Harmonic Oscillator With Amplitude {0}", a));

            for (int i = 0; i < 4; ++i)
            {
                var n = i;
                var exact = 2 * n + 1;
                var q = solutions[i].Item2;
                var u = new double[N];

                q.CopyTo(u, 1);

                var p = plot.AddSignalXY(x, u);
                p.Label = "n = " + n;
                Console.WriteLine("Energy level {2} Measured: {0:0.000} Exact: {1:0.000}", solutions[i].Item1, exact, n);
            }

            plot.SaveFig("plot.png");
            Process.Start("explorer.exe", "plot.png");
        }
    }
}
