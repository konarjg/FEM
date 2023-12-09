using Accord.Math;
using FEM.Symbolics;
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Threading.Tasks.Sources;

namespace FEM
{
    public static class PerturbationTheory
    {
        public static double IntegrateDiscrete(double[] x, double[] values)
        {
            var N = x.Length;
            var dx = x[1] - x[0];
            var y = 0d;

            for (int i = 0; i < N; ++i)
                y += values[i] * dx;

            return y;
        }

        public static List<(double, double[])> Perturb1D((double, double[])[] H, double[] domain, string perturbation, int order)
        {
            var potential = new Expression(perturbation);
            var V = new double[H[0].Item2.Length];
            var x = Generate.LinearSpaced(H[0].Item2.Length, domain[0], domain[1]);

            for (int i = 0; i < V.Length; ++i)
                V[i] = potential.Calculate(x[i]);

            var solution = new List<(double, double[])>();

            for (int n = 0; n < H.Length; ++n)
            {
                var u = H[n].Item2;
                var energy = H[n].Item1;

                var dE = IntegrateDiscrete(x, u.Multiply(u).Multiply(V));
                var du = new double[u.Length];

                for (int k = 0; k < H.Length; ++k)
                {
                    if (k == n)
                        continue;

                    var inner = IntegrateDiscrete(x, H[k].Item2.Multiply(V).Multiply(u));
                    du.Add(H[k].Item2.Multiply(inner / (H[n].Item1 - H[k].Item1)));
                }

                solution.Add((energy + dE, u.Add(du)));
            }

            if (order > 1)
                return Perturb1D(solution.ToArray(), domain, perturbation, order - 1);

            return solution;
        }
    }
}
