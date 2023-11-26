using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Complex = System.Numerics.Complex;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using MathNet.Numerics;
using MathNet.Numerics.Differentiation;
using System.Net.Http.Headers;

namespace FEM
{
    public static class SpectralBasis
    {
        public static double Laguerre(int n, double x)
        {
            if (n == 0)
                return 1;
            else if (n == 1)
                return 1 - x;

            var L2 = 1d;
            var L1 = 1 - x;

            for (int i = 2; i <= n; ++i)
            {
                var L = ((2 * i - 1 - x) * L1 - (i - 1) * L2) / i;
                L2 = L1;
                L1 = L;
            }

            return L1;
        }

        //Legendre polynomials of order l
        public static double Legendre(int l, double x)
        {
            var a = 1d / (Math.Pow(2, l) * SpecialFunctions.Factorial(l));
            var f = new Func<double, double>(x => Math.Pow(x * x - 1, l));

            if (l == 0)
                return a * f(x);

            var df = Differentiate.Derivative(f, x, l);

            return a * df;
        }

        //Associated Legendre polynomials of orders l, m
        public static double AssociatedLegendre(int l, int m, double x)
        {
            var a = Math.Pow(-1d, m) / (Math.Pow(2d, l) * SpecialFunctions.Factorial(l)) * Math.Pow(1 - x * x, m / 2d);
            var f = new Func<double, double>(x => Math.Pow(x * x - 1, l));

            if (l == 0 && m == 0)
                return a * f(x);

            var df = Differentiate.Derivative(f, x, l + m);
            return a * df;
        }

        public static double Chebyshev(int n, double x)
        {
            if (n == 0)
                return 1;
            else if (n == 1)
                return x;

            var T2 = 1d;
            var T1 = x;

            for (int i = 2; i <= n; ++i)
            {
                var Ti = 2 * x * T1 - T2;
                T2 = T1;
                T1 = Ti;
            }

            return T1;
        }

        public static double ChebyshevDerivative(int n, double x, int order)
        {
            return Differentiate.Derivative(new Func<double, double>(x => Chebyshev(n, x)), x, order);
        }

        public static Complex SphericalHarmonic(int l, int m, double phi, double theta)
        {
            return Complex.Exp(Complex.ImaginaryOne * m * phi) * AssociatedLegendre(l, m, Math.Cos(theta));
        }

        public static double Hermite(int n, double x, bool normalized = false)
        {
            if (!normalized)
            {
                var a = 1 / Math.Sqrt(Math.Pow(2d, n) * SpecialFunctions.Factorial(n) * Math.Sqrt(Math.PI));
                return a * Hermite(n, x, true);
            }

            if (n == 0)
                return 1;
            else if (n == 1)
                return 2 * x;

            var H2 = 1d;
            var H1 = 2 * x;

            for (int i = 2; i <= n; ++i)
            {
                var H = 2 * x * H1 - 2 * (i - 1) * H2;
                H2 = H1;
                H1 = H;
            }

            return H1;
        }

        public static double HermiteDerivative(int n, double x, int order)
        {
            switch (order)
            {
                case 0:
                    return Hermite(n, x);

                case 1:
                    if (n == 0)
                        return 0;
                    else if (n == 1)
                        return 2;

                    return 2 * n * Hermite(n - 1, x);

                case 2:
                    if (n == 0 || n == 1)
                        return 0;

                    return 4 * n * (n - 1) * Hermite(n - 2, x);
            }

            return 0;
        }
    }
}
