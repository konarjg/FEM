using Accord.Math;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    public static class MathUtils
    {
        public static int[] Modulo(this int[] x, int a)
        {
            return x.Subtract(x.DivideInteger(a).Multiply(a));
        }

        public static int[] DivideInteger(this int[] x, int a)
        {
            var y = new int[x.Length];

            for (int i = 0; i < x.Length; ++i)
                y[i] = x[i] / a;

            return y;
        }

        public static int[] Subtract(this int[] x, int[] b)
        {
            var y = new int[x.Length];

            for (int i = 0; i < x.Length; ++i)
                y[i] = x[i] - b[i];

            return y;
        }
    }
}
