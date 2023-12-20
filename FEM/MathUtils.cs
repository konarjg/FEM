using Accord.Math;
using Accord.Math.Geometry;
using DelaunatorSharp;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    public static class MathUtils
    {
        public static List<(IPoint, IPoint)> GetEdges(this ITriangle trig)
        {
            var points = trig.Points.ToList();
            var edges = new List<(IPoint, IPoint)>();

            var p0 = points[0];
            var p1 = points[1];
            var p2 = points[2];

            edges.Add((p0, p1));
            edges.Add((p1, p2));
            edges.Add((p2, p0));

            return edges;
        }

        public static IPoint GetMiddle(this (IPoint, IPoint) edge)
        {
            var p0 = edge.Item1;
            var p1 = edge.Item2;

            var x = (p0.X + p1.X) / 2;
            var y = (p0.Y + p1.Y) / 2;
            return new Point(x, y);
        }

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
