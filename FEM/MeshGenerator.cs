using ScottPlot;
using ScottPlot.Drawing;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Point = DelaunatorSharp.Point;

namespace FEM
{
    public static class MeshGenerator
    {
        public static void TestRectangular(double[,] domain, int N)
        {
            var a = domain[0, 1] - domain[0, 0];
            var b = domain[1, 1] - domain[1, 0];

            var c = 2 * a / N;
            var d = 2 * b / N;

            var plot = new Plot();
            plot.SetAxisLimits(domain[0, 0] - 0.5, domain[0, 1] + 0.5, domain[1, 0] - 0.5, domain[1, 1] + 0.5);

            var mesh = new List<Point>();
            var rectangles = new List<Point[]>();

            var n = 0;

            for (int j = 0; j < N/2; ++j)
            {
                for (int i = 0; i < N/2; ++i)
                {
                    var x1 = domain[0, 0] + i * c;
                    var x2 = x1 + c;
                    var y1 = domain[1, 0] + j * d;
                    var y2 = y1 + d;

                    var points = new Point[4] { new(x1, y1), new(x2, y1), new(x2, y2), new(x1, y2) };
                    rectangles.Add(points);

                    for (int k = 0; k < points.Length; ++k)
                    {
                        if (!mesh.Contains(points[k]))
                            mesh.Add(points[k]);
                    }

                    var rect = plot.AddRectangle(x1, x2, y1, y2);
                    rect.Color = Color.Transparent;
                    rect.HatchStyle = HatchStyle.None;
                    rect.BorderColor = Color.Black;
                    var text = plot.AddText("" + n, (x1 + x2) / 2 - 0.1, (y1 + y2) / 2 + 0.05);
                    text.Color = Color.Black;
                    text.FontSize = 20;

                    ++n;
                }
            }

            var Nx = N / 2 + 1;

            Console.WriteLine("{0}x{1}", Nx, Nx);
            plot.SaveFig("plot1.png");
            Process.Start("explorer.exe", "plot1.png");
        }
    }
}
