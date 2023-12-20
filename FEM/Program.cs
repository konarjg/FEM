using Accord.IO;
using FEM;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ScottPlot;
using System.Diagnostics;
using FEM.Symbolics;
using Generate = MathNet.Numerics.Generate;
using Accord;
using System.Drawing;
using Accord.Math;
using ScottPlot.Plottable;
using DelaunatorSharp;
using System.Xml.Schema;
using Point = DelaunatorSharp.Point;
using Meta.Numerics.Functions;
using Accord.Math.Decompositions;
using MathNet.Numerics.Integration;

var solution = FEMSolver2DRectangular.Solve(new string[] { "-1/2", "-1/2", "0", "0", "x^2 + y^2" }, new double[,] { { -6, 6 }, { -6, 6 } }, 30, 25);
var x = Generate.LinearSpaced(50, -6, 6);
var u = new double[50, 50];

var exact = new List<((int, int), double)>();

for (int i = 0; i <= 4; ++i)
{
    for (int j = 0; j <= 4; ++j)
    {
        var Ex = 2 * i + 1;
        var Ey = 2 * j + 1;

        exact.Add(((i, j), (Ex + Ey) / 2));
    }    
}

exact = exact.OrderBy(p => p.Item2).ToList();

for (int i = 0; i < 10; ++i)
    Console.WriteLine("E{0} = {1} Exact: {2}", exact[i].Item1, solution[i].Item1 / solution[0].Item1, exact[i].Item2);

for (int i = 0; i < 50; ++i)
{
    for (int j = 0; j < 50; ++j)
        u[i, j] = solution[1].Item2(x[i], x[j]);
}

var plot = new Plot();
plot.SetAxisLimits(-6, 6, -6, 6);

var map = plot.AddHeatmap(u);
map.OffsetX = -6;
map.OffsetY = -6;
map.CellWidth = x[1] - x[0];
map.CellHeight = x[1] - x[0];

plot.AddColorbar(map);
plot.SaveFig("plot.png");
Process.Start("explorer.exe", "plot.png");