using Accord.IO;
using FEM;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ScottPlot;
using System.Diagnostics;
using FEM.Symbolics;
using Generate = MathNet.Numerics.Generate;

var equation = new string[] { "-1", "0", "-2/x" };
FEMSolver1D.Solve(equation, 0, 40, 600, 1);