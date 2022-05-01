using System;
using System.Collections.Generic;
using System.Windows;
using System.Text;

namespace FiniteElementMethod
{
    public class TriangulationMatrixes
    {

        public Point[] Points { get; }

        public int[,] Triangles { get; }

        public int[][,] Segments { get; }



        public TriangulationMatrixes(Point[] points, int[,] triangles, int[][,] segments)
        {
            this.Points = points;
            this.Triangles = triangles;
            this.Segments = segments;
        }
    }
}
