using System;
using System.Collections.Generic;
using System.Linq;
//using PointW = System.Windows.Point;
using Polygon = TriangleNet.Geometry.Polygon;
using Mesh = TriangleNet.Mesh;
using WPoint = System.Windows.Point;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using System.Windows.Media;

//using Microsoft.Win32;

namespace Triangulation
{
    public class Translator
    {
        const double spaceratio = 1.5;
        public Polygon PolygonIn { get; set; }
        public IMesh Mesh { get; set; }
        public List<WPoint> InputPoints = new List<WPoint>();
        public int ScaleTranform { get; set; }
        public int CanvasPosition { get; set; }

        public Translator(List<WPoint> points, QualityOptions quality)
        {
            ScaleTranform = 1;
            InputPoints = points;

            PolygonIn = new Polygon();

            for (int i = 0; i < points.Count; i++)
            {
                PolygonIn.Add(new Vertex(points[i].X, points[i].Y));
                if (i == points.Count - 1)
                {
                    PolygonIn.Add(new Segment(new Vertex(points[i].X, points[i].Y),
                        new Vertex(points[0].X, points[0].Y)));
                }
                else
                {
                    PolygonIn.Add(new Segment(new Vertex(points[i].X, points[i].Y),
                        new Vertex(points[i + 1].X, points[i + 1].Y)));
                }
            }

            ConstraintOptions options = new ConstraintOptions();

            options.ConformingDelaunay = true;
            try
            {
                Mesh = PolygonIn.Triangulate(options, quality);
            }
            catch
            {
                Console.WriteLine("помилка");
            }
        }

        public void FileOutput()
        {
            List<TriangleModel> triangles = new List<TriangleModel>();
            List<SegmentModel> segments = new List<SegmentModel>();
            List<VertexModel> vertices = new List<VertexModel>();

            try
            {
                using (System.IO.StreamWriter file = new System.IO.StreamWriter(
                           @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/treangulation/triangles.txt"))
                {
                    foreach (var t in Mesh.Triangles)
                    {
                        var info = new TriangleModel(t);
                        triangles.Add(info);
                        file.WriteLine(info.ID + " " + info.Point1ID + " " + info.Point2ID + " " + info.Point3ID);
                    }
                }

                using (System.IO.StreamWriter file =
                       new System.IO.StreamWriter(
                           @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/treangulation/segments.txt"))
                {
                    foreach (var segment in Mesh.Segments)
                    {
                        if (segment.Label == 1)
                        {
                            var info = new SegmentModel(segment);
                            segments.Add(info);
                            file.WriteLine(info.Border + " " + info.P0 + " " + info.P1);
                        }
                    }
                }

                using (System.IO.StreamWriter file =
                       new System.IO.StreamWriter(
                           @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/treangulation/vertices.txt"))
                {
                    foreach (var vertex in Mesh.Vertices)
                    {
                        var info = new VertexModel(vertex);
                        vertices.Add(info);
                        file.WriteLine(info.ID + " " + info.X + " " + info.Y + " " + info.Type);
                    }
                }
            }
            catch
            {
                Console.WriteLine("erorr");
            }
        }
    }
}