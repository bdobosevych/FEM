using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TriangleNet.Geometry;

namespace Triangulation
{
    class VertexModel
    {
        public int ID { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public string Type { get; set; }
        public VertexModel(Vertex vertex)
        {
            X = vertex.X;
            Y = vertex.Y;
            Type = vertex.Type.ToString();
            ID = vertex.ID;
        }
    }
}
