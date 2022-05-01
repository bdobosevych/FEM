using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TriangleNet.Topology;

namespace Triangulation
{
    class TriangleModel
    {
        public int ID { get; set; }
        public int Point3ID { get; set; }

        public int Point2ID { get; set; }

        public int Point1ID { get; set; }

        public TriangleModel(Triangle triangle)
        {
            ID = triangle.ID;
            Point1ID = triangle.GetVertex(0).ID;
            Point2ID = triangle.GetVertex(1).ID;
            Point3ID = triangle.GetVertex(2).ID;
            
        }
    }
}
