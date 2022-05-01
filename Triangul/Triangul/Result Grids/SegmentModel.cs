using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TriangleNet.Geometry;
using TriangleNet.Topology;

namespace Triangulation
{
    class SegmentModel
    {
        public string Border { get; set; }
        public double P0 { get; set; }
        public double P1 { get; set; }
        public SegmentModel(SubSegment segment)
        {
            //Point1ID = segment.GetVertex(2).ID;
            //Point2ID = segment.GetVertex(3).ID;
            Border = segment.GetVertex(2).ID.ToString()+ "-" + segment.GetVertex(3).ID; 
            P0 = segment.P0;
            P1 = segment.P1;

            //Point1 = "(" + segment.GetVertex(0).X + ";" + segment.GetVertex(0).Y + ")";
            //Point2 = "(" + segment.GetVertex(1).X + ";" + segment.GetVertex(1).Y + ")";
        }
    }
}
