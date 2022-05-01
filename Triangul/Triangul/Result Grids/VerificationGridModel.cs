using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Triangulation.Result_Grids
{
    public class VerificationGridModel
    {
        //public double MaximumArea { get; set; }
        public int N { get; set; }
        public double L2Norm { get; set; }
        public double W2Norm { get; set; }
        public double pL2Exact { get; set; }
        public double pL2Eitken { get; set; }
        public double pW2Exact { get; set; }
        public double pW2Eitken { get; set; }
        public double AbsoluteError { get; set; }

        //public VerificationGridModel(double l2, double w2, double area, double error)
        //{
        //    MaximumArea = area; L2Norm = l2; W2Norm = w2; AbsoluteError = error;
        //}
    }
}
