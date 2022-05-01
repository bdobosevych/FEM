using System;
using System.Collections.Generic;
using System.Text;

namespace FiniteElementMethod
{
    public struct BoundaryCondition
    {
        public double beta { get; set; }
        public double sigma { get; set; }
        public double uc { get; set; }

        public BoundaryCondition(double beta, double sigma, double uc)
        {
            this.beta = beta;
            this.sigma = sigma;
            this.uc = uc;
        }

        public override string ToString()
        {
            return beta + "Nu" + ((sigma > 0) ? "+" : "") + sigma + "(u" + ((uc > 0) ? "-" : "+") + Math.Abs(uc) + ")=0";
        }
    }
}
