using System;
using System.Collections.Generic;
using System.Text;

namespace FiniteElementMethod
{
    public struct Equation
    {
        public double a11 { get; }
        public double a22 { get; }
        public double d { get; }
        public Function function { get; }
        public BoundaryCondition[] BoundaryConditions { get; }

        public Equation(double a11, double a22, double d, Function function, BoundaryCondition[] BoundaryConditions)
        {
            this.a11 = a11;
            this.a22 = a22;
            this.d = d;
            this.function = function;
            this.BoundaryConditions = BoundaryConditions;
        }

    }
}
