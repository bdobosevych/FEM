﻿using System;
using System.Collections.Generic;
using FiniteElementMethod;
using Verification;
using WPoint = System.Windows.Point;

namespace Triangul
{
    internal class Program
    {
        public static void Main(string[] args)
        {
            List<WPoint> points = new List<WPoint>();
            points.Add(new WPoint(0,0));
            points.Add(new WPoint(1,0));
            points.Add(new WPoint(1,1));
            points.Add(new WPoint(0,1));
            
            BoundaryCondition[] boundaryConditions = new BoundaryCondition[4];
            boundaryConditions[0] = new BoundaryCondition(Math.Pow(10,5), Math.Pow(10,-5), 0);
            boundaryConditions[1] = new BoundaryCondition(Math.Pow(10,-5), Math.Pow(10,5), 0);
            boundaryConditions[2] = new BoundaryCondition(Math.Pow(10,5), Math.Pow(10,-5), 0);
            boundaryConditions[3] = new BoundaryCondition(Math.Pow(10,-5), Math.Pow(10,5), 5.29);

            double a11 = 1;
            double a22 = 1;
            double d = 0;
            string f = "-7*sin(x)";

            ConvergenceRate rateData = new ConvergenceRate(points,
                new Equation(a11, a22, d, new Function(f), boundaryConditions));
            rateData.get_list_normaL2();
            rateData.get_list_normaW2();
            rateData.EitkenL2();
            rateData.ExactL2();
            rateData.EitkenW2();
            rateData.ExactW2();
        }
    }
}