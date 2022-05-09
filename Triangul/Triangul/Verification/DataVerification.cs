using FiniteElementMethod;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using TriangleNet.Topology;
using Point = System.Windows.Point;

namespace Triangulation
{
    class DataVerification
    {
        double x2;
        int n;
        Point[] pointsOnLine;
        Triangle[] trianglesOnPoints;
        double[] Uh;
        double[] U0;

        double[] U0c;
        double[] Uhc;

        double[] ErrorInCenterPoint;

        public DataAdapter FEMData;
        double a = 0, b = 1;

        public double L2;
        public double W2;
        public DataVerification(DataAdapter data)
        {

            FEMData = data;
            n = FEMData.Mesh.Triangles.Count;
            double step = (b - a) / n;
            trianglesOnPoints = new Triangle[n + 1];
            Uh = new double[n + 1];
            U0 = new double[n + 1];
            U0c = new double[FEMData.Mesh.Triangles.Count];
            Uhc = new double[FEMData.Mesh.Triangles.Count];
            ErrorInCenterPoint = new double[FEMData.Mesh.Triangles.Count];

        }
        public DataVerification(double _x2, int _n, DataAdapter data)
        {
            
            x2 = _x2; n = _n; FEMData = data; double step = (b - a) / n;
            trianglesOnPoints = new Triangle[n + 1];
            Uh = new double[n + 1];
            U0 = new double[n + 1];
            U0c = new double[FEMData.Mesh.Triangles.Count];
            Uhc = new double[FEMData.Mesh.Triangles.Count];
            ErrorInCenterPoint = new double[FEMData.Mesh.Triangles.Count];
            pointsOnLine = new Point[n+1];
            for(int i=0; i<=n;i++)
            {
                pointsOnLine[i] = new Point(step * i, x2);
            
            }
            GetTriangles();
        }
        public double MyFunction(Point point)
        {
            double result;
            
            //  3 4 0 2
             // result = point.X / 3.0 - Math.Pow(point.X, 2) / 3.0;
            
             // 1 1 1 3
            result = (-3 * Math.Pow(Math.E, (1 - point.X)) - 3 * Math.Pow(Math.E, point.X) + 3 + 3 * Math.E) / (1 + Math.E);
            
            // 1 1 0 -7sin(x)
            // result = -10 * Math.Pow(point.X, 4);
            return result;
        }
        public double MyFunctionDerivative(Point point)
        {
            double result;
            //
             // result = (1.0 / 3.0) - point.X * (2.0 / 3.0);
            // 1 1 1 3
             result = -(3 * Math.Pow(Math.E, -point.X) * (Math.Pow(Math.E, 2 * point.X) - Math.E) * Math.Log(Math.E)) / (Math.E + 1); 
            
            // 1 1 0 -7sin(x)
            // result = -10 * 4 * Math.Pow(point.X, 3);
            return result;
        }
        public double[] GetExact()
        {
            for (int i = 0; i <= n; i++)
            {
                U0[i] = MyFunction(pointsOnLine[i]) ;
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/U0.txt"))
            {
                
                for (int i = 0; i <= n; i++)
                {
                    file.WriteLine(U0[i].ToString());
                }
            }
            return U0;
        }
        public double GetError()
        {
            double max=0;

            for (int i = 0; i < ErrorInCenterPoint.Length; i++)
            {
                double error = ErrorInCenterPoint[i];
                if (error >= max) max = error;
            }
            return max;
        }
        private double CalculateArea(Triangle e)
        {
            double result;

            Vertex i = e.GetVertex(0);
            Vertex j = e.GetVertex(1);
            Vertex m = e.GetVertex(2);
            result = (new Matrix(new double[3, 3] { { 1, i.X, i.Y }, { 1, j.X, j.Y }, { 1, m.X, m.Y } })).Determinant3x3() / 2;
        
            return result;
        }
        public double L2Norm(bool flag)
        {
            double result = 0;
            double sum = 0;
            int k = 0;
            foreach (Triangle e in FEMData.Mesh.Triangles)
            {
                Point center = CenterOfTriangle(e);
                double U = MyFunction(center);
                double Uh = UhOnPoint(center, e);

                double add_com;
                if (flag)
                {
                    ErrorInCenterPoint[k] = Math.Abs(U - Uh);
                    add_com=Math.Pow(U-Uh,2.0)* CalculateArea(e);
                }
                else
                {
                    add_com=Math.Pow(Uh,2.0)* CalculateArea(e);
                }

                sum += add_com;


                k++;
            }
            L2 = sum; 
            result = Math.Sqrt(sum);
            
            return result;
        }
        public double W2Norm(bool flag)
        {
            
            L2Norm(flag);
            double result = 0;

            //result = L2Norm();
            double sum = 0;
            foreach(Triangle e in FEMData.Mesh.Triangles)
            {
                Point center = CenterOfTriangle(e);
                double U = MyFunction(center);
                double Uh = UhOnPoint(center, e);
                double add_com;
                if (flag)
                {
                    add_com=Math.Pow(U-Uh,2.0)* CalculateArea(e) + Math.Pow((MyFunctionDerivative(center) -  UhDerivative(e)),2) * CalculateArea(e);
                }
                else
                {
                    add_com=Math.Pow(Uh,2.0)* CalculateArea(e) + Math.Pow(UhDerivative(e),2) * CalculateArea(e);
                }

                sum += add_com;
            }
            result += sum;
            W2 = result;
            
            
            return Math.Sqrt(result);
        }
        public double UhDerivative(Triangle triangle)
        {
            double UhDer;
            double bI = triangle.GetVertex(1).Y - triangle.GetVertex(2).Y;
            double bJ = triangle.GetVertex(2).Y - triangle.GetVertex(0).Y;
            double bM = triangle.GetVertex(0).Y - triangle.GetVertex(1).Y;

            double Ui = FEMData.Un[FEMData.Mesh.Vertices.FirstOrDefault(p => p == triangle.GetVertex(0)).ID];
            double Uj = FEMData.Un[FEMData.Mesh.Vertices.FirstOrDefault(p => p == triangle.GetVertex(1)).ID];
            double Um = FEMData.Un[FEMData.Mesh.Vertices.FirstOrDefault(p => p == triangle.GetVertex(2)).ID];
           
            double sigma = GetArea(triangle) * 2;


            UhDer = (Ui*bI+Uj*bJ+Um*bM)/sigma;

            return UhDer;
        }
        public Point CenterOfTriangle(Triangle triangle)
        {
            Point center = new Point();
            
            center.X = (triangle.GetVertex(0).X + triangle.GetVertex(1).X + triangle.GetVertex(2).X) / 3;
            center.Y = (triangle.GetVertex(0).Y + triangle.GetVertex(1).Y + triangle.GetVertex(2).Y) / 3;
            return center;

        }
        public double[] GetAproximates()
        {
            for(int i=0;i<=n;i++)
            {
                Uh[i] = UhOnPoint(pointsOnLine[i], trianglesOnPoints[i]);
                //MessageBox.Show(Uh[i].ToString());
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/Uh.txt"))
            {

                for (int i = 0; i <= n; i++)
                {
                    file.WriteLine(Uh[i]);
                }
            }
            return Uh;
        }
        double GetArea(Triangle e)
        {
            Vertex i = e.GetVertex(0);
            Vertex j = e.GetVertex(1);
            Vertex m = e.GetVertex(2);
           return (new Matrix(new double[3, 3] { { 1, i.X, i.Y }, { 1, j.X, j.Y }, { 1, m.X, m.Y } })).Determinant3x3() / 2;
        }
        double UhOnPoint(Point point, Triangle triangle)
        {
            double sigma = GetArea(triangle) * 2;

            double aI = triangle.GetVertex(1).X * triangle.GetVertex(2).Y - triangle.GetVertex(2).X * triangle.GetVertex(1).Y;
            double bI = triangle.GetVertex(1).Y - triangle.GetVertex(2).Y;
            double cI = triangle.GetVertex(2).X - triangle.GetVertex(1).X;
            double phiI = (aI + bI * point.X + cI * point.Y)/sigma;

            double aJ = triangle.GetVertex(2).X * triangle.GetVertex(0).Y - triangle.GetVertex(0).X * triangle.GetVertex(2).Y;
            double bJ = triangle.GetVertex(2).Y - triangle.GetVertex(0).Y;
            double cJ = triangle.GetVertex(0).X - triangle.GetVertex(2).X;
            double phiJ = (aJ + bJ * point.X + cJ * point.Y) / sigma;

            double aM = triangle.GetVertex(0).X * triangle.GetVertex(1).Y - triangle.GetVertex(1).X * triangle.GetVertex(0).Y;
            double bM = triangle.GetVertex(0).Y - triangle.GetVertex(1).Y;
            double cM = triangle.GetVertex(1).X - triangle.GetVertex(0).X;
            double phiM = (aM + bM * point.X + cM * point.Y) / sigma;

            double Ui = FEMData.Un[FEMData.Mesh.Vertices.FirstOrDefault(p => p == triangle.GetVertex(0)).ID];
            double Uj = FEMData.Un[FEMData.Mesh.Vertices.FirstOrDefault(p => p == triangle.GetVertex(1)).ID];
            double Um = FEMData.Un[FEMData.Mesh.Vertices.FirstOrDefault(p => p == triangle.GetVertex(2)).ID];
           

            return Ui * phiI + Uj * phiJ + Um * phiM;

        }

        public double GetMaxTriangleArea()
        {
            double max =0;
            for(int i=0;i<FEMData.Mesh.Triangles.Count;i++)
            {
                double area = FEMData.FiniteElementData.AreaTriangles[i];
                if (area > max) max = area;
            }
            return max;
        }
        void GetTriangles()
        {
            for(int i=0;i<=n;i++)
            {
                foreach(Triangle triangle in FEMData.Mesh.Triangles)
                {
                    if (TriangleOnPoint(pointsOnLine[i], triangle)) trianglesOnPoints[i] = triangle;
                }
                //MessageBox.Show(trianglesOnPoints[i].ID.ToString());
            }
        }
        public void GetGraphInfo()
        {
         
            using (StreamWriter file = new StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/GRAPHINFO.txt",true))
            {
                file.WriteLine(FEMData.Mesh.Triangles.Count.ToString() + " " + L2.ToString() + " " + W2.ToString());
                
            }
        }
        bool TriangleOnPoint(Point point, Triangle triangle)
        {
            double a = (triangle.GetVertex(0).X - point.X) * (triangle.GetVertex(1).Y - triangle.GetVertex(0).Y) - (triangle.GetVertex(1).X - triangle.GetVertex(0).X) * (triangle.GetVertex(0).Y - point.Y);
            //int a = (P1.X - PTest.X) * (P2.Y - P1.Y) - (P2.X - P1.X) * (P1.Y - PTest.Y);
            double b = (triangle.GetVertex(1).X - point.X) * (triangle.GetVertex(2).Y - triangle.GetVertex(1).Y) - (triangle.GetVertex(2).X - triangle.GetVertex(1).X) * (triangle.GetVertex(1).Y - point.Y);
            //int b = (P2.X - PTest.X) * (P3.Y - P2.Y) - (P3.X - P2.X) * (P2.Y - PTest.Y);
            double c = (triangle.GetVertex(2).X - point.X) * (triangle.GetVertex(0).Y - triangle.GetVertex(2).Y) - (triangle.GetVertex(0).X - triangle.GetVertex(2).X) * (triangle.GetVertex(2).Y - point.Y);
            //int c = (P3.X - PTest.X) * (P1.Y - P3.Y) - (P1.X - P3.X) * (P3.Y - PTest.Y);

            if ((a >= 0 && b >= 0 && c >= 0) || (a <= 0 && b <= 0 && c <= 0))
            {
                return true;
            }
                
            else
            {
                return false;
            }
                
        }
    }
}
