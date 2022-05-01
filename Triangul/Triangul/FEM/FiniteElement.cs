using System;
using System.Collections.Generic;
using System.Text;
using System.Windows;

namespace FiniteElementMethod
{
    public class FiniteElement
    {
        Equation equation;
        TriangulationMatrixes data;

        public double[] AreaTriangles { get; private set; }

        public double[][] LengthSegments { get; private set; }

        public FiniteElement(Equation equation, TriangulationMatrixes data)
        {
            this.equation = equation;
            this.data = data;
            CalculateArea();
            CalculateLength();
        }
        private void CalculateArea()
        {
            this.AreaTriangles = new double[data.Triangles.GetLength(0)];
            for (int k = 0; k < data.Triangles.GetLength(0); k++)
            {
                Point i = data.Points[data.Triangles[k, 0]];
                Point j = data.Points[data.Triangles[k, 1]];
                Point m = data.Points[data.Triangles[k, 2]];
                AreaTriangles[k] = (new Matrix(new double[3, 3] { { 1, i.X, i.Y }, { 1, j.X, j.Y }, { 1, m.X, m.Y } })).Determinant3x3() / 2;
            }
        }
        private void CalculateLength()
        {
            this.LengthSegments = new double[data.Segments.Length][];
            for (int border = 0; border < data.Segments.Length; border++)
            {
                LengthSegments[border] = new double[data.Segments[border].GetLength(0)];
                for (int k = 0; k < data.Segments[border].GetLength(0); k++)
                {
                    Point p1 = data.Points[data.Segments[border][k, 0]];
                    Point p2 = data.Points[data.Segments[border][k, 1]];
                    LengthSegments[border][k] = Length(p1, p2);
                }
            }
        }
        public Matrix[] Me()
        {
            int N = data.Triangles.GetLength(0); 

            Matrix[] Me = new Matrix[N]; 

            Matrix mat = new Matrix(new double[3, 3] { { 2, 1, 1 }, { 1, 2, 1 }, { 1, 1, 2 } });
            for (int k = 0; k < N; k++)
            {
                double Se = this.AreaTriangles[k];
                Me[k] = ((equation.d * Se) / 12) * mat;
            }
            
            return Me;
        }

        public Matrix[] Qe()
        {
            int N = data.Triangles.GetLength(0); 

            Matrix[] Qe = new Matrix[N]; 

            Matrix mat = new Matrix(new double[3, 3] { { 2, 1, 1 }, { 1, 2, 1 }, { 1, 1, 2 } });
            for (int k = 0; k < N; k++)
            {
                Point i = data.Points[data.Triangles[k, 0]];
                Point j = data.Points[data.Triangles[k, 1]];
                Point m = data.Points[data.Triangles[k, 2]];

                double Se = this.AreaTriangles[k];
                Matrix Fe = new Matrix(new double[3, 1] { { equation.function.Value(i.X, i.Y) }, { equation.function.Value(j.X, j.Y) }, { equation.function.Value(m.X, m.Y) } });
                Qe[k] = (Se / 12) * mat * Fe;
            }
            return Qe;
        }

        public Matrix[] Ke()
        {
            int N = data.Triangles.GetLength(0); 

            Matrix[] Ke = new Matrix[N]; 

            for (int k = 0; k < N; k++)
            {
                Point i = data.Points[data.Triangles[k, 0]];
                Point j = data.Points[data.Triangles[k, 1]];
                Point m = data.Points[data.Triangles[k, 2]];

                Matrix b = new Matrix(new double[1, 3] { { j.Y - m.Y, m.Y - i.Y, i.Y - j.Y } });

                Matrix c = new Matrix(new double[1, 3] { { m.X - j.X, i.X - m.X, j.X - i.X } });

                double Se = this.AreaTriangles[k]; ;

                Matrix Fe = new Matrix(new double[3, 1] { { equation.function.Value(i.X, i.Y) }, { equation.function.Value(j.X, j.Y) }, { equation.function.Value(m.X, m.Y) } });
                Ke[k] = (1.0/(4.0 * Se)) *(equation.a11*b.Transpose()*b+equation.a22*c.Transpose()*c);
            }
            return Ke;
        }

        public Matrix[] Re()
        {
            int N = data.Segments.GetLength(0);

            List<Matrix> Re = new List<Matrix>();

            Matrix mat = new Matrix(new double[2, 2] { { 1.0 / 3.0, 1.0 / 6.0 }, { 1.0 / 6.0, 1.0 / 3.0 } });

            for (int k = 0; k < N; k++)
            {
                for (int s = 0; s < data.Segments[k].GetLength(0); s++)
                {
                    double sigma = equation.BoundaryConditions[k].sigma;
                    double beta = equation.BoundaryConditions[k].beta;
                    Re.Add((sigma/beta) * mat * LengthSegments[k][s]);
                }
            }
            return Re.ToArray();
        }

        public Matrix[] Pe()
        {
            int N = data.Segments.GetLength(0);

            List<Matrix> Pe = new List<Matrix>();

            Matrix mat = new Matrix(new double[2, 1] { { 1.0 / 2.0 }, { 1.0 / 2.0 } });

            for (int k = 0; k < N; k++)
            {
                for (int s = 0; s < data.Segments[k].GetLength(0); s++)
                {
                    double sigma = equation.BoundaryConditions[k].sigma;
                    double beta = equation.BoundaryConditions[k].beta;
                    double uc = equation.BoundaryConditions[k].uc;
                    Pe.Add((uc*sigma/beta) * mat * LengthSegments[k][s]);
                }
            }
            return Pe.ToArray();
        }

        private static double Length(Point p1, Point p2)
        {
            return Math.Sqrt(Math.Pow(p1.X - p2.X, 2) + Math.Pow(p1.Y - p2.Y, 2));
        }
    }
}
