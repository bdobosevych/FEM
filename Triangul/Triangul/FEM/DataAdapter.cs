using FiniteElementMethod;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using TriangleNet.Topology;
using Point = System.Windows.Point;
using MathNet.Numerics.LinearAlgebra.Double;
using Matrix = FiniteElementMethod.Matrix;
namespace Triangulation
{
    class DataAdapter
    {
        //IMesh mesh { get; }
        List<List<SubSegment>> SubSegments = new List<List<SubSegment>>();
        //List<SubSegment> SubsInArray = new List<SubSegment>();
        TriangulationMatrixes data { get; }
        public FiniteElement FiniteElementData { get; }
        public IMesh Mesh { get; }
        Matrix[] Re { get; }
        Matrix[] Qe { get; }
        Matrix[] Pe { get; }
        Matrix[] Me { get; }
        Matrix[] Ke { get; }
        Matrix A { get; }
        Matrix b { get; }
        public Matrix Un { get; set; }
        public DataAdapter(IMesh mesh, IPolygon polygon, Equation equation)
        {
            Mesh = mesh;

            foreach (Segment segment in polygon.Segments)
            {
                List<SubSegment> subs = new List<SubSegment>();
                foreach(SubSegment subSegment in mesh.Segments)
                {
                    if(subSegment.Label==1)
                    {
                        if((subSegment.GetVertex(2)==segment.GetVertex(0) && subSegment.GetVertex(3) == segment.GetVertex(1))|| (subSegment.GetVertex(3) == segment.GetVertex(0) && subSegment.GetVertex(2) == segment.GetVertex(1)))
                        {
                            //MessageBox.Show(" " + subSegment.P0 + " " + subSegment.P1);

                            subs.Add(subSegment);
                        }
                    }
                }
                SubSegments.Add(subs);
            }
            
            int[][,] MatrixSegments = new int[SubSegments.Count][,];
            for(int i=0; i<SubSegments.Count;i++)
            {
                MatrixSegments[i] = new int[SubSegments[i].Count,2];
                for(int j=0;j< SubSegments[i].Count;j++)
                {
                    MatrixSegments[i][j, 0] = SubSegments[i][j].P0;
                    MatrixSegments[i][j, 1] = SubSegments[i][j].P1;
                    //MessageBox.Show(MatrixSegments[i][j, 0].ToString() + MatrixSegments[i][j, 1].ToString());
                }
            }
            int[,] MatrixTriangles = new int[mesh.Triangles.Count, 3];
            int k = 0;
            foreach(Triangle triangle in mesh.Triangles)
            {
                MatrixTriangles[k, 0] = triangle.GetVertex(0).ID;
                MatrixTriangles[k, 1] = triangle.GetVertex(1).ID;
                MatrixTriangles[k, 2] = triangle.GetVertex(2).ID;
                //MessageBox.Show(MatrixTriangles[k, 0].ToString() + MatrixTriangles[k, 1].ToString() + MatrixTriangles[k, 2].ToString());
                k++;
            }
            
            Point[] MatrixPoints = new Point[mesh.Vertices.Count];
            foreach(Vertex vertex in mesh.Vertices)
            {
                MatrixPoints[vertex.ID] = new Point(vertex.X, vertex.Y);
            }
            data = new TriangulationMatrixes(MatrixPoints, MatrixTriangles, MatrixSegments);

            //oOOLLOLOLOLO
            //BoundaryCondition[] boundaryConditions = new BoundaryCondition[SubSegments.Count];
            //Function f = new Function("1");
            //for (int i = 0; i < boundaryConditions.Length; i++)
            //{
            //    boundaryConditions[i] = new BoundaryCondition(1, 1, 1);
            //}
            //Equation equation = new Equation(1, 1, 1, f, boundaryConditions);
            FiniteElementData = new FiniteElement(equation, data);
            Re = FiniteElementData.Re();
            Qe= FiniteElementData.Qe();
            Pe = FiniteElementData.Pe();
            Me = FiniteElementData.Me();
            Ke = FiniteElementData.Ke();
            A = new Matrix(data.Points.Length) ;
            b = new Matrix(data.Points.Length, 1);
            //СЛАР
            for (int e=0; e<Me.GetLength(0);e++)
            {
                //MessageBox.Show(e.ToString() + "e");
                for(int i=0; i<=2;i++)
                {
                    for(int j=0;j<=2;j++)
                    {
                        A[data.Triangles[e, i], data.Triangles[e, j]] += Ke[e][i, j] + Me[e][i, j];
                        //MessageBox.Show(data.Triangles[e, i].ToString() + data.Triangles[e, j].ToString());
                    }
                    b[data.Triangles[e, i],0] += Qe[e][i,0];
                }
            }
             k = 0;

            for (int e=0; e<data.Segments.GetLength(0);e++)
            {
                for(int sub=0;sub<data.Segments[e].GetLength(0);sub++)
                {
                    for(int i=0;i<=1;i++)
                    {
                        for(int j=0;j<=1;j++)
                        {
                            A[data.Segments[e][sub, i], data.Segments[e][sub, j]] += Re[k][i, j];
                        }
                        b[data.Segments[e][sub, i], 0] += Pe[k][i, 0];
                    }
                    k++;
                }
                //for(int i=0;i<data.Segments[e].GetLength(0);i++)
                //{

                //}
            }
            Un = new Matrix(DenseMatrix.OfArray(A.GetValues()).Solve(DenseVector.OfArray(b.GetColumn(0))).ToArray(), 1);
            //Un = new GausMethod(A, b);
            //Un.SolveMatrix();
        }
        public void FileOutput()
        {
            using (System.IO.StreamWriter file =
           new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/CLAR/Res.txt"))
            {
                for(int border =0; border< data.Segments.GetLength(0);border++)
                {
                    for(int subs =0; subs<data.Segments[border].GetLength(0);subs++)
                    {
                        file.WriteLine("Segment: " + data.Segments[border][subs, 0] + "-" + data.Segments[border][subs, 1]+ " with length of: " + FiniteElementData.LengthSegments[border][subs]);

                        for (int i = 0; i < Re[border+subs].RowCount; i++)
                        {

                            string str = "";

                            file.WriteLine(str);
                            for (int j = 0; j < Re[border + subs].ColumnCount; j++)
                            {
                                str += Re[border + subs][i, j].ToString() + " ";
                            }

                            file.WriteLine(str);

                        }
                        file.WriteLine("___________________");
                    }
                }
            }

            
            //MessageBox.Show(Qe.Length.ToString());
            using (System.IO.StreamWriter file=new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/CLAR/Qes.txt"))
            {
                for(int triangle=0; triangle < Qe.GetLength(0); triangle++)
                {
                    file.WriteLine("Triangle: (" + data.Triangles[triangle,0] +"," + data.Triangles[triangle, 1] + "," + data.Triangles[triangle, 2] + ") of area:"+FiniteElementData.AreaTriangles[triangle]);
                    for (int i=0;i< Qe[triangle].RowCount; i++)
                    {
                        string str = "";
                        file.WriteLine(str);
                        for (int j = 0; j < Qe[triangle].ColumnCount; j++)
                        {
                            str += Qe[triangle][i, j].ToString() + " ";
                        }
                        file.WriteLine(str);
                        
                    }
                    file.WriteLine("___________________");
                }
            }


            
            using (System.IO.StreamWriter file =
         new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/CLAR/Pes.txt"))
            {
                for (int border = 0; border < data.Segments.GetLength(0); border++)
                {
                    for (int subs = 0; subs < data.Segments[border].GetLength(0); subs++)
                    {
                        file.WriteLine("Segment: " + data.Segments[border][subs, 0] + "-" + data.Segments[border][subs, 1] + " with length of: " + FiniteElementData.LengthSegments[border][subs]);

                        for (int i = 0; i < Pe[border + subs].RowCount; i++)
                        {

                            string str = "";

                            file.WriteLine(str);
                            for (int j = 0; j < Pe[border + subs].ColumnCount; j++)
                            {
                                str += Pe[border + subs][i, j].ToString() + " ";
                            }

                            file.WriteLine(str);

                        }
                        file.WriteLine("___________________");
                    }
                }
            }

            
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/CLAR/Mes.txt"))
            {
                for (int triangle = 0; triangle < Qe.GetLength(0); triangle++)
                {
                    file.WriteLine("Triangle: (" + data.Triangles[triangle, 0] + "," + data.Triangles[triangle, 1] + "," + data.Triangles[triangle, 2] + ") of area:" + FiniteElementData.AreaTriangles[triangle]);
                    for (int i = 0; i < Me[triangle].RowCount; i++)
                    {
                        string str = "";
                        file.WriteLine(str);
                        for (int j = 0; j < Me[triangle].ColumnCount; j++)
                        {
                            str += Me[triangle][i, j].ToString() + " ";
                        }
                        file.WriteLine(str);

                    }
                    file.WriteLine("___________________");
                }
            }


            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/CLAR/Kes.txt"))
            {
                for (int triangle = 0; triangle < Qe.GetLength(0); triangle++)
                {
                    file.WriteLine("Triangle: (" + data.Triangles[triangle, 0] + "," + data.Triangles[triangle, 1] + "," + data.Triangles[triangle, 2] + ") of area:" + FiniteElementData.AreaTriangles[triangle]);
                    for (int i = 0; i < Ke[triangle].RowCount; i++)
                    {
                        string str = "";
                        file.WriteLine(str);
                        for (int j = 0; j < Ke[triangle].ColumnCount; j++)
                        {
                            str += Ke[triangle][i, j].ToString() + " ";
                        }
                        file.WriteLine(str);

                    }
                    file.WriteLine("___________________");
                }
            }

            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/CLAR/Au=b.txt"))
            {
                string format = "{0,6:00.00}";
                for (int i=0;i<A.RowCount;i++)
                {
                    string str = "\t";
                    for(int j=0;j<A.ColumnCount;j++)
                    {
                        if (j != A.ColumnCount - 1)
                        {
                            string operand;
                            if (A[i, j + 1] >= 0) operand = "   +";
                            else operand = "   ";
                            
                            str += String.Format(format, A[i, j]) + "*U" + j.ToString() + operand;
                        }
                        else str += String.Format(format, A[i, j]) + "*U" + j.ToString();

                    }
                    str += " = " + b[i,0];
                    file.WriteLine(str);
                }
            }
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/CLAR/Un.txt"))
            {
                string str = "";
                int k = 0;
                foreach(Vertex vertex in Mesh.Vertices)
                {
                    str = vertex.X+" "+vertex.Y+" " + Un[k];
                    k++;
                    file.WriteLine(str);
                }
            }
        }
    }
}
