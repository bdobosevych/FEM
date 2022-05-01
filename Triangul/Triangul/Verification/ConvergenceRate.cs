using FiniteElementMethod;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using Triangulation;

namespace Verification
{
    class ConvergenceRate
    {
        List<double> min_areas;
        List<double> list_normaL2;
        List<double> list_normaW2;

        List<DataVerification> verifications;

        public ConvergenceRate(List<System.Windows.Point> polygonPoints, Equation equation)
        {
            min_areas = new List<double>();
            double area = 0.1;
            min_areas.Add(area);
            min_areas.Add(area/2);
            min_areas.Add(area/4);
            min_areas.Add(area/6);
            min_areas.Add(area/8);
            min_areas.Add(area/10);
            min_areas.Add(area/12);


            list_normaL2 = new List<double>();
            list_normaW2 = new List<double>();
            verifications = new List<DataVerification>();

            for (int i = 0; i < min_areas.Count; i++)
            {
                QualityOptions options = new QualityOptions();
                options.MinimumAngle = 25;
                options.MaximumArea = min_areas[i];
                Translator translator = new Translator(polygonPoints, options);
                // translator.FileOutput();

                DataAdapter FEMadapter = new DataAdapter(translator.Mesh, translator.PolygonIn, equation);
                // FEMadapter.FileOutput();
                DataVerification verification = new DataVerification(FEMadapter);
                verifications.Add(verification);

                list_normaW2.Add(verification.W2Norm());
                list_normaL2.Add(verification.L2Norm());
            }
        }

        public List<double> get_list_normaL2()
        {
            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/NormaL2.txt"))
            {
                for (int i = 0; i < min_areas.Count; i++)
                {
                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + list_normaL2[i]);
                }
            }

            return list_normaL2;
        }

        public List<double> get_list_normaW2()
        {
            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/NormaW2.txt"))
            {
                for (int i = 0; i < min_areas.Count; i++)
                {
                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + list_normaW2[i]);
                }
            }

            return list_normaW2;
        }

        public List<double> EitkenL2()
        {
            List<double> p = new List<double>();
            for (int i = 0; i < list_normaL2.Count - 2; i++)
            {
                //p.Add(Abs(Log(list_normaL2[i]) - Log(list_normaL2[i + 1])) / Log(2));
                p.Add(Math.Abs(Math.Log(Math.Abs(list_normaL2[i + 1] - list_normaL2[i])) -
                               Math.Log(Math.Abs(list_normaL2[i + 2] - list_normaL2[i + 1]))) / Math.Log(2));
            }

            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/EitkenL2.txt"))
            {
                for (int i = 0; i < min_areas.Count; i++)
                {
                    string text = " - ";
                    if (i > 1)
                    {
                        text = p[i - 2].ToString();
                    }

                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + text);
                }
            }

            return p;
        }

        public List<double> ExactL2()
        {
            List<double> p = new List<double>();
            for (int i = 0; i < list_normaL2.Count - 1; i++)
            {
                p.Add(Math.Abs(Math.Log(list_normaL2[i]) - Math.Log(list_normaL2[i + 1])) / Math.Log(2));
            }

            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/ExactL2.txt"))
            {
                for (int i = 0; i < min_areas.Count; i++)
                {
                    string text = " - ";
                    if (i != 0)
                    {
                        text = p[i - 1].ToString();
                    }

                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + text);
                }
            }

            return p;
        }

        public List<double> EitkenW2()
        {
            List<double> p = new List<double>();
            for (int i = 0; i < list_normaW2.Count - 2; i++)
            {
                //p.Add(Abs(Log(list_normaL2[i]) - Log(list_normaL2[i + 1])) / Log(2));
                p.Add(Math.Abs(Math.Log(Math.Abs(list_normaW2[i + 1] - list_normaW2[i])) -
                               Math.Log(Math.Abs(list_normaW2[i + 2] - list_normaW2[i + 1]))) / Math.Log(2));
            }

            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/EitkenW2.txt"))
            {
                for (int i = 0; i < min_areas.Count; i++)
                {
                    string text = " - ";
                    if (i > 1)
                    {
                        text = p[i - 2].ToString();
                    }

                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + text);
                }
            }

            return p;
        }

        public List<double> ExactW2()
        {
            List<double> p = new List<double>();
            for (int i = 0; i < list_normaW2.Count - 1; i++)
            {
                p.Add(Math.Abs(Math.Log(list_normaW2[i]) - Math.Log(list_normaW2[i + 1])) / Math.Log(2));
            }

            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/ExactW2.txt"))
            {
                for (int i = 0; i < min_areas.Count; i++)
                {
                    string text = " - ";
                    if (i != 0)
                    {
                        text = p[i - 1].ToString();
                    }

                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + text);
                }
            }

            return p;
        }
    }
}