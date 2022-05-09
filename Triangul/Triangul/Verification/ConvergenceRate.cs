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
        List<double> list_normaW2_er;
        List<double> list_normaL2_er;
        List<double> list_e;

        List<DataVerification> verifications;

        public ConvergenceRate(List<System.Windows.Point> polygonPoints, Equation equation)
        {
            min_areas = new List<double>();
            min_areas.Add(0.035);
            // min_areas.Add(0.019);
            min_areas.Add(0.01);
            // min_areas.Add(0.005);
            min_areas.Add(0.00245);
            // min_areas.Add(0.001265);
            min_areas.Add(0.00063);
            
            


            list_normaL2 = new List<double>();
            list_normaW2 = new List<double>();
            list_normaL2_er = new List<double>();
            list_normaW2_er = new List<double>();
                
            list_e = new List<double>();
            verifications = new List<DataVerification>();

            for (int i = 0; i < min_areas.Count; i++)
            {;
                QualityOptions options = new QualityOptions();
                options.MinimumAngle = 25;
                options.MaximumArea = min_areas[i];
                Translator translator = new Translator(polygonPoints, options);
                // translator.FileOutput();
                
                DataAdapter FEMadapter = new DataAdapter(translator.Mesh, translator.PolygonIn, equation);
                // FEMadapter.FileOutput();
                DataVerification verification = new DataVerification(FEMadapter);
                verifications.Add(verification);
                double l2 = verification.L2Norm(false);
                double w2 = verification.W2Norm(false);
                list_normaL2.Add(l2);
                list_normaW2.Add(w2);
                double l2_er = verification.L2Norm(true);
                double w2_er = verification.W2Norm(true);
                list_normaL2_er.Add( l2_er);
                list_normaW2_er.Add( w2_er);
                
                double area = verification.GetMaxTriangleArea();
                double error = verification.GetError();
               
                list_e.Add(verification.GetError());
                


            }
        }
        
        public void NormaL()
        {
            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/Norma.txt"))
            {
                file.WriteLine("Count" + ";\t" + "L2" + ";\t" + "W2" );
                for (int i = 0; i < min_areas.Count; i++)
                {
                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + list_normaL2[i] + ";\t" + list_normaW2[i] );
                }
            }
        }
        
        public void NormaLEr()
        {
            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/NormaEr.txt"))
            {
                file.WriteLine("Count" + ";\t" + "L2_er" + ";\t" + "W2_er" );
                for (int i = 0; i < min_areas.Count; i++)
                {
                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + list_normaL2_er[i] + ";\t" + list_normaW2_er[i] );
                }
            }
        }
        
        public void Exact()
        {
            List<double> p_L2 = new List<double>();
            for (int i = 0; i < list_normaL2_er.Count - 1; i++)
            {
                p_L2.Add(Math.Abs(Math.Log(list_normaL2_er[i]) - Math.Log(list_normaL2_er[i + 1])) / Math.Log(2));
            }
            List<double> p_W2 = new List<double>();
            for (int i = 0; i < list_normaW2_er.Count - 1; i++)
            {
                p_W2.Add(Math.Abs(Math.Log(list_normaW2_er[i]) - Math.Log(list_normaW2_er[i + 1])) / Math.Log(2));
            }


            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/Exact.txt"))
            {
                file.WriteLine("Count" + ";\t" + "p_L2" + ";\t" + "p_W2" );
                for (int i = 0; i < min_areas.Count; i++)
                {
                    string text = "-;\t-; ";
                    if (i != 0)
                    {
                        text =  p_L2[i - 1].ToString() + ";\t" + p_W2[i - 1].ToString() ;
                    }

                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + text);
                }
            }
            
        }
        public void Eitken()
        {
            List<double> p_L2 = new List<double>();
            for (int i = 0; i < list_normaL2.Count - 2; i++)
            {
                //p.Add(Abs(Log(list_normaL2[i]) - Log(list_normaL2[i + 1])) / Log(2));
                p_L2.Add(
                    (Math.Log(list_normaL2[i + 1] - list_normaL2[i]) - Math.Log(list_normaL2[i + 2] - list_normaL2[i + 1])) 
                         / Math.Log(2)
                    );
            }
            
            List<double> p_W2 = new List<double>();
            for (int i = 0; i < list_normaW2.Count - 2; i++)
            {
                //p.Add(Abs(Log(list_normaL2[i]) - Log(list_normaL2[i + 1])) / Log(2));
                p_W2.Add((Math.Log(list_normaW2[i + 1] - list_normaW2[i]) -
                                  Math.Log(list_normaW2[i + 2] - list_normaW2[i + 1])) / Math.Log(4));
            }
            

            using (System.IO.StreamWriter file =
                   new System.IO.StreamWriter(
                       @"/home/bdobosevych/RiderProjects/Triangul/Triangul/file/verification/Eitken.txt"))
            {
                file.WriteLine("Count" + ";\t" + "p_L2" + ";\t" + "p_W2" );
                for (int i = 0; i < min_areas.Count; i++)
                {
                    string text = "-;\t-; ";
                    if (i > 1)
                    {
                        text = p_L2[i - 2].ToString() + ";\t" + p_W2[i - 2].ToString() ;
                    }

                    file.WriteLine(verifications[i].FEMData.Mesh.Triangles.Count + ";\t" + text);
                }
            }
        }

       

     
    }
}