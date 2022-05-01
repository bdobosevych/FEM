using System;
using System.Collections.Generic;
using System.Text;
using System.Collections;
using System.Windows;


namespace FiniteElementMethod
{
    public class Matrix
    {

        /// <summary>
        /// Contains the rows of the matrix as elements, which
        /// are ArrayLists as well.
        /// </summary>
        private List<List<double>> Values;

        /// <summary>
        /// Number of rows of the matrix.
        /// </summary>
        public int RowCount
        {
            get { return rowCount; }
        }

        /// <summary>
        /// Number of columns of the matrix.
        /// </summary>
        public int ColumnCount
        {
            get { return columnCount; }
        }

        /// <summary>
        /// Number of rows of the matrix.
        /// </summary>
        private int rowCount;

        /// <summary>
        /// Number of columns of the matrix.
        /// </summary>
        private int columnCount;
        public double[,] GetValues()
        {
            double[,] v = new double[RowCount, ColumnCount];
            for (int i = 0; i < RowCount; i++)
            {
                for (int j = 0; j < ColumnCount; j++)
                {
                    v[i, j] = this[i, j];
                }
            }
            return v;
        }
        public double[] GetRows(int i)
        {
            return Values[i].ToArray();
        }
        public double[] GetColumn(int j)
        {
            double[] r = new double[rowCount];
            for (int i = 0; i < rowCount; i++)
            {
                r[i] = this[i, j];
            }
            return r;
        }
        #region Constructors

        /// <summary>
        /// Inits empty matrix 
        /// </summary>
        public Matrix()
        {
            Values = new List<List<double>>();
            rowCount = 0;
            columnCount = 0;
        }

        /// <summary>
        /// Creates m by n matrix filled with zeros; same as Zeros(m, n).
        /// </summary>
        /// <param name="n">Number of rows</param>
        /// <param name="m">Number of columns</param>
        public Matrix(int n, int m, double x = 0.0)
        {
            rowCount = n;
            columnCount = m;

            Values = new List<List<double>>(n);

            for (int i = 0; i < n; i++)
            {
                Values.Add(new List<double>(m));

                for (int j = 0; j < m; j++)
                {
                    Values[i].Add(x);
                }
            }
        }

        /// <summary>
        /// Inits square matrix
        /// </summary>
        /// <param name="n"></param>
        public Matrix(int n, double x = 0.0)
        {
            rowCount = n;
            columnCount = n;

            Values = new List<List<double>>(n);

            for (int i = 0; i < n; i++)
            {
                Values.Add(new List<double>(n));

                for (int j = 0; j < n; j++)
                {
                    Values[i].Add(x);
                }
            }
        }

        /// <summary>
        /// Creates matrix from 2-d double array.
        /// </summary>
        /// <param name="values"></param>
        public Matrix(double[,] values)
        {
            if (values == null)
            {
                Values = new List<List<double>>();
                columnCount = 0;
                rowCount = 0;
            }
            else
            {

                rowCount = (int)values.GetLongLength(0);
                columnCount = (int)values.GetLongLength(1);

                Values = new List<List<double>>(rowCount);

                for (int i = 0; i < rowCount; i++)
                {
                    Values.Add(new List<double>(columnCount));

                    for (int j = 0; j < columnCount; j++)
                    {
                        Values[i].Add(values[i, j]);
                    }
                }
            }
        }

        /// <summary>
        /// Creates matrix from double array.
        /// Type = 0 creates row vector(Default)
        /// Type = 1 creates column vector
        /// </summary>
        /// <param name="values"></param>
        public Matrix(double[] values, int type = 0)
        {
            if (values == null)
            {
                Values = new List<List<double>>();
                columnCount = 0;
                rowCount = 0;
            }
            else
            {
                if (type == 0)
                {
                    rowCount = 1;
                    columnCount = values.Length;

                    Values = new List<List<double>>(1);


                    Values.Add(new List<double>(columnCount));

                    for (int i = 0; i < columnCount; i++)
                    {
                        Values[0].Add(values[i]);
                    }
                }
                else if (type == 1)
                {
                    rowCount = values.Length;
                    columnCount = 1;

                    Values = new List<List<double>>(rowCount);

                    for (int i = 0; i < rowCount; i++)
                    {
                        Values.Add(new List<double>(1));
                        Values[i].Add(values[i]);
                    }
                }
            }

        }

        #endregion

        #region Static func

        /// <summary>
        /// Returns 1 if i = j, and 0 else.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <returns></returns>
        public static double KroneckerDelta(int i, int j)
        {
            return Math.Min(Math.Abs(i - j), 1);
        }

        /// <summary>
        /// Creates m by n chessboard matrix with interchangнng ones and zeros.
        /// 
        /// </summary>
        /// <param name="n">Number of rows.</param>
        /// <param name="m">Number of columns.</param>
        /// <param name="even">Indicates, if matrix entry (1,1) equals zero.</param>
        /// <returns></returns>
        public static Matrix ChessboardMatrix(int n, int m, bool even)
        {
            Matrix M = new Matrix(n, m);

            if (even)
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < m; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 0);
            else
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < m; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 1);

            return M;
        }

        /// <summary>
        /// Creates m by n chessboard matrix with interchangнng ones and zeros.
        /// 
        /// </summary>        
        /// <param name="m">Number of columns.</param>
        /// <param name="even">Indicates, if matrix entry (1,1) equals zero.</param>
        /// <returns></returns>
        public static Matrix ChessboardMatrix(int m, bool even)
        {
            Matrix M = new Matrix(m);

            if (even)
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < m; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 0);
            else
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < m; j++)
                        M[i, j] = KroneckerDelta((i + j) % 2, 1);

            return M;
        }

        /// <summary>
        /// Creates m by n matrix filled with zeros.
        /// </summary>
        /// <param name="n">Number of rows.</param>
        /// <param name="m">Number of columns.</param>
        /// <returns>m by n matrix filled with zeros.</returns>
        public static Matrix Zeros(int n, int m)
        {
            return new Matrix(n, m);
        }

        /// <summary>
        /// Creates n by n matrix filled with zeros.
        /// </summary>       
        /// <param name="n">Number of rows and columns, resp.</param>
        /// <returns>n by n matrix filled with zeros.</returns>
        public static Matrix Zeros(int n)
        {
            return new Matrix(n);
        }

        /// <summary>
        /// Creates m by n matrix filled with ones.
        /// </summary>
        /// <param name="n">Number of rows.</param>
        /// <param name="m">Number of columns.</param>
        /// <returns>m by n matrix filled with ones.</returns>        
        public static Matrix Ones(int n, int m)
        {
            return new Matrix(n, m, 1);
        }

        /// <summary>
        /// Creates n by n matrix filled with ones.
        /// </summary>        
        /// <param name="n">Number of columns.</param>
        /// <returns>n by n matrix filled with ones.</returns>        
        public static Matrix Ones(int n)
        {
            return new Matrix(n, 1);
        }

        /// <summary>
        /// Creates n by n identity matrix.
        /// </summary>
        /// <param name="n">Number of rows and columns respectively.</param>
        /// <returns>n by n identity matrix.</returns>
        public static Matrix Identity(int n)
        {
            return Diag(Ones(n, 1));
        }

        /// <summary>
        /// Creates teh n by n identity matrix.
        /// </summary>
        /// <param name="n">Number of rows and columns, resp.</param>
        /// <returns></returns>
        public static Matrix Eye(int n)
        {
            return Identity(n);
        }


        /// <summary>
        /// Generates diagonal matrix
        /// </summary>
        /// <param name="diag_vector">column vector containing the diag elements</param>
        /// <returns></returns>
        public static Matrix Diag(Matrix diag_vector)
        {
            int dim = diag_vector.VectorLength();

            if (dim == 0)
                throw new ArgumentException("diag_vector must be 1xN or Nx1");

            Matrix M = new Matrix(dim, dim);

            for (int i = 0; i < dim; i++)
            {
                M[i, i] = diag_vector[i];
            }

            return M;

        }


        /// <summary>
        /// Implements the dot product of two vectors.
        /// </summary>
        /// <param name="v">Row or column vector.</param>
        /// <param name="w">Row or column vector.</param>
        /// <returns>Dot product.</returns>
        public static double Dot(Matrix v, Matrix w)
        {
            int n = v.VectorLength();
            int m = w.VectorLength();

            if (m == 0 || n == 0)
                throw new ArgumentException("Arguments need to be vectors.");
            else if (m != n)
                throw new ArgumentException("Vectors must be of the same length.");

            double buf = 0.0;

            for (int i = 0; i < n; i++)
            {
                buf += v[i] * w[i];
            }

            return buf;
        }

        public void BackwardInsertion(Matrix b)
        {
            if (!this.IsUpperTriangular())
                throw new InvalidOperationException("Cannot perform backward insertion for matrix not being upper triangular.");

            if (/*this.Determinant*/this.DiagProd() == 0)
                throw new InvalidOperationException("Warning: Matrix is nearly singular.");

            int n = rowCount;

            if (b.VectorLength() != n)
                throw new ArgumentException("Parameter must vector of the same height as matrix.");

            for (int j = n - 1; j >= 1; j--)
            {
                b[j] /= this[j, j];

                for (int i = 0; i <= j - 1; i++)
                    b[i] -= b[j] * this[i, j];
            }

            b[0] /= this[0, 0];
        }

        public void ForwardInsertion(Matrix b)
        {
            if (!this.IsLowerTriangular())
                throw new InvalidOperationException("Cannot perform forward insertion for matrix not being lower triangular.");

            if (/*this.Determinant*/this.DiagProd() == 0)
                throw new InvalidOperationException("Warning: Matrix is nearly singular.");

            int n = rowCount;

            if (b.VectorLength() != n)
                throw new ArgumentException("Parameter must vector of the same height as matrix.");

            for (int j = 0; j < n - 1; j++)
            {
                b[j] /= this[j, j];

                for (int i = 1; i < n - j; i++)
                    b[j + i] -= b[j] * this[j + i, j];
            }

            b[n - 1] /= this[n - 1, n - 1];
        }

        public Matrix DiagVector()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot get diagonal of non-square matrix.");

            Matrix v = new Matrix(this.columnCount, 1);

            for (int i = 0; i < this.columnCount; i++)
            {
                v[i] = this[i, i];
            }

            return v;
        }

        public Matrix ExtractUpperTrapeze()
        {
            Matrix buf = new Matrix(rowCount, columnCount);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = i; j < columnCount; j++)
                {
                    buf[i, j] = this[i, j];
                }
            }

            return buf;
        }

        public Matrix ExtractLowerTrapeze()
        {
            Matrix buf = new Matrix(rowCount, columnCount);

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    buf[i, j] = this[i, j];
                }
            }

            return buf;
        }

        public static Matrix Solve(Matrix A, Matrix b)
        {
            Matrix A2 = A.Clone();
            Matrix b2 = b.Clone();


            if (!A2.IsSquare())
                throw new InvalidOperationException("Cannot uniquely solve non-square equation system.");

            int n = A2.RowCount;

            Matrix P = A2.LUSafe();



            // We know: PA = LU => [ Ax = b <=> P'LUx = b <=> L(Ux) = (Pb)] since P is orthogonal
            // set y := Ux, solve Ly = Pb by forward insertion
            // and Ux = y by backward insertion

            b2 = P * b2;
            // this solves Ly = Pb
            (A2.ExtractLowerTrapeze() - Diag(A2.DiagVector()) + Identity(n)).ForwardInsertion(b2);

            // this solves Ux = y
            (A2.ExtractUpperTrapeze()).BackwardInsertion(b2);

            return b2;


        }

        #endregion


        #region Matrix manipulations, extractions and decompositions

        /// <summary>
        /// Swaps columns at specified indices. The latter do not have to be ordered.
        /// When equal, nothing is done.
        /// </summary>
        /// <param name="j1">One-based index of first col.</param>
        /// <param name="j2">One-based index of second col.</param>       
        public void SwapColumns(int j1, int j2)
        {
            if (j1 <= 0 || j1 > columnCount || j2 <= 0 || j2 > columnCount)
                throw new ArgumentException("Indices must be positive and <= number of cols.");

            if (j1 == j2)
                return;

            // ArrayList indices are zero-based
            j1--;
            j2--;
            double buf;

            for (int i = 0; i < rowCount; i++)
            {
                buf = Values[i][j1];
                Values[i][j1] = Values[i][j2];
                Values[i][j2] = buf;
            }
        }

        /// <summary>
        /// Swaps rows at specified indices. The latter do not have to be ordered.
        /// When equal, nothing is done.
        /// </summary>
        /// <param name="i1">One-based index of first row.</param>
        /// <param name="i2">One-based index of second row.</param>        
        public void SwapRows(int i1, int i2)
        {
            if (i1 < 0 || i1 > rowCount - 1 || i2 < 0 || i2 > rowCount - 1)
                throw new ArgumentException("Indices must be positive and <= number of rows.");

            if (i1 == i2)
                return;

            List<double> buf = Values[i1];
            Values[i1] = Values[i2];
            Values[i2] = buf;
        }

        /// <summary>
        /// Deletes row at specifies index.
        /// </summary>
        /// <param name="i">One-based index at which to delete.</param>
        public void DeleteRow(int i)
        {
            if (i < 0 || i > rowCount - 1)
                throw new ArgumentException("Index must be positive and <= number of rows.");

            Values.RemoveAt(i);
            rowCount--;
        }

        /// <summary>
        /// Deletes column at specifies index.
        /// </summary>
        /// <param name="j">One-based index at which to delete.</param>
        public void DeleteColumn(int j)
        {
            if (j < 0 || j > columnCount - 1)
                throw new ArgumentException("Index must be positive and <= number of cols.");

            for (int i = 0; i < rowCount; i++)
            {
                Values[i].RemoveAt(j);
            }

            columnCount--;
        }

        /// <summary>
        /// Retrieves row vector at specfifed index and deletes it from matrix.
        /// </summary>
        /// <param name="i">One-based index at which to extract.</param>
        /// <returns>Row vector.</returns>
        public Matrix ExtractRow(int i)
        {
            Matrix buf = this.Row(i);
            this.DeleteRow(i);

            return buf;
        }

        /// <summary>
        /// Retrieves column vector at specfifed index and deletes it from matrix.
        /// </summary>
        /// <param name="j">One-based index at which to extract.</param>
        /// <returns>Row vector.</returns>
        public Matrix ExtractColumn(int j)
        {
            if (j <= 0 || j > columnCount)
                throw new ArgumentException("Index must be positive and <= number of cols.");

            Matrix buf = this.Column(j);
            this.DeleteColumn(j);

            return buf;
        }

        /// <summary>
        /// Inverts square matrix as long as det != 0.
        /// </summary>
        /// <returns>Inverse of matrix.</returns>
        public Matrix Inverse()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot invert non-square matrix.");

            double det = this.Determinant();

            if (det == 0.0)
                throw new InvalidOperationException("Cannot invert (nearly) singular matrix.");

            int n = this.columnCount;

            if (n == 1) return new Matrix(1, 1, 1 / det);

            double[,] buf = new double[n, n];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    buf[i, j] = Math.Pow(-1, i + j) * this.Minor(j, i).Determinant();
                }
            }

            return (new Matrix(buf) / det);
        }

        /// <summary>
        /// Alternative matrix inversion using Leverrier's formula
        /// </summary>
        /// <returns>Inverse of matrix.</returns>
        public Matrix InverseLeverrier()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot invert non-square matrix.");
            //else if (this.Determinant() == 0)
            //    throw new InvalidOperationException("Cannot invert (nearly) singular matrix.");

            int n = this.rowCount;
            Matrix Id = Identity(n);
            Matrix B = Id;
            double alpha;

            for (int k = 1; k < n; k++)
            {
                Matrix buf = (this * B); // DEBUG                
                double buf2 = buf.Trace(); // DEBUG
                alpha = ((double)1 / k) * buf.Trace();
                B = alpha * Id - buf;
            }

            Matrix buf3 = (this * B); // DEBUG                
            double buf4 = buf3.Trace(); // DEBUG
            alpha = (this * B).Trace() / n;
            if (alpha != 0)
                return B / alpha;
            else
                throw new InvalidOperationException("WARNING: Matrix nearly singular or badly scaled.");
        }

        /// <summary>
        /// Calcs the matrix that results in the clearing of a
        /// specified row and a specified column
        /// </summary>        
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <returns></returns>
        public Matrix Minor(int i, int j)
        {
            Matrix A = this.Clone();

            A.DeleteRow(i);
            A.DeleteColumn(j);

            return A;
        }

        /// <summary>
        /// Provides a shallow copy of this matrix in O(m).
        /// </summary>
        /// <returns></returns>
        public Matrix Clone()
        {
            Matrix A = new Matrix(rowCount, columnCount);
            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    A[i, j] = this[i, j];
                }
            }
            return A;
        }

        /// <summary>
        /// Retrieves column with one-based index j.
        /// </summary>
        /// <param name="j"></param>
        /// <returns>j-th column...</returns>
        public Matrix Column(int j)
        {
            if (j < 0 || j > columnCount - 1)
                throw new ArgumentException("Index exceed matrix dimension.");

            Matrix buf = new Matrix(this.rowCount, 1);

            for (int i = 0; i < this.rowCount; i++)
            {
                buf[i, 0] = this[i, j];
            }

            return buf;
        }

        /// <summary>
        /// Retrieves row with one-based index i.
        /// </summary>
        /// <param name="i"></param>
        /// <returns>i-th row...</returns>
        public Matrix Row(int i)
        {
            if (i < 0 || i > rowCount - 1)
                throw new ArgumentException("Index exceed matrix dimension.");

            Matrix buf = new Matrix(columnCount, 1);

            for (int j = 0; j < this.columnCount; j++)
            {
                buf[j, 0] = this[i, j];
            }

            return buf;
        }

        /// <summary>
        /// Swaps each matrix entry A[i, j] with A[j, i].
        /// </summary>
        /// <returns>A transposed matrix.</returns>
        public Matrix Transpose()
        {
            Matrix M = new Matrix(columnCount, rowCount);

            for (int i = 0; i < columnCount; i++)
            {
                for (int j = 0; j < rowCount; j++)
                {
                    M[i, j] = this[j, i];
                }
            }

            return M;
        }

        /// <summary>
        /// Performs LU-decomposition of this instance and saves L and U
        /// within, where the diagonal elements belong to U
        /// (the ones of L are ones...)
        /// </summary>
        public void LU()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot perform LU-decomposition of non-square matrix.");

            int n = this.columnCount;

            for (int j = 0; j < n; j++)
            {
                if (this[j, j] == 0)
                    throw new DivideByZeroException("Warning: Matrix badly scaled or close to singular. Try LUSafe() instead. Check if det != 0.");

                for (int k = 0; k < j; k++)
                {
                    for (int i = k + 1; i < n; i++)
                    {
                        this[i, j] = this[i, j] - this[i, k] * this[k, j];
                    }
                }

                for (int i = j + 1; i < n; i++)
                {
                    this[i, j] = this[i, j] / this[j, j];
                }
            }
        }

        /// <summary>
        /// Performs safe LU-decomposition of this instance with column pivoting 
        /// and saves L and U
        /// within, where the diagonal elements belong to U
        /// (the ones of L are ones...)
        /// </summary>
        /// <returns>Permutation matrix P with P*this = L*U</returns>
        /// <remarks>This needs additional time O(n^2).</remarks>
        public Matrix LUSafe()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot perform LU-decomposition of non-square matrix.");

            int m = this.columnCount;

            Matrix P = Identity(m); // permutation matrix
            int n;

            for (int j = 0; j < m; j++)
            {

                #region Column pivoting

                if (j < m - 1)
                {

                    n = j;

                    for (int i = j + 1; i < m; i++)
                        if (Math.Abs(this[i, j]) > Math.Abs(this[n, j]))
                            n = i;

                    if (n > j) // <=> j2 != j
                    {
                        P.SwapRows(j, n);
                        this.SwapRows(j, n);
                    }

                    if (this[j, j] == 0)
                        throw new DivideByZeroException("Warning: Matrix close to singular.");
                }

                #endregion              

                for (int k = 0; k < j; k++)
                {
                    for (int i = k + 1; i < m; i++)
                    {
                        this[i, j] = this[i, j] - this[i, k] * this[k, j];
                    }
                }

                for (int i = j + 1; i < m; i++)
                {
                    this[i, j] = this[i, j] / this[j, j];
                }
            }

            return P;
        }
        #endregion


        #region Numbers

        public double Determinant3x3()
        {
            if (rowCount != 3 || columnCount != 3)
            {
                throw new InvalidOperationException("Cannot calc determinant of non-3x3 matrix.");
            }
            else
            {
                return (Values[0][0] * Values[1][1] * Values[2][2] + Values[0][1] * Values[1][2] * Values[2][0] + Values[0][2] * Values[1][0] * Values[2][1]) -
                    (Values[0][2] * Values[1][1] * Values[2][0] + Values[1][2] * Values[2][1] * Values[0][0] + Values[2][2] * Values[0][1] * Values[1][0]);
            }
        }

        /// <summary>
        /// Calcs determinant of square matrix
        /// </summary>
        /// <returns></returns>

        public double Determinant()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot calc determinant of non-square matrix.");
            else if (this.columnCount == 1)
                return this[0, 0];
            else if (this.IsTrapeze()) // is square, therefore triangular
            {
                return this.DiagProd();
            }
            else
            {
                // perform LU-decomposition & return product of diagonal elements of U
                Matrix X = this.Clone();


                Matrix P = X.LUSafe();
                //MessageBox.Show(X.DiagProd().ToString());

                return (double)P.Signum() * X.DiagProd();
            }
        }

        /// <summary>
        /// Computes signum of a permutation matrix, which is 1 for an even
        /// number of swaps and -1 for an odd number of swaps. WARNING: 
        /// if *this is not a permutation matrix (e.i. a permutation of Id),
        /// garbage is returned.
        /// </summary>
        /// <returns></returns>
        public double Signum()
        {
            double buf = 1;

            int n = rowCount;
            double fi, fj;

            for (int i = 0; i < n - 1; i++)
            {
                for (fi = 0; fi < n - 1 && this[i, (int)fi] != 1.0; fi++) ;

                for (int j = i + 1; j < n; j++)
                {
                    for (fj = 0; fj < n && this[j, (int)fj] != 1.0; fj++) ;

                    buf *= (fi - fj) / (i - j);
                }
            }

            return buf;
        }

        /// <summary>
        /// Computes product of main diagonal entries.
        /// </summary>
        /// <returns>Product of diagonal elements</returns>
        public double DiagProd()
        {
            double buf = 1.0;
            int dim = Math.Min(this.rowCount, this.columnCount);

            for (int i = 0; i < dim; i++)
            {
                buf *= this[i, i];
            }

            return buf;
        }

        /// <summary>
        /// Calcs trace of the matrix.
        /// </summary>
        /// <returns>Sum of diagonal elements.</returns>
        public double Trace()
        {
            if (!this.IsSquare())
                throw new InvalidOperationException("Cannot calc trace of non-square matrix.");

            double buf = 0.0;

            for (int i = 0; i < this.rowCount; i++)
            {
                buf += this[i, i];
            }

            return buf;
        }

        #endregion

        #region Checks
        /// <summary>
        /// Checks if matrix consists only of zeros and ones.
        /// </summary>
        /// <returns></returns>
        public bool IsZeroOneMatrix()
        {
            for (int i = 0; i < rowCount; i++)
                for (int j = 0; j < columnCount; j++)
                    if (this[i, j] != 0.0 && this[i, j] != 0.0)
                        return false;

            return true;
        }

        /// <summary>
        /// Checks if matrix is permutation of the identity matrix.
        /// </summary>
        /// <returns>True iff matrix is permutation matrix.</returns>
        public bool IsPermutation()
        {
            return (!this.IsSquare() && this.IsZeroOneMatrix() && this.IsInvolutary());

        }

        /// <summary>
        /// Checks if matrix is diagonal matrix.
        /// </summary>
        /// <returns>True iff matrix is diagonal.</returns>


        /// <summary>
        /// Checks if matrix is n by one or one by n.
        /// </summary>
        /// <returns>Length, if vector; zero else.</returns>
        public int VectorLength()
        {
            if (columnCount > 1 && rowCount > 1)
                return 0;
            else return Math.Max(columnCount, rowCount);
        }

        /// <summary>
        /// Checks if number of rows equals number of columns.
        /// </summary>
        /// <returns>True iff matrix is n by n.</returns>
        public bool IsSquare()
        {
            return (this.columnCount == this.rowCount);
        }

        /// <summary>
        /// Checks if matrix is involutary, e.i. if A*A = id.
        /// </summary>
        /// <returns>True iff matrix is involutary.</returns>
        public bool IsInvolutary()
        {
            return (this * this == Identity(rowCount));
        }

        /// <summary>
        /// Checks if A[i, j] == A[j, i].
        /// </summary>
        /// <returns>True iff matrix is symmetric.</returns>
        public bool IsSymmetric()
        {
            for (int i = 1; i <= this.rowCount; i++)
            {
                for (int j = 1; j <= this.columnCount; j++)
                {
                    if (this[i, j] != this[j, i])
                        return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Checks for orthogonality by testing if A*A' == id.
        /// </summary>
        /// <returns>True iff matrix is orthogonal.</returns>
        public bool IsOrthogonal()
        {
            return (this.IsSquare() && this * this.Transpose() == Identity(this.rowCount));
        }

        /// <summary>
        /// Checks if matrix is lower or upper trapeze.
        /// </summary>
        /// <returns>True iff matrix is trapeze.</returns>
        public bool IsTrapeze()
        {
            return (this.IsUpperTrapeze() || this.IsLowerTrapeze());
        }

        /// <summary>
        /// Checks if matrix is trapeze and square.
        /// </summary>
        /// <returns>True iff matrix is triangular.</returns>
        public bool IsTriangular()
        {
            return (this.IsLowerTriangular() || this.IsUpperTriangular());
        }

        /// <summary>
        /// Checks if matrix is square and upper trapeze.
        /// </summary>
        /// <returns>True iff matrix is upper triangular.</returns>
        public bool IsUpperTriangular()
        {
            return (this.IsSquare() && this.IsUpperTrapeze());
        }

        /// <summary>
        /// Checks if matrix is square and lower trapeze.
        /// </summary>
        /// <returns>True iff matrix is lower triangular.</returns>
        public bool IsLowerTriangular()
        {
            return (this.IsSquare() && this.IsLowerTrapeze());
        }

        /// <summary>
        /// Checks if A[i, j] == 0 for i < j.
        /// </summary>
        /// <returns>True iff matrix is upper trapeze.</returns>
        public bool IsUpperTrapeze()
        {
            for (int j = 0; j < columnCount; j++)
                for (int i = j + 1; i < rowCount; i++)
                    if (this[i, j] != 0) return false;

            return true;
        }

        /// <summary>
        /// Checks if A[i, j] == 0 for i > j.
        /// </summary>
        /// <returns>True iff matrix is lower trapeze.</returns>
        public bool IsLowerTrapeze()
        {
            for (int i = 0; i < rowCount; i++)
                for (int j = i + 1; j < columnCount; j++)
                    if (this[i, j] != 0) return false;

            return true;
        }

        #endregion                    

        #region Overrides & Operators

        public static bool operator ==(Matrix A, Matrix B)
        {

            if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
                return false;

            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    if (A[i, j] != B[i, j]) return false;
                }
            }

            return true;
        }

        public static bool operator !=(Matrix A, Matrix B)
        {
            return !(A == B);
        }

        public static Matrix operator +(Matrix A, Matrix B)
        {

            if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
                throw new ArgumentException("Matrices must be of the same dimension.");

            Matrix Result = new Matrix(A.RowCount, A.ColumnCount);

            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    Result[i, j] = A[i, j] + B[i, j];
                }
            }

            return Result;
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {

            if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
                throw new ArgumentException("Matrices must be of the same dimension.");

            Matrix Result = new Matrix(A.RowCount, A.ColumnCount);
            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    Result[i, j] = A[i, j] - B[i, j];
                }
            }

            return Result;
        }

        public static Matrix operator -(Matrix A)
        {
            Matrix Result = new Matrix(A.RowCount, A.ColumnCount);
            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    Result[i, j] = -A[i, j];
                }
            }

            return Result;
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {

            if (A.ColumnCount != B.RowCount)
                throw new ArgumentException("Inner matrix dimensions must agree.");

            Matrix C = new Matrix(A.RowCount, B.ColumnCount);

            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < B.ColumnCount; j++)
                {
                    C[i, j] = Dot(A.Row(i), B.Column(j));
                }
            }

            return C;

        }

        public static Matrix operator *(Matrix A, double x)
        {

            Matrix B = new Matrix(A.rowCount, A.columnCount);

            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    B[i, j] = A[i, j] * x;
                }
            }

            return B;
        }

        public static Matrix operator *(double x, Matrix A)
        {

            Matrix B = new Matrix(A.RowCount, A.ColumnCount);

            for (int i = 0; i < A.RowCount; i++)
            {
                for (int j = 0; j < A.ColumnCount; j++)
                {
                    B[i, j] = A[i, j] * x;
                }
            }

            return B;
        }

        public override int GetHashCode()
        {
            return -1;
        }

        public static Matrix operator /(Matrix A, double x)
        {
            return (1 / x) * A;
        }

        public static Matrix operator ^(Matrix A, int k)
        {
            if (k < 0)
                if (A.IsSquare())
                    return A.InverseLeverrier() ^ (-k);
                else throw new InvalidOperationException("Cannot take non-square matrix to the power of zero.");
            else if (k == 0)
                if (A.IsSquare())
                    return Matrix.Identity(A.RowCount);
                else throw new InvalidOperationException("Cannot take non-square matrix to the power of zero.");
            else if (k == 1)
                if (A.IsSquare())
                    return A;
                else throw new InvalidOperationException("Cannot take non-square matrix to the power of one.");
            else
            {
                Matrix M = A;
                for (int i = 1; i < k; i++)
                {
                    M *= A;
                }

                return M;
            }
        }


        #endregion

        #region Virtuals
        public virtual double this[int i, int j]
        {
            set
            {
                if (i < 0 || j < 0 || i > rowCount - 1 || j > columnCount - 1)
                    throw new ArgumentOutOfRangeException("Indices must be real positive.");

                Values[i][j] = value;
                //this.Values[i - 1, j - 1] = value; 
            }
            get
            {
                if (i >= 0 && i < rowCount && j >= 0 && j < columnCount)
                {
                    return Values[i][j];
                }
                else
                    throw new ArgumentOutOfRangeException("Indices must not exceed size of matrix.");
            }
        }

        public virtual double this[int i]
        {
            set
            {
                if (rowCount == 1)
                {
                    Values[0][i] = value;
                }
                else if (columnCount == 1)
                {
                    Values[i][0] = value;
                }
                else
                    throw new InvalidOperationException("Cannot access multidimensional matrix via single index.");
            }
            get
            {
                if (this.RowCount == 1) // row vector
                    return Values[0][i];
                else if (this.ColumnCount == 1) // column vector
                    return Values[i][0];
                else // neither
                    throw new InvalidOperationException("General matrix acces requires double indexing.");
            }
        }

        public override string ToString()
        {
            string res = "";
            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
                    res += String.Format("{0,15}", Math.Round(this[i, j], 7));
                }
                res += "\n";
            }
            return res;
        }
        #endregion
    }
}



