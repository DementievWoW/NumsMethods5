using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NumsMethods5
{
    internal static class HilbertMatrixSolver
    {
        public static double[,] GilbertMatrix(int size = 8)
        {
            double[,] matrix = new double[size, size];
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    matrix[i, j] = 1.0 / (i + j + 1);
                }
            }
            return matrix;
        }

        public static double[] Gauss(double[,] A, double[] b)
        {
            int size = A.GetLength(0);
            double[,] A_copy = new double[size, size];
            double[] b_copy = new double[size];

            // Копирование массивов, чтобы не изменять исходные данные
            Array.Copy(b, b_copy, size);
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    A_copy[i, j] = A[i, j];
                }
            }


            for (int i = 0; i < size; i++)
            {
                double pivot = A_copy[i, i];
                for (int j = 0; j < size; j++)
                {
                    A_copy[i, j] /= pivot;
                }
                b_copy[i] /= pivot;

                for (int j = i + 1; j < size; j++)
                {
                    double factor = A_copy[j, i];
                    for (int k = 0; k < size; k++)
                    {
                        A_copy[j, k] -= factor * A_copy[i, k];
                    }
                    b_copy[j] -= factor * b_copy[i];
                }
            }

            double[] x = new double[size];
            for (int i = size - 1; i >= 0; i--)
            {
                double sum = 0;
                for (int j = i + 1; j < size; j++)
                {
                    sum += A_copy[i, j] * x[j];
                }
                x[i] = b_copy[i] - sum;
            }

            return x;
        }

        public static double[] MatrixVectorMultiply(double[,] A, double[] x)
        {
            int rows = A.GetLength(0);
            int cols = A.GetLength(1);
            if (cols != x.Length)
            {
                throw new ArgumentException("Размеры матрицы и вектора не совпадают");
            }

            double[] result = new double[rows];
            for (int i = 0; i < rows; i++)
            {
                double sum = 0;
                for (int j = 0; j < cols; j++)
                {
                    sum += A[i, j] * x[j];
                }
                result[i] = sum;
            }
            return result;
        }

        public static double VectorNorm(double[] vector)
        {
            double sumOfSquares = 0;
            foreach (double val in vector)
            {
                sumOfSquares += val * val;
            }
            return Math.Sqrt(sumOfSquares);
        }

        public static double[,] MatrixTranspose(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            double[,] transposedMatrix = new double[cols, rows];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    transposedMatrix[j, i] = matrix[i, j];
                }
            }
            return transposedMatrix;
        }

        public static double[,] MatrixMatrixMultiply(double[,] matrixA, double[,] matrixB)
        {
            int rowsA = matrixA.GetLength(0);
            int colsA = matrixA.GetLength(1);
            int rowsB = matrixB.GetLength(0);
            int colsB = matrixB.GetLength(1);

            if (colsA != rowsB)
            {
                throw new ArgumentException("Недопустимые размеры матриц для умножения.");
            }

            double[,] result = new double[rowsA, colsB];

            for (int i = 0; i < rowsA; i++)
            {
                for (int j = 0; j < colsB; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < colsA; k++)
                    {
                        sum += matrixA[i, k] * matrixB[k, j];
                    }
                    result[i, j] = sum;
                }
            }

            return result;
        }

        public static double[,] MatrixAdd(double[,] matrixA, double[,] matrixB)
        {
            int rowsA = matrixA.GetLength(0);
            int colsA = matrixA.GetLength(1);
            int rowsB = matrixB.GetLength(0);
            int colsB = matrixB.GetLength(1);

            if (rowsA != rowsB || colsA != colsB)
            {
                throw new ArgumentException("Матрицы должны иметь одинаковые размеры для сложения.");
            }

            double[,] result = new double[rowsA, colsA];
            for (int i = 0; i < rowsA; i++)
            {
                for (int j = 0; j < colsA; j++)
                {
                    result[i, j] = matrixA[i, j] + matrixB[i, j];
                }
            }
            return result;
        }
    }
}
