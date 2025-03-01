using static NumsMethods5.HilbertMatrixSolver;
int size = 8;

double[,] A = GilbertMatrix(size);
Console.WriteLine("Матрица Гильберта:");
for (int i = 0; i < size; i++)
{
    for (int j = 0; j < size; j++)
    {
        Console.Write(A[i, j] + " ");
    }
    Console.WriteLine();
}




double[] x_star = Enumerable.Range(1, size).Select(x => (double)x).ToArray();
double[] b = MatrixVectorMultiply(A, x_star);
Console.WriteLine("\nВектор b:");
foreach (double val in b)
{
    Console.WriteLine(val + " ");
}
Console.WriteLine();

// Формирование возмущённой правой части b̃, таким образом, чтобы случайные
// величины вектора возмущения δb = b − b̃ были по модулю порядка 10−7 – 10−8

Random rand = new Random();
double[] delta_b = new double[size];
for (int i = 0; i < size; i++)
{
    // Генерация случайного числа в диапазоне 10^-8 до 10^-6
    double rangeMin = Math.Pow(10, -8);
    double rangeMax = Math.Pow(10, -6);
    delta_b[i] = rangeMin + (rangeMax - rangeMin) * rand.NextDouble();
    delta_b[i] *= (rand.Next(2) == 0 ? -1 : 1); // Случайный знак
}


double[] b_tilde = new double[size];
for (int i = 0; i < size; i++)
{
    b_tilde[i] = b[i] - delta_b[i];
}


// Абсолютная и относительная нормы возмущения
double abs_norm_b = VectorNorm(delta_b);
double rel_norm_b = abs_norm_b / VectorNorm(b);
Console.WriteLine("Абсолютная норма возмущения правой части: " + abs_norm_b);
Console.WriteLine("Относительная норма возмущения правой части: " + rel_norm_b);

// Найдем решение возмущённой системы Ax̃ = b̃, невязку для найденного решения
// r = Ax̃ − b̃ и норму невязки ∥r∥
double[] x_tilde = Gauss(A, b_tilde); // Использовать свой Gauss, а не MathNet.Numerics

double[] r = MatrixVectorMultiply(A, x_tilde);
for (int i = 0; i < size; i++)
{
    r[i] -= b_tilde[i];
}


double r_norm = VectorNorm(r);
Console.WriteLine("Невязка для решения СЛАУ r = Ax - b: ");
foreach (double val in r)
{
    Console.WriteLine(val + " ");
}
Console.WriteLine();
Console.WriteLine("Норма невязки ||r||: " + r_norm);

//возмущение решения δx = x − x̃, абсолютную и относительную нормы возмущения решения
double[] delta_x = new double[size];
for (int i = 0; i < size; i++)
{
    delta_x[i] = x_star[i] - x_tilde[i];
}

double abs_norm_x = VectorNorm(delta_x);
double rel_norm_x = abs_norm_x / VectorNorm(x_star);
Console.WriteLine("Абсолютная норма возмущения решения: " + abs_norm_x);
Console.WriteLine("Относительная норма возмущения решения: " + rel_norm_x);

Console.WriteLine("eps_x / eps_b = " + rel_norm_x / rel_norm_b);


double lambda_reg = 0.0000000000001;

// Создание единичной матрицы
double[,] identity = new double[size, size];
for (int i = 0; i < size; i++)
{
    identity[i, i] = 1.0;
}

double[,] AT = MatrixTranspose(A);
double[,] ATA = MatrixMatrixMultiply(AT, A);
double[,] lambdaI = new double[size, size];

for (int i = 0; i < size; i++)
{
    lambdaI[i, i] = lambda_reg;
}

double[,] leftSide = MatrixAdd(ATA, lambdaI);

double[] rightSide = MatrixVectorMultiply(AT, b);

double[] result = Gauss(leftSide, rightSide);

Console.WriteLine("Решение регуляризацией: ");
foreach (double val in result)
{
    Console.WriteLine(val + " ");
}
Console.WriteLine();

double totalDeviationExact = 0;
for (int i = 0; i < size; i++)
{
    totalDeviationExact += Math.Abs(result[i] - x_star[i]);
}
Console.WriteLine("Отклонение от точного решения: " + totalDeviationExact);


