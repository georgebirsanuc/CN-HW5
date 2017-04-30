
public class Homework5 {
	static int n = 5, kMax = 1000;
	static double epsilon = 0.0000000001d;
	static double[][] V0, V1, A = { { 1, 2, 0, 0, 0 }, { 0, 1, 2, 0, 0 },
			{ 0, 0, 1, 2, 0 }, { 0, 0, 0, 1, 2 }, { 0, 0, 0, 0, 1 } };
	static double[][] B;
	static double[][] A2, A3, A4, A5, A6, A7;

	public static void main(String[] args) {
		printMatrix(Schultz(A), true);
		printMatrix(LiLi1(A), true);
		printMatrix(LiLi2(A), true);

		System.out.println("\nDeducerea formei generale a matricii A cu elemente de forma"
				+ " a[i][i] = 1, a[i][i+1] = 2, pentru i < n si a[n][n] = 1\n");
		A2 = generateA(2);
		A3 = generateA(3);
		A4 = generateA(4);
		A5 = generateA(5);
		A6 = generateA(6);
		A7 = generateA(7);
		printMatrix(LiLi1(A2), true);
		System.out.println();
		printMatrix(LiLi1(A3), true);
		System.out.println();
		printMatrix(LiLi1(A4), true);
		System.out.println();
		printMatrix(LiLi1(A5), true);
		System.out.println();
		printMatrix(LiLi1(A6), true);
		System.out.println();
		printMatrix(LiLi1(A7), true);
	}

	public static double[][] Schultz(double[][] A) {
		int k = 0;
		double deltaV;
		B = negateMatrix(A);
		V1 = divideMatrixByScalar(transposeMatrix(A), calcNorm1(A) * calcNormInf(A));

		do {
			V0 = V1.clone();
			V1 = multiplyMatrices(V0, addXtoDiagonal(multiplyMatrices(B, V0), 2));
			deltaV = calcNorm1(substractMatrices(V1, V0));
			k++;
		} while (deltaV >= epsilon && k <= kMax && deltaV <= Math.pow(10, 10));

		if (deltaV < epsilon) {
			System.out.println("Schultz - Convergenta! \nIteriatii: " + k);
			System.out.println(
					"Norma: " + calcNorm1(addXtoDiagonal(multiplyMatrices(A, V1), -1)));
		} else
			System.out.println("Schultz - Divergenta!");

		return V1;
	}

	public static double[][] LiLi1(double[][] A) {
		int k = 0;
		double deltaV;
		B = negateMatrix(A);
		V1 = divideMatrixByScalar(transposeMatrix(A), calcNorm1(A) * calcNormInf(A));

		// addXtoDiagonal(multiplyMatrices(B, V0), 3); - ultima
		// addXtoDiagonal(multiplyMatrices(multiplyMatrices(B, V0),
		// addXtoDiagonal(multiplyMatrices(B, V0), 3)), 3);

		do {
			V0 = V1.clone();
			V1 = multiplyMatrices(V0,
					addXtoDiagonal(multiplyMatrices(multiplyMatrices(B, V0),
							addXtoDiagonal(multiplyMatrices(B, V0), 3)), 3));
			deltaV = calcNorm1(substractMatrices(V1, V0));
			k++;
		} while (deltaV >= epsilon && k <= kMax && deltaV <= Math.pow(10, 10));

		if (deltaV < epsilon) {
			System.out.println("LiLi1 - Convergenta! \nIteriatii: " + k);
			System.out.println(
					"Norma: " + calcNorm1(addXtoDiagonal(multiplyMatrices(A, V1), -1)));
		} else
			System.out.println("LiLi1 - Divergenta!");

		return V1;
	}

	public static double[][] LiLi2(double[][] A) {
		int k = 0;
		double deltaV;
		B = negateMatrix(A);
		V1 = divideMatrixByScalar(transposeMatrix(A), calcNorm1(A) * calcNormInf(A));
		V0 = V1.clone();

		do {
			V0 = V1.clone();
			V1 = multiplyMatrices(
					addXtoDiagonal(
							multiplyMatrixByScalar(
									multiplyMatrices(
											addXtoDiagonal(
													multiplyMatrices(negateMatrix(V0), A),
													1),
											multiplyMatrices(
													addXtoDiagonal(
															multiplyMatrices(
																	negateMatrix(V0), A),
															3),
													addXtoDiagonal(
															multiplyMatrices(
																	negateMatrix(V0), A),
															3))),
									1d / 4d),
							1),
					V0);

			deltaV = calcNorm1(substractMatrices(V1, V0));
			k++;
		} while (deltaV >= epsilon && k <= kMax && deltaV <= Math.pow(10, 10));

		if (deltaV < epsilon) {
			System.out.println("LiLi2 - Convergenta! \nIteriatii: " + k);
			System.out.println(
					"Norma: " + calcNorm1(addXtoDiagonal(multiplyMatrices(A, V1), -1)));
		} else
			System.out.println("LiLi2 - Divergenta!");

		return V1;
	}

	public static double calcNorm1(double[][] A) {
		double rez = -Double.MAX_VALUE;
		double sum = 0;
		for (int j = 0; j < A.length; ++j) {
			for (int i = 0; i < A[0].length; ++i)
				sum += Math.abs(A[i][j]);
			if (rez < sum)
				rez = sum;
		}
		return rez;
	}

	public static double calcNormInf(double[][] A) {
		double rez = -Double.MAX_VALUE;
		double sum = 0;
		for (int i = 0; i < A.length; ++i) {
			for (int j = 0; j < A[0].length; ++j)
				sum += Math.abs(A[i][j]);
			if (rez < sum)
				rez = sum;
		}
		return rez;
	}

	public static double[][] transposeMatrix(double[][] a) {
		double[][] r = new double[a.length][a.length];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a.length; j++) {
				r[j][i] = a[i][j];
			}
		}
		return r;
	}

	public static double[][] divideMatrixByScalar(double[][] A, double x) {
		double[][] rez = new double[A.length][A[0].length];
		for (int i = 0; i < A.length; ++i)
			for (int j = 0; j < A[0].length; ++j)
				rez[i][j] = A[i][j] / x;
		return rez;
	}

	public static double[][] multiplyMatrixByScalar(double[][] A, double x) {
		double[][] rez = new double[A.length][A[0].length];
		for (int i = 0; i < A.length; ++i)
			for (int j = 0; j < A[0].length; ++j)
				rez[i][j] = A[i][j] * x;
		return rez;
	}

	public static double[][] substractMatrices(double[][] A, double[][] B) {
		double[][] rez = new double[A.length][A[0].length];
		for (int i = 0; i < A.length; ++i)
			for (int j = 0; j < A[0].length; ++j)
				rez[i][j] = A[i][j] - B[i][j];
		return rez;
	}

	public static double[][] addXtoDiagonal(double[][] matrix, int x) {
		for (int i = 0; i < matrix.length; ++i)
			matrix[i][i] += x;
		return matrix;
	}

	public static double[][] multiplyMatrices(double[][] A, double[][] B) {
		double[][] r = new double[A.length][A.length];
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A.length; j++) {
				double s = 0;
				for (int k = 0; k < A.length; k++)
					s += A[i][k] * B[k][j];
				r[i][j] = s;
			}
		}
		return r;
	}

	public static double[][] negateMatrix(double[][] matrix) {
		double[][] rez = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < matrix.length; ++i)
			for (int j = 0; j < matrix[0].length; ++j)
				if (matrix[i][j] != 0)
					rez[i][j] = -matrix[i][j];
		return rez;
	}

	public static double[][] generateA(int n) {
		double[][] rez = new double[n][n];
		for (int i = 0; i < n - 1; ++i) {
			rez[i][i] = 1;
			rez[i][i + 1] = 2;
		}
		rez[n - 1][n - 1] = 1;
		return rez;
	}

	public static void printMatrix(double[][] matrix, boolean scientificNotation) {
		if (scientificNotation)
			for (int i = 0; i < matrix.length; ++i) {
				for (int j = 0; j < matrix[0].length; ++j) {
					System.out.printf("%.20E  ", matrix[i][j]);
				}
				System.out.println();
			}
		else
			for (int i = 0; i < matrix.length; ++i) {
				for (int j = 0; j < matrix[0].length; ++j) {
					System.out.print(matrix[i][j] + " ");
				}
				System.out.println();
			}
		System.out.println();
	}

}
