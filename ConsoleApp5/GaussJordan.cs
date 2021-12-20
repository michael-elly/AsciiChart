using System;
using System.Text;

public static class GaussJordan {
	// Gauss-Jordan Elimination - Solving system of linear equation

	internal static void Demo() {
		/*
		Input Matrix:
				x1  x2  b
				--- --- ---
			1   2   1   4
			2   4   3  -1

		Inverse Matrix becomes
		   1.5   -0.5
		  -2      1

		Solution is:
		  x1 = 6.5
		  x2 = -9

		*/

		double[,] a = new double[,] { { 2, 1 }, { 4, 3 } };
		double[] b = new double[] { 4, -1 };

		int ier = GaussJordan.Eval(ref a, ref b);

		System.Diagnostics.Debug.WriteLine($"Error Flag: {ier}");
		System.Diagnostics.Debug.WriteLine($"Inverse Coeeficients Matrix:\r\n{MatrixToString(a)}");
		System.Diagnostics.Debug.WriteLine($"Solution: \r\n{VectorToString(b)}");

	}

	private static string MatrixToString(double[,] a) {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < a.GetLength(0); i++) {
			for (int j = 0; j < a.GetLength(1); j++) {
				s.Append($"{a[i, j]}  \t");
			}
			s.AppendLine();
		}
		return s.ToString();
	}

	private static string VectorToString(double[] a) {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < a.GetLength(0); i++) {
			s.Append($"{a[i]}  \t");
		}
		return s.ToString();
	}

	public static int Eval(ref double[,] A, ref double[] B) {
		// Solution of a system of linear equations 

		// Solves [ A ] * { X } = { B } using the Gauss-Jordan method

		// Input        
		//  A() = matrix of coefficients   (N rows by N columns) [0..N-1, 0..N-1]
		//  B() = right hand column vector (N rows) [0..N-1]

		// Output

		//  A() = inverse of incoming matrix [ A ] (N rows by N columns)
		//  B() = solution vector of linear system (N rows)
		// 
		// Returns erorr flag	
		//    0 = no error
		//    1 = singular matrix

		// All Arrays are zero based (0..N-1)

		int N = B.Length;

		int L = 0;
		int IR = 0;
		int j = 0;
		int i = 0;
		int K = 0;
		int IC = 0;
		int ll = 0;
		double tmp = 0;
		double PMAX = 0;
		double pivinv = 0;

		int[] IPIVOT = new int[N];
		int[] INDEXR = new int[N];
		int[] INDEXC = new int[N];

		for (j = 0; j < N; j++) {
			IPIVOT[j] = 0;
		}

		for (i = 0; i < N; i++) {
			PMAX = 0.0;

			for (j = 0; j < N; j++) {
				if ((IPIVOT[j] != 1)) {
					for (K = 0; K < N; K++) {
						if ((IPIVOT[K] == 0)) {
							if ((System.Math.Abs(A[j, K]) >= PMAX)) {
								PMAX = System.Math.Abs(A[j, K]);
								IR = j;
								IC = K;
							}
						} else if ((IPIVOT[K] > 1)) {
							return 1;
						}
					}
				}
			}

			IPIVOT[IC] = IPIVOT[IC] + 1;

			if ((IR != IC)) {
				for (L = 0; L < N; L++) {
					tmp = A[IR, L];
					A[IR, L] = A[IC, L];
					A[IC, L] = tmp;
				}

				tmp = B[IR];
				B[IR] = B[IC];
				B[IC] = tmp;
			}

			INDEXR[i] = IR;
			INDEXC[i] = IC;

			if ((A[IC, IC] == 0.0)) {
				return 1;
			}

			pivinv = 1.0 / A[IC, IC];
			A[IC, IC] = 1.0;

			for (L = 0; L < N; L++) {
				A[IC, L] = A[IC, L] * pivinv;
			}

			B[IC] = B[IC] * pivinv;

			for (ll = 0; ll < N; ll++) {
				if ((ll != IC)) {
					tmp = A[ll, IC];
					A[ll, IC] = 0.0;
					for (L = 0; L < N; L++) {
						A[ll, L] = A[ll, L] - A[IC, L] * tmp;
					}
					B[ll] = B[ll] - B[IC] * tmp;
				}
			}
		}

		for (L = N - 1; L >= 0; L += -1) {
			if ((INDEXR[L] != INDEXC[L])) {
				for (K = 0; K < N; K++) {
					tmp = A[K, INDEXR[L]];
					A[K, INDEXR[L]] = A[K, INDEXC[L]];
					A[K, INDEXC[L]] = tmp;
				}
			}
		}

		return 0;

	}
}
