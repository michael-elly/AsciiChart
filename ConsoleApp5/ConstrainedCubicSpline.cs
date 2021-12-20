using System;
using System.Collections.Generic;

public class ConstrainedCubicSpline {
	// https://pages.uoregon.edu/dgavin/software/spline.pdf
	// Constrained Cubic Spline Interpolation for Chemical Engineering Applications by CJC Kruger

	internal static void Demo() {
		// check the constraint cubic spline
		double[] x = new double[] { 1, 2, 3, 4 };
		double[] y = new double[] { 2, 4, 6, 8 };
		ConstrainedCubicSpline ccs = new ConstrainedCubicSpline(x, y);
		System.Diagnostics.Debug.WriteLine("Interpolated values:");
		System.Diagnostics.Debug.WriteLine(ccs.GetValue(2.5));
		System.Diagnostics.Debug.WriteLine(ccs.GetValue(1));
		System.Diagnostics.Debug.WriteLine(ccs.GetValue(3.4));
		System.Diagnostics.Debug.WriteLine(ccs.GetValue(4));

		/* Intrepolated values should be:
			5
			2
			6.8
			8
		*/
	}

	public readonly int n;
	public readonly double[] x;
	public readonly double[] y;
	public ConstrainedCubicSpline(double[] xvals, double[] yvals) {
		x = xvals;
		y = yvals;
		n = x.Length;

		// Input data must be sorted in ascending x.​
		Array.Sort(x, y);
	}

	public ConstrainedCubicSpline(Dictionary<double, double> xy) {
		n = xy.Count;
		x = new double[n];
		y = new double[n];
		int i = 0;
		foreach (KeyValuePair<double, double> kvp in xy) {
			x[i] = kvp.Key;
			y[i] = kvp.Value;
			i++;
		}

		// Input data must be sorted in ascending x.​
		Array.Sort(x, y);
	}

	public double GetValue(double xx) {
		// Constrained cubic spline interpolation. 
		// Returns yy of interpolant at position xx using data vectors x,y

		// edges and extrapolation (not supported)
		if (xx >= x[n - 1]) return y[n - 1]; else if (xx <= x[0]) return y[0];

		// derivative at each data point
		double[] d = new double[n];

		// calculate the derivative at each internal node, and then the end points​
		for (int i = 1; i < n - 1; i++) {
			d[i] = 2 / (((x[i + 1] - x[i]) / (y[i + 1] - y[i])) + ((x[i] - x[i - 1]) / (y[i] - y[i - 1])));
		}
		// calculate the derivative at the edge nodes
		d[0] = 3 * (y[1] - y[0]) / (2 * (x[1] - x[0])) - d[1] / 2;         // end point beginning​
		d[n - 1] = 3 * (y[n - 1] - y[n - 2]) / (2 * (x[n - 1] - x[n - 2])) - d[n - 2] / 2;   // end point end​

		// calculate the polynomial coefficients vector for the segment ​
		// where the interpolated data should be evaluated​		
		int s = 0; // the segment at which data point xx resides
		while (s < n - 1) {
			if (xx < x[s + 1]) break;
			s++;
		}

		double[,] m = new double[,] {
			{1, x[s],    Math.Pow(x[s],2),    Math.Pow(x[s],     3)},
			{1, x[s+1],  Math.Pow(x[s+1],2),  Math.Pow(x[s+1],   3)},
			{0, 1,       2*x[s],              3*Math.Pow(x[s],   2)},
			{0, 1,       2*x[s+1],            3*Math.Pow(x[s+1], 2)} };

		double[] v = new double[] { y[s], y[s + 1], d[s], d[s + 1] };
		int ier = GaussJordan.Eval(ref m, ref v);
		if (ier == 0) {
			// the interpolated data point
			return v[0] + v[1] * xx + v[2] * Math.Pow(xx, 2) + v[3] * Math.Pow(xx, 3);
		} else {
			// error evaluating the Gauss Jordan for solving coefficients (singular matric)
			return double.NaN;
		}
	}
}
