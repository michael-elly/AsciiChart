using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp5 {
	class Program {
		static void Main(string[] args) {
			Console.OutputEncoding = System.Text.Encoding.UTF8;

			Dictionary<double, double> d = new Dictionary<double, double>();
			double n; int c; int r; int ticks;
			string plot;

			d = new Dictionary<double, double>();
			(n, c, r, ticks) = (200, 20, 17, 5);
			for (double i = 0; i < n; i++) d[i + 1] = 300;
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, XLabelformat = "0", Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			(n, c, r, ticks) = (600, 20, 16, 5);
			for (double i = 0; i < n; i++) d[i + 1] = 300;
			d[180] = 2969;
			d[39] = 3;
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, XLabelformat = "0", Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			(n, c, r, ticks) = (16, 16, 16, 4);
			for (double i = 0; i < n; i++) d[i + 1] = i + 1;
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, XLabelformat = "0", Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			(n, c, r, ticks) = (16, 16, 16, 16);
			for (double i = 0; i < n; i++) d[i + 1] = i + 1;
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, XLabelformat = "0", Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			(n, c, r, ticks) = (16, 16, 16, 8);
			for (double i = 0; i < n; i++) d[i + 1] = i + 1;
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, XLabelformat = "0", Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			(n, c, r, ticks) = (16, 16, 16, 8);
			for (double i = 0; i < n; i++) d[i + 1] = i + 1;
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, XLabelformat = "0.0", Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			(n, c, r, ticks) = (16, 16, 16, 8);
			for (double i = 0; i < n; i++) d[i + 1] = i + 1;
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, XLabelformat = "0.00", Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			n = 39; c = 30; r = 26; ticks = 4;
			for (double i = 0; i <= n; i++) d[i] = Math.Sin(i / n * Math.PI * 2);
			plot = Plot(d, new Options() { YLabelformat = "0.0000", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, Title = $"{n + 1} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			n = 49; c = 60; r = 30; ticks = 8;
			for (double i = 0; i <= n; i++) d[i] = Math.Sin(i / n * Math.PI * 2);
			plot = Plot(d, new Options() { YLabelformat = "0.0000", YMargin = 5, Rows = r, Columns = c, XTicks = ticks, Title = $"{n + 1} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			n = 149; c = 80; r = 26; ticks = 8;
			for (double i = 0; i <= n; i++) d[i] = Math.Sin(i / n * Math.PI * 2);
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 0, Rows = r, Columns = c, XTicks = ticks, Title = $"{n + 1} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			n = 149; c = 50; r = 26; ticks = 8;
			for (double i = 0; i <= n; i++) d[i] = Math.Sin(i / n * Math.PI * 2);
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 0, Rows = r, Columns = c, XTicks = ticks, Title = $"{n + 1} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			n = 20; c = 50; r = 26; ticks = 10;
			for (double i = 0; i <= n; i++) d[i] = Math.Sin(i / n * Math.PI * 2);
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 0, Rows = r, Columns = c, XTicks = ticks, Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			n = 20; c = 50; r = 26; ticks = 10;
			for (double i = 0; i <= n; i++) d[i] = Math.Sin(i / n * Math.PI * 2);
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 0, Rows = r, Columns = c, Connected=false, XTicks = ticks, Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			d = new Dictionary<double, double>();
			n = 160; c = 50; r = 26; ticks = 10;
			for (double i = 0; i <= n; i++) d[i] = Math.Sin(i / n * Math.PI * 2);
			plot = Plot(d, new Options() { YLabelformat = "0.0", YMargin = 0, Rows = r, Columns = c, XTicks = ticks, Title = $"{n} data points over [{r},{c}] with {ticks} ticks" });
			Console.WriteLine(plot + "\r\n");

			Console.Write("Press any key to exit..."); Console.ReadKey();
		}

		public static string Plot(Dictionary<double, double> d, Options o) {
			if (o.XTicks >= o.Columns) o.XTicks = o.Columns / 2;
			char[,] cells = new char[o.Rows, o.Columns];
			for (int c = 0; c < o.Columns; c++) {
				for (int r = 0; r < o.Rows; r++) {
					cells[r, c] = ' ';
				}
			}

			double min_x = d.Keys.ToArray().Min();
			double max_x = d.Keys.ToArray().Max();
			double min_y = d.Values.ToArray().Min();
			double max_y = d.Values.ToArray().Max();
			if (min_y == max_y) { min_y--; max_y++; }
			double dx = (max_x - min_x) / (o.Columns - 1);
			double dy = (max_y - min_y) / (o.Rows - 1);

			// calculate it						
			SplineInterpolator sp = new SplineInterpolator(d);			
			//ConstrainedCubicSpline sp = new ConstrainedCubicSpline(d);
			(int i, int j) last = (0, 0);
			double y;
			for (int i = 0; i < o.Columns; i++) {
				// calculate by interpolation
				y = sp.GetValue(min_x + i * dx);
				int j = (int)Math.Round((y - min_y) / dy, MidpointRounding.AwayFromZero);

				// optionally fill in the dots in between last dot for the 2nd dot onwards
				if (o.Connected) {
					if (i > 0) {
						if (j > last.j)
							for (int jj = last.j + 1; jj < j; jj++) SetCell(cells, jj, i, o, '+');
						else
							for (int jj = last.j - 1; jj > j && jj > 0; jj--) SetCell(cells, jj, i, o, '-');
					}
					last = (i, j);
				}

				// set this dot
				SetCell(cells, j, i, o);
			}

			// ============================================================
			// Print the chart
			// ============================================================
			// calculate left pad
			int left_pad = o.YLabelformat.Length + o.YMargin;
			for (int r = 0; r < o.Rows; r++) {
				int pl = (max_y - dy * r).ToString(o.YLabelformat).PadLeft(o.YLabelformat.Length + 2).Length;
				if (pl > left_pad) left_pad = pl;
			}

			// Title
			StringBuilder s = new StringBuilder();
			if (!string.IsNullOrWhiteSpace(o.Title)) {
				s.Append("".PadLeft(left_pad + 2));
				if (o.Title.Length < o.Columns / 2) {
					s.Append("".PadLeft((int)(o.Columns + 2 - o.Title.Length) / 2) + o.Title);
				} else {
					s.Append(o.Title);
				}
				s.AppendLine();
			}

			// top x axis
			s.Append("".PadLeft(left_pad + 2));
			for (int c = 0; c <= o.Columns; c++) {
				s.Append($"─");
			}
			s.Append("".PadLeft(3, '─'));
			s.AppendLine();

			// y margine, labels, axis, and chart points			
			for (int r = 0; r < o.Rows; r++) {
				s.Append($"{Math.Round(max_y - dy * r, 5).ToString(o.YLabelformat).PadLeft(left_pad)} | ");
				for (int c = 0; c < o.Columns; c++) {
					s.Append(cells[r, c]);
				}
				s.Append("   |");
				s.AppendLine();
			}

			// x axis and ticks
			int delta_ticks = o.Columns / o.XTicks;
			s.Append("".PadRight(left_pad + 2));
			s.Append("─");
			for (int c = 0; c <= o.Columns; c += delta_ticks) {
				string label = "|";
				bool last_label = o.Columns - delta_ticks >= c;
				if (label.Length < delta_ticks || last_label) {
					s.Append(label);
					int dd = delta_ticks - label.Length;
					if (c + dd < o.Columns) {
						s.Append("".PadRight(dd, '─'));
					} else {
						s.Append("".PadRight(o.Columns - c + 2, '─'));
					}
				}
			}
			s.AppendLine();
			// x labels
			s.Append("".PadRight(left_pad + 3));
			int x_axis_pos = 0;
			for (int c = 0; c <= o.Columns; c += delta_ticks) {
				double x = min_x + dx * c;
				if (x_axis_pos == c) {
					string x_tick_label = (x.ToString(o.XLabelformat) + " ").PadRight(delta_ticks);
					s.Append(x_tick_label);
					x_axis_pos += x_tick_label.Length;
				} else {
					if (delta_ticks - (x_axis_pos - c) > 0) {
						string x_tick_label = "".PadRight(delta_ticks - (x_axis_pos - c));
						s.Append(x_tick_label);
						x_axis_pos += x_tick_label.Length;
					}
				}
			}
			return s.ToString();

		}

		private static void SetCell(char[,] cells, int y, int x, Options o, char c = '•') {
			//cells[o.Rows - y - 1, x] = '•';   // '─';
			if (y < o.Rows && y >= 0) {
				cells[o.Rows - y - 1, x] = c;
			} else if (y < 0) {
				cells[o.Rows - 1, x] = c;
			} else {
				cells[0, x] = c;
			}
		}

		public class Options {
			public int Rows = 40;
			public int Columns = 80;
			public string YLabelformat = "0.0";
			public string XLabelformat = "0.0";
			public int XTicks = 3;
			public int YMargin = 2;
			public string Title = "Chart Title";
			public bool Connected = true;
		}

		// Akima Spline Interpolation
		public class SplineInterpolator {
			private readonly double[] _keys;

			private readonly double[] _values;

			private readonly IDictionary<double, double> _nodes;

			private readonly double[] _h;

			private readonly double[] _a;

			/// <summary>
			/// Class constructor.
			/// </summary>
			/// <param name="nodes">Collection of known points for further interpolation.
			/// Should contain at least two items.</param>
			public SplineInterpolator(IDictionary<double, double> nodes) {
				if (nodes == null) {
					throw new ArgumentNullException("nodes");
				}

				var n = nodes.Count;

				if (n < 2) {
					throw new ArgumentException("At least two point required for interpolation.");
				}

				_nodes = new Dictionary<double, double>(nodes);
				_keys = nodes.Keys.ToArray();
				_values = nodes.Values.ToArray();
				_a = new double[n];
				_h = new double[n];

				for (int i = 1; i < n; i++) {
					_h[i] = _keys[i] - _keys[i - 1];
				}

				if (n > 2) {
					var sub = new double[n - 1];
					var diag = new double[n - 1];
					var sup = new double[n - 1];

					for (int i = 1; i <= n - 2; i++) {
						diag[i] = (_h[i] + _h[i + 1]) / 3;
						sup[i] = _h[i + 1] / 6;
						sub[i] = _h[i] / 6;
						_a[i] = (_values[i + 1] - _values[i]) / _h[i + 1] - (_values[i] - _values[i - 1]) / _h[i];
					}

					SolveTridiag(sub, diag, sup, ref _a, n - 2);
				}
			}

			/// <summary>
			/// Gets interpolated value for specified argument.
			/// </summary>
			/// <param name="key">Argument value for interpolation. Must be within 
			/// the interval bounded by lowest ang highest <see cref="_keys"/> values.</param>
			public double GetValue(double key) {
				if (_nodes.ContainsKey(key)) return _nodes[key];

				int gap = 0;
				var previous = double.MinValue;

				// At the end of this iteration, "gap" will contain the index of the interval
				// between two known values, which contains the unknown z, and "previous" will
				// contain the biggest z value among the known samples, left of the unknown z
				for (int i = 0; i < _keys.Length; i++) {
					if (_keys[i] < key && _keys[i] > previous) {
						previous = _keys[i];
						gap = i + 1;
					}
				}

				var x1 = key - previous;
				var x2 = _h[gap] - x1;

				return ((-_a[gap - 1] / 6 * (x2 + _h[gap]) * x1 + _values[gap - 1]) * x2 +
					(-_a[gap] / 6 * (x1 + _h[gap]) * x2 + _values[gap]) * x1) / _h[gap];
			}


			/// <summary>
			/// Solve linear system with tridiagonal n*n matrix "a"
			/// using Gaussian elimination without pivoting.
			/// </summary>
			private static void SolveTridiag(double[] sub, double[] diag, double[] sup, ref double[] b, int n) {
				int i;

				for (i = 2; i <= n; i++) {
					sub[i] = sub[i] / diag[i - 1];
					diag[i] = diag[i] - sub[i] * sup[i - 1];
					b[i] = b[i] - sub[i] * b[i - 1];
				}

				b[n] = b[n] / diag[n];

				for (i = n - 1; i >= 1; i--) {
					b[i] = (b[i] - sup[i] * b[i + 1]) / diag[i];
				}
			}
		}

		public class LinearInterpolation {
			private readonly double[] x;
			private readonly IDictionary<double, double> xy;
			private readonly double dx;
			private readonly int n;

			public LinearInterpolation(IDictionary<double, double> nodes, double dx) {
				if (nodes == null) {
					throw new ArgumentNullException("nodes");
				}

				n = nodes.Count;

				if (n < 2) {
					throw new ArgumentException("At least two point required for interpolation.");
				}

				xy = new Dictionary<double, double>(nodes);
				x = nodes.Keys.ToArray();
				Array.Sort(x);
				this.dx = dx;
			}

			public double GetValue(double key) {
				if (xy.ContainsKey(key)) return xy[key];

				// find the range of i, in which the x value to be found at, resides. 
				int i0 = 0, i1 = n - 1;
				while (x[i1] - x[i0] > dx && i1 - i0 > 1) {
					int ic = i0 + (i1 - i0) / 2;
					if (key < ic) i0 = ic;
					if (key > ic) i1 = ic;
				}

				if (i1 - i0 == 1) {
					// perform linear interpolaiton to find find y(key)
					return xy[x[i0]] + (key - x[i0]) * (xy[x[i1]] - xy[x[i0]]) / (x[i1] - x[i0]);
				} else {
					// perform regression which takes into account all the xy values in the dx range, and then finds y(key)
					// or return the average y
					double sum = 0;
					for (int i = i0; i <= i1; i++) sum += xy[x[i]];
					return sum / (i1 - i0 + 1);
				}
			}
		}
	}
}
