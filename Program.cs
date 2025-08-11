using System;
using static System.Console;
using System.Threading;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Globalization;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.RootFinding;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Random;
using MathNet.Numerics.Interpolation;


namespace DoublePrinter
{
    /// <summary>
    /// A class to allow the conversion of doubles to string representations of
    /// their exact decimal values. The implementation aims for readability over
    /// efficiency.
    /// </summary>
    public class DoubleConverter
    {
        /// <summary>
        /// Converts the given double to a string representation of its
        /// exact decimal value.
        /// </summary>
        /// <param name="d">The double to convert.</param>
        /// <returns>A string representation of the double's exact decimal value.</return>
        public static string ToExactString(double d)
        {
            if (double.IsPositiveInfinity(d))
                return "+Infinity";
            if (double.IsNegativeInfinity(d))
                return "-Infinity";
            if (double.IsNaN(d))
                return "NaN";

            // Translate the double into sign, exponent and mantissa.
            long bits = BitConverter.DoubleToInt64Bits(d);
            // Note that the shift is sign-extended, hence the test against -1 not 1
            bool negative = (bits < 0);
            int exponent = (int)((bits >> 52) & 0x7ffL);
            long mantissa = bits & 0xfffffffffffffL;

            // Subnormal numbers; exponent is effectively one higher,
            // but there's no extra normalisation bit in the mantissa
            if (exponent == 0)
            {
                exponent++;
            }
            // Normal numbers; leave exponent as it is but add extra
            // bit to the front of the mantissa
            else
            {
                mantissa = mantissa | (1L << 52);
            }

            // Bias the exponent. It's actually biased by 1023, but we're
            // treating the mantissa as m.0 rather than 0.m, so we need
            // to subtract another 52 from it.
            exponent -= 1075;

            if (mantissa == 0)
            {
                return "0";
            }

            /* Normalize */
            while ((mantissa & 1) == 0)
            {    /*  i.e., Mantissa is even */
                mantissa >>= 1;
                exponent++;
            }

            /// Construct a new decimal expansion with the mantissa
            ArbitraryDecimal ad = new ArbitraryDecimal(mantissa);

            // If the exponent is less than 0, we need to repeatedly
            // divide by 2 - which is the equivalent of multiplying
            // by 5 and dividing by 10.
            if (exponent < 0)
            {
                for (int i = 0; i < -exponent; i++)
                    ad.MultiplyBy(5);
                ad.Shift(-exponent);
            }
            // Otherwise, we need to repeatedly multiply by 2
            else
            {
                for (int i = 0; i < exponent; i++)
                    ad.MultiplyBy(2);
            }

            // Finally, return the string with an appropriate sign
            if (negative)
                return "-" + ad.ToString();
            else
                return ad.ToString();
        }

        /// <summary>Private class used for manipulating
        class ArbitraryDecimal
        {
            /// <summary>Digits in the decimal expansion, one byte per digit
            byte[] digits;
            /// <summary> 
            /// How many digits are *after* the decimal point
            /// </summary>
            int decimalPoint = 0;

            /// <summary> 
            /// Constructs an arbitrary decimal expansion from the given long.
            /// The long must not be negative.
            /// </summary>
            internal ArbitraryDecimal(long x)
            {
                string tmp = x.ToString(CultureInfo.InvariantCulture);
                digits = new byte[tmp.Length];
                for (int i = 0; i < tmp.Length; i++)
                    digits[i] = (byte)(tmp[i] - '0');
                Normalize();
            }

            /// <summary>
            /// Multiplies the current expansion by the given amount, which should
            /// only be 2 or 5.
            /// </summary>
            internal void MultiplyBy(int amount)
            {
                byte[] result = new byte[digits.Length + 1];
                for (int i = digits.Length - 1; i >= 0; i--)
                {
                    int resultDigit = digits[i] * amount + result[i + 1];
                    result[i] = (byte)(resultDigit / 10);
                    result[i + 1] = (byte)(resultDigit % 10);
                }
                if (result[0] != 0)
                {
                    digits = result;
                }
                else
                {
                    Array.Copy(result, 1, digits, 0, digits.Length);
                }
                Normalize();
            }

            /// <summary>
            /// Shifts the decimal point; a negative value makes
            /// the decimal expansion bigger (as fewer digits come after the
            /// decimal place) and a positive value makes the decimal
            /// expansion smaller.
            /// </summary>
            internal void Shift(int amount)
            {
                decimalPoint += amount;
            }

            /// <summary>
            /// Removes leading/trailing zeroes from the expansion.
            /// </summary>
            internal void Normalize()
            {
                int first;
                for (first = 0; first < digits.Length; first++)
                    if (digits[first] != 0)
                        break;
                int last;
                for (last = digits.Length - 1; last >= 0; last--)
                    if (digits[last] != 0)
                        break;

                if (first == 0 && last == digits.Length - 1)
                    return;

                byte[] tmp = new byte[last - first + 1];
                for (int i = 0; i < tmp.Length; i++)
                    tmp[i] = digits[i + first];

                decimalPoint -= digits.Length - (last + 1);
                digits = tmp;
            }

            /// <summary>
            /// Converts the value to a proper decimal string representation.
            /// </summary>
            public override String ToString()
            {
                char[] digitString = new char[digits.Length];
                for (int i = 0; i < digits.Length; i++)
                    digitString[i] = (char)(digits[i] + '0');

                // Simplest case - nothing after the decimal point,
                // and last real digit is non-zero, eg value=35
                if (decimalPoint == 0)
                {
                    return new string(digitString);
                }

                // Fairly simple case - nothing after the decimal
                // point, but some 0s to add, eg value=350
                if (decimalPoint < 0)
                {
                    return new string(digitString) +
                           new string('0', -decimalPoint);
                }

                // Nothing before the decimal point, eg 0.035
                if (decimalPoint >= digitString.Length)
                {
                    return "0." +
                        new string('0', (decimalPoint - digitString.Length)) +
                        new string(digitString);
                }

                // Most complicated case - part of the string comes
                // before the decimal point, part comes after it,
                // eg 3.5
                return new string(digitString, 0,
                                   digitString.Length - decimalPoint) +
                    "." +
                    new string(digitString,
                                digitString.Length - decimalPoint,
                                decimalPoint);
            }
        }
    }
}

namespace CSVHandler
{
    // voir https://www.codeproject.com/Articles/415732/Reading-and-Writing-CSV-Files-in-Csharp
    /// <summary>
    /// Class to store one CSV row
    /// </summary>
    public class CsvRow : List<string>
    {
        /// <summary>
        /// Class to write data to a CSV file
        /// </summary>
        public string LineText { get; set; }
    }

    /// <summary>
    /// Class to write data to a CSV file
    /// </summary>
    public class CsvFileWriter : StreamWriter
    {
        /// <summary>
        /// Class to write data to a CSV file
        /// </summary>
        public CsvFileWriter(Stream stream)
            : base(stream)
        {
        }
        /// <summary>
        /// Class to write data to a CSV file
        /// </summary>
        public CsvFileWriter(string filename)
            : base(filename)
        {
        }

        /// <summary>
        /// Writes a single row to a CSV file.
        /// </summary>
        /// <param name="row">The row to be written</param>
        public void WriteRow(CsvRow row)
        {
            StringBuilder builder = new StringBuilder();
            bool firstColumn = true;
            foreach (string value in row)
            {
                // Add separator if this isn't the first value
                if (!firstColumn)
                    builder.Append(';');
                // Implement special handling for values that contain comma or quote
                // Enclose in quotes and double up any double quotes
                if (value.IndexOfAny(new char[] { '"', ';' }) != -1)
                    builder.AppendFormat("\"{0}\"", value.Replace("\"", "\"\""));
                else
                    builder.Append(value);
                firstColumn = false;
            }
            row.LineText = builder.ToString();
            WriteLine(row.LineText);
        }
    }


    /// <summary>
    /// Class to read data from a CSV file
    /// </summary>
    public class CsvFileReader : StreamReader
    {
        /// <summary>
        /// Class to read data from a CSV file
        /// </summary>
        public CsvFileReader(Stream stream)
            : base(stream)
        {
        }
        /// <summary>
        /// Class to read data from a CSV file
        /// </summary>
        public CsvFileReader(string filename)
            : base(filename)
        {
        }

        /// <summary>
        /// Reads a row of data from a CSV file
        /// </summary>
        /// <param name="row"></param>
        /// <returns></returns>
        public bool ReadRow(CsvRow row)
        {
            row.LineText = ReadLine();
            if (String.IsNullOrEmpty(row.LineText))
                return false;

            int pos = 0;
            int rows = 0;

            while (pos < row.LineText.Length)
            {
                string value;

                // Special handling for quoted field
                if (row.LineText[pos] == '"')
                {
                    // Skip initial quote
                    pos++;

                    // Parse quoted value
                    int start = pos;
                    while (pos < row.LineText.Length)
                    {
                        // Test for quote character
                        if (row.LineText[pos] == '"')
                        {
                            // Found one
                            pos++;

                            // If two quotes together, keep one
                            // Otherwise, indicates end of value
                            if (pos >= row.LineText.Length || row.LineText[pos] != '"')
                            {
                                pos--;
                                break;
                            }
                        }
                        pos++;
                    }
                    value = row.LineText.Substring(start, pos - start);
                    value = value.Replace("\"\"", "\"");
                }
                else
                {
                    // Parse unquoted value
                    int start = pos;
                    while (pos < row.LineText.Length && row.LineText[pos] != ';')
                        pos++;
                    value = row.LineText.Substring(start, pos - start);
                }

                // Add field to list
                if (rows < row.Count)
                    row[rows] = value;
                else
                    row.Add(value);
                rows++;

                // Eat up to and including next comma
                while (pos < row.LineText.Length && row.LineText[pos] != ';')
                    pos++;
                if (pos < row.LineText.Length)
                    pos++;
            }
            // Delete any unused items
            while (row.Count > rows)
                row.RemoveAt(rows);

            // Return true if any columns read
            return (row.Count > 0);
        }
    }
}


namespace ridounet
{
    class SABRFormula
    {
        private static double LOGNORMAL_EPS = 1e-4;

        private double m_maturity;
        private double m_modified_forward;
        private double m_alpha; //Lognormal
        private double m_beta;
        private double m_nu;
        private double m_rho;
        private DiffusionType m_diffusion_type;
        private double[] m_cache;
        public static double gaussianDensity(double x)
        {
            return Math.Exp(-0.5 * x * x) / Math.Sqrt(2 * Math.PI);
        }

        public static double ConvertToAtmLognormalVolatility(double forward, double maturity, double atm_bachelier_volatility)
        {
            double bachelier_time_value = atm_bachelier_volatility * Math.Sqrt(maturity) * gaussianDensity(0.0);
            return ((bachelier_time_value < forward) || (double.Epsilon <= maturity)) ? 2.0 * InverseGaussianCumulative(0.5 * bachelier_time_value / forward + 0.5) / Math.Sqrt(maturity) : atm_bachelier_volatility / forward;
        }

        public static double InverseGaussianCumulative(double x)
        {
            if (x > 1 || x < 0)
            {
                throw new NumericalCalculationException("Argument p must be included in [0,1] in gaussian_inverse_cumulative.");
            }
            else if (x == 0)
            {
                return double.NegativeInfinity;
            }
            else if (x == 1)
            {
                return double.PositiveInfinity; ;
            }
            else
            {
                double guess = Inverse_gaussian_rational_approx(x);
                /*- 
                 * At this point, the relative error of the approximation has Absolute value less than 1.15e-9.
                 * A few iterations (usually only one) of third order Halley's rational method yields 
                 * full machine precision... 
                 */
                double diff = GaussianCumulative(guess) - x;
                double Abs_error = Math.Abs(diff / x);
                double eps = 4 * double.Epsilon;
                int nb_Max_iter = 20;
                for (int iter = 0; Abs_error > eps && iter < nb_Max_iter; ++iter)
                {
                    double dens = gaussianDensity(guess);
                    double t = diff * dens;
                    guess = guess - t / (dens * dens + guess * t * (0.5));              /* Halley's method */
                    diff = GaussianCumulative(guess) - x;
                    Abs_error = Math.Abs(diff / x);
                }
                return guess;
            }
        }


        public class NumericalCalculationException : Exception
        {
            public int code;
            public NumericalCalculationException(String message, int code = 0) : base(message)
            {
                this.code = code;
            }
        }


        public static double Inverse_gaussian_rational_approx(double p)
        {
            double[] a = new double[] { -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00 };
            double[] b = new double[] { -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01 };
            double[] c = new double[] { -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00 };
            double[] d = new double[] { 7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00 };
            double q = Math.Min(p, 1 - p);
            double t, u;
            if (q > 0.02425)
            {
                /* Rational approximation for central region. */
                u = q - 0.5;
                t = u * u;
                u *= (((((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4]) * t + a[5]) / (((((b[0] * t + b[1]) * t + b[2]) * t + b[3]) * t + b[4]) * t + 1);
            }
            else
            {
                /* Rational approximation for tail region. */
                t = Math.Sqrt(-2 * Math.Log(q));
                u = (((((c[0] * t + c[1]) * t + c[2]) * t + c[3]) * t + c[4]) * t + c[5]) / ((((d[0] * t + d[1]) * t + d[2]) * t + d[3]) * t + 1);
            }
            return (p > 0.5 ? -u : u);
        }




        public static double GaussianCumulative(double z)
        {

            IValue f = new CumulativeNormalDistribution();
            return f.Value(z);
        }


        public interface IValue
        {
            double Value(double v);
        }


        public struct Const
        {
            public const double QL_Epsilon = 2.2204460492503131e-016;

            public const double M_SQRT_2 = 0.7071067811865475244008443621048490392848359376887;
            public const double M_1_SQRTPI = 0.564189583547756286948;

            public const double M_LN2 = 0.693147180559945309417;
            public const double M_PI = 3.141592653589793238462643383280;
            public const double M_PI_2 = 1.57079632679489661923;
        }

        public class NormalDistribution : IValue
        {
            private double average_, sigma_, normalizationFactor_, denominator_, derNormalizationFactor_;

            public NormalDistribution() : this(0.0, 1.0) { }
            public NormalDistribution(double average, double sigma)
            {
                average_ = average;
                sigma_ = sigma;

                if (!(sigma_ > 0.0))
                    throw new ApplicationException("sigma must be greater than 0.0 (" + sigma_ + " not allowed)");

                normalizationFactor_ = Const.M_SQRT_2 * Const.M_1_SQRTPI / sigma_;
                derNormalizationFactor_ = sigma_ * sigma_;
                denominator_ = 2.0 * derNormalizationFactor_;
            }

            // function
            public double Value(double x)
            {
                double deltax = x - average_;
                double exponent = -(deltax * deltax) / denominator_;
                // debian alpha had some strange problem in the very-low range
                return exponent <= -690.0 ? 0.0 :  // exp(x) < 1.0e-300 anyway
                    normalizationFactor_ * Math.Exp(exponent);
            }

            public double derivative(double x)
            {
                return (Value(x) * (average_ - x)) / derNormalizationFactor_;
            }
        }

        public class CumulativeNormalDistribution : IValue
        {
            private double average_, sigma_;
            private NormalDistribution gaussian_ = new NormalDistribution();

            public CumulativeNormalDistribution() : this(0.0, 1.0) { }
            public CumulativeNormalDistribution(double average, double sigma)
            {
                average_ = average;
                sigma_ = sigma;

                if (!(sigma_ > 0.0))
                    throw new ApplicationException("sigma must be greater than 0.0 (" + sigma_ + " not allowed)");
            }

            // function
            public double Value(double z)
            {
                //QL_REQUIRE(!(z >= average_ && 2.0*average_-z > average_),
                //           "not a real number. ");
                z = (z - average_) / sigma_;

                double result = 0.5 * (1.0 + erf(z * Const.M_SQRT_2));
                if (result <= 1e-8)
                { //todo: investigate the threshold level
                  // Asymptotic expansion for very negative z following (26.2.12)
                  // on page 408 in M. Abramowitz and A. Stegun,
                  // Pocketbook of Mathematical Functions, ISBN 3-87144818-4.
                    double sum = 1.0, zsqr = z * z, i = 1.0, g = 1.0, x, y,
                         a = double.MaxValue, lasta;
                    do
                    {
                        lasta = a;
                        x = (4.0 * i - 3.0) / zsqr;
                        y = x * ((4.0 * i - 1) / zsqr);
                        a = g * (x - y);
                        sum -= a;
                        g *= y;
                        ++i;
                        a = Math.Abs(a);
                    } while (lasta > a && a >= Math.Abs(sum * Const.QL_Epsilon));
                    result = -gaussian_.Value(z) / z * sum;
                }
                return result;
            }

            public double Derivative(double x)
            {
                double xn = (x - average_) / sigma_;
                return gaussian_.Value(xn) / sigma_;
            }

            #region Sun Microsystems method
            /*
            * ====================================================
            * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
            *
            * Developed at SunPro, a Sun Microsystems, Inc. business.
            * Permission to use, copy, modify, and distribute this
            * software is freely granted, provided that this notice 
            * is preserved.
            * ====================================================
            */

            /* double erf(double x)
            * double erfc(double x)
            *                           x
            *                    2      |\
            *     erf(x)  =  ---------  | exp(-t*t)dt
            *                 sqrt(pi) \| 
            *                           0
            *
            *     erfc(x) =  1-erf(x)
            *  Note that 
            *              erf(-x) = -erf(x)
            *              erfc(-x) = 2 - erfc(x)
            *
            * Method:
            *      1. For |x| in [0, 0.84375]
            *          erf(x)  = x + x*R(x^2)
            *          erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
            *                  = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
            *         where R = P/Q where P is an odd poly of degree 8 and
            *         Q is an odd poly of degree 10.
            *                                               -57.90
            *                      | R - (erf(x)-x)/x | <= 2
            *      
            *
            *         Remark. The formula is derived by noting
            *          erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
            *         and that
            *          2/sqrt(pi) = 1.128379167095512573896158903121545171688
            *         is close to one. The interval is chosen because the fix
            *         point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is
            *         near 0.6174), and by some experiment, 0.84375 is chosen to
            *         guarantee the error is less than one ulp for erf.
            *
            *      2. For |x| in [0.84375,1.25], let s = |x| - 1, and
            *         c = 0.84506291151 rounded to single (24 bits)
            *              erf(x)  = sign(x) * (c  + P1(s)/Q1(s))
            *              erfc(x) = (1-c)  - P1(s)/Q1(s) if x > 0
            *                        1+(c+P1(s)/Q1(s))    if x < 0
            *              |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
            *         Remark: here we use the taylor series expansion at x=1.
            *              erf(1+s) = erf(1) + s*Poly(s)
            *                       = 0.845.. + P1(s)/Q1(s)
            *         That is, we use rational approximation to approximate
            *                      erf(1+s) - (c = (single)0.84506291151)
            *         Note that |P1/Q1|< 0.078 for x in [0.84375,1.25]
            *         where 
            *              P1(s) = degree 6 poly in s
            *              Q1(s) = degree 6 poly in s
            *
            *      3. For x in [1.25,1/0.35(~2.857143)], 
            *              erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
            *              erf(x)  = 1 - erfc(x)
            *         where 
            *              R1(z) = degree 7 poly in z, (z=1/x^2)
            *              S1(z) = degree 8 poly in z
            *
            *      4. For x in [1/0.35,28]
            *              erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2) if x > 0
            *                      = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6<x<0
            *                      = 2.0 - tiny            (if x <= -6)
            *              erf(x)  = sign(x)*(1.0 - erfc(x)) if x < 6, else
            *              erf(x)  = sign(x)*(1.0 - tiny)
            *         where
            *              R2(z) = degree 6 poly in z, (z=1/x^2)
            *              S2(z) = degree 7 poly in z
            *
            *      Note1:
            *         To compute exp(-x*x-0.5625+R/S), let s be a single
            *         precision number and s := x; then
            *              -x*x = -s*s + (s-x)*(s+x)
            *              exp(-x*x-0.5626+R/S) = 
            *                      exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S);
            *      Note2:
            *         Here 4 and 5 make use of the asymptotic series
            *                        exp(-x*x)
            *              erfc(x) ~ ---------- * ( 1 + Poly(1/x^2) )
            *                        x*sqrt(pi)
            *         We use rational approximation to approximate
            *              g(s)=f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
            *         Here is the error bound for R1/S1 and R2/S2
            *              |R1/S1 - f(x)|  < 2**(-62.57)
            *              |R2/S2 - f(x)|  < 2**(-61.52)
            *
            *      5. For inf > x >= 28
            *              erf(x)  = sign(x) *(1 - tiny)  (raise inexact)
            *              erfc(x) = tiny*tiny (raise underflow) if x > 0
            *                      = 2 - tiny if x<0
            *
            *      7. Special case:
            *              erf(0)  = 0, erf(inf)  = 1, erf(-inf) = -1,
            *              erfc(0) = 1, erfc(inf) = 0, erfc(-inf) = 2, 
            *              erfc/erf(NaN) is NaN
            */

            const double tiny = QL_Epsilon,
            one = 1.00000000000000000000e+00, /* 0x3FF00000, 0x00000000 */
                                              /* c = (float)0.84506291151 */
            erx = 8.45062911510467529297e-01, /* 0x3FEB0AC1, 0x60000000 */
                                              //
                                              // Coefficients for approximation to  erf on [0,0.84375]
                                              //
            efx = 1.28379167095512586316e-01, /* 0x3FC06EBA, 0x8214DB69 */
            efx8 = 1.02703333676410069053e+00, /* 0x3FF06EBA, 0x8214DB69 */
            pp0 = 1.28379167095512558561e-01, /* 0x3FC06EBA, 0x8214DB68 */
            pp1 = -3.25042107247001499370e-01, /* 0xBFD4CD7D, 0x691CB913 */
            pp2 = -2.84817495755985104766e-02, /* 0xBF9D2A51, 0xDBD7194F */
            pp3 = -5.77027029648944159157e-03, /* 0xBF77A291, 0x236668E4 */
            pp4 = -2.37630166566501626084e-05, /* 0xBEF8EAD6, 0x120016AC */
            qq1 = 3.97917223959155352819e-01, /* 0x3FD97779, 0xCDDADC09 */
            qq2 = 6.50222499887672944485e-02, /* 0x3FB0A54C, 0x5536CEBA */
            qq3 = 5.08130628187576562776e-03, /* 0x3F74D022, 0xC4D36B0F */
            qq4 = 1.32494738004321644526e-04, /* 0x3F215DC9, 0x221C1A10 */
            qq5 = -3.96022827877536812320e-06, /* 0xBED09C43, 0x42A26120 */
                                               //
                                               // Coefficients for approximation to  erf  in [0.84375,1.25]
                                               //
            pa0 = -2.36211856075265944077e-03, /* 0xBF6359B8, 0xBEF77538 */
            pa1 = 4.14856118683748331666e-01, /* 0x3FDA8D00, 0xAD92B34D */
            pa2 = -3.72207876035701323847e-01, /* 0xBFD7D240, 0xFBB8C3F1 */
            pa3 = 3.18346619901161753674e-01, /* 0x3FD45FCA, 0x805120E4 */
            pa4 = -1.10894694282396677476e-01, /* 0xBFBC6398, 0x3D3E28EC */
            pa5 = 3.54783043256182359371e-02, /* 0x3FA22A36, 0x599795EB */
            pa6 = -2.16637559486879084300e-03, /* 0xBF61BF38, 0x0A96073F */
            qa1 = 1.06420880400844228286e-01, /* 0x3FBB3E66, 0x18EEE323 */
            qa2 = 5.40397917702171048937e-01, /* 0x3FE14AF0, 0x92EB6F33 */
            qa3 = 7.18286544141962662868e-02, /* 0x3FB2635C, 0xD99FE9A7 */
            qa4 = 1.26171219808761642112e-01, /* 0x3FC02660, 0xE763351F */
            qa5 = 1.36370839120290507362e-02, /* 0x3F8BEDC2, 0x6B51DD1C */
            qa6 = 1.19844998467991074170e-02, /* 0x3F888B54, 0x5735151D */
                                              //
                                              // Coefficients for approximation to  erfc in [1.25,1/0.35]
                                              //
            ra0 = -9.86494403484714822705e-03, /* 0xBF843412, 0x600D6435 */
            ra1 = -6.93858572707181764372e-01, /* 0xBFE63416, 0xE4BA7360 */
            ra2 = -1.05586262253232909814e+01, /* 0xC0251E04, 0x41B0E726 */
            ra3 = -6.23753324503260060396e+01, /* 0xC04F300A, 0xE4CBA38D */
            ra4 = -1.62396669462573470355e+02, /* 0xC0644CB1, 0x84282266 */
            ra5 = -1.84605092906711035994e+02, /* 0xC067135C, 0xEBCCABB2 */
            ra6 = -8.12874355063065934246e+01, /* 0xC0545265, 0x57E4D2F2 */
            ra7 = -9.81432934416914548592e+00, /* 0xC023A0EF, 0xC69AC25C */
            sa1 = 1.96512716674392571292e+01, /* 0x4033A6B9, 0xBD707687 */
            sa2 = 1.37657754143519042600e+02, /* 0x4061350C, 0x526AE721 */
            sa3 = 4.34565877475229228821e+02, /* 0x407B290D, 0xD58A1A71 */
            sa4 = 6.45387271733267880336e+02, /* 0x40842B19, 0x21EC2868 */
            sa5 = 4.29008140027567833386e+02, /* 0x407AD021, 0x57700314 */
            sa6 = 1.08635005541779435134e+02, /* 0x405B28A3, 0xEE48AE2C */
            sa7 = 6.57024977031928170135e+00, /* 0x401A47EF, 0x8E484A93 */
            sa8 = -6.04244152148580987438e-02, /* 0xBFAEEFF2, 0xEE749A62 */
                                               //
                                               // Coefficients for approximation to  erfc in [1/.35,28]
                                               //
            rb0 = -9.86494292470009928597e-03, /* 0xBF843412, 0x39E86F4A */
            rb1 = -7.99283237680523006574e-01, /* 0xBFE993BA, 0x70C285DE */
            rb2 = -1.77579549177547519889e+01, /* 0xC031C209, 0x555F995A */
            rb3 = -1.60636384855821916062e+02, /* 0xC064145D, 0x43C5ED98 */
            rb4 = -6.37566443368389627722e+02, /* 0xC083EC88, 0x1375F228 */
            rb5 = -1.02509513161107724954e+03, /* 0xC0900461, 0x6A2E5992 */
            rb6 = -4.83519191608651397019e+02, /* 0xC07E384E, 0x9BDC383F */
            sb1 = 3.03380607434824582924e+01, /* 0x403E568B, 0x261D5190 */
            sb2 = 3.25792512996573918826e+02, /* 0x40745CAE, 0x221B9F0A */
            sb3 = 1.53672958608443695994e+03, /* 0x409802EB, 0x189D5118 */
            sb4 = 3.19985821950859553908e+03, /* 0x40A8FFB7, 0x688C246A */
            sb5 = 2.55305040643316442583e+03, /* 0x40A3F219, 0xCEDF3BE6 */
            sb6 = 4.74528541206955367215e+02, /* 0x407DA874, 0xE79FE763 */
            sb7 = -2.24409524465858183362e+01; /* 0xC03670E2, 0x42712D62 */

            private double erf(double x)
            {
                double R, S, P, Q, s, y, z, r, ax;

                ax = Math.Abs(x);

                if (ax < 0.84375)
                {      /* |x|<0.84375 */
                    if (ax < 3.7252902984e-09)
                    { /* |x|<2**-28 */
                        if (ax < double.MinValue * 16)
                            return 0.125 * (8.0 * x + efx8 * x);  /*avoid underflow */
                        return x + efx * x;
                    }
                    z = x * x;
                    r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
                    s = one + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
                    y = r / s;
                    return x + x * y;
                }
                if (ax < 1.25)
                {      /* 0.84375 <= |x| < 1.25 */
                    s = ax - one;
                    P = pa0 + s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
                    Q = one + s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
                    if (x >= 0) return erx + P / Q; else return -erx - P / Q;
                }
                if (ax >= 6)
                {      /* inf>|x|>=6 */
                    if (x >= 0) return one - tiny; else return tiny - one;
                }

                /* Starts to lose accuracy when ax~5 */
                s = one / (ax * ax);

                if (ax < 2.85714285714285)
                { /* |x| < 1/0.35 */
                    R = ra0 + s * (ra1 + s * (ra2 + s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
                    S = one + s * (sa1 + s * (sa2 + s * (sa3 + s * (sa4 + s * (sa5 + s * (sa6 + s * (sa7 + s * sa8)))))));
                }
                else
                {    /* |x| >= 1/0.35 */
                    R = rb0 + s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
                    S = one + s * (sb1 + s * (sb2 + s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
                }
                r = Math.Exp(-ax * ax - 0.5625 + R / S);
                if (x >= 0) return one - r / ax; else return r / ax - one;
            }
            #endregion
        }


        public const double QL_Epsilon = 2.2204460492503131e-016;


        public static int quadratic_roots(double a, double b, double c, ref double x1, ref double x2)
        {
            double delta = b * b - 4.0 * a * c;
            int res = 2;
            double precision = QL_Epsilon;
            if (delta < 0)
            {
                res = 0;
            }
            else if (Math.Abs(a) < precision)
            {
                double res2 = Math.Abs(c) > precision ? -c / b : 0.0;
                x1 = x2 = res2;
            }
            else if (delta == 0)
            {
                x1 = x2 = -b;
            }
            else
            {
                double inva = 0.5 / a;
                double sqrt_delta = Math.Sqrt(delta);
                x1 = (-b - sqrt_delta) * inva;
                x2 = (-b + sqrt_delta) * inva;
            }
            return res;
        }


        public static int cubic_roots(double a, double b, double c, double d, double[] roots)
        {
            int res = 0;
            if (a != 0)
            {
                double inva = 1 / a;
                res = cubic_roots(b * inva, c * inva, d * inva, roots);
            }
            else
            {
                res = quadratic_roots(b, c, d, ref roots[0], ref roots[1]);
            }
            return res;
        }


        public static int cubic_roots(double p, double q, double r, double[] roots)
        {
            int res = 0;

            // Trivial cases
            if (r == 0)
            {
                res = quadratic_roots(1.0, p, q, ref roots[0], ref roots[1]);
                roots[res++] = 0.0;
            }
            else if (p == 0 && q == 0)
            {
                roots[0] = cubic_root(-r);
                roots[1] = roots[0];
                roots[2] = roots[0];
                res = 3;
            }
            else
            {

                // Transformation into a "X^3 + a X + b" form aka Cardan form with X=x+p/3
                double inv3 = 1.0 / 3.0;
                double po3 = p * inv3;
                // double po3q = po3*q;
                double po3sq = po3 * po3;

                double v = r + po3 * (po3sq + po3sq - q);
                double vsq = v * v;
                double qo3 = q * inv3;
                double uo3 = qo3 - po3sq;
                double u2o3 = uo3 + uo3;

                double uo3sq4 = u2o3 * u2o3;
                double uo3cu4 = uo3sq4 * uo3;
                double wsq = uo3cu4 + vsq;

                // Positive discriminant: one real root
                if (wsq > 0)
                {
                    double mcube = v <= 0 ? (-v + Math.Sqrt(wsq)) * 0.5 : (-v - Math.Sqrt(wsq)) * 0.5;
                    double m1 = cubic_root(mcube);
                    double m2 = m1 != 0 ? -uo3 / m1 : 0;
                    roots[0] = m1 + m2 - po3;
                    res = 1;
                }
                // Negative discriminant : three real roots
                else
                {
                    res = 3;
                    if (uo3 < 0)
                    {
                        double s = p > 0 ? -Math.Sqrt(-uo3) : Math.Sqrt(-uo3);
                        double scube = -s * uo3;
                        double t = -v / (scube + scube);
                        double cosk = Math.Cos(Math.Acos(t) / 3.0);
                        roots[0] = (s + s) * cosk - po3;
                        res = 1;
                        double sinsqk = 1.0 - cosk * cosk;
                        if (sinsqk >= 0)
                        {
                            double rt3sink = Math.Sqrt(3 * sinsqk);
                            roots[1] = s * (-cosk + rt3sink) - po3;
                            roots[1] = s * (-cosk - rt3sink) - po3;
                            res = 3;
                        }
                    }
                    else
                    {
                        roots[0] = cubic_root(v) - po3;
                        roots[1] = roots[0];
                        roots[2] = roots[0];
                        res = 3;
                    }
                }
            }
            newton_solve_cubic_roots(p, q, r, res, roots);
            return res;
        }

        public static double cubic_root(double x)
        {
            double res = 0;
            if (x != 0)
            {
                double absx = x;
                double signx = 1;
                if (x < 0)
                {
                    absx = -x;
                    signx = -1;
                }
                res = Math.Exp(Math.Log(absx) / 3.0) * signx;
            }
            return res;
        }

        public static void newton_solve_cubic_roots(double p, double q, double r, int n, double[] roots)
        {
            int MAX_STEP = 4;
            for (int i = 0; i < n; ++i)
            {
                for (int k = 0; k < MAX_STEP; ++k)
                {
                    double root = roots[i];
                    double err = ((root + p) * root + q) * root + r;
                    double deriv = (3 * root + 2 * p) * root + q;
                    double corr = deriv != 0 ? err / deriv : 0;
                    roots[i] -= corr;
                }
            }
        }



        public enum DiffusionType
        {
            LogNormal,
            Normal
        };




        public SABRFormula(double maturity, double modified_forward, double sigma_atm_normal, double beta, double nu, double rho, DiffusionType diffusion_type = DiffusionType.LogNormal, double is_atm_normal = 1.0, double is_input_alpha = 0.0)
        {
            this.m_maturity = maturity;
            this.m_modified_forward = modified_forward;
            this.m_beta = beta;
            this.m_nu = nu;
            this.m_rho = rho;
            this.m_diffusion_type = diffusion_type;
            if (is_input_alpha == 1.0)
            {
                m_alpha = sigma_atm_normal;
            }
            else
            {
                double vatm = (is_atm_normal == 1.0) ? ConvertToAtmLognormalVolatility(m_modified_forward, maturity, sigma_atm_normal) : sigma_atm_normal;
                compute_alpha(vatm);
            }

            m_cache = new double[1];
        }


        public void get_volatility(double[] modified_strikes,
                                 double[] volatilities,
                              DiffusionType diffusion_type = DiffusionType.LogNormal)
        {
            if (m_maturity < Double.Epsilon)
            {
                for (int i = 0; i < volatilities.Length; i++)
                {
                    volatilities[i] = 0.0;
                }
            }
            else
            {
                if (diffusion_type == DiffusionType.LogNormal)
                {
                    lognormal_smile(modified_strikes, volatilities);
                }
                else
                {
                    normal_smile(modified_strikes, volatilities);
                }
            }
        }

        public void normal_smile(double[] modified_strikes,
                                     double[] volatilities)
        {
            // m_cache will hold m_alpha*a[k] *(1 + b[k]*m_maturiy)
            Array.Resize(ref m_cache, modified_strikes.Length);
            Array.Resize(ref volatilities, modified_strikes.Length);

            double multiplier = m_nu / m_alpha;

            double b1 = (2 - 3 * m_rho * m_rho) * m_nu * m_nu / 24;

            if (m_beta == 0)
            {
                double b = b1;
                double tmp = m_alpha * (1 + b * m_maturity);
                for (int i = 0; i < modified_strikes.Length; ++i)
                {
                    volatilities[i] = (m_modified_forward - modified_strikes[i]) * multiplier;
                    m_cache[i] = tmp;
                }
            }
            else if (m_beta == 1)
            {
                double b = (-m_alpha * m_alpha + 6 * m_rho * m_alpha * m_nu) / 24 + b1;
                double a_multiplier = m_alpha * (1 + b * m_maturity);

                for (int i = 0; i < modified_strikes.Length; ++i)
                {
                    if (Double.Epsilon < modified_strikes[i])
                    {
                        double logq = Math.Log(m_modified_forward / modified_strikes[i]);
                        volatilities[i] = logq * multiplier;

                        // double a = (m_forward-strikes[i])/logq;
                        m_cache[i] = a_multiplier * (1 + b * m_maturity);
                    }
                }
            }
            else if (m_beta > 0 && m_beta < 1)
            {
                double b2part = 0.25 * m_rho * m_nu * m_beta * m_alpha;
                double b3part = -m_beta * (2 - m_beta) * m_alpha * m_alpha / 24;

                for (int i = 0; i < modified_strikes.Length; ++i)
                {
                    if (Double.Epsilon < modified_strikes[i])
                    {
                        double str = modified_strikes[i];
                        double fav = 0.5 * (m_modified_forward + str);
                        double ppow1_betaon2 = Math.Pow(fav, 1 - m_beta);
                        double ppowbetaon2 = Math.Pow(fav, m_beta);
                        double inv_ppowbetaon2 = 1 / ppowbetaon2;

                        double a = 0.0;
                        if (Math.Abs(m_modified_forward - str) < 1e-6)
                        {
                            double logq = Math.Log(m_modified_forward / str);
                            volatilities[i] = ppow1_betaon2 * logq * multiplier;

                            double sqr_logq = logq * logq;
                            double tmp = (1 - m_beta) * logq;
                            tmp *= tmp;
                            a = ppowbetaon2 * (1 + sqr_logq * (1 + sqr_logq / 80) / 24) / (1 + tmp * (1 + tmp / 80) / 24);
                        }
                        else
                        {
                            volatilities[i] = multiplier * (m_modified_forward - str) * inv_ppowbetaon2;
                            a = (1 - m_beta) * (m_modified_forward - str) / (Math.Pow(m_modified_forward, 1 - m_beta) - Math.Pow(str, 1 - m_beta));
                        }

                        double b = b1 + b2part * inv_ppowbetaon2 + b3part * inv_ppowbetaon2 * inv_ppowbetaon2;
                        double tmp2 = 1 + b * m_maturity;
                        if (tmp2 <= 0)
                        {
                            tmp2 = 0;
                        }
                        m_cache[i] = m_alpha * a * tmp2;
                    }
                }
            }

            double inv1_rho = 1 / (1 - m_rho);
            for (int i = 0; i < modified_strikes.Length; ++i)
            {
                if (Double.Epsilon < modified_strikes[i])
                {
                    double z = volatilities[i];
                    if (Math.Abs(z) < 1e-9) // Taylor expansion of z/x for small z
                    {
                        z = 1 + z * (-0.5 * m_rho + z * (1 / 6 - 0.5 * m_rho * m_rho));
                    }
                    else
                    {
                        double x = Math.Log((Math.Sqrt(1 - 2 * m_rho * z + z * z) + z - m_rho) * inv1_rho);
                        z /= x;
                    }
                    volatilities[i] = m_cache[i] * z;
                }
                else
                {
                    volatilities[i] = 0.0;
                }
            }
        }


        public void lognormal_smile(double[] modified_strikes,
                                       double[] volatilities)
        {
            // m_cache will hold m_alpha/a[k] *(1 + b[k]*m_maturiy)
            Array.Resize(ref m_cache, modified_strikes.Length);
            Array.Resize(ref volatilities, modified_strikes.Length);

            double b1 = (2 - 3 * m_rho * m_rho) * m_nu * m_nu / 24;

            if (m_beta <= 1 - LOGNORMAL_EPS)
            {
                double inv_24 = 1.0 / 24.0;
                double inv_1920 = 1.0 / 1920.0;
                double multiplier = m_nu / m_alpha;
                double p1_betaon2 = 0.5 * (1 - m_beta);

                double b2part = 0.25 * m_rho * m_nu * m_beta * m_alpha;
                double b3part = m_alpha * (1 - m_beta);
                b3part *= b3part * inv_24;

                for (int i = 0; i < modified_strikes.Length; ++i)
                {
                    if (Double.Epsilon < modified_strikes[i])
                    {
                        double logq = Math.Log(m_modified_forward / modified_strikes[i]);
                        double ppow1_betaon2 = Math.Pow(m_modified_forward * modified_strikes[i], p1_betaon2);
                        volatilities[i] = multiplier * ppow1_betaon2 * logq;

                        double inv_ppow1_betaon2 = 1.0 / ppow1_betaon2;
                        double b3 = b3part * inv_ppow1_betaon2 * inv_ppow1_betaon2;
                        double b2 = b2part * inv_ppow1_betaon2;
                        double b = b1 + b2 + b3;

                        double tmp = (1 - m_beta) * logq;
                        tmp *= tmp;
                        double a = ppow1_betaon2 * (1 + tmp * (inv_24 + tmp * inv_1920));

                        m_cache[i] = m_alpha / a * (1 + b * m_maturity);
                    }
                }
            }
            else
            {
                // a = 1 and b3 = 0 here.
                double b2 = 0.25 * m_rho * m_nu * m_alpha;
                double b = b1 + b2;
                double tmp = m_alpha * (1 + b * m_maturity);
                double multiplier = m_nu / m_alpha;

                for (int i = 0; i < modified_strikes.Length; ++i)
                {
                    if (Double.Epsilon < modified_strikes[i])
                    {
                        volatilities[i] = Math.Log(m_modified_forward / modified_strikes[i]) * multiplier;
                        m_cache[i] = tmp;
                    }
                }
            }

            double inv1_rho = 1.0 / (1 - m_rho);
            for (int i = 0; i < modified_strikes.Length; ++i)
            {
                if (Double.Epsilon < modified_strikes[i])
                {
                    double z = volatilities[i];
                    if (Math.Abs(z) < 1e-9)
                    {
                        double res = 1 + z * (-0.5 * m_rho + z * (1.0 / 6.0 * m_rho * m_rho));
                        res *= m_cache[i];
                        volatilities[i] = res;
                    }
                    else
                    {
                        double x = Math.Log((Math.Sqrt(1 - 2 * m_rho * z + z * z) + z - m_rho) * inv1_rho);
                        double res = m_cache[i] * z / x;
                        volatilities[i] = res;
                    }
                }
                else
                {
                    volatilities[i] = 0.0;
                }
            }
        }


        private void compute_alpha(double vatm)
        {
            double c0 = 0, c1 = 0, c2 = 0, c3 = 0;

            if (m_diffusion_type == DiffusionType.LogNormal)
            {
                double un = 1 - LOGNORMAL_EPS;
                bool cond = m_beta > un;
                c1 = 1 + (2 - 3 * m_rho * m_rho) * m_nu * m_nu * m_maturity / 24;
                c2 = 0.25 * m_rho * m_beta * m_nu * m_maturity;
                if (cond)
                {
                    c0 = -vatm;
                    c3 = 0;
                }
                else
                {
                    c0 = Math.Pow(m_modified_forward, (1 - m_beta));
                    c2 /= c0;
                    c3 = (1 - m_beta) * (1 - m_beta) / (c0 * c0) * m_maturity / 24;
                    c0 *= -vatm;
                }
            }
            else
            {
                bool cond = m_beta == 0;
                c1 = 1 + (2 - 3 * m_rho * m_rho) * m_nu * m_nu * m_maturity / 24;
                if (cond)
                {
                    c0 = -vatm;
                    c2 = 0;
                    c3 = 0;
                }
                else
                {
                    c0 = Math.Pow(m_modified_forward, (1 - m_beta));
                    c2 = 0.25 * m_rho * m_beta * m_nu * m_maturity / c0;
                    c3 = -m_beta * (2 - m_beta) / (c0 * c0) * m_maturity / 24;
                    c0 *= (-vatm / m_modified_forward);
                }
            }

            double[] roots = { 0.0, 0.0, 0.0 };

            int nb_roots = cubic_roots(c3, c2, c1, c0, roots);

            m_alpha = Double.MaxValue;
            for (int i = 0; i < nb_roots; ++i)
            {
                if (roots[i] > 1e-9 && roots[i] < m_alpha)
                {
                    m_alpha = roots[i];
                }
            }
            if (m_alpha == Double.MaxValue)
            {
                m_alpha = 0;
                throw new ApplicationException("SABRFormula : alpha < 0!!");
            }
        }
    }
}

namespace auxiliary
{
    public static class Math2Constants
    {

        public const double SqrtTwo = 1.41421356237309504880168872421;
        public const double SqrtPi = 1.77245385090551602729816748334;
        public const double SqrtTwoPi = 2.50662827463100050241576528481;
        public const double InvPi = 0.318309886183790671537767526745;
        public const double Inv2Pi = 0.159154943091895335768883763372;
        public const double InvSqrtTwo = 0.707106781186547524400844362105;
        public const double InvSqrtPi = 0.564189583547756286948079451561;
        public const double InvSqrtTwoPi = 0.398942280401432677939946059934;
        public const double Tiny = 1e-13;
        public const double Petit = 1e-10;
        public const double Small = 1e-6;
        public const double Bips = 1e-2;
        public const double Large = 1e6;
        public const double Big = 1e10;
        public const double Hermite_eps = 0.00000000000001;
        public const double sqrt2 = 1.4142135623731;
        public const double pim4 = 0.751125544464943;
        public const int MAXIT = 10;
    }

    public class GaussLegendre2
    {
        public double X(int i)
        {
            return rule.GetAbscissa(i);
        }
        public double W(int i)
        {
            return rule.GetWeight(i);
        }
        MathNet.Numerics.Integration.GaussLegendreRule rule;
        public GaussLegendre2(int NbPoints, double begin, double end)
        {
            rule = new MathNet.Numerics.Integration.GaussLegendreRule(begin, end, NbPoints);
        }

    }

    public static class GaussianDensity
    {

        public static class MultiNormal
        {
            public static double BiNormal(double a, double b, double r)
            {
                double r1, r2, dPb;

                if (r > 1.0 || r < -1.0)
                {
                    throw new Exception("BiNormal::Correlation must be between +1 and -1");
                }

                if (r == 1) return (GaussianDensity.Cdf(Math.Min(a, b)));
                if (r == -1) return (GaussianDensity.Cdf(a) - GaussianDensity.Cdf(-b));
                if (r == 0) return (GaussianDensity.Cdf(a) * GaussianDensity.Cdf(b));

                if (a * b * r <= 0)
                {
                    if (b <= 0) if (a <= 0) if (r < 0) dPb = BinormalDrezner(a, b, r);
                            else if (a == 0) dPb = GaussianDensity.Cdf(b) - BinormalDrezner(0, b, -r);
                            else dPb = GaussianDensity.Cdf(a) - BinormalDrezner(a, 0, -r);
                        else if (r > 0) dPb = GaussianDensity.Cdf(b) - BinormalDrezner(-a, b, -r);
                        else dPb = GaussianDensity.Cdf(a) - 0.5 + BinormalDrezner(-a, 0, r);
                    else if (a <= 0) if (r > 0) dPb = GaussianDensity.Cdf(a) - BinormalDrezner(a, -b, -r);
                        else dPb = GaussianDensity.Cdf(b) - 0.5 + BinormalDrezner(0, -b, r);
                    else dPb = GaussianDensity.Cdf(a) + GaussianDensity.Cdf(b) - 1 + BinormalDrezner(-a, -b, r);
                }
                else // so necessarily a and b and r are <> 0
                {
                    //r2 = sqrt(a * a - 2 * r * a * b + b * b);
                    //r1 = (r * a - b) * Sgn(a) / r2;
                    //r2 = (r * b - a) * Sgn(b) / r2;

                    // More "precise" way
                    double x = r - b / a;
                    if (x != 0)
                    {
                        r1 = Math.Sign(x) * Math.Sqrt(1.0 / (1.0 + (1.0 - r * r) / (x * x)));
                    }
                    else
                    {
                        r1 = 0.0;
                    }

                    double y = r - a / b;
                    if (y != 0)
                    {
                        r2 = Math.Sign(y) * Math.Sqrt(1.0 / (1.0 + (1.0 - r * r) / (y * y)));
                    }
                    else
                    {
                        r2 = 0.0;
                    }

                    dPb = Math.Sign(a) * Math.Sign(b) / 4 - 0.25;

                    if (r1 <= 0) if (a <= 0) dPb += BinormalDrezner(a, 0, r1);
                        else dPb += GaussianDensity.Cdf(a) - 0.5 + BinormalDrezner(-a, 0, r1);
                    else if (a <= 0) dPb += GaussianDensity.Cdf(a) - BinormalDrezner(a, 0, -r1);
                    else dPb += 0.5 - BinormalDrezner(-a, 0, -r1);
                    if (r2 <= 0) if (b <= 0) dPb += BinormalDrezner(b, 0, r2);
                        else dPb += GaussianDensity.Cdf(b) - 0.5 + BinormalDrezner(-b, 0, r2);
                    else if (b <= 0) dPb += GaussianDensity.Cdf(b) - BinormalDrezner(b, 0, -r2);
                    else dPb += 0.5 - BinormalDrezner(-b, 0, -r2);
                }

                return (dPb);
            }

            public static double BiNormal2(double a, double b, double r)
            {
                if (r > 1.0 || r < -1.0)
                {
                    throw new Exception("BiNormal2::Correlation must be between +1 and -1");
                }

                double h1, h2;
                double LH = 0.0; double h12;
                double h3, h5, h6, h7, h8;
                double r1, r2, r3, rr;
                double AA, ab;
                double Bivarcumnorm = 0.0;

                double[] x = { 0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992 };
                double[] W = { 0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042 };

                h1 = a;
                h2 = b;
                h12 = (h1 * h1 + h2 * h2) / 2;

                if (Math.Abs(r) >= 0.7)
                {
                    r2 = 1.0 - r * r;

                    r3 = Math.Sqrt(r2);

                    if (r < 0.0)
                    {
                        h2 = -h2;
                    }
                    h3 = h1 * h2;
                    h7 = GaussianDensity.SafeExp(-h3 / 2);
                    if (Math.Abs(r) < 1.0)
                    {
                        h6 = Math.Abs(h1 - h2);
                        h5 = h6 * h6 / 2.0;
                        h6 = h6 / r3;
                        AA = 0.5 - h3 / 8.0;
                        ab = 3.0 - 2.0 * AA * h5;
                        LH = 0.13298076 * h6 * ab * (1.0 - GaussianDensity.Cdf(h6)) - GaussianDensity.SafeExp(-h5 / r2) * (ab + AA * r2) * 0.053051647;
                        for (int i = 0; i < 5; ++i)
                        {
                            r1 = r3 * x[i];
                            rr = r1 * r1;
                            r2 = Math.Sqrt(1.0 - rr);
                            if (h7 == 0.0)
                            {
                                h8 = 0.0;
                            }
                            else
                            {
                                h8 = GaussianDensity.SafeExp(-h3 / (1.0 + r2)) / r2 / h7;
                            }
                            LH = LH - W[i] * GaussianDensity.SafeExp(-h5 / rr) * (h8 - 1.0 - AA * rr);
                        }
                    }
                    Bivarcumnorm = LH * r3 * h7 + GaussianDensity.Cdf(Math.Min(h1, h2));

                    if (r < 0.0)
                    {
                        Bivarcumnorm = GaussianDensity.Cdf(h1) - Bivarcumnorm;
                    }
                }
                else
                {
                    h3 = h1 * h2;
                    if (r != 0.0)
                    {
                        for (int i = 0; i < 5; ++i)
                        {
                            r1 = r * x[i];
                            r2 = 1.0 - r1 * r1;
                            LH = LH + W[i] * GaussianDensity.SafeExp((r1 * h3 - h12) / r2) / Math.Sqrt(r2);
                        }
                    }
                    Bivarcumnorm = GaussianDensity.Cdf(h1) * GaussianDensity.Cdf(h2) + r * LH;
                }
                return Bivarcumnorm;
            }

            public static double BinormalDrezner(double a, double b, double r)
            {
                double[] pdA = { 0.3253030, 0.4211071, 0.1334425, 0.006374323 };
                double[] pdB = { 0.1337764, 0.6243247, 1.3425378, 2.2626645 };
                double dSum, c, d;
                int i, j;

                // make sure the constraints are satisfied
                if (r > 0.0 || r < -1.0)
                {
                    throw new Exception("Correlation must be between -1 and 0 for the BiNormalDrezner algorithm");
                }

                if (a > 0.0 || b > 0.0)
                {
                    throw new Exception("a and b must be less than or equal to zero for the BiNormalDrezner algorithm");
                }

                dSum = 0;
                d = Math.Sqrt(2.0 * (1.0 - r * r));
                c = a / d;
                d = b / d;

                for (i = 0; i <= 3; i++)
                    for (j = 0; j <= 3; j++)
                        dSum += pdA[i] * pdA[j] * Math.Exp(c * (2 * pdB[i] - c)
                                                     + d * (2 * pdB[j] - d)
                                                     + 2 * r * (pdB[i] - c) * (pdB[j] - d));

                return (Math.Sqrt(1.0 - r * r) / Math.PI * dSum);
            }

            //*************************************************************************
            //This function estimates the joint probability for a multivariate normal
            //law, given an upper bound vector Hi, a correlation matrix Rho.
            //It is based on quadrature estimation, and was described by Zvi Drezner in
            //"Computation of the multivariate normal integral", ACM Transactions on
            //Mathematical Software, Vol. 18, N4, pp 470-480, 1992.
            //*************************************************************************

            public static double MultiNormalS(ref double[] UpperBounds, ref double[][] Correlations, int K)
            {
                int i, j;
                int size = (int)UpperBounds.Length;

                double[] oldHi;
                double[][] oldCorr;
                double oldresu;

                oldHi = new double[size];
                oldCorr = new double[size][];
                for (i = 0; i < size; ++i)
                {
                    oldCorr[i] = new double[size];
                }

                for (i = 0; i < size; ++i)
                {
                    for (j = 0; j < size; ++j)
                    {
                        oldCorr[i][j] = Correlations[i][j];
                        if (i == j)
                        {
                            if (Correlations[i][j] != 1.0) throw new Exception("Invalid Correlation Matrix");
                        }
                        else
                        {
                            if (Correlations[i][j] > 0.99 || Correlations[i][j] < -0.99) throw new Exception("Correlations must be between -0.99 and +0.99 for this approximation");
                            if (Correlations[i][j] != Correlations[j][i]) throw new Exception("Unsymmetric Correlation Matrix");
                        }
                    }
                    oldHi[i] = UpperBounds[i];
                }

                MultiNormQuad_p(oldHi, oldCorr, size, K, out oldresu);

                return oldresu;
            }

            //*******************************************************************************
            public static void MultiNormQuad_p(double[] Hi, double[][] Corr, int size, int K, out double resu)
            {
                if (K > 0)
                {
                    resu = MultiNormQuad(Hi, Corr, size, K);
                }
                else
                {
                    double pk7 = MultiNormQuad(Hi, Corr, size, 9);
                    double pk6 = MultiNormQuad(Hi, Corr, size, 8);
                    double pk5 = MultiNormQuad(Hi, Corr, size, 7);
                    resu = pk7 + (pk7 - pk6) * (pk7 - pk6) / (pk6 - pk5 - pk7 + pk6);
                }
            }

            //*************************************************************************
            //This function estimates the joint probability for a multivariate normal
            //law, given an upper bound vector Hi, a correlation matrix Rho.
            //It is based on quadrature estimation, and was described by Zvi Drezner in
            //"Computation of the multivariate normal integral", ACM Transactions on
            //Mathematical Software, Vol. 18, N4, pp 470-480, 1992.
            //*************************************************************************

            public static double MultiNormQuad(double[] Hi, double[][] Corr, int size, int K)
            {
                int i, j;
                double[] H1, H2;
                double resu;

                H1 = new double[size];
                H2 = new double[size];

                double[][] R1, R2;

                R1 = new double[size][];
                R2 = new double[size][];

                for (i = 0; i < size; i++)
                {
                    R1[i] = new double[size];
                    R2[i] = new double[size];
                }

                for (i = 0; i < size; i++)
                {
                    H1[i] = Hi[i];
                    H2[i] = Hi[i];
                    for (j = 0; j < size; j++)
                    {
                        R1[i][j] = Corr[i][j];
                        R2[i][j] = Corr[i][j];
                    }
                }

                for (i = 0; i < size; i++)
                {
                    if (Hi[i] > 0 && Hi[i] != 1e50)
                    {
                        H1[i] = 1e50;
                        H2[i] = -Hi[i];

                        for (j = 0; j < size; j++)
                        {
                            R1[i][j] = 1e50;
                            R1[j][i] = 1e50;
                            if (j != i)
                            {
                                R2[i][j] = -Corr[i][j];
                                R2[j][i] = -Corr[j][i];
                            }
                        }
                        return (MultiNormQuad(H1, R1, size, K) - MultiNormQuad(H2, R2, size, K));
                    }
                }

                Quad_spec(Hi, Corr, size, K, out resu);
                return resu;
            }

            //This function computes the integral for the case where all the h_i are
            //negative or equal to infinity. It is called in the main program recursively so that
            //it should apply to vectors h_i with negative or/and infinite values (noted as 1e50).

            //First step is to perform a change of variable, which will allow the use of the
            //quadrature with Hermite polynomials to estimate the pdfs
            //First, invert the correlation matrix and compute its determinant. For this,
            //the program uses a matrix inversionthe function Mat_Invert (int size, double **Mat, double *pivot, double det),
            //which call the Nag function performing the LU decomposition f03afc

            //WARNING!!!!!!!!!!!!!!!! The weights and zeros used are DIFFERENT from the usual 
            //values used in Gauss-Hermite quadrature, because the bounds of integration
            //range from zero to infinity. See "Gaussian quadrature for the integrals []",
            //N.M. Steen, G.D. Byrne and E.M. Gelbard, Mathematical Computing, vol. 23, pp 661-671, 1968.

            public static void Quad_spec(double[] Hi, double[][] Corr, int size, int K, out double resu)
            {

                double[][] H_points = new double[][] {
        new double[] {3.00193931060839e-1, 1.25242104533372} ,
        new double[] {1.90554149798192e-1, 8.48251867544577e-1, 1.79977657841573 } ,
        new double[] {1.33776446996068e-1, 6.24324690187190e-1, 1.34253782564499,2.26266447701036} ,
        new double[] {1.00242151968216e-1, 4.82813966046201e-1, 1.06094982152572,1.77972941852026,2.66976035608766} ,
        new double[] {7.86006594130979e-2, 3.86739410270631e-1, 8.66429471682044e-1,1.46569804966352,2.17270779693900,3.03682016932287} ,
        new double[] {6.37164846067008e-2, 3.18192018888619e-1, 7.24198989258373e-1,1.23803559921509,1.83852822027095,2.53148815132768,3.37345643012458} ,
        new double[] {0.0529786439318514, 0.267398372167767, 0.616302884182402,1.06424631211623,1.58885586227006,2.18392115309586,2.86313388370808,3.6860071627244} ,
        new double[] {0.0449390308011934, 0.228605305560535, 0.532195844331646,0.927280745338081,1.39292385519588,1.91884309919743,2.50624783400574,3.17269213348124,3.97889886978978} ,
        new double[] {0.0387385243257289, 0.198233304013083, 0.465201111814767,0.816861885592273,1.23454132402818,1.70679814968913,2.22994008892494,2.80910374689875,3.46387241949586,4.25536180636608} };

                double[][] H_weights = new double[][]  {
        new double[] {6.40529179684379e-1,2.45697745768379e-1},
        new double[] {4.46029770466658e-1,3.96468266998335e-1,4.37288879877644e-2},
        new double[] {3.25302999756919e-1,4.21107101852062e-1,1.33442500357520e-1,6.37432348625728e-3},
        new double[] {2.48406152028443e-1,3.92331066652399e-1,2.11418193076057e-1,3.32466603513439e-2,8.24853344515628e-4},
        new double[] {1.96849675488598e-1,3.49154201525395e-1,2.57259520584421e-1,7.60131375840057e-2,6.85191862513596e-3,9.84716452019267e-5},
        new double[] {1.60609965149261e-1,3.06319808158099e-1,2.75527141784905e-1,1.20630193130784e-1,2.18922863438067e-2,1.23644672831056e-3,1.10841575911059e-5},
        new double[] {0.13410918845336,0.26833075447264,0.275953397988422,0.15744828261879,0.0448141099174625,0.00536793575602526,2.02063649132407e-4,1.19259692659532e-6},
        new double[] {0.114088970242118,0.235940791223685,0.266425473630253,0.183251679101663,0.0713440493066916,0.0139814184155604,0.00116385272078519,0.305670214897831e-4,1.2379051133749e-7},
        new double[] {0.0985520975191087,0.208678066608185,0.252051688403761,0.198684340038387,0.097198422760062,0.0270244164355446,0.00380464962249537,2.28886243044656e-4,4.34534479844469e-6,1.24773714817825e-8}};
                //{5.54433663102343e-2,1.24027738987730e-1,1.75290943892075e-1,1.91488340747342e-1,1.63473797144070e-1,1.05937637278492e-1,5.00270211534535e-2,1.64429690052673e-2,3.57320421428311e-3,4.82896509305201e-4,3.74908650266318e-5,1.49368411589636e-6,2.55270496934465e-8,1.34217679136316e-10,9.56227446736465e-14};
                double determinant, sum;
                int i, j, co = 0, k = 0, dim; // l = 0;
                //int l;
                for (i = 0; i < size; i++)
                {
                    if (Hi[i] == 1e50) co++;
                }
                double[] h_p;
                int[] indice;
                dim = size - co;
                indice = new int[dim];
                if (dim != 0)
                {
                    h_p = new double[dim];
                    k = 0;
                    for (i = 0; i < size; i++)
                    {
                        if (Hi[i] != 1e50)
                        {
                            h_p[k] = Hi[i];
                            indice[k] = i;
                            k++;
                        }
                    }
                    double[][] R_p;
                    R_p = new double[dim][];
                    k = 0;
                    //l = 0;
                    for (i = 0; i < dim; i++) R_p[i] = new double[dim];
                    for (i = 0; i < dim; i++)
                    {
                        for (j = 0; j < dim; j++) R_p[i][j] = Corr[indice[i]][indice[j]];
                    }
                    double[][] R_i;
                    R_i = new double[dim][];
                    for (i = 0; i < dim; i++) R_i[i] = new double[dim];
                    var myMat = CreateMatrix.Dense<double>(dim, dim, 0.0);
                    for (i = 0; i < dim; ++i)
                    {
                        for (j = 0; j < dim; ++j)
                        {
                            myMat[i, j] = R_p[i][j];
                        }
                    }
                    //Invert the Matrix from myMat to myInverse and compute the determinant
                    determinant = myMat.Determinant();
                    var myInverse = myMat.Inverse();
                    for (i = 0; i < dim; ++i)
                    {
                        for (j = 0; j < dim; ++j)
                        {
                            R_i[i][j] = myInverse[i, j];
                        }
                    }
                    double epsilon = 0.0001;  // O.C.
                    double[] h_pp;
                    h_pp = new double[dim];
                    for (i = 0; i < dim; i++)
                        // h_pp[i] = h_p[i] * (Math.Sqrt(R_i[i][i] / 2));
                        h_pp[i] = h_p[i] * (Math.Sqrt(Math.Max(R_i[i][i], epsilon) / 2));  // O.C.
                    double[][] R_pp;
                    R_pp = new double[dim][];
                    for (i = 0; i < dim; i++) R_pp[i] = new double[dim];
                    for (i = 0; i < dim; i++)
                    {
                        // for (j = 0; j < dim; j++) R_pp[i][j] = R_i[i][j] / (Math.Sqrt(R_i[i][i] * R_i[j][j]));
                        for (j = 0; j < dim; j++) R_pp[i][j] = R_i[i][j] / (Math.Sqrt(Math.Max(R_i[i][i] * R_i[j][j], epsilon)));  //O.C.
                    }
                    //Computation of the constant C^2
                    double C_2 = 1;
                    for (i = 0; i < dim; i++) C_2 *= (2 / R_i[i][i]);
                    C_2 /= (Math.Pow(2 * Math.PI, dim) * determinant);
                    double[] int_1;
                    int_1 = new double[dim];
                    double[] B;
                    B = new double[dim];
                    double[] Y;
                    Y = new double[dim];
                    i = 0;
                    int[] kk;
                    kk = new int[dim];
                    for (i = 0; i < dim; i++) kk[i] = 0;
                    int pp = 0;
                    double cut = 15.0;
                    i = 0;
                    double solve = 0;
                    Flag_1: sum = 0;
                    for (i = 0; i < dim; i++)
                    {
                        pp = kk[i];
                        sum += H_points[K - 2][pp] * H_points[K - 2][pp];
                        Y[i] = H_points[K - 2][pp] - h_pp[i];
                    }
                    for (i = 0; i < dim; i++)
                    {
                        for (j = 0; j < dim; j++) sum -= Y[i] * Y[j] * R_pp[i][j];
                    }
                    sum = Math.Exp(sum);
                    for (i = 0; i < dim; i++) sum *= H_weights[K - 2][kk[i]];
                    if (sum < -cut) goto Flag_2;
                    solve += sum;
                    Flag_2: for (i = 0; i < dim; i++)
                    {
                        kk[i] = kk[i] + 1;
                        if (kk[i] < K) goto Flag_1;
                        kk[i] = 0;
                    }
                    solve *= Math.Sqrt(C_2);
                    resu = solve;
                }
                else resu = 1.0;
                return;
            }

            //-----------------------------------------------------------
            //This function is used in the Drezner's algorithm
            //-----------------------------------------------------------

            public static double F_mult(ref double[] inp, int size)
            {
                double res = 1;
                int i;
                for (i = 0; i < size; i++) res *= inp[i];
                return res;
            }

            //Computation of the function f
            public static double ff(double[] inp, double[] h_pp, double[][] R_pp, int size)
            {
                int i, j;
                double outRes = 0;
                double[] vect1 = new double[size];
                double[] vect2 = new double[size];

                for (i = 0; i < size; i++) vect2[i] = 0;
                for (i = 0; i < size; i++) vect1[i] = (h_pp[i] - inp[i]);

                for (i = 0; i < size; i++)
                {
                    for (j = 0; j < size; j++) vect2[i] += R_pp[i][j] * vect1[j];
                }

                for (i = 0; i < size; i++) outRes += inp[i] * inp[i] - vect1[i] * vect2[i];

                outRes = Math.Exp(outRes);

                return outRes;
            }
        }

        #region Normal density function (1)

        /// <summary>
        /// Computes the Gaussian density function.
        /// </summary>
        /// <param name="z">The z.</param>
        /// <returns></returns>
        public static double Pdf(double z)
        {
            return Math2Constants.InvSqrtTwoPi * Math.Exp(-0.5 * z * z);
        }

        #endregion

        #region Normal cumulative function (1)

        /// <summary>
        /// Computes the Gaussian cumulative density function.
        /// </summary>
        /// <param name="z">The z.</param>
        /// <returns></returns>
        public static double Cdf(double z)
        {
            //
            // The computation algorithm is inspired from the following implementation:
            //    http://mpa.itc.it/markus/grass63progman/normp_8c-source.html
            //
            // The function returns normal distribution probabilities accurate to 1.e-15.
            //
            // This algorithlm is based upon algorithm 5666 for the error function, from
            // Hart, J.F. et al, 'Computer Approximations', Wiley 1968.
            //
            double cutoff = 7.071;
            double zabs = 0.0;
            double expntl = 0.0;
            double p = 0.0;
            double pdf = 0.0;

            zabs = Math.Abs(z);

            // Specific case where we are out of the bound of the algorithm precision.
            if (zabs > 37.0)
            {
                if (z > 0.0)
                    return 1.0;
                else
                    return 0.0;
            }

            // General case computation.
            expntl = Math.Exp(-0.5 * zabs * zabs); ;
            pdf = expntl * Math2Constants.InvSqrtTwoPi;

            if (zabs < cutoff)
                p = expntl * ((((((_p[6] * zabs + _p[5]) * zabs + _p[4])
                * zabs + _p[3]) * zabs + _p[2]) * zabs + _p[1]) * zabs + _p[0])
                / (((((((_q[7] * zabs + _q[6]) * zabs + _q[5]) * zabs + _q[4])
                * zabs + _q[3]) * zabs + _q[2]) * zabs + _q[1]) * zabs + _q[0]);
            else
                p = pdf / (zabs + 1.0 / (zabs + 2.0 / (zabs + 3.0 / (zabs + 4.0 /
                (zabs + 0.65)))));

            if (z > 0.0)
                p = 1.0 - p;

            return p;
        }

        #endregion

        #region Inverse normal cumulative function (1)

        public static double SafeExp(double x)
        {
            double y = Math.Min(x, 300.0);
            return Math.Exp(y);
        }

        public static double InvCdfSafe(double p)
        {
            double q = Math.Max(p, Math2Constants.Tiny);
            q = Math.Min(q, 1.0 - Math2Constants.Tiny);
            return InvCdf(q);
        }

        /// <summary>
        /// Computes the inverse of the normla cumulative function.
        /// </summary>
        /// <param name="p">The p.</param>
        /// <returns></returns>
        public static double InvCdf(double p)
        {
            //
            // The computation algorithm is inspired from the following imlpementation
            //    http://home.online.no/~pjacklam/notes/invnorm/impl/sprouse/ltqnorm.c
            //
            // Lower tail quantile for standard normal distribution function.
            // This function returns an approximation of the inverse cumulative
            // standard normal distribution function.  I.e., given P, it returns
            // an approximation to the Xnorm satisfying P = Pr{Z <= Xnorm} where Z is a
            // random variable from the standard normal distribution.
            // The algorithm uses a minimax approximation by rational functions
            // and the result has a relative error whose absolute value is less
            // than 1.15e-9.
            //
            double low = 0.02425;
            double high = 0.97575;

            double q = 0.0;
            double r = 0.0;
            double e = 0.0;
            double u = 0.0;
            double x = 0.0;

            // Standard case.
            if (p > 0.0 && p < 1.0)
            {
                // Rational approximation for lower region
                if (p < low)
                {
                    q = Math.Sqrt(-2.0 * Math.Log(p));
                    x = (((((_c[0] * q + _c[1]) * q + _c[2]) * q + _c[3]) * q + _c[4]) * q + _c[5]) /
                        ((((_d[0] * q + _d[1]) * q + _d[2]) * q + _d[3]) * q + 1.0);
                }

                // Rational approximation for upper region
                else if (p > high)
                {
                    q = Math.Sqrt(-2.0 * Math.Log(1 - p));
                    x = -(((((_c[0] * q + _c[1]) * q + _c[2]) * q + _c[3]) * q + _c[4]) * q + _c[5]) /
                        ((((_d[0] * q + _d[1]) * q + _d[2]) * q + _d[3]) * q + 1.0);
                }
                // Rational approximation for central region
                else
                {
                    q = p - 0.5;
                    r = q * q;
                    x = (((((_a[0] * r + _a[1]) * r + _a[2]) * r + _a[3]) * r + _a[4]) * r + _a[5]) * q /
                        (((((_b[0] * r + _b[1]) * r + _b[2]) * r + _b[3]) * r + _b[4]) * r + 1.0);
                }

                // Halley's rational precision improvement.
                e = GaussianDensity.Cdf(x) - p;
                u = e * Math2Constants.SqrtTwoPi * Math.Exp(0.5 * x * x);
                x = x - u / (1 + 0.5 * u * x);

                return x;
            }

            // Error case. Should never reach this point.
            throw new Exception("Invalid input probability for inverse normal cumulative computation");
        }

        #endregion

        #region Arrays for the direct normal cumulative

        // Coefficients for the direct cumulative rational approximation.

        private static readonly double[] _p =
        {
            220.2068679123761,
            221.2135961699311,
            112.079291497870,
            33.91286607838300,
            6.37396220353165,
            0.7003830644436881,
            0.352624965998910e-1
        };

        private static readonly double[] _q =
        {
            440.4137358247522,
            793.8265125199484,
            637.3336333788311,
            296.5642487796737,
            86.78073220294608,
            16.06417757920695,
            1.755667163182642,
            0.8838834764831844e-1
        };

        #endregion

        #region Arrays for the inverse cumulative

        // Coefficients for the inverse cumulative rational approximation.

        private static readonly double[] _a =
        {
       -3.969683028665376e+01,
       2.209460984245205e+02,
       -2.759285104469687e+02,
       1.383577518672690e+02,
       -3.066479806614716e+01,
       2.506628277459239e+00
        };

        private static readonly double[] _b =
        {
       -5.447609879822406e+01,
       1.615858368580409e+02,
       -1.556989798598866e+02,
       6.680131188771972e+01,
       -1.328068155288572e+01
        };

        private static readonly double[] _c =
        {
        -7.784894002430293e-03,
       -3.223964580411365e-01,
       -2.400758277161838e+00,
       -2.549732539343734e+00,
       4.374664141464968e+00,
       2.938163982698783e+00
        };

        private static readonly double[] _d =
        {
       7.784695709041462e-03,
       3.224671290700398e-01,
       2.445134137142996e+00,
       3.754408661907416e+00
        };

        #endregion

    }

    public class MultiNormal2
    {
        static double NormalCDFfromMathNet(double x)
        {
            return MathNet.Numerics.Distributions.Normal.CDF(0, 1, x);
        }


        public static double BiNormalApprox(double X, double Y, double rho, int nbpoints, double Integrationstart)
        {
            GaussLegendre2 gl = new GaussLegendre2(nbpoints, Integrationstart, Y);
            double res = 0.0;
            double denom = Math.Sqrt(1.0 - rho * rho);
            for (int i = 0; i < nbpoints; i++)
            {
                double y = gl.X(i);
                res += Math.Exp(-y * y / 2.0) * GaussianDensity.Cdf((X - y * rho) / denom) * gl.W(i);
            }
            return res * Math2Constants.InvSqrtTwoPi;

        }

        public static double DL3Drho12(double h1, double h2, double h3, double rho12, double rho13, double rho23)
        {
            double t = 1.0 - rho12 * rho12 - rho13 * rho13 - rho23 * rho23 + 2.0 * rho12 * rho13 * rho23;
            double q12 = (1.0 - rho12 * rho12);
            if (t > 0)
            {
                return Math.Exp(-(h1 * h1 + h2 * h2 - 2.0 * rho12 * h1 * h2) / (2.0 * q12)) * GaussianDensity.Cdf((h3 * q12 - h1 * (rho13 - rho23 * rho12) - h2 * (rho23 - rho13 * rho12)) / Math.Sqrt(q12 * t)) / Math.Sqrt(q12) * Math2Constants.Inv2Pi;

            }
            else
            {
                return Math.Exp(-(h1 * h1 + h2 * h2 - 2.0 * rho12 * h1 * h2) / (2.0 * q12)) / Math.Sqrt(q12) * Math2Constants.Inv2Pi;
            }
        }

        public static double DL3Drho13(double h1, double h2, double h3, double rho12, double rho13, double rho23)
        {
            return DL3Drho12(h1, h3, h2, rho13, rho12, rho23);
        }

        public static double TriNormalApprox(double h1, double h2, double h3, double rho12, double rho13, double rho23, int nbpoints, double Integrationstart)
        {
            double res = GaussianDensity.Cdf(h1) * BiNormalApprox(h2, h3, rho23, nbpoints, Integrationstart);
            GaussLegendre2 gl = new GaussLegendre2(nbpoints, 0.0, 1.0);
            for (int i = 0; i < nbpoints; i++)
            {
                double t = gl.X(i);
                res += (rho12 * DL3Drho12(h1, h2, h3, t * rho12, t * rho13, rho23) + rho13 * DL3Drho13(h1, h2, h3, t * rho12, t * rho13, rho23)) * gl.W(i);
            }
            return res;
        }
        public static double TriNormalApprox2(double h1, double h2, double h3, double rho12, double rho13, double rho23, int nbpoints)
        {
            double res = GaussianDensity.Cdf(h1) * GaussianDensity.MultiNormal.BiNormal2(h2, h3, rho23);
            GaussLegendre2 gl = new GaussLegendre2(nbpoints, 0.0, 1.0);
            for (int i = 0; i < nbpoints; i++)
            {
                double t = gl.X(i);
                res += (rho12 * DL3Drho12(h1, h2, h3, t * rho12, t * rho13, rho23) + rho13 * DL3Drho13(h1, h2, h3, t * rho12, t * rho13, rho23)) * gl.W(i);
            }
            return res;
        }
    }

    public class BlackScholes
    {
        public static double Call(double f, double k, double t, double v)
        {
            if (k == 0)
            {
                return f;
            }
            else
            {
                return f * GaussianDensity.Cdf((Math.Log(f / k) + v * v * t / 2.0) / (v * Math.Sqrt(t))) -
                            k * GaussianDensity.Cdf((Math.Log(f / k) - v * v * t / 2.0) / (v * Math.Sqrt(t)));
            }
        }
        public static double Call2(double f, double k, double t, double v)
        {
            if (k == 0)
            {
                return f;
            }
            else
            {
                return f * GaussianDensity.Cdf((Math.Log(f / k) + v * v * t / 2.0) / (v * Math.Sqrt(t))) -
                            k * GaussianDensity.Cdf((Math.Log(f / k) - v * v * t / 2.0) / (v * Math.Sqrt(t)));
            }
        }
        public static double Put(double f, double k, double t, double v)
        {
            return Call(f, k, t, v) - (f - k);
        }

        public static double PhiVol(double m, double u, double target)
        {
            return Math.Exp(m) * GaussianDensity.Cdf(m / u + u / 2.0) - GaussianDensity.Cdf(m / u - u / 2.0) - target;
        }

        public static double PhiVolAsymth(double m, double u, double target)
        {
            return Math.Exp(m) * GaussianDensity.Cdf(m / u + u / 2.0) - GaussianDensity.Cdf(m / u - u / 2.0) - target;
        }

        public static double ImpVolBS(double f, double k, double t, double opt)
        {
            double o = opt / k; double m = Math.Log(f / k);
            double phivol(double u) => PhiVol(m, u, o);
            double solu = Brent.FindRoot(phivol, 0.0001, 7, 0.0000000001, 5000000);
            return solu / Math.Sqrt(t);
        }


        public static double ImpliedVolatilityX(double forward, double strike, double maturity, double raw_target_price, int type)
        {
            double IntrinseqVal;
            if (type == 0 || type == 1)
            {
                IntrinseqVal = Math.Max(forward - strike, 0.0);  // Call
            }
            else
            {
                if (type == -1)
                {
                    IntrinseqVal = Math.Max(strike - forward, 0.0); //Put
                }
                else
                {
                    if (strike >= forward)
                    {
                        IntrinseqVal = Math.Max(forward - strike, 0.0);  // Call
                    }
                    else
                    {
                        IntrinseqVal = Math.Max(strike - forward, 0.0); //Put
                    }
                }
            }
            double target_price = raw_target_price - IntrinseqVal;
            return ImpVolBS(forward, strike, maturity, target_price);
        }



    }
    public static class XLMultiNormal
    {
        public static double calculateMultiNormalS(
        double[] UpperBounds,
        double[,] Correlations)
        {
            int n = UpperBounds.Length;
            double[][] _Correlations = new double[n][];
            for (int i = 0; i < n; i++)
            {
                _Correlations[i] = new double[n];
                for (int j = 0; j < n; j++)
                {
                    _Correlations[i][j] = Correlations[i, j];
                }
            }
            int K = 0;
            double result = GaussianDensity.MultiNormal.MultiNormalS(ref UpperBounds, ref _Correlations, K);
            return result;
        }
    }

    public static class WorstOf3
    {

        public static double N3Dis(double x, double y, double z, double rho12, double rho13, double rho23)
        {

            double[] UpperBounds = new double[3] { x, y, z };
            double[][] _Correlations = new double[3][];
            for (int i = 0; i < 3; i++)
            {
                _Correlations[i] = new double[3];
                _Correlations[i][i] = 1.0;
            }

            _Correlations[0][1] = _Correlations[1][0] = Math.Min(rho12, 0.98);

            _Correlations[0][2] = _Correlations[2][0] = Math.Min(rho13, 0.98);

            _Correlations[1][2] = _Correlations[2][1] = Math.Min(rho23, 0.98);

            int K = 9;
            double result = GaussianDensity.MultiNormal.MultiNormalS(ref UpperBounds, ref _Correlations, K);
            return result;
        }

        public static double N3Dis2(double x, double y, double z, double rho12, double rho13, double rho23, int N3Flag, int nb = 35, double start = -7)
        {
            if (N3Flag == 1) return MultiNormal2.TriNormalApprox(x, y, z, rho12, rho13, rho23, nb, start);
            if (N3Flag == 2) return MultiNormal2.TriNormalApprox2(x, y, z, rho12, rho13, rho23, nb);
            return N3Dis(x, y, z, rho12, rho13, rho23);
        }


        public static double projectedCorrelation(double[] sigma, double[,] rho, int i, int j, int k)
        {
            double cap = 0.99 - Math2Constants.Small;
            double V2 = (sigma[i] * sigma[i] + sigma[k] * sigma[k] - 2.0 * rho[i, k] * sigma[i] * sigma[k]) * (sigma[j] * sigma[j] + sigma[k] * sigma[k] - 2.0 * rho[j, k] * sigma[j] * sigma[k]);
            double V1 = (rho[i, j] * sigma[i] * sigma[j] - rho[i, k] * sigma[i] * sigma[k] - rho[k, j] * sigma[k] * sigma[j] + sigma[k] * sigma[k]);
            double result;
            if (V2 > 0.0)
            {
                result = V1 / Math.Sqrt(V2);
            }
            else
            {
                if (V1 >= 0.0)
                {
                    result = 1.0;
                }
                else
                {
                    result = -1.0;
                }
            }
            result = Math.Max(result, -cap);
            result = Math.Min(result, cap);

            return result;

        }

        public static double dplusij(double[] s, double T, double[] sigma, double[,] rho, int i, int j)
        {
            double sigmaij2 = sigma[i] * sigma[i] + sigma[j] * sigma[j] - 2.0 * rho[i, j] * sigma[i] * sigma[j];
            if (sigmaij2 < 0.000001) return 10.0 * Math.Sign(Math.Log(s[i] / s[j]));
            return (Math.Log(s[i] / s[j]) + 0.5 * sigmaij2 * T) / Math.Sqrt(sigmaij2 * T);
        }

        public static double dmoinsij(double[] s, double T, double[] sigma, double[,] rho, int i, int j)
        {
            double sigmaij2 = sigma[i] * sigma[i] + sigma[j] * sigma[j] - 2.0 * rho[i, j] * sigma[i] * sigma[j];
            if (sigmaij2 < 0.000001) return 10.0 * Math.Sign(Math.Log(s[i] / s[j]));
            double r = (Math.Log(s[i] / s[j]) - 0.5 * sigmaij2 * T) / Math.Sqrt(sigmaij2 * T);
            return r;
        }

        public static double dplusi(double[] s, double K, double T, double[] sigma, int i)
        {
            if (sigma[i] < 0.000001) return 10.0 * Math.Sign(Math.Log(s[i] / K));
            return (Math.Log(s[i] / K) + 0.5 * sigma[i] * sigma[i] * T) / (sigma[i] * Math.Sqrt(T));
        }

        public static double dmoinsi(double[] s, double K, double T, double[] sigma, int i)
        {
            if (sigma[i] < 0.000001) return 10.0 * Math.Sign(Math.Log(s[i] / K));
            return (Math.Log(s[i] / K) - 0.5 * sigma[i] * sigma[i] * T) / (sigma[i] * Math.Sqrt(T));

        }

        //Max{S1,S2,S3,S4}
        public static double Best4Underlyings(double S1, double S2, double S3, double S4, double T, double sigma1, double sigma2, double sigma3, double sigma4, double rho12, double rho13, double rho14, double rho23, double rho24, double rho34, int N3Flag, int nb, double start) // N3Flag=0 -> N3Dis usuelle, N3Flag=1 -> N3Dis2
        {
            double[] s = new double[] { S1, S2, S3, S4 };
            double[] s1 = new double[] { S2, S3, S4 };
            double[] s2 = new double[] { S1, S3, S4 };
            double[] s3 = new double[] { S1, S2, S4 };
            double[] sig = new double[] { sigma1, sigma2, sigma3, sigma4 };
            double[] sig1 = new double[] { sigma2, sigma3, sigma4 };
            double[] sig2 = new double[] { sigma1, sigma3, sigma4 };
            double[] sig3 = new double[] { sigma1, sigma2, sigma3 };
            double[,] rho = new double[,] { { 1.0, rho12, rho13, rho14 }, { rho12, 1.0, rho23, rho24 }, { rho13, rho23, 1.0, rho34 }, { rho14, rho24, rho34, 1.0 } };
            int lag = 1;
            double rho23_1 = projectedCorrelation(sig, rho, 2 - lag, 3 - lag, 1 - lag);
            double rho24_1 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 1 - lag);
            double rho34_1 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 1 - lag);

            double rho13_2 = projectedCorrelation(sig, rho, 1 - lag, 3 - lag, 2 - lag);
            double rho14_2 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 2 - lag);
            double rho34_2 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 2 - lag);

            double rho12_3 = projectedCorrelation(sig, rho, 1 - lag, 2 - lag, 3 - lag);
            double rho14_3 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 3 - lag);
            double rho24_3 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 3 - lag);

            double rho12_4 = projectedCorrelation(sig, rho, 1 - lag, 2 - lag, 4 - lag);
            double rho13_4 = projectedCorrelation(sig, rho, 1 - lag, 3 - lag, 4 - lag);
            double rho23_4 = projectedCorrelation(sig, rho, 2 - lag, 3 - lag, 4 - lag);

            double d21m = dmoinsij(s, T, sig, rho, 2 - lag, 1 - lag); double d31m = dmoinsij(s, T, sig, rho, 3 - lag, 1 - lag); double d41m = dmoinsij(s, T, sig, rho, 4 - lag, 1 - lag);
            double d12m = dmoinsij(s, T, sig, rho, 1 - lag, 2 - lag); double d32m = dmoinsij(s, T, sig, rho, 3 - lag, 2 - lag); double d42m = dmoinsij(s, T, sig, rho, 4 - lag, 2 - lag);
            double d13m = dmoinsij(s, T, sig, rho, 1 - lag, 3 - lag); double d23m = dmoinsij(s, T, sig, rho, 2 - lag, 3 - lag); double d43m = dmoinsij(s, T, sig, rho, 4 - lag, 3 - lag);
            double d14m = dmoinsij(s, T, sig, rho, 1 - lag, 4 - lag); double d24m = dmoinsij(s, T, sig, rho, 2 - lag, 4 - lag); double d34m = dmoinsij(s, T, sig, rho, 3 - lag, 4 - lag);
            double n1, n2, n3, n4;
            n1 = N3Dis2(-d21m, -d31m, -d41m, rho23_1, rho24_1, rho34_1, N3Flag, nb, start);
            n2 = N3Dis2(-d12m, -d32m, -d42m, rho13_2, rho14_2, rho34_2, N3Flag, nb, start);
            n3 = N3Dis2(-d13m, -d23m, -d43m, rho12_3, rho14_3, rho24_3, N3Flag, nb, start);
            n4 = N3Dis2(-d14m, -d24m, -d34m, rho12_4, rho13_4, rho23_4, N3Flag, nb, start);

            return S1 * n1 + S2 * n2 + S3 * n3 + S4 * n4;

        }

        // Max{Min{S1,S2,S3},S4}
        public static double MaxWorst3Underlyings(double S1, double S2, double S3, double S4, double T, double sigma1, double sigma2, double sigma3, double sigma4, double rho12, double rho13, double rho14, double rho23, double rho24, double rho34, int N3Flag, int nb, double start) // N3Flag=0 -> N3Dis usuelle, N3Flag=1 -> N3Dis2
        {
            double[] s = new double[] { S1, S2, S3, S4 };
            double[] s1 = new double[] { S2, S3, S4 };
            double[] s2 = new double[] { S1, S3, S4 };
            double[] s3 = new double[] { S1, S2, S4 };
            double[] sig = new double[] { sigma1, sigma2, sigma3, sigma4 };
            double[,] rho = new double[,] { { 1.0, rho12, rho13, rho14 }, { rho12, 1.0, rho23, rho24 }, { rho13, rho23, 1.0, rho34 }, { rho14, rho24, rho34, 1.0 } };

            int lag = 1;
            double rho23_1 = projectedCorrelation(sig, rho, 2 - lag, 3 - lag, 1 - lag);
            double rho24_1 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 1 - lag);
            double rho34_1 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 1 - lag);

            double rho13_2 = projectedCorrelation(sig, rho, 1 - lag, 3 - lag, 2 - lag);
            double rho14_2 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 2 - lag);
            double rho34_2 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 2 - lag);

            double rho12_3 = projectedCorrelation(sig, rho, 1 - lag, 2 - lag, 3 - lag);
            double rho14_3 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 3 - lag);
            double rho24_3 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 3 - lag);

            double rho12_4 = projectedCorrelation(sig, rho, 1 - lag, 2 - lag, 4 - lag);
            double rho13_4 = projectedCorrelation(sig, rho, 1 - lag, 3 - lag, 4 - lag);
            double rho23_4 = projectedCorrelation(sig, rho, 2 - lag, 3 - lag, 4 - lag);

            double d21m = dmoinsij(s, T, sig, rho, 2 - lag, 1 - lag); double d31m = dmoinsij(s, T, sig, rho, 3 - lag, 1 - lag); double d41m = dmoinsij(s, T, sig, rho, 4 - lag, 1 - lag);
            double d12m = dmoinsij(s, T, sig, rho, 1 - lag, 2 - lag); double d32m = dmoinsij(s, T, sig, rho, 3 - lag, 2 - lag); double d42m = dmoinsij(s, T, sig, rho, 4 - lag, 2 - lag);
            double d13m = dmoinsij(s, T, sig, rho, 1 - lag, 3 - lag); double d23m = dmoinsij(s, T, sig, rho, 2 - lag, 3 - lag); double d43m = dmoinsij(s, T, sig, rho, 4 - lag, 3 - lag);
            double d14m = dmoinsij(s, T, sig, rho, 1 - lag, 4 - lag); double d24m = dmoinsij(s, T, sig, rho, 2 - lag, 4 - lag); double d34m = dmoinsij(s, T, sig, rho, 3 - lag, 4 - lag);
            double n1, n2, n3, n4;

            n1 = N3Dis2(d21m, d31m, -d41m, rho23_1, -rho24_1, -rho34_1, N3Flag, nb, start);
            n2 = N3Dis2(-d12m, d32m, -d42m, rho13_2, -rho14_2, -rho34_2, N3Flag, nb, start);
            n3 = N3Dis2(d13m, d23m, -d43m, rho12_3, -rho14_3, -rho24_3, N3Flag, nb, start);
            n4 = (1.0 - N3Dis2(d14m, d24m, d34m, rho12_4, rho13_4, rho23_4, N3Flag, nb, start));

            return S1 * n1 + S2 * n2 + S3 * n3 + S4 * n4;

        }

        // Max{Min{S1,S2,S3}-K,0}
        public static double Worst3UnderlyingsCall(double S1, double S2, double S3, double K, double T, double sigma1, double sigma2, double sigma3, double rho12, double rho13, double rho23, int N3Flag, int nb, double start) // N3Flag=0 -> N3Dis usuelle, N3Flag=1 -> N3Dis2
        {
            double[] s = new double[] { S1, S2, S3, 0 };
            double[] s1 = new double[] { S2, S3, 0 };
            double[] s2 = new double[] { S1, S3, 0 };
            double[] s3 = new double[] { S1, S2, 0 };
            double[] sig = new double[] { sigma1, sigma2, sigma3, 0 };
            double[,] rho = new double[,] { { 1.0, rho12, rho13, 0 }, { rho12, 1.0, rho23, 0 }, { rho13, rho23, 1.0, 0 }, { 0, 0, 0, 1.0 } };

            int lag = 1;

            double rho23_1 = projectedCorrelation(sig, rho, 2 - lag, 3 - lag, 1 - lag);
            double rho24_1 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 1 - lag);
            double rho34_1 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 1 - lag);

            double rho13_2 = projectedCorrelation(sig, rho, 1 - lag, 3 - lag, 2 - lag);
            double rho14_2 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 2 - lag);
            double rho34_2 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 2 - lag);

            double rho12_3 = projectedCorrelation(sig, rho, 1 - lag, 2 - lag, 3 - lag);
            double rho14_3 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 3 - lag);
            double rho24_3 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 3 - lag);


            double d21m = dmoinsij(s, T, sig, rho, 2 - lag, 1 - lag); double d31m = dmoinsij(s, T, sig, rho, 3 - lag, 1 - lag);
            double d12m = dmoinsij(s, T, sig, rho, 1 - lag, 2 - lag); double d32m = dmoinsij(s, T, sig, rho, 3 - lag, 2 - lag);
            double d13m = dmoinsij(s, T, sig, rho, 1 - lag, 3 - lag); double d23m = dmoinsij(s, T, sig, rho, 2 - lag, 3 - lag);

            double d1p = dplusi(s, K, T, sig, 1 - lag);
            double d2p = dplusi(s, K, T, sig, 2 - lag);
            double d3p = dplusi(s, K, T, sig, 3 - lag);
            double d1m = dmoinsi(s, K, T, sig, 1 - lag);
            double d2m = dmoinsi(s, K, T, sig, 2 - lag);
            double d3m = dmoinsi(s, K, T, sig, 3 - lag);
            double n1, n2, n3, n4;
            n1 = N3Dis2(d21m, d31m, d1p, rho23_1, -rho24_1, -rho34_1, N3Flag, nb, start);
            n2 = N3Dis2(d12m, d32m, d2p, rho13_2, -rho14_2, -rho34_2, N3Flag, nb, start);
            n3 = N3Dis2(d13m, d23m, d3p, rho12_3, -rho14_3, -rho24_3, N3Flag, nb, start);
            n4 = N3Dis2(d1m, d2m, d3m, rho12, rho13, rho23, N3Flag, nb, start);
            double opt = S1 * n1 + S2 * n2 + S3 * n3 - K * n4;
            return opt;
        }

        public static double ForwardWorst3Underlyings(double S1, double S2, double S3, double T, double sigma1, double sigma2, double sigma3, double rho12, double rho13, double rho23, int N3Flag, int nb, double start) // N3Flag=0 -> N3Dis usuelle, N3Flag=1 -> N3Dis2
        {
            double[] s = new double[] { S1, S2, S3, 0 };
            double[] s1 = new double[] { S2, S3, 0 };
            double[] s2 = new double[] { S1, S3, 0 };
            double[] s3 = new double[] { S1, S2, 0 };
            double[] sig = new double[] { sigma1, sigma2, sigma3, 0 };
            double[] sig1 = new double[] { sigma2, sigma3, 0 };
            double[] sig2 = new double[] { sigma1, sigma3, 0 };
            double[] sig3 = new double[] { sigma1, sigma2, sigma3 };
            double[,] rho = new double[,] { { 1.0, rho12, rho13, 0 }, { rho12, 1.0, rho23, 0 }, { rho13, rho23, 1.0, 0 }, { 0, 0, 0, 1.0 } };
            int lag = 1;
            double rho23_1 = projectedCorrelation(sig, rho, 2 - lag, 3 - lag, 1 - lag);
            double rho24_1 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 1 - lag);
            double rho34_1 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 1 - lag);
            double rho13_2 = projectedCorrelation(sig, rho, 1 - lag, 3 - lag, 2 - lag);
            double rho14_2 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 2 - lag);
            double rho34_2 = projectedCorrelation(sig, rho, 3 - lag, 4 - lag, 2 - lag);
            double rho12_3 = projectedCorrelation(sig, rho, 1 - lag, 2 - lag, 3 - lag);
            double rho14_3 = projectedCorrelation(sig, rho, 1 - lag, 4 - lag, 3 - lag);
            double rho24_3 = projectedCorrelation(sig, rho, 2 - lag, 4 - lag, 3 - lag);
            double d21m = dmoinsij(s, T, sig, rho, 2 - lag, 1 - lag);
            double d31m = dmoinsij(s, T, sig, rho, 3 - lag, 1 - lag);
            double d12m = dmoinsij(s, T, sig, rho, 1 - lag, 2 - lag);
            double d32m = dmoinsij(s, T, sig, rho, 3 - lag, 2 - lag);
            double d13m = dmoinsij(s, T, sig, rho, 1 - lag, 3 - lag);
            double d23m = dmoinsij(s, T, sig, rho, 2 - lag, 3 - lag);
            double epsilon = 0.0001;
            double d1p = dplusi(s, epsilon, T, sig, 1 - lag);
            double d2p = dplusi(s, epsilon, T, sig, 2 - lag);
            double d3p = dplusi(s, epsilon, T, sig, 3 - lag);
            double d1m = dmoinsi(s, epsilon, T, sig, 1 - lag);
            double d2m = dmoinsi(s, epsilon, T, sig, 2 - lag);
            double d3m = dmoinsi(s, epsilon, T, sig, 3 - lag);
            double n1, n2, n3;
            n1 = N3Dis2(d21m, d31m, d1p, rho23_1, -rho24_1, -rho34_1, N3Flag, nb, start);
            n2 = N3Dis2(d12m, d32m, d2p, rho13_2, -rho14_2, -rho34_2, N3Flag, nb, start);
            n3 = N3Dis2(d13m, d23m, d3p, rho12_3, -rho14_3, -rho24_3, N3Flag, nb, start);
            return S1 * n1 + S2 * n2 + S3 * n3;
        }

        public static double Worst3Forward(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, int N3Flag, int nb, double start)
        {
            double ForwardCF;
            ForwardCF = ForwardWorst3Underlyings(Svect[0] * Math.Exp(muvect[0] * T), Svect[1] * Math.Exp(muvect[1] * T), Svect[2] * Math.Exp(muvect[2] * T),
                     T, sigmavect[0], sigmavect[1], sigmavect[2], Rhomatrix[0, 1], Rhomatrix[0, 2], Rhomatrix[1, 2], N3Flag, nb, start);
            return ForwardCF;
        }
        public static double[] Worst3ForwardBatch(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double[,] basketDescriptor, int N3Flag, int nb, double start)
        {
            int nbBaskets = basketDescriptor.GetLength(0);
            int nbUnder = basketDescriptor.GetLength(1);
            double[] ForwardCF = new double[nbBaskets];
            int[] firstLiveIndex = new int[nbBaskets];
            int[] secondLiveIndex = new int[nbBaskets];
            int[] thirdLiveIndex = new int[nbBaskets];
            for (int k = 0; k < nbBaskets; k++)    // recherche du premier stock du panier
            {
                firstLiveIndex[k] = nbUnder - 1;
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < firstLiveIndex[k])) firstLiveIndex[k] = iunder; }
            }
            for (int k = 0; k < nbBaskets; k++)    // recherche du second stock du panier
            {
                secondLiveIndex[k] = nbUnder - 1;
                for (int iunder = firstLiveIndex[k] + 1; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < secondLiveIndex[k])) secondLiveIndex[k] = iunder; }
            }
            for (int k = 0; k < nbBaskets; k++)    // recherche du troisieme stock du panier
            {
                thirdLiveIndex[k] = nbUnder - 1;
                for (int iunder = secondLiveIndex[k] + 1; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < thirdLiveIndex[k])) thirdLiveIndex[k] = iunder; }
            }
            for (int ipan = 0; ipan < nbBaskets; ipan++)
            {
                int i = firstLiveIndex[ipan];
                int j = secondLiveIndex[ipan];
                int k = thirdLiveIndex[ipan];
                ForwardCF[ipan] = ForwardWorst3Underlyings(Svect[i] * Math.Exp(muvect[i] * T), Svect[j] * Math.Exp(muvect[j] * T), Svect[k] * Math.Exp(muvect[k] * T),
                     T, sigmavect[i], sigmavect[j], sigmavect[k], Rhomatrix[i, j], Rhomatrix[i, k], Rhomatrix[j, k], N3Flag, nb, start);

            }
            return ForwardCF;
        }

        public static double Worst3Call(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double strike, int N3Flag, int nb, double start)
        {
            double CallCF;
            CallCF = Worst3UnderlyingsCall(Svect[0] * Math.Exp(muvect[0] * T), Svect[1] * Math.Exp(muvect[1] * T), Svect[2] * Math.Exp(muvect[2] * T),
                    strike, T, sigmavect[0], sigmavect[1], sigmavect[2], Rhomatrix[0, 1], Rhomatrix[0, 2], Rhomatrix[1, 2], N3Flag, nb, start);
            return CallCF;
        }

        public static double[] Worst3CallBatch(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double[] strike, double[,] basketDescriptor, int N3Flag, int nb, double start)
        {
            int nbBaskets = basketDescriptor.GetLength(0);
            int nbUnder = basketDescriptor.GetLength(1);
            double[] CallCF = new double[nbBaskets];
            int[] firstLiveIndex = new int[nbBaskets];
            int[] secondLiveIndex = new int[nbBaskets];
            int[] thirdLiveIndex = new int[nbBaskets];
            for (int k = 0; k < nbBaskets; k++)    // recherche du premier stock du panier
            {
                firstLiveIndex[k] = nbUnder - 1;
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < firstLiveIndex[k])) firstLiveIndex[k] = iunder; }
            }
            for (int k = 0; k < nbBaskets; k++)    // recherche du second stock du panier
            {
                secondLiveIndex[k] = nbUnder - 1;
                for (int iunder = firstLiveIndex[k] + 1; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < secondLiveIndex[k])) secondLiveIndex[k] = iunder; }
            }
            for (int k = 0; k < nbBaskets; k++)    // recherche du troisieme stock du panier
            {
                thirdLiveIndex[k] = nbUnder - 1;
                for (int iunder = secondLiveIndex[k] + 1; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < thirdLiveIndex[k])) thirdLiveIndex[k] = iunder; }
            }
            for (int ipan = 0; ipan < nbBaskets; ipan++)
            {
                int i = firstLiveIndex[ipan];
                int j = secondLiveIndex[ipan];
                int k = thirdLiveIndex[ipan];
                CallCF[ipan] = Worst3UnderlyingsCall(Svect[i] * Math.Exp(muvect[i] * T), Svect[j] * Math.Exp(muvect[j] * T), Svect[k] * Math.Exp(muvect[k] * T),
                    strike[ipan], T, sigmavect[i], sigmavect[j], sigmavect[k], Rhomatrix[i, j], Rhomatrix[i, k], Rhomatrix[j, k], N3Flag, nb, start);

            }
            return CallCF;
        }
    }


}

namespace pricing
{

    public class Export_ClosedForms
    {
        //  indispensable
        public static double Export_ImpliedCallVolatility(double forward, double strike, double maturity, double target_price)
        {
            return auxiliary.BlackScholes.ImpVolBS(forward, strike, maturity, target_price);
        }

        public static double Export_Worst3CFRoughCall(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double strike, int N3Flag = 0, int nb = 35, double start = -7)
        {
            return auxiliary.WorstOf3.Worst3Call(Svect, muvect, sigmavect, Rhomatrix, T, strike, N3Flag, nb, start);
        }

        public static double Export_Worst3CFRoughForward(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, int N3Flag = 0, int nb = 35, double start = -7)
        {
            return auxiliary.WorstOf3.Worst3Forward(Svect, muvect, sigmavect, Rhomatrix, T, N3Flag, nb, start);
        }

        public static double Export_BasketAvgBSCall(double[] Svect, double[] sigmavect, double[,] Rhomatrix, double T, double strike)
        {
            double f = (Svect[1] + Svect[2] + Svect[0]) / 3;
            double sigma = Math.Sqrt(sigmavect[0] * sigmavect[0] + sigmavect[1] * sigmavect[1] + sigmavect[2] * sigmavect[2] +
                   2 * (Rhomatrix[0, 1] * sigmavect[1] * sigmavect[0] + Rhomatrix[0, 2] * sigmavect[0] * sigmavect[2] + Rhomatrix[1, 2] * sigmavect[1] * sigmavect[2]));

            return auxiliary.BlackScholes.Call(f, strike, T, sigma);
        }


    }

    public static class MonteCarlo
    {
        public static Matrix<double> ConvertToMatrix(double[,] matrix)
        {
            int dim = matrix.GetLength(0);
            var myMat = CreateMatrix.Dense<double>(dim, dim, 0.0);
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    myMat[i, j] = matrix[i, j];
                }
            }
            return myMat;
        }


        public static double[] EigenValues(double[,] rhomatrix)
        {
            int nbUnder = rhomatrix.GetLength(0);
            Matrix<double> CovMatrix = ConvertToMatrix(rhomatrix);
            MathNet.Numerics.LinearAlgebra.Factorization.Svd<double> SvdmatrixA = CovMatrix.Svd();
            Matrix<double> U = SvdmatrixA.U;
            Vector<double> S = SvdmatrixA.S;
            double[] res = new double[nbUnder];
            for (int i = 0; i < nbUnder; i++)
            {
                res[i] = S[i];
            }
            return res;
        }

        public static double[,] GenerateTrajectoire(Random random, double[] Svect, double[] muvect, double[] sigmavect, double[,] rhomatrix, double Tmax, int nbdate)
        {
            int nbUnder = Svect.Length;
            double[,] Trajec = new double[nbdate, nbUnder];
            for (int idate = 0; idate < nbdate; idate++)
            {
                for (int iunder = 0; iunder < nbUnder; iunder++) { Trajec[idate, iunder] = 0.0; }
            }
            double deltaT = Tmax / (nbdate - 1.0);
            for (int iunder = 0; iunder < nbUnder; iunder++) { Trajec[0, iunder] = Svect[iunder]; }
            double[] Muvect = new double[nbUnder];
            for (int iunder = 0; iunder < nbUnder; iunder++) { Muvect[iunder] = (muvect[iunder] - sigmavect[iunder] * sigmavect[iunder] / 2.0) * deltaT; }
            double[,] Covmatrix = new double[nbUnder, nbUnder];
            for (int i1 = 0; i1 < nbUnder; i1++)
            {
                for (int i2 = 0; i2 < nbUnder; i2++) { Covmatrix[i1, i2] = sigmavect[i1] * sigmavect[i2] * rhomatrix[i1, i2] * deltaT; }
            }
            Matrix<double> CovMatrix = ConvertToMatrix(Covmatrix);
            MathNet.Numerics.LinearAlgebra.Factorization.Svd<double> SvdmatrixA = CovMatrix.Svd();
            Matrix<double> U = SvdmatrixA.U;
            Vector<double> S = SvdmatrixA.S;
            Matrix<double> VT = SvdmatrixA.VT;
            Vector<double> S1 = S.Clone();
            double epsilon = 0.00000000001; // regularisation des valeur singuliere de la matrice de covariance
            for (int i = 0; i < nbUnder; i++)
            {
                S1[i] = Math.Sqrt(Math.Max(S1[i], epsilon));
            }

            // Matrix<double> Choleskymatrix = MathNet.Numerics.LinearAlgebra.CreateMatrix.DiagonalOfDiagonalVector<double>(S1);
            // Choleskymatrix.Multiply(VT);
            Matrix<double> DS = MathNet.Numerics.LinearAlgebra.CreateMatrix.DiagonalOfDiagonalVector<double>(S1);
            Matrix<double> Choleskymatrix = U;
            Choleskymatrix = Choleskymatrix.Multiply(DS);

            //MathNet.Numerics.LinearAlgebra.Factorization.Cholesky<double> CholeskymatrixA = NewCovMatrix.Cholesky();
            //Matrix<double> Choleskymatrix = CholeskymatrixA.Factor;
            double[] randnoisevect = new double[nbUnder];
            Vector<double> randnoisevect1 = CreateVector.Dense<double>(randnoisevect);
            Vector<double> randnoisevect2 = CreateVector.Dense<double>(nbUnder);
            //  Random random = new MersenneTwister(42); // seed 42
            for (int idate = 1; idate < nbdate; idate++)
            {
                Normal.Samples(random, randnoisevect, 0.0, 1.0);
                randnoisevect1.SetValues(randnoisevect);
                // Choleskymatrix.Transpose().Multiply(randnoisevect1, randnoisevect2);
                Choleskymatrix.Multiply(randnoisevect1, randnoisevect2);
                for (int iunder = 0; iunder < nbUnder; iunder++)
                {
                    double lambda = randnoisevect2[iunder] + Muvect[iunder];
                    double x = Trajec[idate - 1, iunder];
                    Trajec[idate, iunder] = Trajec[idate - 1, iunder] * Math.Exp(lambda);
                }
            }
            return Trajec;
        }

        public static double[] GenerateTrajectoire1dim(Random random, double S, double mu, double sigma, double Tmax, int nbdate)
        {

            double[] Trajec = new double[nbdate];
            for (int idate = 0; idate < nbdate; idate++)
            {
                Trajec[idate] = 0.0;
            }
            double deltaT = Tmax / (nbdate - 1.0);
            Trajec[0] = S;
            double Mu;
            Mu = (mu - sigma * sigma / 2.0) * deltaT;

            //  Random random = new MersenneTwister(42); // seed 42
            for (int idate = 1; idate < nbdate; idate++)
            {
                double randnoise = Normal.Sample(random, 0.0, 1.0);
                double randnoise2 = sigma * randnoise * Math.Sqrt(deltaT) + Mu;
                Trajec[idate] = Trajec[idate - 1] * Math.Exp(randnoise2);
            }
            return Trajec;
        }

        public static double[,] GenerateTrajectoireVolSto1dim(Random random, double F0, double alpha0,
                double beta, double beta2, double d, double gamma, double nu, double lambda, double omega, double rho, double Tmax, int nbdate, int milsteinflag, int expmartingaleflag)
        {
            double[,] TrajecLog = new double[nbdate, 2];
            for (int idate = 0; idate < nbdate; idate++)
            {
                TrajecLog[idate, 0] = 0.0;
            }
            double deltaT = Tmax / (nbdate - 1.0);
            TrajecLog[0, 0] = Math.Log(F0);
            if (expmartingaleflag > 0)
            {
                TrajecLog[0, 1] = Math.Log(alpha0);
            }
            else
            {
                TrajecLog[0, 1] = alpha0;
            }
            //  Random random = new MersenneTwister(42); // seed 42
            double fbeta1;
            double deltaf = F0 * beta2 * Math.Exp(-d) / beta;
            double M0 = Math.Pow(F0, beta) / Math.Pow(deltaf, beta2);
            double M1 = F0 - deltaf;
            for (int idate = 1; idate < nbdate; idate++)
            {
                double f = TrajecLog[idate - 1, 0];
                double sigma = TrajecLog[idate - 1, 1];
                double Frandnoise = Normal.Sample(random, 0.0, 1.0);
                double Vrandnoise = Frandnoise * rho + Normal.Sample(random, 0.0, 1.0) * Math.Sqrt(1.0 - rho * rho);
                if (expmartingaleflag > 0)
                {
                    if (f <= F0)
                    {
                        fbeta1 = Math.Exp(f * (beta - 1) + sigma);
                    }
                    else
                    {
                        fbeta1 = M0 * (Math.Exp(f) - M1) * Math.Exp(sigma - f);
                    }

                }
                else
                {
                    fbeta1 = Math.Exp(f * (beta - 1)) * sigma;
                }
                //Euler shema
                double df = fbeta1 * (Frandnoise * Math.Sqrt(deltaT) - fbeta1 * deltaT / 2.0);
                // Milstein correction
                if (milsteinflag > 0)
                {
                    df += (beta - 1.0) * fbeta1 * fbeta1 * Math.Exp(sigma) * (Frandnoise * Frandnoise - 1) * deltaT / 2;
                }
                TrajecLog[idate, 0] = TrajecLog[idate - 1, 0] + df;
                if (expmartingaleflag > 0)
                {
                    double sigmagamma1 = Math.Exp(sigma * (gamma - 1)) * nu;
                    double dsigma = sigmagamma1 * (Vrandnoise * Math.Sqrt(deltaT) - sigmagamma1 * deltaT / 2.0) + (omega * Math.Exp(-sigma) - 1.0) * lambda * deltaT;
                    // Milstein correction
                    if (milsteinflag > 0)
                    {
                        dsigma += (gamma - 1.0) * sigmagamma1 * sigmagamma1 * nu * (Vrandnoise * Vrandnoise - 1) * deltaT / 2;
                    }
                    TrajecLog[idate, 1] = TrajecLog[idate - 1, 1] + dsigma;
                }
                else
                {
                    double dAlpha = sigma * nu * Math.Sqrt(deltaT) * Vrandnoise;
                    TrajecLog[idate, 1] = TrajecLog[idate - 1, 1] + dAlpha;
                }


            }
            return TrajecLog; // une trajectoire en Log  pour S et V
        }

        public static double CallVanilleMC(double S, double mu, double sigma, double T, double strike, int nbTraj)
        {
            Random random = new MersenneTwister(42); // seed 42
            double sum = 0.0; double sumSQ = 0.0;
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[] traj = GenerateTrajectoire1dim(random, S, mu, sigma, T, 2);
                double underlying = traj[1];
                double tirage;
                if (underlying >= strike)
                {
                    tirage = underlying - strike;
                }
                else
                {
                    tirage = 0.0;
                }
                sum += tirage;
                sumSQ += tirage * tirage;
            }
            double CallMC = sum / nbTraj;
            double varianceMC = (sumSQ - CallMC * CallMC * nbTraj) / (nbTraj - 1.0);
            return CallMC;
        }


        public static double CallVanilleVolStoMC(double F0, double alpha0,
                double beta, double beta2, double d, double gamma, double nu, double lambda, double omega, double rho,
                double Tmax, int nbdate, double strike, int nbTraj, int milsteinflag, int expmartingaleflag)
        {
            Random random = new MersenneTwister(42); // seed 42
            double sum = 0.0; double sumSQ = 0.0;
            for (int itraj = 0; itraj < nbTraj; itraj++)
            {
                double[,] traj = GenerateTrajectoireVolSto1dim(random, F0, alpha0,
                                beta, beta2, d, gamma, nu, lambda, omega, rho, Tmax, nbdate, milsteinflag, expmartingaleflag);
                double underlying = Math.Exp(traj[nbdate - 1, 0]);
                double tirage;
                if (underlying >= strike)
                {
                    tirage = underlying - strike;
                }
                else
                {
                    tirage = 0.0;
                }
                sum += tirage;
                sumSQ += tirage * tirage;
            }
            double CallMC = sum / nbTraj;
            double varianceMC = (sumSQ - CallMC * CallMC * nbTraj) / (nbTraj - 1.0);
            return CallMC;
        }

        public static double[] CallVanilleSetVolStoMC(double F0, double alpha0,
               double beta, double beta2, double d, double gamma, double nu, double lambda, double omega, double rho,
               double Tmax, int nbdate, double[] strikes, int nbTraj, int milsteinflag, int expmartingaleflag)
        {
            Random random = new MersenneTwister(42); // seed 42          
            int nbstrikes = strikes.Length;
            double[] tirages = new double[nbstrikes];
            double[] sum = new double[nbstrikes];
            for (int k = 0; k < nbstrikes; k++)
            {
                sum[k] = 0.0;
            }
            for (int itraj = 0; itraj < nbTraj; itraj++)
            {
                double[,] traj = GenerateTrajectoireVolSto1dim(random, F0, alpha0,
                                beta, beta2, d, gamma, nu, lambda, omega, rho, Tmax, nbdate, milsteinflag, expmartingaleflag);
                double underlying = Math.Exp(traj[nbdate - 1, 0]);

                for (int k = 0; k < nbstrikes; k++)
                {
                    if (underlying >= strikes[k])
                    {
                        tirages[k] = underlying - strikes[k];
                    }
                    else
                    {
                        tirages[k] = 0.0;
                    }
                    sum[k] += tirages[k];
                }

            }
            double[] CallMC = new double[nbstrikes];
            for (int k = 0; k < nbstrikes; k++)
            {
                CallMC[k] = sum[k] / nbTraj;
            }
           ;
            return CallMC;
        }

        public static double[,] VolVanilleSetVolStoMC(double F0, double alpha0,
               double beta, double beta2, double d, double gamma, double nu, double lambda, double omega, double rho,
               double Tmax, int nbdate, double[] strikes, int nbTraj, int milsteinflag, int expmartingaleflag)
        {
            int nbstrikes = strikes.Length;
            double[,] Vols = new double[2, nbstrikes];
            double[] opt = CallVanilleSetVolStoMC(F0, alpha0, beta, beta2, d, gamma, nu, lambda, omega, rho, Tmax, nbdate, strikes, nbTraj, milsteinflag, expmartingaleflag);
            for (int k = 0; k < nbstrikes; k++)
            {
                double timevalue = opt[k] - F0 + strikes[k];
                if (timevalue <= 0)
                {
                    Vols[1, k] = 0.0;
                }

                else
                {
                    if (opt[k] >= F0)
                    {
                        Vols[1, k] = 0.0;
                    }
                    else
                    {
                        try
                        {
                            Vols[1, k] = auxiliary.BlackScholes.ImpVolBS(F0, strikes[k], Tmax, opt[k]);
                        }
                        catch (Exception)
                        {
                            Vols[1, k] = 0.0;
                        }
                    }

                }
                Vols[0, k] = opt[k];


            }
            return Vols; // Vols;

        }

        public static double CallWorstMC(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double strike, int nbTraj)
        {
            Random random = new MersenneTwister(42); // seed 42
            double sum = 0.0; double sumSQ = 0.0; int nbUnder = Svect.Length;
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] traj = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, 2); ;
                double underlying = traj[1, 0];
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                {
                    underlying = Math.Min(underlying, traj[1, iunder]);
                }
                double tirage;
                if (underlying >= strike)
                {
                    tirage = underlying - strike;
                }
                else
                {
                    tirage = 0.0;
                }
                sum += tirage;
                sumSQ += tirage * tirage;
            }
            double CallMC = sum / nbTraj;
            double varianceMC = (sumSQ - CallMC * CallMC * nbTraj) / (nbTraj - 1.0);
            return CallMC;
        }

        public static void CallAndForwardWorstMC(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix,
            double T, double strike, int nbTraj, ref double forward, ref double option)
        {
            Random random = new MersenneTwister(42); // seed 42
            double sum = 0.0; double sumU = 0.0; int nbUnder = Svect.Length;
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] traj = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, 2); ;
                double underlying = traj[1, 0];
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                {
                    underlying = Math.Min(underlying, traj[1, iunder]);
                }
                double tirage;
                if (underlying >= strike)
                {
                    tirage = underlying - strike;
                }
                else
                {
                    tirage = 0.0;
                }
                sum += tirage;
                sumU += underlying;


            }
            option = sum / nbTraj;
            forward = sumU / nbTraj;

            return;
        }

        public static int[,] generate_table2power(int n)
        {
            int dimtotale = (int)Math.Pow(2, n);
            int[] dimpartielle = new int[n];
            for (int i = 0; i < n; i++) dimpartielle[i] = (int)Math.Pow(2, i);
            int[,] res = new int[dimtotale, n];
            for (int col = 0; col < n; col++)
            {
                int basis = dimpartielle[col]; int count = 1; int state = 0;
                for (int row = 0; row < dimtotale; row++)
                {
                    res[row, col] = state;
                    count++;
                    if (count > basis)
                    {
                        count = 1;
                        state = 1 - state;
                    }

                }
            }
            return res;
        }

        public static double[,] generate_tableCnp(int n, int p)
        {
            int dimreelle = 1;
            for (int i = 0; i < p; i++)
            {
                dimreelle = dimreelle * (n - i) / (1 + i);
            }
            double[,] res = new double[dimreelle, n];
            int[,] res_aux = generate_table2power(n);
            int sumdigits; int ires = 0;
            for (int i = 0; i < res_aux.GetLength(0); i++)
            {
                sumdigits = 0; for (int j = 0; j < n; j++) sumdigits += res_aux[i, j];
                if (sumdigits == p)
                {
                    for (int j = 0; j < n; j++) res[ires, j] = res_aux[i, j];
                    ires++;
                }
            }
            return res;
        }

        public static double[] oneslice_generate_tableCnp(int n, int p, int j)
        {
            double[,] tab = generate_tableCnp(n, p);
            double[] res = new double[n];
            for (int i = 0; i < n; i++)
            {
                res[i] = tab[j, i];
            }
            return res;
        }


        public static double testSVD(double a)
        {
            double[,] Covmatrix = new double[2, 2];
            Covmatrix[0, 0] = 1;
            Covmatrix[0, 1] = 2;
            Covmatrix[1, 0] = 1;
            Covmatrix[1, 1] = 2;
            Matrix<double> CovMatrix = ConvertToMatrix(Covmatrix);
            MathNet.Numerics.LinearAlgebra.Factorization.Svd<double> SvdmatrixA = CovMatrix.Svd();
            Matrix<double> U = SvdmatrixA.U;
            Vector<double> S = SvdmatrixA.S;
            Matrix<double> VT = SvdmatrixA.VT;
            Vector<double> S1 = S.Clone();
            double epsilon = 0.00000000001; // regularisation des valeur singuliere de la matrice de covariance
            for (int i = 0; i < 2; i++)
            {
                S1[i] = Math.Sqrt(Math.Max(S1[i], epsilon));
            }

            // Matrix<double> Choleskymatrix = MathNet.Numerics.LinearAlgebra.CreateMatrix.DiagonalOfDiagonalVector<double>(S1);
            // Choleskymatrix.Multiply(VT);
            Matrix<double> DS = MathNet.Numerics.LinearAlgebra.CreateMatrix.DiagonalOfDiagonalVector<double>(S1);
            Matrix<double> Choleskymatrix = U;
            Choleskymatrix.Multiply(DS);
            return a + 1;
        }

        public static double CallSumMC(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double strike, int nbTraj)
        {
            Random random = new MersenneTwister(42); // seed 42
            double sum = 0.0; double sumSQ = 0.0; int nbUnder = Svect.Length;
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] traj = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, 2); ;
                double underlying = 0;
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                {
                    underlying += traj[1, iunder];
                }
                underlying /= nbUnder;
                double tirage;
                if (underlying >= strike)
                {
                    tirage = underlying - strike;
                }
                else
                {
                    tirage = 0.0;
                }
                sum += tirage;
                sumSQ += tirage * tirage;
            }
            double CallMC = sum / nbTraj;
            double varianceMC = (sumSQ - CallMC * CallMC * nbTraj) / (nbTraj - 1.0);
            return CallMC;
        }

        public static double ImpliedVolCallSumMC(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double strike, int nbTraj)
        {
            double option = CallSumMC(Svect, muvect, sigmavect, Rhomatrix, T, strike, nbTraj);
            double forward = 0.0;
            for (int i = 0; i < Svect.Length; i++) forward += Svect[i];
            forward /= Svect.Length;
            double vol = auxiliary.BlackScholes.ImpVolBS(forward, strike, T, option);

            return vol;
        }

        public static double[] CallWorstMCForwardBatch(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, int nbTraj, double[,] basketDescriptor)
        {
            Random random = new MersenneTwister(42); // seed 42
            int nbBaskets = basketDescriptor.GetLength(0);
            int nbUnder = basketDescriptor.GetLength(1);
            double[] payoff = new double[nbBaskets]; double[] PDIFlag = new double[nbBaskets];
            double[] worst = new double[nbBaskets];
            int[] firstLiveIndex = new int[nbBaskets];
            double[] CallMC = new double[nbBaskets];
            for (int k = 0; k < nbBaskets; k++) { PDIFlag[k] = 0; payoff[k] = 0; worst[k] = 0.0; }
            for (int k = 0; k < nbBaskets; k++)    // recherche du premier stock du panier
            {
                firstLiveIndex[k] = nbUnder - 1;
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < firstLiveIndex[k])) firstLiveIndex[k] = iunder; }
                CallMC[k] = 0.0;
            }


            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] traj = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, 2); ;
                for (int k = 0; k < nbBaskets; k++)
                {
                    worst[k] = traj[1, firstLiveIndex[k]];
                }
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                {
                    for (int k = 0; k < nbBaskets; k++)
                    { if (basketDescriptor[k, iunder] > 0) worst[k] = Math.Min(worst[k], traj[1, iunder]); }
                }
                for (int k = 0; k < nbBaskets; k++)
                {
                    CallMC[k] += worst[k];
                }
            }
            for (int k = 0; k < nbBaskets; k++) { CallMC[k] /= nbTraj; }
            return CallMC;
        }

        public static double[] CallWorstMCBatch(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double[] strike, int nbTraj, double[,] basketDescriptor)
        {
            Random random = new MersenneTwister(42); // seed 42
            int nbBaskets = basketDescriptor.GetLength(0);
            int nbUnder = basketDescriptor.GetLength(1);
            double[] payoff = new double[nbBaskets]; double[] PDIFlag = new double[nbBaskets];
            double[] worst = new double[nbBaskets];
            int[] firstLiveIndex = new int[nbBaskets];
            for (int k = 0; k < nbBaskets; k++) { PDIFlag[k] = 0; payoff[k] = 0; worst[k] = 0.0; }
            for (int k = 0; k < nbBaskets; k++)    // recherche du premier stock du panier
            {
                firstLiveIndex[k] = nbUnder - 1;
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                { if ((basketDescriptor[k, iunder] > 0) && (iunder < firstLiveIndex[k])) firstLiveIndex[k] = iunder; }
            }

            double[] CallMC = new double[nbBaskets];
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] traj = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, 2); ;
                for (int k = 0; k < nbBaskets; k++) { worst[k] = traj[1, firstLiveIndex[k]]; }
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                {
                    for (int k = 0; k < nbBaskets; k++)
                    { if (basketDescriptor[k, iunder] > 0) worst[k] = Math.Min(worst[k], traj[1, iunder]); }
                }
                double[] tirage = new double[nbBaskets];
                for (int k = 0; k < nbBaskets; k++)
                {
                    if (worst[k] >= strike[k])
                    { tirage[k] = worst[k] - strike[k]; }
                    else
                    { tirage[k] = 0.0; }
                    CallMC[k] += tirage[k];
                }
            }
            for (int k = 0; k < nbBaskets; k++) { CallMC[k] /= nbTraj; }
            return CallMC;
        }

        public static double YetiPhenix3DimWorstPayoff(double[,] trajec,
            double Bonus, double YetiBarrier, double YetiCoupon, double PhoenixBarrier, double PhoenixCoupon,
            double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, out double[] forwards)
        {
            double payoff = 0.0; double PDIFlag = 0;
            int nbunder = trajec.GetLength(1);
            int nbdate = trajec.GetLength(0);
            forwards = new double[nbdate]; ;
            for (int idate = 1; idate < nbdate; idate++)
            {
                double w = trajec[idate, 0];
                for (int iunder = 0; iunder < nbunder; ++iunder)
                {
                    w = Math.Min(w, trajec[idate, iunder]);
                }
                forwards[idate] = w;
            }
            for (int idate = 1; idate < nbdate; idate++)
            {
                if (forwards[idate] <= YetiBarrier)
                {
                    payoff += Bonus + YetiCoupon;
                }
                else
                {
                    payoff += Bonus;
                }
                if (forwards[idate] >= PhoenixBarrier)
                {
                    payoff += PhoenixCoupon; break;
                }
                if (forwards[idate] <= PDIbarrier)
                {
                    PDIFlag = 1.0;
                }

            }
            return payoff + PDIFlag * PDIGearing * Math.Max(0.0, PDItype * (forwards[nbdate - 1] - PDIStrike));
        }

        public static double LookbackCall3DimWorstPayoff(double[,] trajec,
           double PDIStrike, double PDItype, out double[] forwards)
        {
            int nbunder = trajec.GetLength(1);
            int nbdate = trajec.GetLength(0);
            forwards = new double[nbdate];
            double w = trajec[0, 0];
            for (int iunder = 1; iunder < nbunder; ++iunder)
            {
                w = Math.Min(w, trajec[0, iunder]);
            }
            forwards[0] = w;

            for (int idate = 1; idate < nbdate; idate++)
            {
                w = trajec[idate, 0];
                for (int iunder = 0; iunder < nbunder; ++iunder)
                {
                    w = Math.Min(w, trajec[idate, iunder]);
                }
                forwards[idate] = Math.Min(w, forwards[idate - 1]);
            }
            if (PDItype >= 0)
            {
                return Math.Max(0.0, (forwards[nbdate - 1] - PDIStrike));
            }
            else
            {
                return Math.Max(0.0, -(forwards[nbdate - 1] - PDIStrike));
            }

        }


        public static double[] YetiPhenix3DimWorstPayoffBatch(double[,] trajec, double PhoenixBarrier, double PhoenixCoupon,
            double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, double[,] basketDescriptor)
        {
            int nbBaskets = basketDescriptor.GetLength(0);
            int nbUnder = basketDescriptor.GetLength(1);
            double[] payoff = new double[nbBaskets]; double[] PDIFlag = new double[nbBaskets];
            int nbdate = trajec.GetLength(0);
            //  trajec.GetLength(1) doit etre egal  a basketDescriptor.GetLength(1)
            // 
            double[] worst = new double[nbBaskets];
            int[] firstLiveIndex = new int[nbBaskets];
            int[] actifFlag = new int[nbBaskets];
            for (int k = 0; k < nbBaskets; k++) { actifFlag[k] = 1; PDIFlag[k] = 0; payoff[k] = 0; worst[k] = 0.0; }
            // recherche du premier stock du panier
            for (int k = 0; k < nbBaskets; k++)
            {
                firstLiveIndex[k] = nbUnder - 1;
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                {
                    if ((basketDescriptor[k, iunder] > 0) && (iunder < firstLiveIndex[k])) firstLiveIndex[k] = iunder;
                }

            }
            for (int idate = 1; idate < nbdate; idate++)
            {
                for (int k = 0; k < nbBaskets; k++)
                {
                    worst[k] = trajec[idate, firstLiveIndex[k]];
                }
                for (int iunder = 0; iunder < nbUnder; ++iunder)
                {
                    for (int k = 0; k < nbBaskets; k++)
                    {
                        if (basketDescriptor[k, iunder] > 0) worst[k] = Math.Min(worst[k], trajec[idate, iunder]);
                    }
                }
                for (int k = 0; k < nbBaskets; k++)
                {
                    if ((worst[k] >= PhoenixBarrier) && (actifFlag[k] > 0))
                    {
                        payoff[k] = PhoenixCoupon; actifFlag[k] = 1;
                    }
                    if ((worst[k] <= PDIbarrier) && (actifFlag[k] > 0))
                    {
                        PDIFlag[k] = 1.0;
                    }
                }
            }
            for (int k = 0; k < nbBaskets; k++) { payoff[k] += PDIFlag[k] * PDIGearing * Math.Max(0.0, PDItype * (worst[k] - PDIStrike)); }
            return payoff;
        }

        public static double YetiPhenix3DimWorstPayoff1dim(double[] trajec,
           double Bonus, double YetiBarrier, double YetiCoupon, double PhoenixBarrier, double PhoenixCoupon,
           double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype)
        {
            double payoff = 0.0; double PDIFlag = 0;
            int nbdate = trajec.GetLength(0);
            double forward = 0.0;
            for (int idate = 1; idate < nbdate; idate++)
            {
                forward = trajec[idate];
                if (forward <= YetiBarrier)
                {
                    payoff += Bonus + YetiCoupon;
                }
                if (forward >= PhoenixBarrier)
                {
                    payoff = PhoenixCoupon; break;
                }
                if (forward <= PDIbarrier)
                {
                    PDIFlag = 1.0;
                }
            }
            return payoff + PDIFlag * PDIGearing * Math.Max(0.0, PDItype * (forward - PDIStrike));
        }


        public static double YetiPhoenixPriceMC(Random random, double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T,
            double Bonus, double YetiBarrier, double YetiCoupon, double PhoenixBarrier, double PhoenixCoupon,
            double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, int nbdate, int nbTraj, out double[] forwards)
        {

            double sum = 0.0; int nbUnder = Svect.Length;
            forwards = new double[nbdate];
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] trajec = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, nbdate + 1); ;
                double tirage = YetiPhenix3DimWorstPayoff(trajec, Bonus, YetiBarrier, YetiCoupon, PhoenixBarrier, PhoenixCoupon, PDIbarrier, 
                    PDIGearing, PDIStrike, PDItype, out double[] forwardsvec);

                sum += tirage;
                for (int j = 0; j < nbdate; j++)
                {
                    forwards[j] += forwardsvec[j];
                }
            }
            double CallMC = sum / nbTraj;
            for (int j = 0; j < nbdate; j++)
            {
                forwards[j] /= nbTraj;
            }

            return CallMC;
        }

        public static double LookbackCallWorstPriceMC(Random random, double[] Svect, double[] muvect,
            double[] sigmavect, double[,] Rhomatrix, double T, double PDIStrike, double PDIType,
             int nbdate, int nbTraj, out double[] forwards)
        {

            double sum = 0.0; int nbUnder = Svect.Length;
            forwards = new double[nbdate];
            double[] forwardsvec;
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] trajec = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, nbdate + 1);
                double tirage = LookbackCall3DimWorstPayoff(trajec, PDIStrike, PDIType, out forwardsvec);

                sum += tirage;
                for (int j = 0; j < nbdate; j++)
                {
                    forwards[j] += forwardsvec[j];
                }
            }
            double CallMC = sum / nbTraj;
            for (int j = 0; j < nbdate; j++)
            {
                forwards[j] /= nbTraj;
            }

            return CallMC;
        }


        public static double[] YetiPhoenixPriceMCBatch(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T, double PhoenixBarrier, double PhoenixCoupon,
           double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, int nbdate, int nbTraj, double[,] basketDescriptor)
        {
            Random random = new MersenneTwister(42); // seed 42
            int nbBasket = basketDescriptor.GetLength(0);
            double[] sum = new double[nbBasket]; double[] CallMC = new double[nbBasket]; int nbUnder = Svect.Length;
            for (int k = 0; k < nbBasket; k++)
            {
                sum[k] = 0.0;
            }
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[,] trajec = GenerateTrajectoire(random, Svect, muvect, sigmavect, Rhomatrix, T, nbdate + 1); ;
                double[] tirage = YetiPhenix3DimWorstPayoffBatch(trajec, PhoenixBarrier, PhoenixCoupon, PDIbarrier, PDIGearing, PDIStrike, PDItype, basketDescriptor);
                for (int k = 0; k < nbBasket; k++)
                {
                    sum[k] += tirage[k];
                }
            }
            for (int k = 0; k < nbBasket; k++)
            {
                CallMC[k] = sum[k] / nbTraj;
            }

            return CallMC;
        }


        public static double YetiPhoenixPriceMC_1dim(Random random, double S, double mu, double sigma, double T,
             double Bonus, double YetiBarrier, double YetiCoupon, double PhoenixBarrier, double PhoenixCoupon,
            double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, int nbdate, int nbTraj)
        {
            // Random random = new MersenneTwister(42); // seed 42
            double sum = 0.0; 
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[] trajec = GenerateTrajectoire1dim(random, S, mu, sigma, T, nbdate + 1); ;
                double tirage = YetiPhenix3DimWorstPayoff1dim(trajec, Bonus, YetiBarrier, YetiCoupon,PhoenixBarrier, PhoenixCoupon, PDIbarrier, PDIGearing, PDIStrike, PDItype);

                sum += tirage;
            
            }
            double CallMC = sum / nbTraj;
            return CallMC;
        }


        public static double[] WorstLocalVolMC(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix,
            double T, double deltaT, double[] Strikevect, int nbTraj)
        {
            double forward = CallWorstMC(Svect, muvect, sigmavect, Rhomatrix, T, 0, nbTraj);
            int nbunder = Svect.Length;
            int nbStrike = Strikevect.Length;
            double[] wVectm = new double[nbStrike]; double[] wVectp = new double[nbStrike];
            double[] wVect = new double[nbStrike]; double[] yVect = new double[nbStrike];
            double[] OptT = new double[nbStrike];
            double opt;
            for (int i = 0; i < nbStrike; i++)
            {
                opt = CallWorstMC(Svect, muvect, sigmavect, Rhomatrix, T + deltaT, Strikevect[i], nbTraj);
                OptT[i] = opt;
                wVectp[i] = auxiliary.BlackScholes.ImpVolBS(forward, Strikevect[i], T + deltaT, opt);
                opt = CallWorstMC(Svect, muvect, sigmavect, Rhomatrix, T - deltaT, Strikevect[i], nbTraj);
                wVectm[i] = auxiliary.BlackScholes.ImpVolBS(forward, Strikevect[i], T - deltaT, opt);
                wVect[i] = (wVectp[i] + wVectm[i]) / 2.0; wVect[i] = wVect[i] * wVect[i] * T;
                yVect[i] = Math.Log(Strikevect[i] / forward);
            }

            CubicSpline smileObject = CubicSpline.InterpolateAkima(yVect, wVect);
            double[] Deryw = new double[nbStrike]; double[] Dery2w = new double[nbStrike]; double[] DerTw = new double[nbStrike];
            for (int i = 0; i < nbStrike; i++)
            {
                Deryw[i] = smileObject.Differentiate(yVect[i]);
                Dery2w[i] = smileObject.Differentiate2(yVect[i]);
                DerTw[i] = (wVectm[i] * wVectm[i] - wVectp[i] * wVectp[i]) * T / (2.0 * deltaT);
            }
            double[] localvol = new double[nbStrike]; double[] yw = new double[nbStrike]; double[] result = new double[nbStrike];
            for (int i = 0; i < nbStrike; i++)
            {
                yw[i] = yVect[i] / wVect[i];
                localvol[i] = DerTw[i] / (1.0 - yw[i] * Deryw[i] + 0.25 * (-0.25 - 1.0 / wVect[i] + yw[i] * yw[i]) * (Deryw[i] * Deryw[i] + 0.5 * Dery2w[i]));
                if (localvol[i] < 0.0001)
                {
                    result[i] = 0.0001;
                }
                else
                {
                    result[i] = Math.Sqrt(localvol[i]);
                }

            }


            return result;

        }


        public static double[] GenerateSmiledTrajectoire1dim(Random random, double S, double mu, CubicSpline[] smileObject, double Tmax, int nbdate)
        {
            double[] Trajec = new double[nbdate];
            for (int idate = 0; idate < nbdate; idate++)
            {
                Trajec[idate] = 0.0;
            }
            double deltaT = Tmax / (nbdate - 1.0);
            Trajec[0] = S;
            double Mu;
            //  Random random = new MersenneTwister(42); // seed 42
            for (int idate = 1; idate < nbdate; idate++)
            {
                double sigma = smileObject[idate-1].Interpolate(Trajec[idate - 1]/S);
                Mu = (mu - sigma * sigma / 2.0) * deltaT;
                double randnoise = Normal.Sample(random, 0.0, 1.0);
                double randnoise2 = sigma * randnoise * Math.Sqrt(deltaT) + Mu;
                Trajec[idate] = Trajec[idate - 1] * Math.Exp(randnoise2);
            }
            return Trajec;
        }



        public static double YetiPhoenixSmiledPriceMC_1dim(Random random, double S, double mu, double[] StrikeList, int NbDate, double[,] localVol,
            double T, double Bonus, double YetiBarrier, double YetiCoupon,
            double PhoenixBarrier, double PhoenixCoupon,
            double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, int nbdate, int nbTraj)
        {
            int nbStrike = localVol.GetLength(1);
            int nbDate = localVol.GetLength(0);
            CubicSpline[] smileObject = new CubicSpline[nbDate];
            
            
            double[] smile = new double[nbStrike];
            for (int iDate = 0; iDate < nbDate; iDate++)
            {
                for (int iStrike = 0; iStrike < nbStrike; iStrike++)
                {
                    smile[iStrike] = localVol[iDate, iStrike];
                }
               
                smileObject[iDate] = CubicSpline.InterpolateAkima(StrikeList, smile);
            }
            

            double sum = 0.0;           
            for (int itraj = 1; itraj < nbTraj; itraj++)
            {
                double[] trajec = GenerateSmiledTrajectoire1dim(random, S, mu, smileObject, T, nbdate + 1); ;
                double tirage = YetiPhenix3DimWorstPayoff1dim(trajec, Bonus, YetiBarrier, YetiCoupon,
                    PhoenixBarrier, PhoenixCoupon, PDIbarrier, 
                    PDIGearing, PDIStrike, PDItype);
                sum += tirage;              
            }
            double CallMC = sum / nbTraj;         
            return CallMC;
        }


    }
    public class Export_MonteCarlo
    {

        public static double Export_YetiPhoenixSmiledPriceMC_1dim(double S, double mu, double[] StrikeList, double[,] localVol, double T,
            double Bonus, double YetiBarrier, double YetiCoupon, double PhoenixBarrier, double PhoenixCoupon,
            double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, int nbdate, int nbTraj)
        {
            Random random = new MersenneTwister(41); // seed 42
            return MonteCarlo.YetiPhoenixSmiledPriceMC_1dim(random,S, mu, StrikeList, nbdate, localVol, T,
                Bonus, YetiBarrier, YetiCoupon, PhoenixBarrier, PhoenixCoupon,
             PDIbarrier, PDIGearing, PDIStrike, PDItype, nbdate, nbTraj);

        }

        public static double Export_YetiPhoenixPriceMC(double[] Svect, double[] muvect, double[] sigmavect, double[,] Rhomatrix, double T,
            double Bonus, double YetiBarrier, double YetiCoupon, double PhoenixBarrier, double PhoenixCoupon,
            double PDIbarrier, double PDIGearing, double PDIStrike, double PDItype, int nbdate, int nbTraj, out double[] forwardVector)
        {
            Random random = new MersenneTwister(41); // seed 42
            return MonteCarlo.YetiPhoenixPriceMC(random, Svect, muvect, sigmavect, Rhomatrix, T,
                 Bonus, YetiBarrier, YetiCoupon, PhoenixBarrier, PhoenixCoupon,
             PDIbarrier, PDIGearing, PDIStrike, PDItype, nbdate, nbTraj, out forwardVector);

        }


        public static double[] Export_generate_tableCnp(int n, int p, int j)
        {
            return MonteCarlo.oneslice_generate_tableCnp(n, p, j);
        }
    }

    public class SABR
    {
        public static double MACHINE_PRECISION = 1e-16;

        public static double[,] vanilla_call(double maturity, double forward, double alpha, double beta, double nu,
             double rho, double[] strikes, int nbdate, int nb_tirages, int milsteinflag, int expmartingaleflag)
        {
            Random generator = new Random(2);
            double[,] result = new double[6, strikes.Length];

            double tmp, tmp3, deltaT;
            int nbstrikes = strikes.Length;
            deltaT = maturity / (nbdate - 1.0);
            Random random = new MersenneTwister(42); // seed 42
            for (int i = 0; i < nb_tirages; i++)
            {
                double[,] TrajecLog = new double[nbdate, 2];
                for (int idate = 0; idate < nbdate; idate++)
                {
                    TrajecLog[idate, 0] = 0.0;
                }
                TrajecLog[0, 0] = Math.Log(forward);
                TrajecLog[0, 1] = alpha;
                #region   trajectoire Rida
                for (int idate = 1; idate < nbdate; idate++)
                {
                    double f = TrajecLog[idate - 1, 0];
                    double sigma = TrajecLog[idate - 1, 1];
                    double Frandnoise = Normal.Sample(random, 0.0, 1.0);
                    double Vrandnoise = Frandnoise * rho + Normal.Sample(random, 0.0, 1.0) * Math.Sqrt(1.0 - rho * rho);
                    double fbeta1 = Math.Exp(f * (beta - 1)) * sigma;
                    double df = fbeta1 * (Frandnoise * Math.Sqrt(deltaT) - fbeta1 * deltaT / 2.0);
                    double dAlpha = sigma * nu * Math.Sqrt(deltaT) * Vrandnoise;
                    TrajecLog[idate, 0] = TrajecLog[idate - 1, 0] + df;
                    TrajecLog[idate, 1] = TrajecLog[idate - 1, 1] + dAlpha;
                }
                #endregion

                #region  trajectoire3
                //Random random, double F0, double alpha0, double beta, double beta2, double d, double gamma, double nu, double lambda, double omega, double rho, double Tmax, int nbdate, int milsteinflag, int expmartingaleflag
                double[,] TrajecLog3 = pricing.MonteCarlo.GenerateTrajectoireVolSto1dim(random, forward, alpha,
                        beta, beta, 0, 1, nu, 0, 0, rho, maturity, nbdate, milsteinflag, expmartingaleflag);
                #endregion
                double underlying = Math.Exp(TrajecLog[nbdate - 1, 0]);
                double underlying3 = Math.Exp(TrajecLog3[nbdate - 1, 0]);
                for (int j = 0; j < strikes.Length; j++)
                {
                    tmp = Math.Max(underlying - strikes[j], 0.0);
                    result[0, j] += tmp / nb_tirages;
                    result[1, j] += tmp * tmp / nb_tirages;
                    tmp3 = Math.Max(underlying3 - strikes[j], 0.0);
                    result[3, j] += tmp3 / nb_tirages;
                    result[4, j] += tmp3 * tmp3 / nb_tirages;
                }
            }
            //calcul des standard devs & implied vols
            for (int j = 0; j < strikes.Length; j++)
            {
                result[1, j] -= result[0, j] * result[0, j];
                result[1, j] = Math.Sqrt(result[1, j] / (nb_tirages - 1));
                result[2, j] = auxiliary.BlackScholes.ImpVolBS(forward, strikes[j], maturity, result[0, j]);
                result[4, j] -= result[3, j] * result[3, j];
                result[4, j] = Math.Sqrt(result[4, j] / (nb_tirages - 1));
                result[5, j] = auxiliary.BlackScholes.ImpVolBS(forward, strikes[j], maturity, result[3, j]);
            }

            return result;
        }


        public static double[] generate_strikes(double forward, double maturite, double vol_atm, int n_gauche, int n_droite, double delta_call = 0.006, double delta_put = -0.006, double epsilon = 0.03)
        {
            double[] resultat = new double[n_gauche + n_droite - 1];

            double stdev = vol_atm * Math.Sqrt(maturite);

            double x_min = 0.0;
            double x_max = 1.0 - epsilon;


            double k_max = forward * Math.Exp(0.5 * stdev * stdev - stdev * auxiliary.GaussianDensity.InvCdfSafe(delta_call));
            double k_min = forward * Math.Exp(0.5 * stdev * stdev - stdev * auxiliary.GaussianDensity.InvCdfSafe(delta_put + 1.0));

            for (int i = 0; i < n_gauche; i++)
            {
                resultat[i] = ((k_min - forward) / fisher_transform(x_max)) * fisher_transform(x_min + i * (x_max - x_min) / (n_gauche - 1)) + forward;
            }

            for (int i = 0; i < n_gauche; i++)
            {
                resultat[n_gauche - i - 1] = ((k_min - forward) / fisher_transform(x_max)) * fisher_transform(x_min + i * (x_max - x_min) / (n_gauche - 1)) + forward;
            }

            for (int i = 1; i < n_droite; i++)
            {
                resultat[n_gauche + i - 1] = ((k_max - forward) / fisher_transform(x_max)) * fisher_transform(x_min + i * (x_max - x_min) / (n_droite - 1)) + forward;
            }

            return resultat;
        }
        public static double fisher_transform(double x)
        {

            double result = 0.5 * Math.Log((1 + x) / (1 - x));

            return result;
        }


        public class SABRVolatilityProcess
        {

            double sigma0;
            double nu;

            public SABRVolatilityProcess(double sigma0, double nu)
            {
                this.sigma0 = sigma0;
                this.nu = nu;
            }
            public double get_drift(double x)
            {
                return 0.0;
            }

            public double get_voatility(double x)
            {
                return nu * x;
            }
        }

        public class Mrg32k3aGiles
        {
            private const ulong mrg32k3a_2pow32 = 4294967296;

            private const ulong mrg32k3a_m1 = 4294967087;
            private const ulong mrg32k3a_m1c = 209;
            private const ulong mrg32k3a_m2 = 4294944443;
            private const ulong mrg32k3a_m2c = 22853;

            private const ulong mrg32k3a_a12 = 1403580;
            private const ulong mrg32k3a_a13 = (4294967087 - 810728);
            private const ulong mrg32k3a_a13n = 810728;

            private const ulong mrg32k3a_a21 = 527612;
            private const ulong mrg32k3a_a23 = (4294944443 - 1370589);
            private const ulong mrg32k3a_a23n = 1370589;


            private const double norm = 2.3283065498378288e-10;

            // State Variables
            private ulong[] v1 = new ulong[3];
            private ulong[] v2 = new ulong[3];

            #region .ctor
            public Mrg32k3aGiles(ulong seed1, ulong seed2)
            {
                if ((seed1 <= 0) || (seed1 >= mrg32k3a_m1))
                {
                    throw new Exception("Mrg32k3aGiles instanciation with bad first seed");
                }

                if ((seed2 <= 0) || (seed2 >= mrg32k3a_m2))
                {
                    throw new Exception("Mrg32k3aGiles instanciation with bad second seed");
                }

                v1[0] = v1[1] = v1[2] = seed1;
                v2[0] = v2[1] = v2[2] = seed2;
            }
            #endregion

            #region Copy Constructor
            private Mrg32k3aGiles(Mrg32k3aGiles copy)
            {
                this.v1 = new ulong[3];
                this.v2 = new ulong[3];
                for (int i = 0; i < 3; ++i)
                {
                    this.v1[i] = copy.v1[i];
                    this.v2[i] = copy.v2[i];
                }
            }
            #endregion


            #region NextUniform
            public double NextUniform()
            {
                ulong ul = NextStep();
                return (double)ul * norm;
            }
            #endregion

            #region NextStep
            private ulong NextStep()
            {
                ulong p = mrg32k3a_a12 * v1[1] + mrg32k3a_a13n * (mrg32k3a_m1 - v1[0]);

                // alternative 1: 64-bit modulus
                p = p % mrg32k3a_m1;

                // alternative 2: avoids 64-bit modulus
                //p = (p & (mrg32k3a_2pow32-1)) + (p>>32)*mrg32k3a_m1c;
                //if (p>=mrg32k3a_m1) p -= mrg32k3a_m1;

                v1[0] = v1[1]; v1[1] = v1[2]; v1[2] = p;

                p = mrg32k3a_a21 * v2[2] + mrg32k3a_a23n * (mrg32k3a_m2 - v2[0]);

                // alternative 1: 64-bit modulus
                p = p % mrg32k3a_m2;

                // alternative 2: avoids 64-bit modulus
                // p = (p & (mrg32k3a_2pow32-1)) + (p>>32)*mrg32k3a_m2c;
                // p = (p & (mrg32k3a_2pow32-1)) + (p>>32)*mrg32k3a_m2c;
                // if (p>=mrg32k3a_m2) p -= mrg32k3a_m2;

                v2[0] = v2[1]; v2[1] = v2[2]; v2[2] = p;

                p = v1[2] - v2[2];
                if (v1[2] < v2[2]) p += mrg32k3a_m1;

                return p;
            }
            #endregion

            #region Skip
            public void Skip(int offset)
            {
                ulong[] A1 = new ulong[9];
                ulong[] A2 = new ulong[9];
                ulong[] A1b = new ulong[9];
                ulong[] A2b = new ulong[9];
                ulong[] v1b = new ulong[3];
                ulong[] v2b = new ulong[3];

                ulong[] CPU_mrg32k3a_v1 = new ulong[3];
                ulong[] CPU_mrg32k3a_v2 = new ulong[3];
                for (int n = 0; n < 3; n++)
                {
                    CPU_mrg32k3a_v1[n] = v1[n];
                    CPU_mrg32k3a_v2[n] = v2[n];
                }

                // initialise skip-ahead matrices
                A1[0] = 0; A1[3] = 1; A1[6] = 0;
                A1[1] = 0; A1[4] = 0; A1[7] = 1;
                A1[2] = mrg32k3a_a13; A1[5] = mrg32k3a_a12; A1[8] = 0;

                A2[0] = 0; A2[3] = 1; A2[6] = 0;
                A2[1] = 0; A2[4] = 0; A2[7] = 1;
                A2[2] = mrg32k3a_a23; A2[5] = 0; A2[8] = mrg32k3a_a21;

                while (offset > 0)
                {
                    // apply offset to seed vectors
                    if (offset % 2 == 1)
                    {
                        for (int i = 0; i < 3; i++)
                        {
                            v1b[i] = ((A1[i + 3 * 0] * CPU_mrg32k3a_v1[0]) % mrg32k3a_m1
                                     + (A1[i + 3 * 1] * CPU_mrg32k3a_v1[1]) % mrg32k3a_m1
                                     + (A1[i + 3 * 2] * CPU_mrg32k3a_v1[2]) % mrg32k3a_m1) % mrg32k3a_m1;
                            v2b[i] = ((A2[i + 3 * 0] * CPU_mrg32k3a_v2[0]) % mrg32k3a_m2
                                     + (A2[i + 3 * 1] * CPU_mrg32k3a_v2[1]) % mrg32k3a_m2
                                     + (A2[i + 3 * 2] * CPU_mrg32k3a_v2[2]) % mrg32k3a_m2) % mrg32k3a_m2;
                        }
                        for (int i = 0; i < 3; i++)
                        {
                            CPU_mrg32k3a_v1[i] = v1b[i];
                            CPU_mrg32k3a_v2[i] = v2b[i];
                        }
                    }
                    offset = offset / 2;

                    // square skip-ahead matrices
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            ulong a1 = 0;
                            ulong a2 = 0;
                            for (int k = 0; k < 3; k++)
                            {
                                a1 += (A1[i + 3 * k] * A1[k + 3 * j]) % mrg32k3a_m1;
                                a2 += (A2[i + 3 * k] * A2[k + 3 * j]) % mrg32k3a_m2;
                            }
                            A1b[i + 3 * j] = a1 % mrg32k3a_m1;
                            A2b[i + 3 * j] = a2 % mrg32k3a_m2;
                        }
                    }
                    for (int i = 0; i < 9; i++)
                    {
                        A1[i] = A1b[i];
                        A2[i] = A2b[i];
                    }
                }

                for (int n = 0; n < 3; n++)
                {
                    v1[n] = CPU_mrg32k3a_v1[n];
                    v2[n] = CPU_mrg32k3a_v2[n];
                }
            }
            #endregion


        }


        public static double[,] compute_sabr_vanille(double maturity, double forward, double alpha, double beta,
            double nu, double rho, double[] strikes, int nb_tirages, int nbdate)
        {
            double deltaT = maturity / (nbdate - 1.0);
            double[,] result = new double[3, strikes.Length];

            double spot_lsv = Math.Log(forward);
            double vol_sto = alpha;

            double spot_lsv_anti = Math.Log(forward);
            double vol_sto_anti = alpha;

            double temps = 0.0;


            //double u1, u2;
            double g1, g2;
            Random random = new MersenneTwister(42); // seed 42
            for (int i = 0; i < nb_tirages; i++)
            {

                spot_lsv = Math.Log(forward);
                vol_sto = alpha;

                spot_lsv_anti = Math.Log(forward);
                vol_sto_anti = alpha;

                temps = 0.0;



                for (int t = 0; t < nbdate - 1; t++)
                {


                    /*
                    u1 = new_gen.NextUniform();
                    u2 = new_gen.NextUniform();
                    u1 = Math.Max(Math.Min(u1, 1 - MACHINE_PRECISION), MACHINE_PRECISION);
                    u2 = Math.Max(Math.Min(u2, 1 - MACHINE_PRECISION), MACHINE_PRECISION);
                    g1 = ridounet.SABRFormula.inverseGaussianCumulative(u1);
                    g2 = ridounet.SABRFormula.inverseGaussianCumulative(u2);
                    */
                    g1 = Normal.Sample(random, 0.0, 1.0);
                    g2 = Normal.Sample(random, 0.0, 1.0);

                    spot_lsv = spot_lsv + Math.Exp((beta - 1) * spot_lsv) * vol_sto * (-0.5 * Math.Exp((beta - 1) * spot_lsv) * vol_sto * deltaT + Math.Sqrt(deltaT) * g1);
                    vol_sto = vol_sto + nu * vol_sto * Math.Sqrt(deltaT) * (rho * g1 + Math.Sqrt(1 - rho * rho) * g2);


                    spot_lsv_anti = spot_lsv_anti + Math.Exp((beta - 1) * spot_lsv_anti) * vol_sto_anti * (-0.5 * Math.Exp((beta - 1) * spot_lsv_anti) * vol_sto_anti * deltaT - Math.Sqrt(deltaT) * g1);
                    vol_sto_anti = vol_sto_anti + vol_sto_anti * nu * Math.Sqrt(deltaT) * (-rho * g1 - Math.Sqrt(1 - rho * rho) * g2);

                    temps += deltaT;
                }


                for (int j = 0; j < strikes.Length; j++)
                {
                    double tmp = Math.Max(Math.Exp(spot_lsv) - strikes[j], 0.0);

                    result[0, j] += tmp / (2.0 * nb_tirages);
                    result[1, j] += tmp * tmp / (2.0 * nb_tirages);

                    double tmp_anti = Math.Max(Math.Exp(spot_lsv_anti) - strikes[j], 0.0);

                    result[0, j] += tmp_anti / (2.0 * nb_tirages);
                    result[1, j] += tmp_anti * tmp_anti / (2.0 * nb_tirages);

                }

            }
            for (int j = 0; j < strikes.Length; j++)
            {
                result[2, j] = auxiliary.BlackScholes.ImpVolBS(forward, strikes[j], maturity, result[0, j]);
            }
            return result;
        }


        public static double Hagan_Formula(double f, double alpha, double beta, double rho, double nu, double K, double T)
        {
            double z = (Math.Pow(f, 1 - beta) - Math.Pow(K, 1 - beta)) / (alpha * (1 - beta));
            double x = Math.Log((-rho + nu * z + Math.Sqrt(1 - 2.0 * nu * rho * z + nu * nu * z * z)) / (1.0 - rho)) / nu;
            double LogfSurK = Math.Log(f / K);
            double EffetMat = (1.0 + (1.0 / 24 * alpha * alpha * (1 - beta) * (1 - beta) * Math.Pow(f * K, -1 + beta) + 1.0 / 4 * alpha * beta * Math.Pow(f * K, (-1 + beta) / 2.0) * nu * rho + 1.0 / 24 * nu * nu * (2.0 - 3.0 * rho * rho)) * T) /
                 (1.0 + 1.0 / 24 * (1 - beta) * (1 - beta) * LogfSurK * LogfSurK + (1 - beta) * (1 - beta) * (1 - beta) * (1 - beta) * LogfSurK * LogfSurK * LogfSurK * LogfSurK / 1920.0);
            if (Math.Abs(z) < 0.00001)
            {
                return alpha / Math.Pow(f * K, (1 - beta) / 2.0) * 1.0 / (1.0 + nu * rho * z / 2.0 + (-nu * nu / 6.0 + nu * nu * rho * rho / 2.0) * z * z) * EffetMat;
            }
            else
            {
                return alpha / Math.Pow(f * K, (1 - beta) / 2.0) * z / x * EffetMat;
            }

        }


    }
}



namespace MachineLearningPriceGenerator
{
    public class DeterministicVolSmile
    {
        double a, b, c;
        int nbstrikes=0;

        public DeterministicVolSmile(double[] smile, double v1, double v2, double v3)
        {
            nbstrikes = smile.Length;
            double x1 = smile[0]; double x2 = 1.0; double x3 = smile[nbstrikes - 1];
           
            a = -((v1 * x3 + v2 * x1 + v3 * x2 - v1 * x2 - v2 * x3 - v3 * x1) / ((x2 - x1) * (x2 - x3) * (x3 - x1)));
            b = -((v1 * (x2 * x2) + v2 * (x3 * x3) + v3 * (x1 * x1) - v1 * (x3 * x3) -
                   v2 * (x1 * x1) - v3 * (x2 * x2)) / ((x1 - x2) * (x1 - x3) * (x2 - x3)));
            c = -((v1 * x2 * (x3 * x3) + v2 * (x1 * x1) * x3 + v3 * x1 * (x2 * x2) - v1 * (x2 * x2) * x3 -
                    v2 * x1 * (x3 * x3) - v3 * (x1 * x1) * x2) / ((x1 - x2) * (x1 - x3) * (x2 - x3)));
        }
        public double get_Volatility(double Strike, double smileApplatissement)
        {
            double vol = a * smileApplatissement * Strike * Strike + b * Strike + c;
            return vol;
        }
        
    }
   

    public class ThreadComputationBloc

    {
        #region variable locales

        public static int nbTraj;
        public static int nbprice;
        public static int nbBlocDisplay;
        public static int NbThread;
        public static int n_gauche;
        public static int n_droite;
        public static int nbstrikes;
        public static int volSurfNbDate;

        public static int N3Flag;
        public static int nbIntegSteps;
        public static double startIntegration;

        public static double smallestNbDate;
        public static double largestNbDate;

        public static double smallestmaturity;
        public static double largestmaturity;
        public static double smallestF;
        public static double largestF;
        public static double smallestmu;
        public static double largestmu;
        public static double smallestSigma;
        public static double largestSigma;
        public static double smallestRho;
        public static double largestRho;

        public static double smallestBonus;
        public static double largestBonus;
        public static double smallestYetiBarrier;
        public static double largestYetiBarrier;
        public static double smallestYetiCoupon;
        public static double largestYetiCoupon;

        public static double smallestPhoenixBarrier;
        public static double largestPhoenixBarrier;
        public static double smallestPhoenixCoupon;
        public static double largestPhoenixCoupon;
        public static double smallestPDIBarrier;
        public static double largestPDIBarrier;
        public static double smallestPDIGearing;
        public static double largestPDIGearing;
        public static double smallestPDIStrike;
        public static double largestPDIStrike;
        public static double smallestPDIType;
        public static double largestPDIType;

        public static double smallestF0;
        public static double largestF0;

        public static double smallestalpha0;
        public static double largestalpha0;

        public static double smallestbeta;
        public static double largestbeta;

        public static double smallestbeta2;
        public static double largestbeta2;

        public static double smallestd;
        public static double largestd;

        public static double smallestgamma;
        public static double largestgamma;

        public static double smallestnu;
        public static double largestnu;

        public static double smallestomega;
        public static double largestomega;

        public static double smallestlambda;
        public static double largestlambda;

        public static double smallestrho;
        public static double largestrho;


        public int threadId;
        public double[,] correlationMatrix;
        public double[] Svect;
        public double[] muvect;
        public double[] sigmavect;
        public static int random_init;

        public static int[,] NbDate { set; get; }

        public static double[,] maturity { set; get; }
        public static double[,] S1 { set; get; }
        public static double[,] S2 { set; get; }
        public static double[,] S3 { set; get; }
        public static double[,] mu1 { set; get; }
        public static double[,] mu2 { set; get; }
        public static double[,] mu3 { set; get; }
        public static double[,] sigma1 { set; get; }
        public static double[,] sigma2 { set; get; }
        public static double[,] sigma3 { set; get; }
        public static double[,] rho12 { set; get; }
        public static double[,] rho13 { set; get; }
        public static double[,] rho23 { set; get; }

        public static double[,,,] localVol { set; get; }
        public static double[] InitialStrikeList { set; get; }
        public static double[,] StrikeList { set; get; }
        public static double[,] VolSurfDateList { set; get; }

        public static double[,] Bonus { set; get; }
        public static double[,] YetiBarrier { set; get; }
        public static double[,] YetiCoupon { set; get; }
        public static double[,] PhoenixBarrier { set; get; }

        public static double[,] PhoenixCoupon { set; get; }
        public static double[,] PDIBarrier { set; get; }
        public static double[,] PDIGearing { set; get; }
        public static double[,] PDIStrike { set; get; }
        public static double[,] PDIType { set; get; }
        public static double[,,] forwardsTable { set; get; }

        public static double[,] F0 { set; get; }
        public static double[,] alpha0 { set; get; }
        public static double[,] beta { set; get; }
        public static double[,] beta2 { set; get; }
        public static double[,] d { set; get; }
        public static double[,] gamma { set; get; }
        public static double[,] nu { set; get; }
        public static double[,] lambda { set; get; }
        public static double[,] omega { set; get; }
        public static double[,] rho { set; get; }

        public static int milsteinflag = 0;
        public static int expmartingaleflag = 0;
        public static double[,] prices { set; get; }
        public static double[,,] StrikedOption { set; get; }
        public static double[,,] StrikedVol { set; get; }
        public static double[,,] strike { set; get; }

        public static double[] StrikeListInitiale;

        public static string filename;
        public static string origin;
        #endregion
        static public string makeFileName(string filename, string origin, int ifile)
        {
            return filename + "-"+origin + ifile + ".CSV";
        }

        public ThreadComputationBloc(int ThreadNb)
        {
            threadId = ThreadNb;
            correlationMatrix = new double[3, 3];
            Svect = new double[3];
            muvect = new double[3];
            sigmavect = new double[3];
        }
        public double generateRandom(double xmin, double xmax, double mean, double sigma)
        {
            double resu = 2 * xmax;
            while ((resu > xmax) || (resu < xmin))
            {
                Random random = new MersenneTwister();
                resu = Normal.Sample(random, mean, sigma);
            }
            return resu;
        }
        public double exemplePricing_YetiPhoenixPriceMC()
        {
            double maturity = 3.0;
            int nbDate = 12;
            double S1 = 100.0;
            double S2 = 100.0;
            double S3 = 100.0;
            double mu1 = 0.01;
            double mu2 = 0.01;
            double mu3 = 0.01;
            double sigma1 = 0.25;
            double sigma2 = 0.26;
            double sigma3 = 0.27;
            double rho12 = 0.87;
            double rho13 = 0.89;
            double rho23 = 0.88;
            double Bonus = 0.01;
            double YetiBarrier = 95;
            double YetiCoupon = 2;
            double PhoenixBarrier = 90.0;
            double PhoenixCoupon = 2.25;
            double PDIBarrier = 50.0;
            double PDIGearing = -1.0;
            double PDIStrike = 40.0;
            double PDIType = -1.0;
            double[] forwards;

            Svect[0] = S1; Svect[1] = S2; Svect[2] = S3;
            muvect[0] = mu1; muvect[1] = mu2; muvect[2] = mu3;
            sigmavect[0] = sigma1; sigmavect[1] = sigma2; sigmavect[2] = sigma3;
            correlationMatrix[0, 0] = 1.0; correlationMatrix[0, 1] = rho12; correlationMatrix[0, 2] = rho13;
            correlationMatrix[1, 0] = rho12; correlationMatrix[1, 1] = 1.0; correlationMatrix[1, 2] = rho23;
            correlationMatrix[2, 0] = rho13; correlationMatrix[2, 1] = rho23; correlationMatrix[2, 2] = 1.0;

            double price = pricing.Export_MonteCarlo.Export_YetiPhoenixPriceMC(Svect, muvect, sigmavect, correlationMatrix, maturity,
                 Bonus, YetiBarrier, YetiCoupon, PhoenixBarrier, PhoenixCoupon,
             PDIBarrier, PDIGearing, PDIStrike, PDIType, nbDate, nbTraj, out forwards);
            Console.WriteLine("price=");
            Console.WriteLine(price);
            return price;

        }

        public void Computeprices_Lookback(int threadid)
        {

            for (int isample = 0; isample < nbprice; isample++)
            {
                NbDate[threadid, isample] = (int)Math.Floor(0.5 + MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestNbDate, largestNbDate));

                maturity[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmaturity, largestmaturity);
                S1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);
                S2[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);
                S3[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);
                mu1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);
                mu2[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);
                mu3[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);
                bool indic = false;
                while (indic == false)
                {
                    rho12[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestRho, largestRho);
                    rho13[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestRho, largestRho);
                    rho23[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestRho, largestRho);
                    indic = ((1.0 - rho12[threadid, isample] * rho12[threadid, isample]) * (1.0 - rho13[threadid, isample] * rho13[threadid, isample]) >
                        (rho23[threadid, isample] - rho12[threadid, isample] * rho13[threadid, isample]) * (rho23[threadid, isample] - rho12[threadid, isample] * rho13[threadid, isample]));
                }
                sigma1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);
                sigma2[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);
                sigma3[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);

                Bonus[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestBonus, largestBonus);
                YetiBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiBarrier, largestYetiBarrier);
                YetiCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiCoupon, largestYetiCoupon);

                PhoenixBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixBarrier, largestPhoenixBarrier);
                PhoenixCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixCoupon, largestPhoenixCoupon);
                PDIBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIBarrier, largestPDIBarrier);
                PDIGearing[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIGearing, largestPDIGearing);
                PDIStrike[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIStrike, largestPDIStrike);
                PDIType[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIType, largestPDIType);

                Svect[0] = S1[threadid, isample]; Svect[1] = S2[threadid, isample]; Svect[2] = S3[threadid, isample];
                muvect[0] = mu1[threadid, isample]; muvect[1] = mu2[threadid, isample]; muvect[2] = mu3[threadid, isample];
                sigmavect[0] = sigma1[threadid, isample]; sigmavect[1] = sigma2[threadid, isample]; sigmavect[2] = sigma3[threadid, isample];
                correlationMatrix[0, 0] = 1.0; correlationMatrix[0, 1] = rho12[threadid, isample]; correlationMatrix[0, 2] = rho13[threadid, isample];
                correlationMatrix[1, 0] = rho12[threadid, isample]; correlationMatrix[1, 1] = 1.0; correlationMatrix[1, 2] = rho23[threadid, isample];
                correlationMatrix[2, 0] = rho13[threadid, isample]; correlationMatrix[2, 1] = rho23[threadid, isample]; correlationMatrix[2, 2] = 1.0;
                double[] forwards;
                Random random = new MersenneTwister(random_init); // seed 42

                /*double price = MonteCarlo.YetiPhoenixPriceMC(random, Svect, muvect, sigmavect, 
                  correlationMatrix, maturity[threadid, isample],
                     Bonus[threadid, isample], YetiBarrier[threadid, isample], YetiCoupon[threadid, isample],
                    PhoenixBarrier[threadid, isample], PhoenixCoupon[threadid, isample],
                    PDIBarrier[threadid, isample],
                    PDIGearing[threadid, isample], PDIStrike[threadid, isample], PDIType[threadid, isample],
                    NbDate[threadid, isample], nbTraj, out forwards);
                */

                double price = pricing.MonteCarlo.LookbackCallWorstPriceMC(random, Svect, muvect, sigmavect,

                    correlationMatrix, maturity[threadid, isample],
                    PDIStrike[threadid, isample], PDIType[threadid, isample],
                    NbDate[threadid, isample], nbTraj, out forwards);

                prices[threadid, isample] = price;
                int nbdate = NbDate[threadid, isample];
                for (int k = 0; k < nbdate; k++)
                {
                    forwardsTable[threadid, isample, k] = forwards[k];
                }
                if (isample % nbBlocDisplay == 0)
                {
                    Console.WriteLine("threadid =" + threadid + "  computing i =" + isample);
                }

            }
        }
        public void Computeprices_SABRVolMC(int threadid)
        {

            for (int isample = 0; isample < nbprice; isample++)
            {
                NbDate[threadid, isample] = (int)Math.Floor(0.5 + MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestNbDate, largestNbDate));

                maturity[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmaturity, largestmaturity);
                F0[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF0, largestF0);
                alpha0[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestalpha0, largestalpha0);
                beta[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestbeta, largestbeta);
                beta2[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestbeta2, largestbeta2);
                d[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestd, largestd);
                gamma[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestgamma, largestgamma);
                nu[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestnu, largestnu);
                lambda[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestlambda, largestlambda);
                omega[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestomega, largestomega);
                rho[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestrho, largestrho);

                Random random = new MersenneTwister(random_init); // seed 42
                double vol_atm = 0.2;
                int milsteinflag = 0;
                int expmartingaleflag = 1;
                double delta_call = 0.006, delta_put = -0.006, epsilon = 0.03;

                double[] strikes = pricing.SABR.generate_strikes(F0[threadid, isample], maturity[threadid, isample], vol_atm, n_gauche, n_droite, delta_call, delta_put, epsilon);
                double[,] CallVols = pricing.MonteCarlo.VolVanilleSetVolStoMC(F0[threadid, isample], alpha0[threadid, isample], beta[threadid, isample], beta2[threadid, isample], d[threadid, isample],
                    gamma[threadid, isample], nu[threadid, isample], lambda[threadid, isample], omega[threadid, isample], rho[threadid, isample], maturity[threadid, isample],
                            NbDate[threadid, isample], strikes, nbTraj, milsteinflag, expmartingaleflag);

                for (int k = 0; k < nbstrikes; k++)
                {
                    strike[threadid, isample, k] = strikes[k];
                    StrikedOption[threadid, isample, k] = CallVols[0, k];
                    StrikedVol[threadid, isample, k] = CallVols[1, k];
                }

                if (isample % nbBlocDisplay == 0)
                {
                    Console.WriteLine("threadid =" + threadid + "  computing i =" + isample);
                }

            }
        }
        public void Computeprices_YetiPhoenixPriceMC(int threadid)
        {

            for (int isample = 0; isample < nbprice; isample++)
            {
                NbDate[threadid, isample] = (int)Math.Floor(0.5 + MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestNbDate, largestNbDate));

                maturity[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmaturity, largestmaturity);
                S1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);
                S2[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);
                S3[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);
                mu1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);
                mu2[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);
                mu3[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);
                bool indic = false;
                while (indic == false)
                {
                    rho12[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestRho, largestRho);
                    rho13[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestRho, largestRho);
                    rho23[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestRho, largestRho);
                    indic = ((1.0 - rho12[threadid, isample] * rho12[threadid, isample]) * (1.0 - rho13[threadid, isample] * rho13[threadid, isample]) >
                        (rho23[threadid, isample] - rho12[threadid, isample] * rho13[threadid, isample]) * (rho23[threadid, isample] - rho12[threadid, isample] * rho13[threadid, isample]));
                }
                sigma1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);
                sigma2[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);
                sigma3[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);

                Bonus[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestBonus, largestBonus);
                YetiBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiBarrier, largestYetiBarrier);
                YetiCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiCoupon, largestYetiCoupon);

                PhoenixBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixBarrier, largestPhoenixBarrier);
                PhoenixCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixCoupon, largestPhoenixCoupon);
                PDIBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIBarrier, largestPDIBarrier);
                PDIGearing[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIGearing, largestPDIGearing);
                PDIStrike[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIStrike, largestPDIStrike);
                PDIType[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIType, largestPDIType);

                Svect[0] = S1[threadid, isample]; Svect[1] = S2[threadid, isample]; Svect[2] = S3[threadid, isample];
                muvect[0] = mu1[threadid, isample]; muvect[1] = mu2[threadid, isample]; muvect[2] = mu3[threadid, isample];
                sigmavect[0] = sigma1[threadid, isample]; sigmavect[1] = sigma2[threadid, isample]; sigmavect[2] = sigma3[threadid, isample];
                correlationMatrix[0, 0] = 1.0; correlationMatrix[0, 1] = rho12[threadid, isample]; correlationMatrix[0, 2] = rho13[threadid, isample];
                correlationMatrix[1, 0] = rho12[threadid, isample]; correlationMatrix[1, 1] = 1.0; correlationMatrix[1, 2] = rho23[threadid, isample];
                correlationMatrix[2, 0] = rho13[threadid, isample]; correlationMatrix[2, 1] = rho23[threadid, isample]; correlationMatrix[2, 2] = 1.0;
                double[] forwards;
                Random random = new MersenneTwister(random_init); // seed 42

                double price = pricing.MonteCarlo.YetiPhoenixPriceMC(random, Svect, muvect, sigmavect, correlationMatrix, maturity[threadid, isample],
                     Bonus[threadid, isample], YetiBarrier[threadid, isample], YetiCoupon[threadid, isample],
                    PhoenixBarrier[threadid, isample], PhoenixCoupon[threadid, isample], PDIBarrier[threadid, isample],
                    PDIGearing[threadid, isample], PDIStrike[threadid, isample], PDIType[threadid, isample],
                    NbDate[threadid, isample], nbTraj, out forwards);
                prices[threadid, isample] = price;
                int nbdate = NbDate[threadid, isample];
                for (int k = 0; k < nbdate; k++)
                {
                    forwardsTable[threadid, isample, k] = forwards[k];
                }
                if (isample % nbBlocDisplay == 0)
                {
                    Console.WriteLine("threadid =" + threadid + "  computing i =" + isample);
                }

            }
        }
        public void Computeprices_YetiPhoenixPriceMC_1dim(int threadid)
        {
            
            for (int isample = 0; isample < nbprice; isample++)
            {
                NbDate[threadid, isample] = (int)Math.Floor(0.5 + MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestNbDate, largestNbDate));

                maturity[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmaturity, largestmaturity);
                S1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);
                
                mu1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);
               
                
                sigma1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);
                

                Bonus[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestBonus, largestBonus);
                YetiBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiBarrier, largestYetiBarrier);
                YetiCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiCoupon, largestYetiCoupon);

                PhoenixBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixBarrier, largestPhoenixBarrier);
                PhoenixCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixCoupon, largestPhoenixCoupon);
                PDIBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIBarrier, largestPDIBarrier);
                PDIGearing[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIGearing, largestPDIGearing);
                PDIStrike[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIStrike, largestPDIStrike);
                PDIType[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIType, largestPDIType);
                
                Random random = new MersenneTwister(random_init); // seed 42
             
                double price = pricing.MonteCarlo.YetiPhoenixPriceMC_1dim(random, S1[threadid, isample], mu1[threadid, isample], sigma1[threadid, isample], maturity[threadid, isample],
                     Bonus[threadid, isample], YetiBarrier[threadid, isample], YetiCoupon[threadid, isample],
                    PhoenixBarrier[threadid, isample], PhoenixCoupon[threadid, isample], PDIBarrier[threadid, isample],
                    PDIGearing[threadid, isample], PDIStrike[threadid, isample], PDIType[threadid, isample],
                    NbDate[threadid, isample], nbTraj);
                prices[threadid, isample] = price;
                int nbdate = NbDate[threadid, isample];
                
                if (isample % nbBlocDisplay == 0)
                {
                    Console.WriteLine("threadid =" + threadid + "  computing i =" + isample);
                }

            }
        }
        public void Computeprices_YetiPhoenixPriceMC_1dim_VolLoc(int threadid)
        {         
            double[,] LocalVol=new double[volSurfNbDate,nbstrikes];
            double voldebutT0, volfinT0, volmoneyT0, smileApplatissement;
            double[] smilestrikelist = new double[nbstrikes];
            for (int isample = 0; isample < nbprice; isample++)
            {
                NbDate[threadid, isample] = volSurfNbDate;
                maturity[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmaturity, largestmaturity);
                for (int iDate = 0; iDate < volSurfNbDate; iDate++)
                {
                    VolSurfDateList[threadid,iDate] = iDate/ NbDate[threadid, isample] * maturity[threadid, isample];
                }
                for (int istrike=0;istrike< ThreadComputationBloc.nbstrikes;istrike++)
                {
                    smilestrikelist[istrike] = ThreadComputationBloc.StrikeList[threadid, istrike];
                }               
                S1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestF, largestF);

                mu1[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestmu, largestmu);

                volfinT0 = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestSigma, largestSigma);
                voldebutT0 = volfinT0 + MathNet.Numerics.Distributions.ContinuousUniform.Sample(3.0* smallestSigma, 3.0 * largestSigma);
                volmoneyT0 = volfinT0 * MathNet.Numerics.Distributions.ContinuousUniform.Sample(0.8, 0.9);
                DeterministicVolSmile volsmile = new DeterministicVolSmile(smilestrikelist, voldebutT0, volmoneyT0, volfinT0);
                for (int istrike = 0; istrike < nbstrikes; istrike++)
                {
                    for (int iDate = 0; iDate < volSurfNbDate; iDate++)
                    {
                        smileApplatissement = 1.0 - (iDate / volSurfNbDate) * MathNet.Numerics.Distributions.ContinuousUniform.Sample(0.1, 0.9);
                        LocalVol[iDate,istrike] = volsmile.get_Volatility(smilestrikelist[istrike],smileApplatissement) +
                                MathNet.Numerics.Distributions.ContinuousUniform.Sample(0.0, 0.1);
                    }
                }
                Bonus[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestBonus, largestBonus);
                YetiBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiBarrier, largestYetiBarrier);
                YetiCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestYetiCoupon, largestYetiCoupon);
                PhoenixBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixBarrier, largestPhoenixBarrier);
                PhoenixCoupon[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPhoenixCoupon, largestPhoenixCoupon);
                PDIBarrier[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIBarrier, largestPDIBarrier);
                PDIGearing[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIGearing, largestPDIGearing);
                PDIStrike[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIStrike, largestPDIStrike);
                PDIType[threadid, isample] = MathNet.Numerics.Distributions.ContinuousUniform.Sample(smallestPDIType, largestPDIType);

                Random random = new MersenneTwister(random_init); // seed 42

                double price = pricing.MonteCarlo.YetiPhoenixSmiledPriceMC_1dim(random, S1[threadid, isample],
                     mu1[threadid, isample], smilestrikelist, NbDate[threadid, isample],LocalVol, maturity[threadid, isample],
                     Bonus[threadid, isample], YetiBarrier[threadid, isample], YetiCoupon[threadid, isample],
                     PhoenixBarrier[threadid, isample], PhoenixCoupon[threadid, isample], PDIBarrier[threadid, isample],
                     PDIGearing[threadid, isample], PDIStrike[threadid, isample], PDIType[threadid, isample],
                     NbDate[threadid, isample], nbTraj);
                prices[threadid, isample] = price;
                int nbdate = NbDate[threadid, isample];
                for (int istrike = 0; istrike < nbstrikes; istrike++)
                {
                    for (int iDate = 0; iDate < volSurfNbDate; iDate++)
                    {
                        localVol[threadid, isample, iDate, istrike] = LocalVol[ iDate, istrike];
                    }
                }
                if (isample % nbBlocDisplay == 0)
                {
                    Console.WriteLine("threadid =" + threadid + "  computing i =" + isample);
                }

            }
        }

        static public void SaveToCVS_SABRVolMC()
        {
            using (CSVHandler.CsvFileWriter writer = new CSVHandler.CsvFileWriter(filename))
            {

                CSVHandler.CsvRow rowtitle = new CSVHandler.CsvRow();
                rowtitle.Add("F0");
                rowtitle.Add("alpha0");
                rowtitle.Add("beta");
                rowtitle.Add("beta2");
                rowtitle.Add("d");
                rowtitle.Add("gamma");
                rowtitle.Add("nu");
                rowtitle.Add("omega");
                rowtitle.Add("lambda");
                rowtitle.Add("rho");
                rowtitle.Add("maturity");
                rowtitle.Add("strike");
                rowtitle.Add("option");
                rowtitle.Add("vol");

                writer.WriteRow(rowtitle);


                for (int threadid = 0; threadid < NbThread; threadid++)
                {
                    for (int i = 0; i < nbprice; i++)
                    {
                        for (int k = 0; k < nbstrikes; k++)
                        {
                            int iligne = i + threadid * nbprice + 2;
                            CSVHandler.CsvRow row = new CSVHandler.CsvRow();
                            row.Add(String.Format("{0}", F0[threadid, i]));
                            row.Add(String.Format("{0}", alpha0[threadid, i]));
                            row.Add(String.Format("{0}", beta[threadid, i]));
                            row.Add(String.Format("{0}", beta2[threadid, i]));
                            row.Add(String.Format("{0}", d[threadid, i]));
                            row.Add(String.Format("{0}", gamma[threadid, i]));
                            row.Add(String.Format("{0}", nu[threadid, i]));
                            row.Add(String.Format("{0}", omega[threadid, i]));
                            row.Add(String.Format("{0}", lambda[threadid, i]));
                            row.Add(String.Format("{0}", rho[threadid, i]));
                            row.Add(String.Format("{0}", maturity[threadid, i]));
                            row.Add(String.Format("{0}", strike[threadid, i, k]));
                            row.Add(String.Format("{0}", StrikedOption[threadid, i, k]));
                            row.Add(String.Format("{0}", StrikedVol[threadid, i, k]));

                            writer.WriteRow(row);
                        }

                    }
                }
            }

        }
        static public void SaveToCVS_YetiPhoenixPriceMC()
        {
            using (CSVHandler.CsvFileWriter writer = new CSVHandler.CsvFileWriter(filename))
            {
                CSVHandler.CsvRow rowtitle = new CSVHandler.CsvRow();
                rowtitle.Add("S1");
                rowtitle.Add("S2");
                rowtitle.Add("S3");
                rowtitle.Add("mu1");
                rowtitle.Add("mu2");
                rowtitle.Add("mu3");
                rowtitle.Add("sigma1");
                rowtitle.Add("sigma2");
                rowtitle.Add("sigma3");
                rowtitle.Add("rho12");
                rowtitle.Add("rho13");
                rowtitle.Add("rho23");
                rowtitle.Add("bonus");
                rowtitle.Add("YetiBarrier");
                rowtitle.Add("YetiCoupon");
                rowtitle.Add("PhoenixBarrier");
                rowtitle.Add("PhoenixCoupon");
                rowtitle.Add("PDIBarrier");
                rowtitle.Add("PDIGearing");
                rowtitle.Add("PDIStrike");
                rowtitle.Add("PDIType");
                rowtitle.Add("maturity");
                rowtitle.Add("nbDates");
                rowtitle.Add("price");
                for (int j = 1; j < ThreadComputationBloc.largestNbDate; j++)
                {
                    rowtitle.Add("forward" + String.Format("{0}", j));
                }
                writer.WriteRow(rowtitle);
                for (int threadid = 0; threadid < NbThread; threadid++)
                {
                    for (int i = 0; i < nbprice; i++)
                    {
                        int iligne = i + threadid * nbprice + 2;
                        CSVHandler.CsvRow row = new CSVHandler.CsvRow();
                        row.Add(String.Format("{0}", S1[threadid, i]));
                        row.Add(String.Format("{0}", S2[threadid, i]));
                        row.Add(String.Format("{0}", S3[threadid, i]));
                        row.Add(String.Format("{0}", mu1[threadid, i]));
                        row.Add(String.Format("{0}", mu2[threadid, i]));
                        row.Add(String.Format("{0}", mu3[threadid, i]));
                        row.Add(String.Format("{0}", sigma1[threadid, i]));
                        row.Add(String.Format("{0}", sigma2[threadid, i]));
                        row.Add(String.Format("{0}", sigma3[threadid, i]));
                        row.Add(String.Format("{0}", rho12[threadid, i]));
                        row.Add(String.Format("{0}", rho13[threadid, i]));
                        row.Add(String.Format("{0}", rho23[threadid, i]));
                        row.Add(String.Format("{0}", Bonus[threadid, i]));
                        row.Add(String.Format("{0}", YetiBarrier[threadid, i]));
                        row.Add(String.Format("{0}", YetiCoupon[threadid, i]));
                        row.Add(String.Format("{0}", PhoenixBarrier[threadid, i]));
                        row.Add(String.Format("{0}", PhoenixCoupon[threadid, i]));
                        row.Add(String.Format("{0}", PDIBarrier[threadid, i]));
                        row.Add(String.Format("{0}", PDIGearing[threadid, i]));
                        row.Add(String.Format("{0}", PDIStrike[threadid, i]));
                        row.Add(String.Format("{0}", PDIType[threadid, i]));
                        row.Add(String.Format("{0}", maturity[threadid, i]));
                        row.Add(String.Format("{0}", NbDate[threadid, i]));
                        row.Add(String.Format("{0}", prices[threadid, i]));

                        for (int j = 1; j < ThreadComputationBloc.NbDate[threadid, i]; j++)
                        {
                            row.Add(String.Format("{0}", forwardsTable[threadid, i, j]));
                        }

                        writer.WriteRow(row);
                    }
                }


            }
        }
        static public void SaveToCVS_YetiPhoenixPriceMC_1dim()
        {
            using (CSVHandler.CsvFileWriter writer = new CSVHandler.CsvFileWriter(filename))
            {
                CSVHandler.CsvRow rowtitle = new CSVHandler.CsvRow();
                rowtitle.Add("S1");              
                rowtitle.Add("mu1");             
                rowtitle.Add("sigma1");              
                rowtitle.Add("bonus");
                rowtitle.Add("YetiBarrier");
                rowtitle.Add("YetiCoupon");
                rowtitle.Add("PhoenixBarrier");
                rowtitle.Add("PhoenixCoupon");
                rowtitle.Add("PDIBarrier");
                rowtitle.Add("PDIGearing");
                rowtitle.Add("PDIStrike");
                rowtitle.Add("PDIType");
                rowtitle.Add("maturity");
                rowtitle.Add("nbDates");
                rowtitle.Add("price");
               
                writer.WriteRow(rowtitle);
                for (int threadid = 0; threadid < NbThread; threadid++)
                {
                    for (int i = 0; i < nbprice; i++)
                    {
                        int iligne = i + threadid * nbprice + 2;
                        CSVHandler.CsvRow row = new CSVHandler.CsvRow();
                        row.Add(String.Format("{0}", S1[threadid, i]));                      
                        row.Add(String.Format("{0}", mu1[threadid, i]));                     
                        row.Add(String.Format("{0}", sigma1[threadid, i]));                      
                        row.Add(String.Format("{0}", Bonus[threadid, i]));
                        row.Add(String.Format("{0}", YetiBarrier[threadid, i]));
                        row.Add(String.Format("{0}", YetiCoupon[threadid, i]));
                        row.Add(String.Format("{0}", PhoenixBarrier[threadid, i]));
                        row.Add(String.Format("{0}", PhoenixCoupon[threadid, i]));
                        row.Add(String.Format("{0}", PDIBarrier[threadid, i]));
                        row.Add(String.Format("{0}", PDIGearing[threadid, i]));
                        row.Add(String.Format("{0}", PDIStrike[threadid, i]));
                        row.Add(String.Format("{0}", PDIType[threadid, i]));
                        row.Add(String.Format("{0}", maturity[threadid, i]));
                        row.Add(String.Format("{0}", NbDate[threadid, i]));
                        row.Add(String.Format("{0}", prices[threadid, i]));

                        writer.WriteRow(row);
                    }
                }


            }
        }

      
        static public void SaveToCVS_YetiPhoenixPriceMC_1dim_VolLoc(int ifile)
        {
            using (CSVHandler.CsvFileWriter writer = new CSVHandler.CsvFileWriter(makeFileName(filename, origin, ifile)))
            {
                CSVHandler.CsvRow rowtitle = new CSVHandler.CsvRow();
                rowtitle.Add("S1");
                rowtitle.Add("mu1");
                rowtitle.Add("bonus");
                rowtitle.Add("YetiBarrier");
                rowtitle.Add("YetiCoupon");
                rowtitle.Add("PhoenixBarrier");
                rowtitle.Add("PhoenixCoupon");
                rowtitle.Add("PDIBarrier");
                rowtitle.Add("PDIGearing");
                rowtitle.Add("PDIStrike");
                rowtitle.Add("PDIType");
                rowtitle.Add("maturity");
                rowtitle.Add("nbDates");
                rowtitle.Add("price");
                for (int idate = 0; idate < ThreadComputationBloc.volSurfNbDate; idate++)
                {
                    for (int istrike = 0; istrike < ThreadComputationBloc.nbstrikes; istrike++)
                    {
                        rowtitle.Add("vol" + "-date" + idate + "-strike" + istrike);
                    }
                }

                writer.WriteRow(rowtitle);
                for (int threadid = 0; threadid < NbThread; threadid++)
                {
                    for (int i = 0; i < nbprice; i++)
                    {
                        int iligne = i + threadid * nbprice + 2;
                        CSVHandler.CsvRow row = new CSVHandler.CsvRow();
                        row.Add(String.Format("{0}", S1[threadid, i]));
                        row.Add(String.Format("{0}", mu1[threadid, i]));
                        row.Add(String.Format("{0}", Bonus[threadid, i]));
                        row.Add(String.Format("{0}", YetiBarrier[threadid, i]));
                        row.Add(String.Format("{0}", YetiCoupon[threadid, i]));
                        row.Add(String.Format("{0}", PhoenixBarrier[threadid, i]));
                        row.Add(String.Format("{0}", PhoenixCoupon[threadid, i]));
                        row.Add(String.Format("{0}", PDIBarrier[threadid, i]));
                        row.Add(String.Format("{0}", PDIGearing[threadid, i]));
                        row.Add(String.Format("{0}", PDIStrike[threadid, i]));
                        row.Add(String.Format("{0}", PDIType[threadid, i]));
                        row.Add(String.Format("{0}", maturity[threadid, i]));
                        row.Add(String.Format("{0}", NbDate[threadid, i]));
                        row.Add(String.Format("{0}", prices[threadid, i]));
                        for (int idate = 0; idate < ThreadComputationBloc.volSurfNbDate; idate++)
                        {
                            for (int istrike = 0; istrike < ThreadComputationBloc.nbstrikes; istrike++)
                            {
                                row.Add(String.Format("{0}", localVol[threadid, i, idate, istrike]));
                            }
                        }

                        writer.WriteRow(row);
                    }
                }


            }
        }

        static public void SetUp_SABR(int iserie)
        {
            #region main params
            nbTraj = 100000;
            nbprice = 1000;
            nbBlocDisplay = 10;
            n_gauche = 11;
            n_droite = 17;
            nbstrikes = n_gauche + n_droite - 1;
            #endregion
            #region les limites de simulation
            ThreadComputationBloc.random_init = 44;
            ThreadComputationBloc.smallestNbDate = 1000;
            ThreadComputationBloc.largestNbDate = 1000;

            ThreadComputationBloc.smallestmaturity = 0.05;
            ThreadComputationBloc.largestmaturity = 10.0;

            ThreadComputationBloc.smallestF0 = 1.0;
            ThreadComputationBloc.largestF0 = 1.0;

            ThreadComputationBloc.smallestalpha0 = 0.0;
            ThreadComputationBloc.largestalpha0 = 0.5;

            ThreadComputationBloc.smallestbeta = 0.01;
            ThreadComputationBloc.largestbeta = 0.9;

            ThreadComputationBloc.smallestbeta2 = 0.01;
            ThreadComputationBloc.largestbeta2 = 0.9;

            ThreadComputationBloc.smallestd = 0.0;
            ThreadComputationBloc.largestd = 1.0;

            ThreadComputationBloc.smallestgamma = 0.01;
            ThreadComputationBloc.largestgamma = 2.0;

            ThreadComputationBloc.smallestnu = 0.00;
            ThreadComputationBloc.largestnu = 3.0;

            ThreadComputationBloc.smallestomega = 0.00;
            ThreadComputationBloc.largestomega = 0.5;

            ThreadComputationBloc.smallestlambda = 0;
            ThreadComputationBloc.largestlambda = 0.5;

            ThreadComputationBloc.smallestrho = -0.95;
            ThreadComputationBloc.largestrho = 0.95;

            #endregion
            #region fichier de sortie

            filename = @"G:\Fcommun\Papiers_Olivier2\Resultat_SABR\ResultatGeneratedprices-extendedSABR2-" +
             nbprice + "-" + ThreadComputationBloc.smallestmaturity + "-" + ThreadComputationBloc.largestmaturity + "-" +
             ThreadComputationBloc.smallestF + "-" + ThreadComputationBloc.largestF + "-" + nbTraj + "-" + iserie + ".CSV";

            #endregion
            #region les containers
            NbDate = new int[NbThread, nbprice];
            maturity = new double[NbThread, nbprice];
            F0 = new double[NbThread, nbprice];
            alpha0 = new double[NbThread, nbprice];
            beta = new double[NbThread, nbprice];
            beta2 = new double[NbThread, nbprice];
            d = new double[NbThread, nbprice];
            gamma = new double[NbThread, nbprice];
            nu = new double[NbThread, nbprice];
            omega = new double[NbThread, nbprice];
            lambda = new double[NbThread, nbprice];
            rho = new double[NbThread, nbprice];

            StrikedOption = new double[NbThread, nbprice, nbstrikes];
            StrikedVol = new double[NbThread, nbprice, nbstrikes];
            strike = new double[NbThread, nbprice, nbstrikes];

            #endregion
            #region general parameters
            N3Flag = 2;
            nbIntegSteps = 35;
            startIntegration = -7;
            #endregion
        }
        static public void SetUp_YetiPhoenix_3dim()
        {
            #region main params
            nbTraj = 10000;
            nbprice = 500;
            nbBlocDisplay = 10;
            #endregion
            #region les limites de simulation
            ThreadComputationBloc.random_init = 43;
            ThreadComputationBloc.smallestNbDate = 12;
            ThreadComputationBloc.largestNbDate = 12;

            ThreadComputationBloc.smallestmaturity = 1;
            ThreadComputationBloc.largestmaturity = 3;
            ThreadComputationBloc.smallestF = 90;
            ThreadComputationBloc.largestF = 110;
            ThreadComputationBloc.largestmu = 0.05;
            ThreadComputationBloc.smallestmu = 0.0;

            ThreadComputationBloc.smallestSigma = 0.15;
            ThreadComputationBloc.largestSigma = 0.3;
            ThreadComputationBloc.smallestRho = 0.4;
            ThreadComputationBloc.largestRho = 1;

            ThreadComputationBloc.smallestBonus = -2;
            ThreadComputationBloc.largestBonus = 2;
            ThreadComputationBloc.smallestYetiBarrier = 90;
            ThreadComputationBloc.largestYetiBarrier = 110;
            ThreadComputationBloc.smallestYetiCoupon = 0;
            ThreadComputationBloc.largestYetiCoupon = 2;

            ThreadComputationBloc.smallestPhoenixBarrier = 80;
            ThreadComputationBloc.largestPhoenixBarrier = 100;
            ThreadComputationBloc.smallestPhoenixCoupon = 0.5;
            ThreadComputationBloc.largestPhoenixCoupon = 2;
            ThreadComputationBloc.smallestPDIBarrier = 40;
            ThreadComputationBloc.largestPDIBarrier = 70;
            ThreadComputationBloc.smallestPDIGearing = -5;
            ThreadComputationBloc.largestPDIGearing = +5;
            ThreadComputationBloc.smallestPDIStrike = 80;
            ThreadComputationBloc.largestPDIStrike = 120;
            ThreadComputationBloc.smallestPDIType = 1;
            ThreadComputationBloc.largestPDIType = 1;
            #endregion
            #region fichier de sortie

            filename = "D:\\DDisk\\AI\\Projet_AI\\Generateprice_YetiPhoenix_1dim\\Generated_Prices\\YetiPhoenix_1dim-" +
             nbprice + "-" + ThreadComputationBloc.smallestmaturity + "-" + ThreadComputationBloc.largestmaturity + "-" +
             ThreadComputationBloc.smallestF + "-" + ThreadComputationBloc.largestF + "-" + nbTraj + ".CSV";


            #endregion
            #region les containers
            NbDate = new int[NbThread, nbprice];
            maturity = new double[NbThread, nbprice];
            S1 = new double[NbThread, nbprice];
            S2 = new double[NbThread, nbprice];
            S3 = new double[NbThread, nbprice];
            mu1 = new double[NbThread, nbprice];
            mu2 = new double[NbThread, nbprice];
            mu3 = new double[NbThread, nbprice];
            sigma1 = new double[NbThread, nbprice];
            sigma2 = new double[NbThread, nbprice];
            sigma3 = new double[NbThread, nbprice];
            rho12 = new double[NbThread, nbprice];
            rho13 = new double[NbThread, nbprice];
            rho23 = new double[NbThread, nbprice];
            Bonus = new double[NbThread, nbprice];
            YetiBarrier = new double[NbThread, nbprice];
            YetiCoupon = new double[NbThread, nbprice];

            PhoenixBarrier = new double[NbThread, nbprice];
            PhoenixCoupon = new double[NbThread, nbprice];
            PDIBarrier = new double[NbThread, nbprice];
            PDIGearing = new double[NbThread, nbprice];
            PDIStrike = new double[NbThread, nbprice];
            PDIType = new double[NbThread, nbprice];
            prices = new double[NbThread, nbprice];
            forwardsTable = new double[NbThread, nbprice, (int)ThreadComputationBloc.largestNbDate];
            #endregion
            #region general parameters
            N3Flag = 2;
            nbIntegSteps = 35;
            startIntegration = -7;
            #endregion
        }
        static public void SetUp_YetiPhoenix_1dim()
        {
            #region main params
            nbTraj = 100000;
            nbprice = 120000;
            nbBlocDisplay = 10;
            #endregion
            #region les limites de simulation
            ThreadComputationBloc.random_init = 43;
            ThreadComputationBloc.smallestNbDate = 12;
            ThreadComputationBloc.largestNbDate = 12;

            ThreadComputationBloc.smallestmaturity = 0.5;
            ThreadComputationBloc.largestmaturity = 7;
            ThreadComputationBloc.smallestF = 90;
            ThreadComputationBloc.largestF = 110;
            ThreadComputationBloc.largestmu = 0.05;
            ThreadComputationBloc.smallestmu = 0.0;

            ThreadComputationBloc.smallestSigma = 0.1;
            ThreadComputationBloc.largestSigma = 0.3;


            ThreadComputationBloc.smallestBonus = -2;
            ThreadComputationBloc.largestBonus = 3;
            ThreadComputationBloc.smallestYetiBarrier = 90;
            ThreadComputationBloc.largestYetiBarrier = 110;
            ThreadComputationBloc.smallestYetiCoupon = 0;
            ThreadComputationBloc.largestYetiCoupon = 2;

            ThreadComputationBloc.smallestPhoenixBarrier = 80;
            ThreadComputationBloc.largestPhoenixBarrier = 100;
            ThreadComputationBloc.smallestPhoenixCoupon = 0.5;
            ThreadComputationBloc.largestPhoenixCoupon = 2;
            ThreadComputationBloc.smallestPDIBarrier = 30;
            ThreadComputationBloc.largestPDIBarrier = 80;
            ThreadComputationBloc.smallestPDIGearing = -3;
            ThreadComputationBloc.largestPDIGearing = +3;
            ThreadComputationBloc.smallestPDIStrike = 80;
            ThreadComputationBloc.largestPDIStrike = 120;
            ThreadComputationBloc.smallestPDIType = 1;
            ThreadComputationBloc.largestPDIType = 1;
            #endregion
            #region fichier de sortie
            filename = "D:\\DDisk\\AI\\Projet_AI\\Generateprice_YetiPhoenix_1dim\\Generated_Prices\\YetiPhoenix_1dim-" +
             nbprice + "-" + ThreadComputationBloc.smallestmaturity + "-" + ThreadComputationBloc.largestmaturity + "-" +
             ThreadComputationBloc.smallestF + "-" + ThreadComputationBloc.largestF + "-" + nbTraj + ".CSV";


            #endregion
            #region les containers
            NbDate = new int[NbThread, nbprice];
            maturity = new double[NbThread, nbprice];
            S1 = new double[NbThread, nbprice];
            mu1 = new double[NbThread, nbprice];
            sigma1 = new double[NbThread, nbprice];

            Bonus = new double[NbThread, nbprice];
            YetiBarrier = new double[NbThread, nbprice];
            YetiCoupon = new double[NbThread, nbprice];

            PhoenixBarrier = new double[NbThread, nbprice];
            PhoenixCoupon = new double[NbThread, nbprice];
            PDIBarrier = new double[NbThread, nbprice];
            PDIGearing = new double[NbThread, nbprice];
            PDIStrike = new double[NbThread, nbprice];
            PDIType = new double[NbThread, nbprice];
            prices = new double[NbThread, nbprice];
            forwardsTable = new double[NbThread, nbprice, (int)ThreadComputationBloc.largestNbDate];
            #endregion
            #region general parameters
            N3Flag = 2;
            nbIntegSteps = 35;
            startIntegration = -7;
            #endregion
        }
        static public void SetUp_YetiPhoenix_1dim_VolLoc()
        {
            #region main params
            nbTraj = 100000;
            nbprice = 30;
            nbBlocDisplay = 10;
            nbstrikes = 15;
            double[] StrikeListInitiale0 = { 0.4, 0.6, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0 };
            volSurfNbDate = 12;
            #endregion
            #region les limites de simulation
            ThreadComputationBloc.random_init = 43;
            ThreadComputationBloc.smallestNbDate = 12;
            ThreadComputationBloc.largestNbDate = 12;

            ThreadComputationBloc.smallestmaturity = 0.5;
            ThreadComputationBloc.largestmaturity = 7;
            ThreadComputationBloc.smallestF = 90;
            ThreadComputationBloc.largestF = 110;
            ThreadComputationBloc.largestmu = 0.03;
            ThreadComputationBloc.smallestmu = 0.0;

            ThreadComputationBloc.smallestSigma = 0.1;
            ThreadComputationBloc.largestSigma = 0.3;

            ThreadComputationBloc.smallestBonus = -2;
            ThreadComputationBloc.largestBonus = 3;
            ThreadComputationBloc.smallestYetiBarrier = 90;
            ThreadComputationBloc.largestYetiBarrier = 110;
            ThreadComputationBloc.smallestYetiCoupon = 0;
            ThreadComputationBloc.largestYetiCoupon = 2;

            ThreadComputationBloc.smallestPhoenixBarrier = 80;
            ThreadComputationBloc.largestPhoenixBarrier = 100;
            ThreadComputationBloc.smallestPhoenixCoupon = 0.5;
            ThreadComputationBloc.largestPhoenixCoupon = 2;
            ThreadComputationBloc.smallestPDIBarrier = 30;
            ThreadComputationBloc.largestPDIBarrier = 80;
            ThreadComputationBloc.smallestPDIGearing = -3;
            ThreadComputationBloc.largestPDIGearing = +3;
            ThreadComputationBloc.smallestPDIStrike = 80;
            ThreadComputationBloc.largestPDIStrike = 120;
            ThreadComputationBloc.smallestPDIType = 1;
            ThreadComputationBloc.largestPDIType = 1;
            #endregion
            #region fichier de sortie
            origin = "Olivier";
            filename = "P:\\Services\\Echanges\\MachineLeaning\\PocNVIDIA\\"+origin+"\\YetiPhoenix_1dim_VolLoc-" +
             nbprice + "-" + ThreadComputationBloc.smallestmaturity + "-" + ThreadComputationBloc.largestmaturity + "-" +
             ThreadComputationBloc.smallestF + "-" + ThreadComputationBloc.largestF + "-" + nbTraj;
           

            #endregion
            #region les containers
            NbDate = new int[NbThread, nbprice];
            maturity = new double[NbThread, nbprice];
            S1 = new double[NbThread, nbprice];
            mu1 = new double[NbThread, nbprice];

            Bonus = new double[NbThread, nbprice];
            YetiBarrier = new double[NbThread, nbprice];
            YetiCoupon = new double[NbThread, nbprice];

            PhoenixBarrier = new double[NbThread, nbprice];
            PhoenixCoupon = new double[NbThread, nbprice];
            PDIBarrier = new double[NbThread, nbprice];
            PDIGearing = new double[NbThread, nbprice];
            PDIStrike = new double[NbThread, nbprice];
            PDIType = new double[NbThread, nbprice];
            prices = new double[NbThread, nbprice];

            InitialStrikeList = new double[nbprice];
            StrikeList = new double[NbThread, nbstrikes];
            
            for(int istrike=0;istrike<nbstrikes;istrike++)
            {
                InitialStrikeList[istrike] = StrikeListInitiale0[istrike];
            }
            
            VolSurfDateList = new double[NbThread, volSurfNbDate];
            localVol = new double[NbThread, nbprice, volSurfNbDate, nbstrikes];
            #endregion
            #region general parameters
            N3Flag = 2;
            nbIntegSteps = 35;
            startIntegration = -7;
            #endregion
        }

    }


    class Program
    {
        static void Submain_YetiPhoenixMC(object threadid)
        {

            Console.WriteLine("threadid=" + (int)threadid + " Hello World!");
            ThreadComputationBloc threadblocInstance = new ThreadComputationBloc((int)threadid);
            threadblocInstance.Computeprices_YetiPhoenixPriceMC((int)threadid);
            Console.WriteLine("threadid=" + (int)threadid + " Calcul des prix Fini");
            Console.WriteLine("threadid=" + (int)threadid + "Ecriture des prix Fini");
        }
        static void Submain_YetiPhoenixMC_1dim(object threadid)
        {

            Console.WriteLine("threadid=" + (int)threadid + " Hello World!");
            ThreadComputationBloc threadblocInstance = new ThreadComputationBloc((int)threadid);
            threadblocInstance.Computeprices_YetiPhoenixPriceMC_1dim((int)threadid);
            Console.WriteLine("threadid=" + (int)threadid + " Calcul des prix Fini");
            Console.WriteLine("threadid=" + (int)threadid + "Ecriture des prix Fini");
        }
        static void Submain_YetiPhoenixMC_1dim_VolLoc(object threadid)
        {
            Console.WriteLine("threadid=" + (int)threadid + " Hello World!");
            ThreadComputationBloc threadblocInstance = new ThreadComputationBloc((int)threadid);
            threadblocInstance.Computeprices_YetiPhoenixPriceMC_1dim_VolLoc((int)threadid);
            Console.WriteLine("threadid=" + (int)threadid + " Calcul des prix Fini");
            Console.WriteLine("threadid=" + (int)threadid + "Ecriture des prix Fini");
        }
        static void Submain_SABRVolMC(object threadid)
        {

            Console.WriteLine("threadid=" + (int)threadid + " Hello World!");
            ThreadComputationBloc threadblocInstance = new ThreadComputationBloc((int)threadid);
            threadblocInstance.Computeprices_SABRVolMC((int)threadid);
            Console.WriteLine("threadid=" + (int)threadid + " Calcul des prix Fini");
            Console.WriteLine("threadid=" + (int)threadid + "Ecriture des prix Fini");
        }

        static void Main_Sabr2(string[] args)
        {
            double F0 = 1.0;
            double alpha0 = 0.5;
            double beta = 0.5;
            double beta2 = 0.6;
            double d = 0.1;
            double gamma = 1;
            double nu = 0.8;
            double lambda = 0;
            double omega = 0;
            double rho = -0.5;
            double Tmax = 3;
            int nbdate = 1000;
            int milsteinflag = 0;
            int expmartingaleflag = 0;
            // double[] strikes= { 0.7,0.8,0.9, 1.0, 1.1,1.2,1.3 };
            int nbTraj = 500000;
            double[,] CallVs;
            double delta_call = 0.01, delta_put = -0.01, epsilon = 0.03;
            double vol_atm = 0.2;
            int n_gauche = 11, n_droite = 17;
            double[] strikes = pricing.SABR.generate_strikes(F0, Tmax, vol_atm, n_gauche, n_droite, delta_call, delta_put, epsilon);
            var watch = System.Diagnostics.Stopwatch.StartNew();
            CallVs = pricing.MonteCarlo.VolVanilleSetVolStoMC(F0, alpha0, beta, beta2, d, gamma, nu, lambda, omega, rho, Tmax, nbdate, strikes, nbTraj, milsteinflag, expmartingaleflag);
            watch.Stop();
            var elapsedMs = watch.ElapsedMilliseconds;

            for (int k = 0; k < strikes.Length; k++)
            {
                Console.WriteLine(String.Concat(Convert.ToString(strikes[k]), " = ", Convert.ToString(CallVs[1, k])));
            }
            Console.WriteLine(String.Concat(" time (ms) = ", Convert.ToString(elapsedMs)));
            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();

        }
        static void Main_YetiPhoenixMC2(string[] args)   //YetiPhoenixPriceMC
        {

            Console.Out.WriteLine("Thread number :");
            string ligne = Console.In.ReadLine();
            ThreadComputationBloc.NbThread = int.Parse(ligne);

            DateTime start = DateTime.Now;
            ThreadComputationBloc.SetUp_YetiPhoenix_3dim();
            List<Thread> th = new List<Thread>();
            for (int i = 0; i < ThreadComputationBloc.NbThread; i++)
            {
                Thread thread = new Thread(new ParameterizedThreadStart(Submain_YetiPhoenixMC));
                th.Add(thread);
                thread.Start(i);
            }

            foreach (var thread in th)
            {
                thread.Join();
            }

            ThreadComputationBloc.SaveToCVS_YetiPhoenixPriceMC();


            DateTime end = DateTime.Now;
            Console.WriteLine("Computation time (sec): " + (end - start).TotalSeconds);
            Console.WriteLine("Computation time  poun un chunk de prix (sec): " + (end - start).TotalSeconds / ThreadComputationBloc.NbThread);

            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();
        }
        static void Main_Sabr1(string[] args)     //ExtendedSABR  iSerie = 0 dans SetUp_SABR(iserie)
        {
            Console.Out.WriteLine("Thread number :");
            string ligne = Console.In.ReadLine();
            ThreadComputationBloc.NbThread = int.Parse(ligne);
            DateTime start = DateTime.Now;
            ThreadComputationBloc.SetUp_SABR(0);
            List<Thread> th = new List<Thread>();
            for (int i = 0; i < ThreadComputationBloc.NbThread; i++)
            {
                Thread thread = new Thread(new ParameterizedThreadStart(Submain_SABRVolMC));
                th.Add(thread);
                thread.Start(i);
            }
            foreach (var thread in th)
            {
                thread.Join();
            }
            ThreadComputationBloc.SaveToCVS_SABRVolMC();
            DateTime end = DateTime.Now;
            Console.WriteLine("Computation time (sec): " + (end - start).TotalSeconds);
            Console.WriteLine("Computation time  poun un chunk de prix (sec): " + (end - start).TotalSeconds / ThreadComputationBloc.NbThread);
            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();
        }
        static void Main_Sabr(string[] args)     //ExtendedSABR
        {
            ThreadComputationBloc.NbThread = 5;
            DateTime start = DateTime.Now;
            for (int iserie = 0; iserie < 100; iserie++)
            {
                ThreadComputationBloc.SetUp_SABR(iserie);
                List<Thread> th = new List<Thread>();
                for (int i = 0; i < ThreadComputationBloc.NbThread; i++)
                {
                    Thread thread = new Thread(new ParameterizedThreadStart(Submain_SABRVolMC));
                    th.Add(thread);
                    thread.Start(i);
                }
                foreach (var thread in th)
                {
                    thread.Join();
                }
                ThreadComputationBloc.SaveToCVS_SABRVolMC();
            }
            DateTime end = DateTime.Now;
            Console.WriteLine("Computation time (sec): " + (end - start).TotalSeconds);
            Console.WriteLine("Computation time  poun un chunk de prix (sec): " + (end - start).TotalSeconds / ThreadComputationBloc.NbThread);
            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();
        }
        static void Main_Lookback(string[] args) // Pas de Setup (directement dans le main)
        {
            double[,] correlationMatrix = new double[3, 3];
            double[] Svect = new double[3];
            double[] muvect = new double[3];
            double[] sigmavect = new double[3];
            double S1 = 71.063;
            double S2 = 70.7741;
            double S3 = 133.3036;
            double mu1 = 0.0074;
            double mu2 = 0.0343;
            double mu3 = 0.0049;
            double sigma1 = 0.1558;
            double sigma2 = 0.2052;
            double sigma3 = 0.2748;
            double rho12 = 0.5566;
            double rho13 = 0.7757;
            double rho23 = 0.4006;
            double PDIStrike = 59.7044;
            double PDIType = 0.2062;
            double maturity = 4.2224;
            double[] forwards;
            int nbDate = 7;
            int nbTraj = 100000;
            S1 = 91.41581521;
            S2 = 109.0606908;
            S3 = 100.8166408;
            mu1 = 0.01539818;
            mu2 = 0.018952227;
            mu3 = 0.041056664;
            sigma1 = 0.212203208;
            sigma2 = 0.15366158;
            sigma3 = 0.1717465;
            rho12 = 0.85864402;
            rho13 = 0.793020666;
            rho23 = 0.991187713;
            PDIStrike = 87.38449682;
            PDIType = 0.969624225;
            maturity = 1.071456235;
            nbDate = 12;

            Svect[0] = S1; Svect[1] = S2; Svect[2] = S3;
            muvect[0] = mu1; muvect[1] = mu2; muvect[2] = mu3;
            sigmavect[0] = sigma1; sigmavect[1] = sigma2; sigmavect[2] = sigma3;
            correlationMatrix[0, 0] = 1.0; correlationMatrix[0, 1] = rho12; correlationMatrix[0, 2] = rho13;
            correlationMatrix[1, 0] = rho12; correlationMatrix[1, 1] = 1.0; correlationMatrix[1, 2] = rho23;
            correlationMatrix[2, 0] = rho13; correlationMatrix[2, 1] = rho23; correlationMatrix[2, 2] = 1.0;

            Random random = new MersenneTwister(41);


            double price = pricing.MonteCarlo.LookbackCallWorstPriceMC(random, Svect, muvect, sigmavect,

                    correlationMatrix, maturity,
                    PDIStrike, PDIType,
                    nbDate, nbTraj, out forwards);

            Console.WriteLine("price=");
            Console.WriteLine(price);
            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();
        }
        static void Main_Ridounet(string[] args)        //Test Hagan - Montevarlo
        {
            double F0 = 1.0;
            double alpha0 = 0.5;
            double beta = 0.5;
            double beta2 = 0.5;
            double d = 0.0;
            double gamma = 1;
            double nu = 0.8;
            double lambda = 0;
            double omega = 0;
            double rho = -0.5;
            double Tmax = 3;
            int nbdate = 1000;
            int milsteinflag = 0;
            int expmartingaleflag = 0;
            int nbTraj = 5000;
            double[,] CallVs;
            double[] strikes = { 0.8, 0.9, 1.0, 1.1, 1.2 };

            //double[] strikes = SABR.generate_strikes(F0, Tmax, vol_atm, n_gauche, n_droite, delta_call, delta_put , epsilon );

            double[] volatilities = new double[strikes.Length];
            var watch = System.Diagnostics.Stopwatch.StartNew();
            CallVs = pricing.MonteCarlo.VolVanilleSetVolStoMC(F0, alpha0, beta, beta2, d, gamma, nu, lambda, omega, rho, Tmax, nbdate, strikes, nbTraj, milsteinflag, expmartingaleflag);
            watch.Stop();
            var elapsedMs = watch.ElapsedMilliseconds;

            ridounet.SABRFormula sa = new ridounet.SABRFormula(Tmax, F0, alpha0, beta, nu, rho, ridounet.SABRFormula.DiffusionType.LogNormal, 1.0, 1.0);

            double option11 = pricing.MonteCarlo.CallVanilleVolStoMC(F0, alpha0, beta, beta2, d, gamma, nu, lambda, omega, rho, Tmax, nbdate, 1.1, nbTraj, milsteinflag, expmartingaleflag);
            double vol11 = auxiliary.BlackScholes.ImpVolBS(F0, 1.1, Tmax, option11);
            sa.get_volatility(strikes, volatilities);

            double[,] new_vol_Ridounet = pricing.SABR.compute_sabr_vanille(Tmax, F0, alpha0, beta, nu, rho, strikes, nbTraj, nbdate);

            double[] volHagan = new double[strikes.Length];
            for (int i = 0; i < strikes.Length; i++)
            {
                volHagan[i] = pricing.SABR.Hagan_Formula(F0, alpha0, beta, rho, nu, strikes[i], Tmax);
            }
            double[,] opt;
            opt = pricing.SABR.vanilla_call(Tmax, F0, alpha0, beta, nu, rho, strikes, nbdate, nbTraj, milsteinflag, expmartingaleflag);
            Console.WriteLine("strike = VolVanilleSetVolStoMC/ ridounet.SABRFormula / volHagan/ MonteCarlo.CallVanilleVolStoMC /  SABR.vanilla_call :Rida/ SABR.vanilla_call:oliv/compute_sabr_vanille");
            for (int k = 0; k < strikes.Length; k++)
            {
                Console.WriteLine(String.Concat(Convert.ToString(strikes[k]), " = ", Convert.ToString(CallVs[1, k]), " / ", Convert.ToString(volatilities[k]), " / ", Convert.ToString(volHagan[k]), " / ",
                    Convert.ToString(vol11), " / ", Convert.ToString(opt[2, k]), " / ", Convert.ToString(opt[5, k]), " / ", Convert.ToString(new_vol_Ridounet[2, k])));
            }
            Console.WriteLine(String.Concat(" time (ms) = ", Convert.ToString(elapsedMs)));
            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();

        }
        static void Main_YetiPhoenixMC_1dim(string[] args)   //YetiPhoenixMC_1dim
        {
            Console.Out.WriteLine("Thread number :");
            string ligne = Console.In.ReadLine();
            ThreadComputationBloc.NbThread = int.Parse(ligne);
            DateTime start = DateTime.Now;
            ThreadComputationBloc.SetUp_YetiPhoenix_1dim();
            List<Thread> th = new List<Thread>();
            for (int i = 0; i < ThreadComputationBloc.NbThread; i++)
            {
                Thread thread = new Thread(new ParameterizedThreadStart(Submain_YetiPhoenixMC_1dim));
                th.Add(thread);
                thread.Start(i);
            }
            foreach (var thread in th)
            {
                thread.Join();
            }
            ThreadComputationBloc.SaveToCVS_YetiPhoenixPriceMC_1dim();


            DateTime end = DateTime.Now;
            Console.WriteLine("Computation time (sec): " + (end - start).TotalSeconds);
            Console.WriteLine("Computation time  poun un chunk de prix (sec): " + (end - start).TotalSeconds / ThreadComputationBloc.NbThread);

            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();
        }
        static void Main_YetiPhoenixMC_1dim_VolLoc(string[] args)   //YetiPhoenixMC_1dim_VolLoc
        {
            int maxbatch = 5;
            Console.Out.WriteLine("Thread number :");
            string ligne = Console.In.ReadLine();
            ThreadComputationBloc.NbThread = int.Parse(ligne);
            for (int ibatch = 0; ibatch < maxbatch; ibatch++)
            {
                
                DateTime start = DateTime.Now;
                ThreadComputationBloc.SetUp_YetiPhoenix_1dim_VolLoc();
                for (int ithread = 0; ithread < ThreadComputationBloc.NbThread; ithread++)
                {
                    for (int istrike = 0; istrike < ThreadComputationBloc.nbstrikes; istrike++)
                    {
                        double x = ThreadComputationBloc.InitialStrikeList[istrike];
                        ThreadComputationBloc.StrikeList[ithread, istrike] = x;
                    }
                }
                List<Thread> th = new List<Thread>();
                for (int i = 0; i < ThreadComputationBloc.NbThread; i++)
                {
                    Thread thread = new Thread(new ParameterizedThreadStart(Submain_YetiPhoenixMC_1dim_VolLoc));
                    th.Add(thread);
                    thread.Start(i);
                }

                foreach (var thread in th)
                {
                    thread.Join();
                }

                ThreadComputationBloc.SaveToCVS_YetiPhoenixPriceMC_1dim_VolLoc(ibatch);
                DateTime end = DateTime.Now;
                Console.WriteLine("Computation time (sec): " + (end - start).TotalSeconds);
                Console.WriteLine("Computation time  poun un chunk de prix (sec): " + (end - start).TotalSeconds / ThreadComputationBloc.NbThread);

            }

            Console.WriteLine("Finished ...Press any key to exit.");
            Console.ReadLine();
        }

        static void Main(string[] args)   //YetiPhoenixMC_1dim
        {
            //Main_YetiPhoenixMC_1dim(args);
            Main_YetiPhoenixMC_1dim_VolLoc(args);
        }

    }
}
