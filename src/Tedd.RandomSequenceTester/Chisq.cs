using System;
using System.Buffers.Text;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;

namespace Tedd
{
    /// <summary>
    /// Compute probability of measured Chi Square value.
    /// </summary>
    internal sealed class Chisq
    {
        // Gary Perlman, Wang Institute, Tyngsboro, MA 01879

        private const double Z_MAX = 6.0;

        /// <summary>
        /// probability of normal z value.
        /// </summary>
        /// <remarks>This routine has six digit accuracy, so it is only useful for absolute z values < 6.  For z values >= to 6.0, poz() returns 0.0.</remarks>
        /// <param name="z"></param>
        /// <returns>returns cumulative probability from -oo to z</returns>
        private static double poz(double z)  /*VAR normal z value */
        {
            // Adapted from a polynomial approximation in:
            // Ibbetson D, Algorithm 209
            // Collected Algorithms of the CACM 1963 p. 616
            double y, x, w;

            if (z == 0.0D)
            {
                x = 0.0D;
            }
            else
            {
                y = 0.5D * Math.Abs(z);
                if (y >= (Z_MAX * 0.5D))
                {
                    x = 1.0D;
                }
                else if (y < 1.0D)
                {
                    w = y * y;
                    x = ((((((((0.000124818987D * w
                        - 0.001075204047D) * w + 0.005198775019D) * w
                        - 0.019198292004D) * w + 0.059054035642D) * w
                        - 0.151968751364D) * w + 0.319152932694D) * w
                        - 0.531923007300D) * w + 0.797884560593D) * y * 2.0D;
                }
                else
                {
                    y -= 2.0;
                    x = (((((((((((((-0.000045255659D * y
                        + 0.000152529290D) * y - 0.000019538132D) * y
                        - 0.000676904986D) * y + 0.001390604284D) * y
                        - 0.000794620820D) * y - 0.002034254874D) * y
                        + 0.006549791214D) * y - 0.010557625006D) * y
                        + 0.011630447319D) * y - 0.009279453341D) * y
                        + 0.005353579108D) * y - 0.002141268741D) * y
                        + 0.000535310849D) * y + 0.999936657524D;
                }
            }

            return (z > 0.0D ? ((x + 1.0D) * 0.5D) : ((1.0D - x) * 0.5D));
        }

        private const double LOG_SQRT_PI = 0.5723649429247000870717135D; // log (sqrt (pi))
        private const double I_SQRT_PI = 0.5641895835477562869480795D;   // 1 / sqrt (pi)
        private const double BIGX = 20.0D;                               // max value to represent exp (x)
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double ex(double x) => (((x) < -BIGX) ? 0.0D : Math.Exp(x));

        /// <summary>
        /// Compute probability of chi square value.
        /// </summary>
        /// <param name="ax">Obtained chi-square value</param>
        /// <param name="df">Degrees of freedom</param>
        /// <returns></returns>
        public static double pochisq(double ax, int df)
        {
            // Adapted from:
            // Hill, I.D.and Pike, M. C.Algorithm 299
            //  Collected Algorithms for the CACM 1967 p. 243
            // Updated for rounding errors based on remark in
            //    ACM TOMS June 1985, page 185

            double x = ax;
            double a, y = 0D, s;
            double e, c, z;
            bool even;               /* true if df is an even number */

            if (x <= 0.0D || df < 1D)
            {
                return 1.0D;
            }

            a = 0.5D * x;
            even = (2 * (df / 2)) == df;
            if (df > 1)
            {
                y = ex(-a);
            }
            s = (even ? y : (2.0D * poz(-Math.Sqrt(x))));
            if (df > 2)
            {
                x = 0.5D * (df - 1.0D);
                z = (even ? 1.0D : 0.5D);
                if (a > BIGX)
                {
                    e = (even ? 0.0D : LOG_SQRT_PI);
                    c = Math.Log(a);
                    while (z <= x)
                    {
                        e = Math.Log(z) + e;
                        s += ex(c * z - a - e);
                        z += 1.0D;
                    }
                    return (s);
                }
                else
                {
                    e = (even ? 1.0D : (I_SQRT_PI / Math.Sqrt(a)));
                    c = 0.0D;
                    while (z <= x)
                    {
                        e = e * (a / z);
                        c = c + e;
                        z += 1.0D;
                    }
                    return (c * y + s);
                }
            }
            else
            {
                return s;
            }
        }

    }
}
