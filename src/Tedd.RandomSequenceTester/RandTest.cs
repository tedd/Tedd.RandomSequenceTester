using System;
using System.Runtime.CompilerServices;

namespace Tedd
{
    internal sealed class RandTest : IDisposable
    {
        private const double log2of10 = 3.32192809488736234787D;

        private long[] ccount = new long[256];          // Bins to count occurrences of values
        private long totalc = 0;                        // Total bytes counted

        /// <summary>
        /// Calculate log to the base 2
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double rt_log2(double x) => log2of10 * Math.Log10(x);

        private const int MONTEN = 6;             // Bytes used as Monte Carlo co-ordinates. This should be no more bits than the mantissa of your "double" floating point type.

        private bool sccfirst = true;
        private int mp;
        private uint[] monte = new uint[MONTEN];
        private long inmont, mcount;

        private double incirc = Math.Pow(Math.Pow(256.0D, (double)(MONTEN / 2)) - 1, 2.0D); // 65535.0 * 65535.0;     // In-circle distance for Monte Carlo
        private double montex, montey, sccun, sccu0, scclast, scct1, scct2, scct3;

        /// <summary>
        /// Add data.
        /// </summary>
        /// <param name="bytes"></param>
        [MethodImpl()]
        public void Add(Span<byte> bytes)
        {
            for (var i = 0; i < bytes.Length; i++)
            {
                var c = bytes[i];
                ccount[(byte)c]++;          /* Update counter for this bin */
                totalc++;

                /* Update inside / outside circle counts for Monte Carlo
                   computation of PI */

                monte[mp++] = c;       /* Save character for Monte Carlo */
                if (mp >= MONTEN)
                {     /* Calculate every MONTEN character */
                    int mj;

                    mp = 0;
                    mcount++;
                    montex = montey = 0;
                    for (mj = 0; mj < MONTEN / 2; mj++)
                    {
                        montex = (montex * 256.0D) + monte[mj];
                        montey = (montey * 256.0D) + monte[(MONTEN / 2) + mj];
                    }
                    if ((montex * montex + montey * montey) <= incirc)
                        inmont++;
                }

                /* Update calculation of serial correlation coefficient */

                sccun = c;
                if (sccfirst)
                {
                    sccfirst = false;
                    scclast = 0;
                    sccu0 = sccun;
                }
                else
                {
                    scct1 = scct1 + scclast * sccun;
                }
                scct2 = scct2 + sccun;
                scct3 = scct3 + (sccun * sccun);
                scclast = sccun;
            }
        }

        public (long totalc, double ent, double chisq, double mean, double montepi, double scc) Finish()
        {
            var mean = CalculateMean();
            var scc = CalculateSerialCorrelationCoefficient();
            var chisq = CalculateChiSquareDistribution(out var prob);
            var ent = CalculateEntropy(prob);
            var montepi = CalculateMonteCarloPI();

            return (totalc, ent, chisq, mean, montepi, scc);
        }

        /// <summary>
        /// Calculate mean.
        /// </summary>
        /// <returns></returns>
        private double CalculateMean()
        {
            double sum = 0D;
            var bmm = 256;
            for (var i = 0; i < bmm; i++)
            {
                sum += ((double)i) * ccount[i];
            }
            return sum / totalc;
        }

        /// <summary>
        /// Calculate serial correlation coefficient
        /// </summary>
        /// <returns></returns>
        private double CalculateSerialCorrelationCoefficient()
        {
            var pscct1 = scct1 + scclast * sccu0;
            var pscct2 = scct2 * scct2;
            var scc = totalc * scct3 - pscct2;
            if (scc == 0.0)
                scc = -100000;
            else
                scc = (totalc * pscct1 - pscct2) / scc;
            return scc;
        }

        /// <summary>
        /// Scan bins and calculate probability for each bin and Chi-Square distribution.
        /// The probability will be reused in the entropy calculation below.  
        /// While we're at it, we sum of all the data which will be used to compute the mean.
        /// </summary>
        /// <param name="prob"></param>
        /// <returns></returns>
        private double CalculateChiSquareDistribution(out double[] prob)
        {
            prob = new double[256];
            double chisq = 0;
            var bmm = 256;
            var cexp = totalc / (double)bmm;  /* Expected count per bin */
            for (var i = 0; i < bmm; i++)
            {
                var a = ccount[i] - cexp; ;

                prob[i] = ((double)ccount[i]) / totalc;
                chisq += (a * a) / cexp;
            }
            return chisq;
        }

        /// <summary>
        /// Calculate entropy.
        /// </summary>
        public double CalculateEntropy(double[] prob)
        {
            double ent = 0D;
            var bmm = 256;
            for (var i = 0; i < bmm; i++)
            {
                if (prob[i] > 0.0)
                {
                    ent += prob[i] * rt_log2(1 / prob[i]);
                }
            }
            return Math.Round(ent, 6);
        }

        /// <summary>
        /// Calculate Monte Carlo value for PI from percentage of hits within the circle.
        /// </summary>
        private double CalculateMonteCarloPI() => 4.0 * (((double)inmont) / mcount);

        public void Dispose()
        {
            throw new NotImplementedException();
        }
    }
}