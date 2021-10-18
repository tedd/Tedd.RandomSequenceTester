using System;
using System.Runtime.CompilerServices;

namespace Tedd
{
    internal sealed class RandTest : IDisposable
    {
        private const double log2of10 = 3.32192809488736234787D;

        private readonly long[] _charCount = new long[256];
        private long _totalBytes = 0;
        private const int MONTEN = 6;                               // Bytes used as Monte Carlo co-ordinates. This should be no more bits than the mantissa of double.
        private bool _sccFirst = true;
        private int _mp;
        private readonly uint[] _monte = new uint[MONTEN];
        private long _inmont, _mcount;
        private readonly double _incirc = Math.Pow(Math.Pow(256.0D, (double)(MONTEN / 2)) - 1, 2.0D); // 65535.0 * 65535.0;     // In-circle distance for Monte Carlo
        private double sccu0, scclast, scct1, scct2, scct3;

        /// <summary>
        /// Calculate log to the base 2
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double rt_log2(double x) => log2of10 * Math.Log10(x);

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
                _charCount[(byte)c]++;
                _totalBytes++;

                // Update inside / outside circle counts for Monte Carlo computation of PI.
                _monte[_mp++] = c;
                if (_mp >= MONTEN)
                {
                    // Calculate every MONTEN character
                    int mj;

                    _mp = 0;
                    _mcount++;
                    var montex = 0D;
                    var montey = 0D;
                    for (mj = 0; mj < MONTEN / 2; mj++)
                    {
                        montex = (montex * 256.0D) + _monte[mj];
                        montey = (montey * 256.0D) + _monte[(MONTEN / 2) + mj];
                    }
                    if ((montex * montex + montey * montey) <= _incirc)
                        _inmont++;
                }

                // Update serial correlation coefficiency
                var sccun = c;
                if (_sccFirst)
                {
                    _sccFirst = false;
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

            return (_totalBytes, ent, chisq, mean, montepi, scc);
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
                sum += ((double)i) * _charCount[i];
            }
            return sum / _totalBytes;
        }

        /// <summary>
        /// Calculate serial correlation coefficient
        /// </summary>
        /// <returns></returns>
        private double CalculateSerialCorrelationCoefficient()
        {
            var pscct1 = scct1 + scclast * sccu0;
            var pscct2 = scct2 * scct2;
            var scc = _totalBytes * scct3 - pscct2;
            if (scc == 0.0)
                scc = -100000;
            else
                scc = (_totalBytes * pscct1 - pscct2) / scc;
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
            var cexp = _totalBytes / (double)bmm;  // Expected count per bin
            for (var i = 0; i < bmm; i++)
            {
                var a = _charCount[i] - cexp; ;

                prob[i] = ((double)_charCount[i]) / _totalBytes;
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
        private double CalculateMonteCarloPI() => 4.0 * (((double)_inmont) / _mcount);

        public void Dispose()
        {

        }
    }
}