using System;
using System.Collections.Generic;
using System.IO;
using System.Security.Cryptography;
using System.Text;

namespace Tedd
{
    public class RandomSequenceTester
    {
        private readonly RandTest _randTest;

        public RandomSequenceTester()
        {
            _randTest = new RandTest();

            //var ccount = new long[256];
            //long totalc = 0;

        }

        /// <summary>
        /// Add data.
        /// </summary>
        /// <param name="bytes"></param>
        public void Add(Span<byte> bytes) => _randTest.Add(bytes);

        /// <summary>
        /// Add data from file in 32KB chunks.
        /// </summary>
        /// <param name="filename">Filename to read</param>
        public void AddFile(string filename)
        {
            using (var source = File.OpenRead(filename))
            {

                byte[] buffer = new byte[1024 * 32];
                int bytesRead;
                while ((bytesRead = source.Read(buffer, 0, buffer.Length)) > 0)
                    Add(((Span<byte>)buffer).Slice(0, bytesRead));
            }
        }

        public RandomSequenceTestResult Finish()
        {
            var sb = new StringBuilder();
            var samp = "byte";
            var result = _randTest.Finish();

            var chip = Chisq.pochisq(result.chisq, 255);
            sb.AppendLine($"Entropy = {result.ent} bits per {samp}.");
            var compPct = (short)((100 * (8 - result.ent) / 8.0D));
            sb.AppendLine($"\nOptimum compression would reduce the size of this {result.totalc} {samp} file by {compPct} percent.");
            sb.AppendLine();
            sb.Append($"Chi square distribution for {result.totalc} samples is {Math.Round(result.chisq, 2)}, and randomly ");
            if (chip < 0.0001)
            {
                sb.AppendLine("would exceed this value less than 0.01 percent of the times.");
            }
            else if (chip > 0.9999)
            {
                sb.AppendLine("would exceed this value more than than 99.99 percent of the times.");
            }
            else
            {
                sb.AppendLine($"would exceed this value {Math.Round(chip * 100, 2)} percent of the times.");
            }
            sb.AppendLine();
            var meanRandom = 127.5;
            sb.AppendLine($"Arithmetic mean value of data {samp} is {Math.Round(result.mean, 4)} ({Math.Round(meanRandom, 2)} = random).");
            var montecarloErrorPct = 100.0 * (Math.Abs(Math.PI - result.montepi) / Math.PI);
            sb.AppendLine($"Monte Carlo value for Pi is {Math.Round(result.montepi, 9)} (error {Math.Round(montecarloErrorPct, 2)} percent).");
            sb.Append("Serial correlation coefficient is ");
            if (result.scc >= -99999)
            {
                sb.AppendLine($"{Math.Round(result.scc, 6)} (totally uncorrelated = 0.0).\n");
            }
            else
            {
                sb.AppendLine("undefined (all values equal!).\n");
            }

            return new RandomSequenceTestResult()
            {
                Total = result.totalc,
                Entropy = result.ent,
                OptimalCompressionPct = compPct,
                Mean = result.mean,
                MeanRandom = meanRandom,
                MonteCarloPI = result.montepi,
                MonteCarloPIErrorPct = montecarloErrorPct,
                ChiSquared = result.chisq,
                ChiSquaredProbability = chip,
                SerialCorrelationCoefficien = result.scc,
                Text = sb.ToString()
            };
        }
    }
}
