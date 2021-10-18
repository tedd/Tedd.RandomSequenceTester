namespace Tedd
{
    public class RandomSequenceTestResult
    {
        public long Total { get; set; }
        public double Entropy { get; set; }
        public double OptimalCompressionPct { get; set; }
        public double Mean { get; set; }
        public double MeanRandom { get; set; }
        public double MonteCarloPI { get; set; }
        public double MonteCarloPIErrorPct { get; set; }
        public double ChiSquared { get; set; }
        public double ChiSquaredProbability { get; set; }
        public double SerialCorrelationCoefficien { get; set; }
        public string Text { get; set; }
    }
}
