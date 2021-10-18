using System;
using System.Runtime.InteropServices;

using Xunit;
using FluentAssertions;
using Xunit.Abstractions;
using Xunit.Sdk;
using System.IO;

namespace Tedd.RandomSequenceTester_.Tests
{
    public class RandomTests
    {
        private readonly ITestOutputHelper _testOutputHelper;

        public RandomTests(ITestOutputHelper testOutputHelper)
        {
            _testOutputHelper = testOutputHelper;
        }


        [Fact]
        public void TestSystemRandom()
        {
            var rnd = new Random();
            var tester = new global::Tedd.RandomSequenceTester();

            //var ints = new uint[1_000_000];
            //var bytes = MemoryMarshal.Cast<uint, byte>(ints);
            var bytes = new byte[1_000_000];

            for (var i = 0; i < 1_000; i++)
            {
                rnd.NextBytes(bytes);
                tester.Add(bytes);
            }

            var result = tester.Finish();

            _testOutputHelper.WriteLine(result.Text);

            result.Entropy.Should().BeGreaterThan(7.99D);
            result.OptimalCompressionPct.Should().BeLessThan(0.001D);
            result.ChiSquaredProbability.Should().BeGreaterThan(0.20D);
            Math.Abs(result.Mean - result.MeanRandom).Should().BeLessThan(0.1D);
            result.MonteCarloPIErrorPct.Should().BeLessThan(0.2D);
            result.SerialCorrelationCoefficien.Should().BeGreaterOrEqualTo(-99999);
        }
        [Fact]
        public void TestFile()
        {
            var bytes = File.ReadAllBytes(@"D:\Temp\Random.bin");
            var tester = new global::Tedd.RandomSequenceTester();

            tester.Add(bytes);

            var result = tester.Finish();

            _testOutputHelper.WriteLine(result.Text);

            result.Entropy.Should().BeGreaterThan(7.99D);
            result.OptimalCompressionPct.Should().BeLessThan(0.001D);
            result.ChiSquaredProbability.Should().BeGreaterThan(0.20D);
            Math.Abs(result.Mean - result.MeanRandom).Should().BeLessThan(0.1D);
            result.MonteCarloPIErrorPct.Should().BeLessThan(0.2D);
            result.SerialCorrelationCoefficien.Should().BeGreaterOrEqualTo(-99999);
        }
    }
}