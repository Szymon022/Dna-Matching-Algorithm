using System;
using System.Linq;
using System.Text;
using ASD;

namespace Lab2
{

    public abstract class DnaMatchingTestCase : TestCase
    {
        protected readonly string seq1;
        protected readonly string seq2;
        protected int matchValue;
        protected int mismatchValue;
        protected int gapStartValue;
        protected int gapContinuationValue;
        private bool checkStrings;
        private readonly int expectedVal;

        protected (string res1, string res2, int val) result;

        protected DnaMatchingTestCase(string seq1, string seq2, int expectedVal, bool checkStrings, double timeLimit, string description) : base(timeLimit, null, description)
        {
            this.seq1 = seq1;
            this.seq2 = seq2;
            this.expectedVal = expectedVal;
            this.checkStrings = checkStrings;
        }

        protected override (Result resultCode, string message) VerifyTestCase(object settings)
        {
            var (code, msg) = checkSolution(result.res1, result.res2, result.val);
            return (code, $"{msg} [{this.Description}]");
        }

        protected (Result resultCode, string message) checkSolution(string res1, string res2, int returnedValue)
        {
            if (expectedVal != returnedValue)
            {
                return (Result.WrongResult, $"Zwrócona wartość ({returnedValue}) jest inna niż oczekiwana ({expectedVal})");
            }

            if (!checkStrings)
            {
                return OkResult("OK");
            }

            if (res1 == null || res2 == null)
            {
                return (Result.WrongResult, "Zwrócono null jako jeden z ciągów dopasowania, w tym etapie należy zwrócić stringi");
            }

            if (res1.Length != res2.Length)
            {
                return (Result.WrongResult, $"Oba zwrócone ciągi powinny mieć taką samą długość (mają {res1.Length} i {res2.Length}");
            }

            var possibleChars = "ATCG-".ToCharArray();

            var res = 0;
            for (int i = 0; i < res1.Length; i++)
            {
                if (!possibleChars.Contains(res1[i]) || !possibleChars.Contains(res2[i]))
                {
                    return (Result.WrongResult, $"Znaleziono niedozwolony znak w zwróconych stringach. Powinny one składać się jedynie ze znaków \"ACGT-\"");
                }
                if (res1[i] == res2[i])
                {
                    res += matchValue;
                }
                else if (res1[i] == '-' || res2[i] == '-')
                {
                    if (i > 0 && (res1[i - 1] == '-' || res2[i - 1] == '-'))
                    {
                        res += gapContinuationValue;
                    }
                    else
                    {
                        res += gapStartValue;
                    }
                }
                else
                {
                    res += mismatchValue;
                }
            }

            if (res != expectedVal)
            {
                return (Result.WrongResult, $"Zwrócone dopasowanie ma inną wartość ({res}) niż zwrócona wartość ({returnedValue})");
            }

            if (res1.Replace("-", "") != seq1 || res2.Replace("-", "") != seq2)
            {
                return (Result.WrongResult, "Zwrócone dopasowanie ma inną kolejność znaków niż wejściowe sekwencje");
            }

            return OkResult("OK");
        }

        public (Result resultCode, string message) OkResult(string message) => (TimeLimit < PerformanceTime ? Result.LowEfficiency : Result.Success, $"{message} {PerformanceTime.ToString("#0.00")}s");
    }


    public class V1TestCase : DnaMatchingTestCase
    {
        public V1TestCase(string seq1, string seq2, int expectedVal, bool checkStrings, double timeLimit, string description) : base(seq1, seq2, expectedVal, checkStrings, timeLimit, description)
        {
            matchValue = 1;
            mismatchValue = -3;
            gapStartValue = -2;
            gapContinuationValue = -2;
        }

        protected override void PerformTestCase(object prototypeObject)
        {
            result = ((DnaMatching)prototypeObject).FindMatchingV1(seq1, seq2);
        }
    }

    public class V2TestCase : DnaMatchingTestCase
    {
        public V2TestCase(string seq1, string seq2, int expectedVal, bool checkStrings, double timeLimit, string description) : base(seq1, seq2, expectedVal, checkStrings, timeLimit, description)
        {
            matchValue = 1;
            mismatchValue = -3;
            gapStartValue = -5;
            gapContinuationValue = -2;
        }

        protected override void PerformTestCase(object prototypeObject)
        {
            result = ((DnaMatching)prototypeObject).FindMatchingV2(seq1, seq2);
        }
    }

    public class DnaMatchingTests : TestModule
    {
        TestSet v1Stage1 = new TestSet(prototypeObject: new DnaMatching(), description: "Wariant I, sama wartość liczbowa - etap 1", settings: true);
        TestSet v1Stage2 = new TestSet(prototypeObject: new DnaMatching(), description: "Wariant I, całe dopasowanie - etap 2", settings: true);
        TestSet v2Stage3 = new TestSet(prototypeObject: new DnaMatching(), description: "Wariant II, sama wartość liczbowa - etap 3", settings: true);
        TestSet v2Stage4 = new TestSet(prototypeObject: new DnaMatching(), description: "Wariant II, całe dopasowanie - etap 4", settings: true);

        public override void PrepareTestSets()
        {
            TestSets["V1Stage1"] = v1Stage1;
            TestSets["V1Stage2"] = v1Stage2;
            TestSets["V2Stage3"] = v2Stage3;
            TestSets["V2Stage4"] = v2Stage4;

            prepareStage1();
            prepareStage2();
            prepareStage3();
            prepareStage4();
        }

        private void addV1Stage1(V1TestCase v1TestCase)
        {
            v1Stage1.TestCases.Add(v1TestCase);
        }

        private void addV1Stage2(V1TestCase v1TestCase)
        {
            v1Stage2.TestCases.Add(v1TestCase);
        }

        private void addV2Stage3(V2TestCase v2TestCase)
        {
            v2Stage3.TestCases.Add(v2TestCase);
        }

        private void addV2Stage4(V2TestCase v2TestCase)
        {
            v2Stage4.TestCases.Add(v2TestCase);
        }

        private string getRandomDna(Random rand, int length)
        {
            var chars = "ACTG".ToCharArray();
            var result = new StringBuilder();

            for (int i = 0; i < length; i++)
            {
                result.Append(chars[rand.Next(0, 4)]);
            }

            return result.ToString();
        }

        private void prepareStage1()
        {
            addV1Stage1(new V1TestCase("GATAC", "GTCCAG", expectedVal: -6, timeLimit: 1, description: "Przykład z zadania", checkStrings: false));
            addV1Stage1(new V1TestCase("AAATTATATGC", "AAATTATATGC", expectedVal: 11, timeLimit: 1, description: "Dwie identyczne sekwencje", checkStrings: false));
            addV1Stage1(new V1TestCase("A", "A", expectedVal: 1, timeLimit: 1, description: "Jeden identyczny znak", checkStrings: false));
            addV1Stage1(new V1TestCase("A", "T", expectedVal: -3, timeLimit: 1, description: "Dwa różne znaki", checkStrings: false));
            addV1Stage1(new V1TestCase("AAGTGC", "ACTAGCTAG", expectedVal: -5, timeLimit: 1, description: "Prosty test 1", checkStrings: false));
            addV1Stage1(new V1TestCase("TTAGCTGA", "AAAGCTTAAAG", expectedVal: -10, timeLimit: 1, description: "Prosty test 2", checkStrings: false));
            addV1Stage1(new V1TestCase("CGAC", "A", expectedVal: -5, timeLimit: 1, description: "Prosty test 3", checkStrings: false));
            addV1Stage1(new V1TestCase(getRandomDna(new Random(23), 2000), getRandomDna(new Random(7), 1000), expectedVal: -1588, timeLimit: 3, description: "Długie losowe sekwencje 1", checkStrings: false));
            addV1Stage1(new V1TestCase(getRandomDna(new Random(24), 999), getRandomDna(new Random(8), 2001), expectedVal: -1554, timeLimit: 3, description: "Długie losowe sekwencje 2", checkStrings: false));
        }

        private void prepareStage2()
        {
            addV1Stage2(new V1TestCase("GATAC", "GTCCAG", expectedVal: -6, timeLimit: 1, description: "Przykład z zadania", checkStrings: true));
            addV1Stage1(new V1TestCase("AAATTATATGC", "AAATTATATGC", expectedVal: 11, timeLimit: 1, description: "Dwie identyczne sekwencje", checkStrings: true));
            addV1Stage2(new V1TestCase("A", "A", expectedVal: 1, timeLimit: 1, description: "Jeden identyczny znak", checkStrings: true));
            addV1Stage2(new V1TestCase("A", "T", expectedVal: -3, timeLimit: 1, description: "Dwa różne znaki", checkStrings: true));
            addV1Stage2(new V1TestCase("AAGTGC", "ACTAGCTAG", expectedVal: -5, timeLimit: 1, description: "Prosty test 1", checkStrings: true));
            addV1Stage2(new V1TestCase("TTAGCTGA", "AAAGCTTAAAG", expectedVal: -10, timeLimit: 1, description: "Prosty test 2", checkStrings: true));
            addV1Stage2(new V1TestCase("CGAC", "A", expectedVal: -5, timeLimit: 1, description: "Prosty test 3", checkStrings: true));
            addV1Stage2(new V1TestCase(getRandomDna(new Random(23), 2000), getRandomDna(new Random(7), 1000), expectedVal: -1588, timeLimit: 3, description: "Długie losowe sekwencje 1", checkStrings: true));
            addV1Stage2(new V1TestCase(getRandomDna(new Random(24), 999), getRandomDna(new Random(8), 2001), expectedVal: -1554, timeLimit: 3, description: "Długie losowe sekwencje 2", checkStrings: true));
        }

        private void prepareStage3()
        {
            addV2Stage3(new V2TestCase("GATAC", "GTCCAG", expectedVal: -12, timeLimit: 1, description: "Przykład z zadania 1", checkStrings: false));
            addV2Stage3(new V2TestCase("AAATTATATGC", "AAATTATATGC", expectedVal: 11, timeLimit: 1, description: "Dwie identyczne sekwencje", checkStrings: false));
            addV2Stage3(new V2TestCase("CAAAAAATTTTTG", "CAAATTTG", expectedVal: -5, timeLimit: 1, description: "Przykład z zadania 2", checkStrings: false));
            addV2Stage3(new V2TestCase("A", "A", expectedVal: 1, timeLimit: 1, description: "Jeden identyczny znak", checkStrings: false));
            addV2Stage3(new V2TestCase("A", "T", expectedVal: -3, timeLimit: 1, description: "Dwa różne znaki", checkStrings: false));
            addV2Stage3(new V2TestCase("AAGTGC", "ACTAGCTAG", expectedVal: -14, timeLimit: 1, description: "Prosty test 1", checkStrings: false));
            addV2Stage3(new V2TestCase("TTAGCTGA", "AAAGCTTAAAG", expectedVal: -13, timeLimit: 1, description: "Prosty test 2", checkStrings: false));
            addV2Stage3(new V2TestCase("CGAC", "A", expectedVal: -11, timeLimit: 1, description: "Prosty test 3", checkStrings: false));
            addV2Stage3(new V2TestCase(getRandomDna(new Random(23), 1400), getRandomDna(new Random(7), 700), expectedVal: -1759, timeLimit: 3, description: "Długie losowe sekwencje 1", checkStrings: false));
            addV2Stage3(new V2TestCase(getRandomDna(new Random(24), 701), getRandomDna(new Random(8), 1401), expectedVal: -1723, timeLimit: 3, description: "Długie losowe sekwencje 2", checkStrings: false));
        }

        private void prepareStage4()
        {
            addV2Stage4(new V2TestCase("GATAC", "GTCCAG", expectedVal: -12, timeLimit: 1, description: "Przykład z zadania 1", checkStrings: true));
            addV2Stage3(new V2TestCase("AAATTATATGC", "AAATTATATGC", expectedVal: 11, timeLimit: 1, description: "Dwie identyczne sekwencje", checkStrings: true));
            addV2Stage4(new V2TestCase("CAAAAAATTTTTG", "CAAATTTG", expectedVal: -5, timeLimit: 1, description: "Przykład z zadania 2", checkStrings: true));
            addV2Stage4(new V2TestCase("A", "A", expectedVal: 1, timeLimit: 1, description: "Jeden identyczny znak", checkStrings: true));
            addV2Stage4(new V2TestCase("A", "T", expectedVal: -3, timeLimit: 1, description: "Dwa różne znaki", checkStrings: true));
            addV2Stage4(new V2TestCase("AAGTGC", "ACTAGCTAG", expectedVal: -14, timeLimit: 1, description: "Prosty test 1", checkStrings: true));
            addV2Stage4(new V2TestCase("TTAGCTGA", "AAAGCTTAAAG", expectedVal: -13, timeLimit: 1, description: "Prosty test 2", checkStrings: true));
            addV2Stage4(new V2TestCase("CGAC", "A", expectedVal: -11, timeLimit: 1, description: "Prosty test 4", checkStrings: true));
            addV2Stage4(new V2TestCase(getRandomDna(new Random(23), 1400), getRandomDna(new Random(7), 700), expectedVal: -1759, timeLimit: 3, description: "Długie losowe sekwencje 1", checkStrings: true));
            addV2Stage4(new V2TestCase(getRandomDna(new Random(24), 701), getRandomDna(new Random(8), 1401), expectedVal: -1723, timeLimit: 3, description: "Długie losowe sekwencje 2", checkStrings: true));
        }
    }


    class Program
    {
        static void Main(string[] args)
        {
            var tests = new DnaMatchingTests();
            tests.PrepareTestSets();
            foreach (var ts in tests.TestSets)
            {
                ts.Value.PerformTests(verbose: true, checkTimeLimit: false);
            }
            Console.WriteLine("");
        }
    }
}
