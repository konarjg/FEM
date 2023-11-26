namespace FEM_Testing
{
    [TestClass]
    public class FEMTesting
    {
        [TestMethod]
        public void TestLTOG()
        {
            var exact = new int[,] { { 0, 1, 2 }, { 1, 2, 3 } };
            var calculated = FEMSolver1D.LTOG(4, 2);

            CollectionAssert.AreEqual(exact, calculated);
        }

        [TestMethod]
        public void TestRPNAddSpaces()
        {
            var expression = "-x+-x+1^2*2";
            var spaces = expression.AddSpaces();

            Assert.AreEqual(spaces, "0 - x + 0 - x + 1 ^ 2 * 2");
        }

        [TestMethod]
        public void TestRPNParse()
        {
            var expression = new Expression("3+4*2/(1-5)^2");
            var parsed = "3 4 2 * 1 5 - 2 ^ / +";

            Assert.AreEqual(expression.ToString(), parsed);
        }
    }
}