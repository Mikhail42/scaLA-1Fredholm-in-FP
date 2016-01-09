
import math._
import scala.Array._
import scala.collection.mutable.ArraySeq
  
import LAfunction._
import LAfunction.mathFun
  
class BigDecimalData(a: BigDecimal, b: BigDecimal, N: Int)
  extends Data[BigDecimal](a, b, N, mathFun.mySqrt, mathFun.cos, mathFun.div)
  
class FloatData(a: Float, b: Float, N: Int)  
  extends Data[Float]     (a, b, N, mathFun.mySqrt, mathFun.cos, mathFun.div)
  
class DoubleData(a: Double, b: Double, N: Int) 
  extends Data[Double]    (a, b, N, mathFun.mySqrt, mathFun.cos, mathFun.div)

class Data[T](
    a: T, b: T, N: Int, mySqrt: T => T, cos: T => T, div: (T, T) => T)
      (implicit num: Numeric[T])  {
  
  import BasicData._
  import num._
  
  val funct = new Function[T](mySqrt, cos, div)
    
  import funct._
  import funct.Matrix._
  import funct.Number._
  import funct.Vector._
  
  val nameSolution = "Решение интегрального уравнения \u222Bsec(ts) \u03D5(s)ds=4t"
  
  def createSupAAAs(n1: Int, n2: Int) = {
    val res = new Array[Array[ArraySeq[T]]](n1)
    for (i <- 0 until n1) res(i) = new Array[ArraySeq[T]](n2)
    res
  }
  
  def createAAS(m: Int, n: Int) = {
    val res = new Array[ArraySeq[T]](m)
    for (i <- 0 until m) res(i) = new ArraySeq[T](n)
    res
  }
  def createAAS(n: Int): Array[ArraySeq[T]] = createAAS(n,n)
  def createAr(n: Int) = new ArraySeq[T](n)
  
  object BasicData{
    val h = div(b - a, num.fromInt(N))
    val size = N+1
    println("size="+size)
    val name = nameSolution +"s:=" +a+":"+h+":"+b
    val xs = {
      val res = new ArraySeq[T](size)  
      res(0) = a
      for (i <- 1 until size) res(i) = res(i-1)+h
      res
    }
    val ts =  {
      val res = new ArraySeq[T](size)  
      res(0) = a
      for (i <- 1 until size) res(i) = res(i-1)+h
      res
    }
    
    val K = {
      val res = createAAS(size,size)
      for (i <- 0 until size; j <- 0 until size) res(i)(j) = getK(i,j)
      res
    }
  
    def f(t: T)         = four(t)
    def getK(i: Int, j: Int) = sec(xs(i)*ts(j))
 
    lazy val as = {
      val res = new ArraySeq[T](size)
      val oneDivH = inverse(h)
      val hDiv3   = div3(h)
      // Simpson's rule
      if (isDiv2(N)) {
        res(0) =  hDiv3; res(N) = hDiv3
        val fourHDiv3 = four(hDiv3)
        for (i <- 1 until N by 2) res(i) = fourHDiv3
        val twoHDiv3 = two(hDiv3)
        for (i <- 2 until N by 2) res(i) = twoHDiv3
      }
      // Trapezoidal rule
      else {
        val hDiv2 = div2(h); res(0) = hDiv2; res(N) = hDiv2
        for (i <- 1 until N) res(i) = h
      }      
      res
    }
 
    val nBestSqrAlphas = 4
    lazy val bestSqrAlphas = {
      val res  = createAr(nBestSqrAlphas)
      val div  = inverse(num.fromInt(10))
      val step = sqr(div)
      res(0) = div; for (i <- 1 until nBestSqrAlphas) res(i) = res(i-1)*step 
      res
    }      

    lazy val namesBeta = List("1", "\u03B1", "sqrt(\u03B1^2+0.5*\u03C3_m^2(A))")
    
    val nBetas = 3     // don't edit!!! see down    
    lazy val betas = 
      for (i <- 0 until nBetas) yield  
      for (j <- 0 until nBestSqrAlphas) yield getBeta(i,j)
    
    import SLAE._
    def getBeta(i: Int, j: Int): T = {
      if (i == 0)    num.one 
      else if (i==1) mySqrt(bestSqrAlphas(j))
      else           mySqrt((div2(sigma_m) + bestSqrAlphas(j))) 
    }
  }
  
  object Eignvalues{
    import BasicData._
    import SLAE._
    
    lazy val matEigns = {
      val res = new Array[ArraySeq[T]](nBestSqrAlphas)
      for (j <- 0 until nBestSqrAlphas) 
        try { res(j) = eigenvaluesA_k(ATAsWithAlpha(j))._1
        } catch { case e: Exception => {} }
      res
    }
          
    lazy val superMatEigns = 
      for (i <- range(0, nBetas))
        yield for (j <- range(0, nBestSqrAlphas))
          yield getEigns(bestSqrAlphas(j), betas(i)(j))
    
    private def getEigns(alpha: T, beta: T) = {
      val elem = -div(alpha, beta)
      def add(i: Int) = if (i<size) beta else elem
      val matr = SLAE.resultMat.clone()
      for (l <- 0 until twoSize) matr(l)(l) += add(l)
      val res = eigenvaluesA_k(matr)._1        
      res
    }
  }
  
  object Cond{
    import SLAE._
    import BasicData._
    def getCond(eigns: ArraySeq[T]) = {
      val maxMin = getMaxAndMin(eigns)
      div(maxMin._1, maxMin._2) 
    }
    
    lazy val condA = cond2(A)
    
    lazy val condATA = cond2(ATA)
    
    /**
     * using Alpha
     */
    lazy val vectCond = for (i <- 0 until nBestSqrAlphas) yield getCond(Eignvalues.matEigns(i))
    
    /**
     * using Alpha, Eignvalues (Alpha), Beta (constants), Cond (Alpha), SLAE (BasicData)
     */
    lazy val matCond =  {
      val res = createAAS(nBetas, nBestSqrAlphas)
      for (i <- 0 until nBetas)
        for (j <- 0 until nBestSqrAlphas)
            res(i)(j) = getCond(Eignvalues.superMatEigns(i)(j))  
      res
    }
         
  }
  
  /**
   * using BasicData, As 
   */
  object SLAE{
    
    import BasicData._
    
    val A = multiByString(K, as) 
    
    private val n = A.length
    val b = {
      val res = createAr(n)
      for (i <- 0 until n) res(i) = f(xs(i))
      res
    }
    
    /** for (A^T A + (alpha^2)*E)x = A^T b */
    lazy val ATA = multiATA(A)
    lazy val ATb = multiATb(A, b)
    lazy val sigma_m = getMaxAndMin(eigenvaluesA_k(SLAE.resultMat)._1)._2
    
    /** for
     * 	(b   A     ) (r/b) = (b)
     * 	(AT  -a^2/b) (x  ) = (0)
     */
    val twoSize = (size<<1)
    
    lazy val resultMat = {
      val res = createAAS(twoSize, twoSize)
      for (k <- 0 until size; l <- size until twoSize) res(k)(l) = A(k)(l-size)
      for (k <- size until twoSize; l <- 0 until size) res(k)(l) = A(l)(k-size)
      res
    }
    
    lazy val ATAsWithAlpha = {
      val res = new Array[Array[ArraySeq[T]]](nBestSqrAlphas)
      for (i <- 0 until nBestSqrAlphas){
        res(i) = ATA.clone() 
        plusDiag(res(i), bestSqrAlphas(i))
      }
      res
    }
    
    lazy val rightPart = b.clone()
    
    lazy val ATAsWithAlphaAndBeta = {      
      val res = new Array[Array[ArraySeq[T]]](nBetas)
      for (i <- 0 until nBetas) {
        res(i) = resultMat.clone
        val alpha = bestSqrAlphas(Solution.betstIndAlpha)
        val betaElem = betas(i)(if (i==0) 0 else Solution.betstIndAlpha)
        val eleme = -div(alpha, betaElem)
        for (l <- 0 until twoSize)       res(i)(l)(l) += betaElem
        for (l <- twoSize until twoSize) res(i)(l)(l) += eleme
      }
      res
    }
  }
  
  object Solution {
    import SLAE._
    import BasicData._
    
    val betstIndAlpha = 3 // 0 <= betstIndAlpha < nBestAlphas
    
    def solveATAsWithAlphaAndBeta(
        f: (Array[ArraySeq[T]], ArraySeq[T]) =>  ArraySeq[T]) = {
      val res = createAAS(nBetas, xs.length)
      for (i <- 0 until nBetas){
        val solve = f(ATAsWithAlphaAndBeta(i), rightPart)
        Array.copy(solve, xs.length-1, res(i), 0, xs.length)
      }
      res
    }
    
    def solveATAsWithAlpha(
        fun: (Array[ArraySeq[T]], ArraySeq[T]) =>  ArraySeq[T])  = {
      val res = createAAS(nBestSqrAlphas, size)
      for (i <- 0 until nBestSqrAlphas) res(i) = fun(ATAsWithAlpha(i), ATb)
      res
    }
    
    object Jacobi{
      lazy val solveATAsWithAlpha        = Solution.solveATAsWithAlpha(solveJacobiAndSeidel) 
      lazy val solveATAsWithAlphaAndBeta = Solution.solveATAsWithAlphaAndBeta(solveJacobiAndSeidel) 
    }
    
    object LU {
      lazy val solveATAsWithAlpha        = Solution.solveATAsWithAlpha(solveLU) 
      lazy val solveATAsWithAlphaAndBeta = Solution.solveATAsWithAlphaAndBeta(solveLU)      
    }
    
    object LUP {
      lazy val solveATAsWithAlpha        = Solution.solveATAsWithAlpha(solveLUP)
      lazy val solveATAsWithAlphaAndBeta = Solution.solveATAsWithAlphaAndBeta(solveLUP) 
    }
    
    object QR {
      lazy val solveATAsWithAlpha        = Solution.solveATAsWithAlpha(solveQR)
      lazy val solveATAsWithAlphaAndBeta = Solution.solveATAsWithAlphaAndBeta(solveQR) 
    }
  }
}
