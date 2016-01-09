package LAfunction

import scala.Array._
import scala.collection.mutable.ArraySeq
import scala.BigDecimal

/** see also comments in of end a file */
/** for bibliography, see my first project: NLA1 */
object mathFun  {
  def mySqrt(x: BigDecimal): BigDecimal = BigDecimal(math.sqrt(x.toDouble)) 
  def mySqrt(x: Float):      Float      = math.sqrt(x.toDouble).toFloat
  def mySqrt(x: Double):     Double     = math.sqrt(x)
  
  def cos(x: Double):     Double        = math.cos(x)  
  def cos(x: Float):      Float         = math.cos(x.toDouble).toFloat
  def cos(x: BigDecimal): BigDecimal    = BigDecimal(math.cos(x.toDouble))
   
  def div(x: BigDecimal, y: BigDecimal) = x/y
  def div(x: Float, y: Float)           = x/y
  def div(x: Double, y: Double)         = x/y 
}

object BigDecimalFunction extends Function[BigDecimal](mathFun.mySqrt, mathFun.cos, mathFun.div)
object FloatFunction      extends Function[Float]     (mathFun.mySqrt, mathFun.cos, mathFun.div)
object DoubleFunction     extends Function[Double]    (mathFun.mySqrt, mathFun.cos, mathFun.div)

class Function[T](mySqrt: T => T, cos: T => T, div: (T, T) => T)(implicit num: Numeric[T]) {
  
  import num._
  import Number._
  
  def zero(n: Int): ArraySeq[T] = new ArraySeq[T](n)
  def zero(m: Int, n: Int): Array[ArraySeq[T]] = {
    val res = new Array[ArraySeq[T]](m)
    for(i <- 0 until m) res(i) = zero(n)
    res
  }
  
  var error = 1e-15
  var QR_ERROR = 1e-4
  var QR_ITER  = 1000
  
  object Number{
    
    def sec(x: T)  = inverse(cos(x))   
    
    def sqr(x: T) = x*x
    
    def two(x: T)  = x+x
    def four(x: T) = two(x+x)
    def tree(x: T) = two(x)+x
    
    def difAbs(x: T, y: T) = (x-y).toDouble.abs
    
    def inverse(x: T): T = div(num.one, x)
    def div2(x: T) = div(x, two(num.one)) 
    private val oneDiv3 = inverse(tree(num.one))
    def div3(x: T): T = x*oneDiv3
    
    def compare(a: T, b: T) = 
      if ((a-b).abs.toDouble < error) 0 else if (a>b) 1 else -1
    
    def specSign(x: T): T = if (compare(x, num.zero)>=0) num.one else -num.one
    
    def isZero(x: T) = x.abs.toDouble < error
    
     /** x mod 2 == 0
     *  @param x is positive number
     */
    def isDiv2(x: Int) = (x&1)==0
    
  }
  
  import Vector._
  
  object Vector {
    
    def plus(a: ArraySeq[T], b: ArraySeq[T]): ArraySeq[T] = a.zip(b).map(f => f._1 + f._2)
    def minus(a: ArraySeq[T], b: ArraySeq[T]): ArraySeq[T] = a zip b map Function.tupled(_-_)
    def multi(a: ArraySeq[T], f: T): ArraySeq[T] = a.map { x => x*f }
        
    def product(a: Array[T]): T = a.product
    
    /** a(i1:i2-1).*b(i1:i2-1) = sum(ai*bi : i in [i1,i2)) */
    def scalarPart(a: ArraySeq[T], b: ArraySeq[T], i1: Int, i2: Int): T = {
      var sum: T = num.zero; for (i <- i1 until i2) sum += a(i)*b(i); sum
    }
    
    /** a(i1:i2-1).*B(i1:i2-1)(j) */
    def scalarPart(a: ArraySeq[T], B: Array[ArraySeq[T]], j: Int, i1: Int, i2: Int): T = {
      var sum: T = num.zero; for (i <- i1 until i2) sum += a(i)*B(i)(j); sum
    }
    
    def scalar(a: ArraySeq[T], b: ArraySeq[T]): T = scalarPart(a, b, 0, a.length)
    
    /** multi(a,B,j,0,a.length) */
    def scalar(a: ArraySeq[T], B: Array[ArraySeq[T]], j: Int): T = scalarPart(a, B, j, 0, a.length)
    
    def scalar(A: Array[ArraySeq[T]], k: Int, B: Array[ArraySeq[T]], l: Int): T =  
      A.zip(B.map { x => x(l) }).map(f => f._1(k)*f._2).sum
    
    def multiByElements(a: ArraySeq[T], b: ArraySeq[T]): ArraySeq[T] = 
      a.zip(b.par.map { x => x }).map(f => f._1*f._2)
 
    /** norm(a.-b)/norm(a) */
    def error(a: ArraySeq[T], b: ArraySeq[T]) =
      norm1(minus(a,b)).toDouble/norm1(b).toDouble
    
    def normCubic(a: ArraySeq[T]): T            = a.map { x => x.abs }.max
    def norm1(a: Array[T]): T                   = a.map { x => x.abs }.sum
    def norm1(a: ArraySeq[T])                   = a.map { x => x.abs }.sum
    def norm1(A: Array[ArraySeq[T]], j: Int): T = A.map { x => x(j).abs }.sum
    def norm2(a: ArraySeq[T]): T                = mySqrt(a.map { x => x*x }.sum)
    def norm2(a: Array[T]): T                   = mySqrt(a.map { x => x*x }.sum)

     /** rearranges a vector xs in accordance with the permutation vector q 
     *  @param x = (1,2)
     *  @param q = (1,0)
     *  @return (2,1)
     *  */
    def swapInverse(xs: ArraySeq[T], q: Array[Int]): Unit = {
      for (i <- 0 until xs.length) {
        val ind = indexOf(q,i)
        if (i!=ind) {
          val x = xs(i); xs(i)=xs(ind); xs(ind)=x;
          val y =  q(i);  q(i)= q(ind);  q(ind)=y;
        }
      }        
    }
  
    /** 
     *  Search index maximal element in the abs(a(beginIndex:a.length-1)). Search begins with the beginIndex. For example: 
     *  @beginIndex = 1
     *  @param a (-11, -10, 9, 2)
     *  @return 1
     *  */ 
    def indexOfMaxAbs(a: ArraySeq[T], beginIndex: Int): Int = {
      var curInd = beginIndex
      for (i <- curInd+1 until a.length) 
        if (a(curInd).abs < a(i).abs) curInd = i
      curInd    
    }  
    
    def indexOfMaxAbs(A: Array[ArraySeq[T]], j: Int, beginIndex: Int) = {
      var curInd = beginIndex
      for (i <- curInd+1 until A.length) 
        if (A(i)(j).abs < A(curInd)(j).abs) curInd = i
      curInd    
    }  
    
    /** 
     *  @return (max(abs(a)),min(abs(a)))
     *  */
    def getMaxAndMin(a: ArraySeq[T]) = {val aud = a.map { x => x.abs }; (aud.max, aud.min)}
    
    def indexOf[U](a: Array[U], elem: U) = a.indexOf(elem)    
  
    /** 
   *  @count the number of digits after the decimal point. 
   *  format = '%'+@count.toString+'f'.  
   *  @return sum(format.format(a(:))+' ')
   **/
  def toStringStr(a: ArraySeq[T], count1: Int, count2: Int): String = 
    toStringStr(a, "%"+count1+"."+count2+"f")
    
  /** 
   *  @count the number of digits after the decimal point. 
   *  format = '%'+@count.toString+'f'.  
   *  @return sum(format.format(a(:))+' ')
   *  */
  def toStringStr(a: ArraySeq[T], format: String): String = {
    val out = new java.lang.StringBuilder()
    a.foreach { x => out.append(format.format(x)).append(' ')}
    out.toString()
  }
    
  def toStringStr(a: Array[Int]): String = {
    val out = new java.lang.StringBuilder()
    a.foreach { x => out.append(x).append(' ')}
    out.toString()
  }  
    
  def toStringStr(a: ArraySeq[T]): String = toStringStr(a, "%."+4+"f")
  
  def toStringStr(a: Array[T]): String = {
    val format = "%."+4+"f"
    val out = new java.lang.StringBuilder()
    a.foreach { x => out.append(format.format(x)).append(' ')}
    out.toString
  }
  
  }
  
  object Matrix {  
    
    /** A^T A */
    def multiATA(A: Array[ArraySeq[T]]) = 
      for (i <- range(0, A.length)) 
        yield for (j <- range(0, A.length)) 
          yield scalar(A, i, A, j)
    
    /** A^T b */
    def multiATb(A: Array[ArraySeq[T]], b: ArraySeq[T]) = 
      for (j <- Array.range(0, A.length)) 
        yield scalar(b, A, j)
    
    /** f*A */
    def multi(A: Array[ArraySeq[T]], f: T) =
      for (x <- A) 
        yield for (y <- x) 
          yield y*f
    
    /** res(:) = A(:).*b  */
    def multiByString(A: Array[ArraySeq[T]], b: ArraySeq[T]) = 
      for (x <- A) 
        yield multiByElements(x, b)
    
    /** A*b */
    def multi(A: Array[ArraySeq[T]], b: ArraySeq[T])  =
      for (x <- A) yield scalar(x, b)
    
    /** A*B */
    def multi(A: Array[ArraySeq[T]], B: Array[ArraySeq[T]]) = {
      val n = A.length
      val res = zero(n,n)
      for (i <- 0 until n; j <- 0 until n; k <- 0 until n)
        res(i)(j) += A(i)(k)*B(k)(j)
      res
    }
    
    /** mat(:)(i1) <=> mat(:)(i2) */
    def swapColumn(mat: Array[ArraySeq[T]], i1: Int, i2: Int) = 
      if (i1 != i2) mat.foreach { x => {val y = x(i1); x(i1) = x(i2); x(i2) = y; } } 
     
    /** mat(i1)(:) <=> mat(i2)(:) */
    def swapString(mat: Array[ArraySeq[T]], i1: Int, i2: Int) = 
      if (i1 != i2) { val x = mat(i1).clone(); mat(i1) = mat(i2).clone(); mat(i2) = x }
      
    /** mat^T */
    def transpose(mat: Array[Array[T]]) = mat.transpose  
 
    def detLUviaU(U: Array[ArraySeq[T]]) = detDiagMatr(U)
    
    def det(A: Array[ArraySeq[T]]) = detLUviaU(getLU(A, false, null)._2)
    
    def cond(A: Array[ArraySeq[T]], L: Array[ArraySeq[T]], U: Array[ArraySeq[T]]) = 
        normMS(A)*normMS(inverseMatrix(L, U))
    
    def cond(A: Array[ArraySeq[T]]) = normMS(A)*normMS(inverseMatrix(A)) 
    
    /** |A|_{infinity} - norm */
    def normMS(A: Array[ArraySeq[T]]) = 
      A.map { x => normCubic(x) }.max
    
    /** |A|_{1} - norm */
    def normMC(A: Array[ArraySeq[T]]) = 
      A.map { x => norm1(x) }.max
    
    /** |A|_{2} - norm */
    def normSpectr(eigenvalues: ArraySeq[T]): T = {
      val maxmin = getMaxAndMin(eigenvalues)
      div(maxmin._1, maxmin._2)
    }
    
    private def beginError = 
        "Последовательность {Ak}, k=0:infinity не сходится по норме на "+
        QR_ITER+" итерациях (с точностью до QR_ERROR = "+QR_ERROR+").\n"+
        "Предполагаемые собственные числа: "
    
    private def getA_k(A: Array[ArraySeq[T]]): (Array[ArraySeq[T]], String)  = {  
      val n = A.length      
      val QR = getQR_HR(A)
      var curA = multi(QR._2, QR._1)
      
      var norm2 = normMS(curA).toDouble()
      var norm1 = norm2.toDouble+QR_ERROR*10+1.0
      var iter = QR_ITER
      while ((norm2-norm1).abs>QR_ERROR && iter>0) {
        val QR = getQR_GR(curA)
        curA = multi(QR._2, QR._1)
        norm1 = norm2
        norm2 = normMS(curA).toDouble()
        iter-=1
      } 
        
      (curA, {
          if (iter == 0) (beginError + toStringStr(for (i <- range(0,n)) yield curA(i)(i)))
          else null
      })
    }
    
    def eigenvaluesA_k(A: Array[ArraySeq[T]]): (ArraySeq[T], String) = { 
      val res = getA_k(A: Array[ArraySeq[T]])       
      val r = res._1
      val diag = zero(r.length)      
      for (i <- 0 until diag.length) diag(i) = r(i)(i)
      (diag, res._2)
    }
    
    def cond2(A: Array[ArraySeq[T]])= {
      val res = eigenvaluesA_k(A)        
      val maxAndMin = getMaxAndMin(res._1)
      (maxAndMin._1.toDouble/maxAndMin._2.toDouble, res._2)
    } 

   def inverseTriangleUp(U: Array[ArraySeq[T]]) = {
      val n = U.length
      val res = zero(n, n)
      for (i <- 0 until n) res(i)(i) = div(num.one, U(i)(i))
      for (i <- 0 until n; j <- i+1 until n) 
        res(i)(j) = -res(j)(j)*scalarPart(res(i), U, j, 0, j)
      res
    }
    
    def inverseTriangleDown(L: Array[ArraySeq[T]]) = {
      val n = L.length
      val res = zero(n, n)
      //[http://alexandr4784.narod.ru/kaplan5/kaplan_5_06.pdf p. 149] 
      for (i <- 0 until n; j <- 0 to i) 
        res(i)(j) = 
          if (i==j) inverse(L(i)(i))
          else      -res(i)(i)*scalarPart(L(i), res, j, 0, i)
        res
    }
    
    /**
    * A = LU 
    * [1, eqs. 1.8-1.10]
    * @return A^(-1)
    */
    def inverseMatrix(L: Array[ArraySeq[T]], U: Array[ArraySeq[T]]) = {
      val n = L.length
      val mat = zero(n, n)
      /** [1, eqs. 1.8-1.10] */
      for (i <- n-1 to 0; j <- n-1 to 0) mat(i)(j) = 
              if (i == j)    div(num.one - scalarPart(U(i), mat, j, j+1, n), U(j)(j))
              else if (i < j)        div(- scalarPart(U(i), mat, j, i+1, n), U(i)(i))
              else                      (- scalarPart(U(i), L,   j, j+1, n))
      mat
    }
    
    def inverseMatrix(A: Array[ArraySeq[T]]):  Array[ArraySeq[T]] = {
      val LU = getLU(A, false, null)
      inverseMatrix(LU._1, LU._2)
    }
    
    /** solution of SLAE Ax = b  */
    def solveLUP(A: Array[ArraySeq[T]], b: ArraySeq[T]) =  { 
      val q = new Array[Int](A.length)
      for (i <- 0 until A.length) q(i) = i
      val LU = getLU(A, true, q)
      val x = solveLU(LU._1, LU._2, b)
      swapInverse(x, q); 
      x
    }
 
    /** [1], p. 7, eq. (1.5) */
    def solveLU(L: Array[ArraySeq[T]], U: Array[ArraySeq[T]], b: ArraySeq[T]) = {
      val n = L.length
      
      val ys = zero(n)
      for (i <- 0 until n) ys(i) = b(i) - scalarPart(L(i), ys, 0, i)
      
      val xs = zero(n)
      for (i <- n-1 to 0 by -1) 
        xs(i) = div(num.one, U(i)(i))*(ys(i) - scalarPart(U(i), xs, i+1, n))  
      xs
    }
    
    /** [1], p. 7, eq. (1.5) */
    def solveLU(A: Array[ArraySeq[T]], b: ArraySeq[T]): ArraySeq[T] = {
      val LU = getLU(A, false, null)
      solveLU(LU._1, LU._2, b)
    }
    
    def solveQR(Q: Array[ArraySeq[T]], R: Array[ArraySeq[T]], b: ArraySeq[T]): ArraySeq[T] = {
      val n = Q.length
      val newB = multiATb(Q, b) 
      val xs = zero(n)
      for (i <- n-1 to 0 by -1)
        xs(i) = div(num.one,R(i)(i))*(newB(i) - scalarPart(R(i),xs,i+1,n)) 
      xs
    }
    
    def solveQR(A: Array[ArraySeq[T]], b: ArraySeq[T]): ArraySeq[T] = {
      val QR = getQR_HR(A)
      solveQR(QR._1, QR._2, b)
    }
    
    // Philipp Ludwig von Seidel
    
    def xsSeidel(A: Array[ArraySeq[T]], b: ArraySeq[T], ERROR: Double) =
      solveJacobiAndSeidel(A, b, "GS", ERROR: Double)
    
    def xsJacobi(A: Array[ArraySeq[T]], b: ArraySeq[T], ERROR: Double) = 
      solveJacobiAndSeidel(A, b, "Jacobi", ERROR: Double)
    
    def solveJacobiAndSeidel(A: Array[ArraySeq[T]], b: ArraySeq[T]): ArraySeq[T] = 
      xsSeidel(A, b, 1e-7)
    /** @return solution of SLAE Ax=b (with ERROR) */
    def solveJacobiAndSeidel(
              A: Array[ArraySeq[T]], b: ArraySeq[T], name: String, ERROR: Double) : ArraySeq[T] = {
      
      val n = A.length
      
      val xs  = b.clone
      val xs2 = if (name == "GS") b.clone else zero(n)
      
      val mInvDiag = new ArraySeq[T](n)
      for (i <- 0 until n) mInvDiag(i) = -inverse(A(i)(i))
      
      /** see [1: p. 21, theor. 2]
      * return (norm(B)<1)
      */
      def converge = {
            
        val B: Array[ArraySeq[T]] = 
          if (name == "GS")     //  Jacobi
          {
            val res = new Array[ArraySeq[T]](n)
            for (i <- 0 until n) res(i) = Vector.multi(A(i),mInvDiag(i))
            for (i <- 0 until n) res(i)(i) += num.one
            res
          } 
          else if (name == "Jacobi") // Seidel 
          {
            val mat = zero(n,n)  
            for (i <- 0 until n; j <- 0 to i) mat(i)(j) = -A(i)(j)
            
            val R = zero(n,n)  
            for (i <- 0 until n; j <- i+1 until n) R(i)(j) = A(i)(j)
                  
            val res = Matrix.multi(inverseTriangleDown(mat),R)
            res
          }
          else null
              
        val q = Matrix.normMC(B).toDouble
        val mInvDiagB = new ArraySeq[T](b.length)
        for (i <- 0 until n) mInvDiagB(i)*=b(i)
        val xs2 = Vector.plus(Matrix.multi(B, b), mInvDiagB)
        val dX = Vector.norm1(Vector.minus(xs2, b)).toDouble
        
        val iter = {
          val m1dX = 1.0/dX
          var k = 0
          var mult = 1.0
          if (q < 1.0) 
          while (mult > ERROR*(1.0 - q)*m1dX) {
            mult *= q
            k += 1
          }
          else if (mult < ERROR*(1.0-q)*m1dX) {k = 1}
          else k = -1;
          k
        }
        
        Array.copy(xs2, 0, xs, 0, xs.length)
        val norm = normMC(B)
      
        (iter, norm < num.one)        
      }
      
      def iter = converge._1
      def isConverge = converge._2
      
      lazy val isStability = {
        var flag  = true 
        for (i <- 0 until n) flag&&=(!isZero(A(i)(i)))
        flag
      }
         
      // [1: 3.2]
      if (!isStability) throw 
        new java.lang.IllegalArgumentException(
            "На диагонали матрицы А есть нуль, решение невозможно")
      if (!isConverge) throw 
        new java.lang.IllegalArgumentException(
            "Решение невозможно -- норма матрицы B не меньше 1.0")
      if (iter == -1) throw
        new java.lang.IllegalArgumentException(
            "Для решения необходимо слишком много итераций")
      
      for (it <- 0 until iter; i <- 0 until n) 
        xs2(i) = (scalar(A(i), xs) - A(i)(i)*xs(i) - b(i))*mInvDiag(i)
      xs2    
    }
    
      /** (c,s) **/
      def getGivensMatrixCS(v: ArraySeq[T], i: Int, j: Int): (T, T) = {
        val m1tay = inverse(mySqrt((v(i)*v(i) + v(j)*v(j))))
        (v(j)*m1tay, v(i)*m1tay)
      }
      
      /**
      * 1 0  0 0    11 12 13 14   11      12      13      14 
      * 0 c  s 0  * 21 22 23 24 = c21+s31 c22+s32 c23+s33 c24+s34 
      * 0 -s c 0    31 32 33 34  -s21+c31-s22+c32
      * 0 0  0 1    41 42 43 44   41      42      43      44
      **/
      def multiGR(B: Array[ArraySeq[T]], c: T, s: T,  i00: Int, j00: Int, n: Int) = {
        val res = zero(n, n) 
        for (i <- 0 until n) 
          if (i!=i00 && i!=j00) 
            for (j <- 0 until n) res(i)(j) = B(i)(j)
        var i0 = math.min(i00,j00); var j0 = math.max(i00,j00); 
        for (j <- 0 until n) res(i0)(j) =  c*B(i0)(j)+s*B(j0)(j)
        for (j <- 0 until n) res(j0)(j) = -s*B(i0)(j)+c*B(j0)(j)
        res
      }
            
      /**  Householder reflection 
      *  [1: 1.3--Метод отражений]
      **/
      def getQR_HR(A: Array[ArraySeq[T]]) = {
        val n = A.length
        val curH2 = zero(n,n)
        var curR = A.clone()
          
        for (k <- 0 until n-1) { 
          val v = {             
            val res = new ArraySeq[T](n-k)
            for (i <- k until n) res(i-k) = (curR(i)(k))
            res
          }
          
          def norm2(list: ArraySeq[T]): T = 
            mySqrt(list.map { x => x*x }.sum)
          
          // p = v; [1, 1.3 -- eq 1.14]
          v(0) += specSign(v(0))*norm2(v)
          
          def multiGM(vec: ArraySeq[T], mul: T) = {
            val n = vec.length
            val res = new Array[ArraySeq[T]](n)
            for (i <- 0 until n) 
              res(i) = new ArraySeq[T](n)
            for (i <- 0 until n; j <- 0 until n) res(i)(j) = vec(i)*vec(j)*mul
            res
          }
          
          /** Householders Matrix */
          val curH1: Array[ArraySeq[T]] = {
            val mat = multiGM(v, div(-two(num.one), v.map { x => x*x }.sum))
            for (i <- 0 until v.length) mat(i)(i) += num.one
            mat
          }
              
          val i2 = n - v.length
          for (i <- 0 until i2; j <- 0 until i2) curH2(i)(j) = if (i == j) num.one else num.zero
          for (i <- 0 until i2; j <- i2 until n) curH2(i)(j) = num.zero
          for (j <- 0 until i2; i <- i2 until n) curH2(i)(j) = num.zero
          for (i <- i2 until n; j <- i2 until n) curH2(i)(j) = curH1(i-i2)(j-i2)
          
          curR = multi(curH2, curR)  
        }
        
        val Q = multi(A, inverseTriangleUp(curR));
        (Q, curR)
      }
      
      /** Givens rotations */
      def getQR_GR(A: Array[ArraySeq[T]] ) = { 
        val n = A.length
        var curR = A.clone()
        
        for (i <- 0 until n; j <- 0 until i) 
          if (!isZero(curR(i)(j))) {
            val colm = for (k <- range(0, n)) yield curR(k)(j)
            val (c,s) = getGivensMatrixCS(colm, i, j)
            curR = multiGR(curR, c, s, i, j, A.length)
          }
        
        val Q = multi(A, inverseTriangleUp(curR))
        (Q, curR)
      }
      
      
    def getLLT(matr: Array[ArraySeq[T]]) = {
      val n = matr.length
      val L = zero(n, n)
      
      for (i <- 0 until n; j <- 0 until i)
        L(i)(j) =
          if (i==j) mySqrt(matr(i)(i) - scalarPart(L(i), L(i), 0, i))
          else div(num.one, L(j)(j))*matr(i)(j) - scalarPart(L(i), L(j), 0, j)
      L
    }
    
    /**
     * LU-decomposition. 
     * [1 : 1.1,1.2], [2: 1.2, 4.4] 
     * @param dwpFlag -- (true => decomposition with pivoting, false => simple decomposition)
     * @q -- the resulting vector permutations (if dwpFlag == true). 
     * @return (L,U). q is UPDATE
     */
    @throws(classOf[java.lang.IllegalArgumentException])
    def getLU(matr: Array[ArraySeq[T]], dwpFlag: Boolean, q: Array[Int]): (Array[ArraySeq[T]], Array[ArraySeq[T]])= {       
      val n = matr.length
      val U = matr.clone()
      val L = matr.clone()
      for (i <- 0 until n; j <- 0 until n) {
        U(i)(j) = num.zero
        L(i)(j) = num.zero
      }
      
      /**
       *  [2]: 
       *   стр. 11: Можно чередовать вычисление строк U(k)(:) и столбцов L(:)(k)
       *   стр. 42: Устойчивость не гарантирована 
       */
      for (k <- 0 until n) {
        for (j <- k until n) 
          U(k)(j) = matr(k)(j) - scalarPart(L(k),U,j,0,k);
        if (dwpFlag) {
          // индекс элемента с наибольшим модулем
          // !!!!!!!!!!!!!!!!!!!!!!!!!!!
          // val ind = LAfunction.Vector.indexOfMaxAbs(matr(k),k) // val ind = LAfunction.Vector.indexOfMaxAbs(U(k),k)
          val ind = indexOfMaxAbs(matr(k),k) 
          swapColumn(matr, k, ind)
          swapColumn(U, k, ind)
          if (k!=ind) { 
            val g = q(k)
            q(k) = q(ind)
            q(ind) = g
          }
        }
        if (isZero(U(k)(k))) throw new java.lang.IllegalArgumentException(
              "U("+k+")("+k+")="+U(k)(k)
              +". \n\t\t\t LU-разложение невозможно ввиду необходимости деления на сверхмалое число. \n "
              +"Вектор перестановок: " + toStringStr(q)+'\n'
              +k+"-я строка: " + toStringStr(U(k))
            )
        
        val m1Ukk = 
          implicitly[Numeric[T]] match {
            case num: Fractional[_] => import num._; num.one/U(k)(k)
            case num: Integral[_] => import num._; num.one/U(k)(k)
            case _ => sys.error("Undivisable numeric!")
          }
        L(k)(k) = num.one
        for (i <- k+1 until n) 
          L(i)(k) = (matr(i)(k) - scalarPart(L(i), U, k, 0, k))*m1Ukk
      }      
      (L,U)
    }
    
    def plusDiag(A: Array[ArraySeq[T]], x: T) = 
      for (i <- 0 until A.length) A(i)(i) += x
      
    def detDiagMatr(A: Array[ArraySeq[T]]) = 
      (for (i <- range(0, A.length)) yield A(i)(i)).product 
    
      // edit
    def countZeros(A: Array[ArraySeq[T]]) = 
      A.map { x => x.map { y => y }.filter {
                z => isZero(z) }.length
      }
    
    def toString(A: Array[ArraySeq[T]]): String = {
      val format = "%."+4+'f'
      val maxLenght = A.map { x => x.map { y =>  format.format(y).length}.max }.max
      toString(A, maxLenght, 4)
    }
    
    def toString(A: Array[ArraySeq[T]], count1: Int, count2: Int): String = {
      val out = new java.lang.StringBuilder()
      for (x <- A) out.append(toStringStr(x, count1, count2)).append('\n')
      out.toString()
    }
    
    def input(strs: Array[String]) = 
      for (x <- strs) yield x.split(' ').map { x => BigDecimal(x)} 
    
    def copy(inp: Array[ArraySeq[T]], out: Array[ArraySeq[T]], size: Int) = 
      for (i <- 0 until size) Array.copy(inp(i), 0, out(i), 0, size)
    
    def copy(inp: Array[ArraySeq[T]], out: Array[ArraySeq[T]]): Unit =
      copy(inp, out, inp.length)
      
  }
}

/*{
      if (x < num.zero) throw new IllegalArgumentException("negative numbers not allowed")
      val threshold = if (x < num.one) x.toDouble * 1e-6 else 1e-4
      
      def sqrtx(x: T, p: T, threshold: Double): T = p match {
        case q if (q == div(x, q)) => q // without this condition, non-termination with 1e50
        case q if ((q * q - x).abs.toDouble < threshold) => {
          def diff1 = (x - p * p).abs
          def diff2 = (x - div(x*x, p*p)).abs
          if (diff1 < diff2) p else div(x, p)
        }
        case _ => sqrtx(x, div2(div(p + x, p)), threshold)
      }
      
      sqrtx(x, div2(x), threshold)
    }
    *    def cos(x: T): T = {
      
      var sum = num.zero
      var iter = num.one
      
      var v1 = num.zero
      var v2 = num.one
  
      var powX = num.one
      var fact = num.one
      var sign = true
          
      while (difAbs(v1,v2)>error){
        val dv = div(powX,fact)
        sum += (if (sign) dv else -dv)
        
        powX*=sqr(x)
        iter+=num.one; fact*=iter
        iter+=num.one; fact*=iter  
        sign = !sign
        
        v1 = v2; v2 = sum
      }
      
      sum
    }
    */
    /*def div(x: T, y: T): T = {
      implicitly[Numeric[T]] match {
        case num2: Fractional[T] => 
          if (x == num.zero && y == num.zero) num.one
          else num2.div(x,y)
        case num2: Integral[T]   =>
          import num2._
          x./(y) 
        case _ => throw new Exception("Undivisable numeric!")
      }
    }*/
