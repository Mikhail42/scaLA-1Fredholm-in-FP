import math._
import LAfunction._

object Main {
  import LAfunction._
  import math._
  
  /** \int_{a}^{b} {K(s,t)\phi(s)} \dd[s] = f(t) */
  def main(args: Array[String]): Unit = {

    val in = new java.util.Scanner(System.in)
    println("Введите a, b, N")
    val a = in.next().toDouble
    val b = in.next().toDouble
    val N = in.nextInt    
    val data = new DoubleData(a, b, N)
    
    import data.SLAE.A
    println(DoubleFunction.Matrix.toString(A))
    println("cond_2(A)     = "+data.Cond.condA._1
      +(if (data.Cond.condA._2==null) "" else ". cond_2(A) не точно"))
    println("cond_2(A^T A) = "+data.Cond.condATA._1
      +(if (data.Cond.condATA._2==null) "" else ". cond_2(A^T A) не точно"))
    println
    for (i <- 0 until data.BasicData.nBestSqrAlphas)
      println("cond_2(A^T A + diag[{"+data.BasicData.bestSqrAlphas(i)+"}]) = "+data.Cond.vectCond(i))
    println
    println("Числа обусловленности матрицы")
    println("    "+"(\u03B2E      A   )")
    println("A = "+"(            )")
    println("    "+"(A^T (\u03B1^2/\u03B2)E)")
    println("при \u03B1^2 = "+
        DoubleFunction.Vector.toStringStr(data.BasicData.bestSqrAlphas, "%.7f")
        + ", ")
    println("при \u03B2 = {1, sqrt(\u03B1), sqrt(\u03B1^2+0.5*\u03C3_m^2(A))}")
    println(DoubleFunction.Matrix.toString(data.Cond.matCond))
    println    
    val xs = data.BasicData.xs
    println("frames")
    val frames = List((new draw.MyFrame, "LU"), (new draw.MyFrame, "LUP"), 
                      (new draw.MyFrame, "QR"), (new draw.MyFrame, "Gauss-Seidel"))
    println("for frames")
    for (frame <- frames) {
      try{
        for (i <- 0 until data.BasicData.nBestSqrAlphas){
          println("for frames " + frame._2 + " " + i)
          val solut =  frame._2 match {
            case "LU" => data.Solution.LU.solveATAsWithAlpha(i)
            case "LUP" => data.Solution.LUP.solveATAsWithAlpha(i)
            case "Gauss-Seidel" => data.Solution.Jacobi.solveATAsWithAlpha(i)
            case _ => data.Solution.QR.solveATAsWithAlpha(i)
          } 
          println("for frames sol")
          val x = new Array[Double](xs.length)
          val y = new Array[Double](solut.length)
          for (i <- 0 until xs.length) {
            y(i) = solut(i).toDouble
            x(i) = xs(i).toDouble
          }
          frame._1.addGraph(x, y, "\u03B1^2="+data.BasicData.bestSqrAlphas(i))
              
          println("for frames add")
        }
        println("for frames data.BasicData.nBetas")
        for (i <- 0 until data.BasicData.nBetas){
          println("for frames data.BasicData.nBetas " + i)
          val solut = frame._2 match {
            case "LU" => data.Solution.LU.solveATAsWithAlphaAndBeta(i)
            case "LUP" => data.Solution.LUP.solveATAsWithAlphaAndBeta(i)
            case "Gauss-Seidel" => data.Solution.Jacobi.solveATAsWithAlphaAndBeta(i)
            case _ => data.Solution.QR.solveATAsWithAlphaAndBeta(i)
          }        
          
          println("for frames data.BasicData.nBetas sol" + i)
          val x = new Array[Double](xs.length)
          val y = new Array[Double](solut.length)
          for (i <- 0 until xs.length) {
            y(i) = solut(i).toDouble
            x(i) = xs(i).toDouble
          }
          frame._1.addGraph(
              x, y, 
              "\u03B1="+data.BasicData.bestSqrAlphas(3) +
              ", \u03B2=" + data.BasicData.namesBeta(i))
        }
        println("for frames building "+frame._2)
        frame._1.building(frame._2+": "+data.BasicData.name)
      } catch {
        case e: Exception => {}
      }
    }
  }
}
