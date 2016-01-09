package draw

import org.jfree.chart.ChartFactory
import org.jfree.chart.ChartPanel
import org.jfree.chart.plot.PlotOrientation
import org.jfree.data.xy.XYSeries
import org.jfree.data.xy.XYSeriesCollection


class MyFrame extends javax.swing.JFrame {
  
  val serialVersionUID = -886569502438160491L
  var xSizeWindow = 900
  var ySizeWindow = 600 
  
  val xySC  = new XYSeriesCollection
  
  def addGraph(xs: Array[Double], ys: Array[Double], methodName: String): Unit = {
    val series = new XYSeries(methodName)
    for(i <- 0 until xs.length) series.add(xs(i), ys(i))
    xySC.addSeries(series)
  }
  
  def building(figureName: String){
    val chart = ChartFactory.createXYLineChart(
                      figureName, 
                      "x", "y",
                      xySC, 
                      PlotOrientation.VERTICAL,
                      true, true, true
                )
    this.getContentPane().add(new ChartPanel(chart))
    this.setSize(xSizeWindow, ySizeWindow)
    if (xySC.getSeriesCount()>0) this.setVisible(true)
  }
}
