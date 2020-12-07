package apps.molarSteps

import java.io.File

import apps.molar.Paths.rawPath
import breeze.linalg.DenseVector
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.mesh.MeshMetrics
import scalismo.ui.api.{ScalismoUI, ScalismoUIHeadless}
import scalismo.ui.util.FileUtil

object Step6_GradientBasedRegistration_crowns {

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val modelFile = new File(rawPath, s"reference/gp_model_2024-components_smooth_aligned_crown.h5")
    val model = StatisticalModelIO.readStatisticalMeshModel(modelFile).get //.truncate(1000)

    val targetMeshes = new File("/Volumes/storage/Dropbox/Workspace/uni-data/albert/surface/trainingsets_teeth/uk6er_extended/").listFiles(_.getName.endsWith(".stl")).sorted

    println(s"Model: ${modelFile}, vertices: ${model.referenceMesh.pointSet.numberOfPoints}, rank: ${model.rank}")

    //    val meshFile = targetMeshes(5)
    targetMeshes.par.foreach { meshFile =>

      val targetName = FileUtil.basename(meshFile)

      println(s"TARGET NAME: ${targetName}")

      val targetMesh = MeshIO.readMesh(meshFile).get

//      val ui = ScalismoUI(s"target: ${targetName}")
      val ui = ScalismoUIHeadless()
      val modelGroup = ui.createGroup("model")
      val targetGroup = ui.createGroup("target")
      val showModel = ui.show(modelGroup, model, "model")
      val showTarget = ui.show(targetGroup, targetMesh, "target")
      ui.show(targetGroup, model.mean, "init")


      val regWeights = Seq(1e-1, 1e-2, 1e-4, 1e-5, 1e-7, 1e-8)
      val initialCoefficients = DenseVector.zeros[Double](model.rank)
      val finalCoefficients = regWeights.foldLeft(initialCoefficients) { (modelCoefficients, regParameters) =>
        println(s"Registration, current regularization weights: ${regParameters}")
        GradientBasedRegistration.fitting(model, targetMesh, None, robustMetric = true, numOfIterations = 100, showModel = Option(showModel), regWeight = regParameters, modelCoefficients = modelCoefficients)
      }

      val registered = model.instance(finalCoefficients)
      val avg = MeshMetrics.avgDistance(registered, targetMesh)
      val haus = MeshMetrics.hausdorffDistance(registered, targetMesh)
      print(s"registration eval, avg: ${avg}, haus: ${haus}")

      MeshIO.writeMesh(registered, new File(rawPath, s"computed/crop/premolar2/meshGradientBasedSmooth_crown/${targetName}.ply"))
    }
  }
}
