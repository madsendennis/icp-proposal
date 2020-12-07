package apps.molar

import java.io.File

import apps.molar.Paths.rawPath
import breeze.linalg.DenseVector
import scalismo.geometry._3D
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.MeshMetrics
import scalismo.ui.api.ScalismoUI
import scalismo.ui.util.FileUtil

object GradientBasedRegistrationIOSonlyCrown {

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val modelFile = new File(rawPath, "reference/gp_model_1503-components_smooth_crown_aligned2IOS.h5")
    val model = StatisticalModelIO.readStatisticalMeshModel(modelFile).get
    val modelLms = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_crown_smooth.json")).get

    val targetMeshes = new File("/Volumes/storage/Dropbox/Workspace/uni-data/albert/surface/trainingsets_teeth/uk6er_extended/").listFiles(_.getName.endsWith(".stl")).sorted

    val meshFile = targetMeshes(5)
    //    targetMeshes.foreach { meshFile =>

    val targetName = FileUtil.basename(meshFile)

    println(s"TARGET NAME: ${targetName}")

    val targetMesh = MeshIO.readMesh(meshFile).get

    val ui = ScalismoUI(s"target: ${targetName}")
    val modelGroup = ui.createGroup("model")
    val targetGroup = ui.createGroup("target")
    val showModel = ui.show(modelGroup, model, "model")
    val showTarget = ui.show(targetGroup, targetMesh, "target")
    ui.show(targetGroup, model.mean, "init")

    val centerLandmark = modelLms.find(_.id == "MiddlePit").get.point
    val lm1 = modelLms.find(_.id == "F1").get.point
    val lm2 = modelLms.find(_.id == "L2").get.point
    val dist = (lm1 - lm2).norm / 1.2
    println(s"Tooth cross distance: ${(lm1 - lm2).norm}")

    val modelDeci = model.decimate(3000)

    val regWeights = Seq(1e-1, 1e-2, 1e-4, 1e-5, 1e-7, 1e-8)
    val initialCoefficients = DenseVector.zeros[Double](model.rank)
    val finalCoefficients = regWeights.foldLeft(initialCoefficients)((modelCoefficients, regParameters) =>
      GradientBasedRegistration.fitting(modelDeci, targetMesh, None, robustMetric = true, numOfIterations = 100, showModel = Option(showModel), regWeight = regParameters, modelCoefficients = modelCoefficients)
    )

    val registered = model.instance(finalCoefficients)
    val avg = MeshMetrics.avgDistance(registered, targetMesh)
    val haus = MeshMetrics.hausdorffDistance(registered, targetMesh)
    print(s"registration eval, avg: ${avg}, haus: ${haus}")

    //      MeshIO.writeMesh(registered, new File(rawPath, s"computed/crop/premolar2/meshGradientBasedDetailed/${targetName}.ply"))
    //    }
  }
}
