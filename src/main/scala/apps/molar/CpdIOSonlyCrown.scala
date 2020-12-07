package apps.molar

import java.io.File

import api.registration.CPDNonRigid
import apps.molar.Paths.rawPath
import breeze.linalg.DenseVector
import scalismo.geometry._3D
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.{MeshMetrics, TriangleMesh}
import scalismo.ui.api.ScalismoUI
import scalismo.ui.util.FileUtil

object CpdIOSonlyCrown {

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val modelFile = new File(rawPath, "reference/gp_model_1503-components_smooth_crown_aligned2IOS.h5")
    val model = StatisticalModelIO.readStatisticalMeshModel(modelFile).get
    val modelLms = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_crown_smooth.json")).get
    val modelRef = model.referenceMesh.operations.decimate(10000)

    val targetMeshes = new File("/Volumes/storage/Dropbox/Workspace/uni-data/albert/surface/trainingsets_teeth/uk6er_extended/").listFiles(_.getName.endsWith(".stl")).sorted

    val meshFile = targetMeshes(5)

    val targetName = FileUtil.basename(meshFile)

    println(s"TARGET NAME: ${targetName}")

    val targetMesh = MeshIO.readMesh(meshFile).get.operations.decimate(2000)

    val ui = ScalismoUI(s"target: ${targetName}")
    val modelGroup = ui.createGroup("model")
    val targetGroup = ui.createGroup("target")
    val showTarget = ui.show(targetGroup, targetMesh, "target")
    ui.show(modelGroup, modelRef, "init")

    val centerLandmark = modelLms.find(_.id == "MiddlePit").get.point
    val lm1 = modelLms.find(_.id == "F1").get.point
    val lm2 = modelLms.find(_.id == "L2").get.point
    val dist = (lm1-lm2).norm/1.2
    println(s"Tooth cross distance: ${(lm1-lm2).norm}")


    val cpd = new CPDNonRigid(targetMesh, modelRef, lamdba = 2, beta = 2, w = 0)
    val registered = cpd.Registration(50)

//    val registered = model.instance(finalCoefficients)
    val avg = MeshMetrics.avgDistance(registered, targetMesh)
    val haus = MeshMetrics.hausdorffDistance(registered, targetMesh)
    print(s"registration eval, avg: ${avg}, haus: ${haus}")
    ui.show(modelGroup, registered, "fit")
//      MeshIO.writeMesh(registered, new File(rawPath, s"computed/crop/premolar2/meshGradientBasedDetailed/${targetName}.ply"))
  }
}
