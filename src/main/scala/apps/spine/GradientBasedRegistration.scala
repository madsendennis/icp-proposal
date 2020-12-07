package apps.spine

import java.io.File

import apps.spine.Paths.generalPath
import breeze.linalg.DenseVector
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.{MeshMetrics, TriangleMesh3D}
import scalismo.numerics.{FixedPointsUniformMeshSampler3D, LBFGSOptimizer}
import scalismo.registration.{GaussianProcessTransformationSpace, L2Regularizer, MeanHuberLossMetric, MeanSquaresMetric, Registration}
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.ui.api.{ScalismoUI, StatisticalMeshModelViewControls}
import scalismo.ui.util.FileUtil
import scalismo.utils.Random.implicits._


object GradientBasedRegistration {
  def fitting(model: StatisticalMeshModel, targetMesh: TriangleMesh3D, robustMetric: Boolean, numOfIterations: Int, showModel: Option[StatisticalMeshModelViewControls], regWeight: Double = 0.1, modelCoefficients: DenseVector[Double]): DenseVector[Double] = {

    val t0 = System.currentTimeMillis()

    val initialCoefficients = modelCoefficients //if(initialParameters) initialParameters.get.shapeParameters.parameters else DenseVector.zeros[Double](model.rank)

    val modelDecimated = model.decimate(2000)


    val referenceMesh = modelDecimated.referenceMesh
    val transformationSpace = GaussianProcessTransformationSpace(modelDecimated.gp.interpolate(TriangleMeshInterpolator3D()))
    val fixedImage = referenceMesh.operations.toDistanceImage
    val movingImage = targetMesh.operations.toDistanceImage
    val sampler = FixedPointsUniformMeshSampler3D(
      referenceMesh,
      referenceMesh.pointSet.numberOfPoints
    )
    val metric = if(robustMetric)
      MeanHuberLossMetric(fixedImage, movingImage, transformationSpace, sampler)
    else
    MeanSquaresMetric(fixedImage, movingImage, transformationSpace, sampler)

    val optimizer = LBFGSOptimizer(numOfIterations)
    val regularizer = L2Regularizer(transformationSpace)
    val registration = Registration(
      metric,
      regularizer,
      regWeight,
      optimizer
    )
    val registrationIterator = registration.iterator(initialCoefficients)
    val visualizingRegistrationIterator = for ((it, itnum) <- registrationIterator.zipWithIndex) yield {
      println(s"object value in iteration $itnum is ${it.value}")
      showModel.get.shapeModelTransformationView.shapeTransformationView.coefficients = it.parameters
      it
    }
    val registrationResult = visualizingRegistrationIterator.toSeq.last

    val t1 = System.currentTimeMillis()
    println(s"Registration-Timing: ${(t1 - t0) / 1000.0} sec")
    registrationResult.parameters
  }


  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val modelFile = new File(generalPath, "reference/gp_model_470-components.h5")
    val modelFull = StatisticalModelIO.readStatisticalMeshModel(modelFile).get
    val modelLms = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "reference/landmarks/01_L1.json")).get

    val targetMeshes = new File(generalPath, "aligned/meshes/").listFiles(!_.getName.contains("microCT")).sorted
    val lmPath = new File(generalPath, "aligned/landmarks/")

    val outputPath = new File(generalPath, "registered/meshes")

//    val meshFile = targetMeshes(1)
    targetMeshes.par.foreach { meshFile =>

      val targetName = FileUtil.basename(meshFile)

      println(s"TARGET NAME: ${targetName}")

      val lmFile = lmPath.listFiles().find(f => FileUtil.basename(f) == targetName).get

      val targetMesh = MeshIO.readMesh(meshFile).get
      val targetLms = LandmarkIO.readLandmarksJson[_3D](lmFile).get

      val commonLmNames = modelLms.map(_.id) intersect targetLms.map(_.id)
      val corrPoints = commonLmNames.map{id =>
        val mLM = modelLms.find(_.id == id).get
        val tLM = targetLms.find(_.id == id).get
        (modelFull.referenceMesh.pointSet.findClosestPoint(mLM.point).id, targetMesh.pointSet.findClosestPoint(tLM.point).point)
      }.toIndexedSeq
      val model = modelFull.posterior(corrPoints, 1.0)

      val ui = ScalismoUI(s"target: ${targetName}")
      val modelGroup = ui.createGroup("model")
      val targetGroup = ui.createGroup("target")
      val showModel = ui.show(modelGroup, model, "model")
      val showTarget = ui.show(targetGroup, targetMesh, "target")

      val regWeights = Seq(1e-1, 1e-2, 1e-4, 1e-6)  // USE MORE OR LESS REGULARIZATION HERE DEPENDING ON THE TARGET...
      val initialCoefficients = DenseVector.zeros[Double](model.rank)
      val finalCoefficients = regWeights.foldLeft(initialCoefficients){
        case (modelCoefficients, regParameters)=>
          println(s"Registering with regularization: ${regParameters}")
          fitting(model, targetMesh, robustMetric = true, numOfIterations = 100, showModel = Option(showModel), regWeight = regParameters, modelCoefficients = modelCoefficients)
      }
      val registered = model.instance(finalCoefficients)
      val avg = MeshMetrics.avgDistance(registered, targetMesh)
      val haus = MeshMetrics.hausdorffDistance(registered, targetMesh)
      print(s"registration eval, avg: ${avg}, haus: ${haus}")

      MeshIO.writeMesh(registered, new File(outputPath, s"${targetName}.ply"))
    }
  }
}
