package apps.molarSteps

import java.io.File

import apps.molarSteps.Paths.rawPath
import apps.util.AlignmentTransforms
import breeze.linalg.DenseVector
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.{MeshMetrics, TriangleMesh, TriangleMesh3D}
import scalismo.numerics.{FixedPointsUniformMeshSampler3D, LBFGSOptimizer, Sampler}
import scalismo.registration._
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.ui.api.{ScalismoUI, StatisticalMeshModelViewControls}
import scalismo.ui.util.FileUtil
import scalismo.utils.Random
import scalismo.utils.Random.implicits._

case class MyLandmarkSampler(mesh: TriangleMesh[_3D], filterPoint: Point3D, filterDistance: Double, numberOfPoints: Int)(implicit rng: Random)
  extends Sampler[_3D] {
  val filteredPointset = mesh.pointSet.points.toIndexedSeq.filter(p => (p - filterPoint).norm > filterDistance)
  println(s"Mesh points: ${mesh.pointSet.numberOfPoints}, filtered points: ${filteredPointset.length}")
  override val volumeOfSampleRegion = mesh.area
  val samplePoints = filteredPointset.map(p => (p, 1.0))

  override def sample() = samplePoints
}


object GradientBasedRegistration {
  def fitting(model: StatisticalMeshModel, targetMesh: TriangleMesh3D, mySampler: Option[Sampler[_3D]], robustMetric: Boolean = false, numOfIterations: Int, showModel: Option[StatisticalMeshModelViewControls], regWeight: Double = 0.1, modelCoefficients: DenseVector[Double]): DenseVector[Double] = {

    val t0 = System.currentTimeMillis()

    val initialCoefficients = modelCoefficients //if(initialParameters) initialParameters.get.shapeParameters.parameters else DenseVector.zeros[Double](model.rank)

    val referenceMesh = model.referenceMesh
    val transformationSpace = GaussianProcessTransformationSpace(model.gp.interpolate(TriangleMeshInterpolator3D()))
    val fixedImage = referenceMesh.operations.toDistanceImage
    val movingImage = targetMesh.operations.toDistanceImage
    val sampler = mySampler.getOrElse(FixedPointsUniformMeshSampler3D(referenceMesh, referenceMesh.pointSet.numberOfPoints))
    val metric = if (robustMetric)
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

    val modelFile = new File(rawPath, "reference/gp_model_186-components_smooth_aligned.h5")
    val modelInit = StatisticalModelIO.readStatisticalMeshModel(modelFile).get
    val modelLmsInit = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_smooth_aligned.json")).get

    val targetMeshes = new File(rawPath, "specified/crop/premolar2/mesh/").listFiles(!_.getName.contains("microCT")).sorted //.filter(_.getName.contains("9.1.81"))
    val lmPath = new File(rawPath, "specified/crop/premolar2/landmarks/")

    //    val meshFile = targetMeshes(0)
    targetMeshes.foreach { meshFile =>

      val targetName = FileUtil.basename(meshFile)

      println(s"TARGET NAME: ${targetName}")

      val lmFile = lmPath.listFiles().find(f => FileUtil.basename(f) == targetName).get

      val targetMesh = MeshIO.readMesh(meshFile).get
      val targetLms = LandmarkIO.readLandmarksJson[_3D](lmFile).get

      val alignTransform = AlignmentTransforms.computeTransform(modelLmsInit, targetLms, Point3D(0, 0, 0))
      val modelFull = modelInit.transform(alignTransform)
      val modelLms = modelLmsInit.map(_.transform(alignTransform))

      val commonLmNames = modelLms.map(_.id) intersect targetLms.map(_.id)
      val corrPoints = commonLmNames.map { id =>
        val mLM = modelLms.find(_.id == id).get
        val tLM = targetLms.find(_.id == id).get
        (modelFull.referenceMesh.pointSet.findClosestPoint(mLM.point).id, targetMesh.pointSet.findClosestPoint(tLM.point).point)
      }.toIndexedSeq
      val model = modelFull.posterior(corrPoints, 0.2)

      val ui = ScalismoUI(s"target: ${targetName}")
      val modelGroup = ui.createGroup("model")
      val targetGroup = ui.createGroup("target")
      val showModel = ui.show(modelGroup, model, "model")
      val showTarget = ui.show(targetGroup, targetMesh, "target")

      val centerLandmark = modelLms.find(_.id == "MiddlePit").get.point
      val lm1 = modelLms.find(_.id == "F1").get.point
      val lm2 = modelLms.find(_.id == "L2").get.point
      val dist = (lm1 - lm2).norm / 1.2
      println(s"Tooth cross distance: ${(lm1 - lm2).norm}")

      val modelDeci = model.decimate(3000)

      val sampler = MyLandmarkSampler(modelDeci.referenceMesh, centerLandmark, dist, 1000)

      val regWeights = Seq(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
      val initialCoefficients = DenseVector.zeros[Double](model.rank)
      val finalCoefficients = regWeights.foldLeft(initialCoefficients)((modelCoefficients, regParameters) =>
        fitting(modelDeci, targetMesh, Option(sampler), robustMetric = true, numOfIterations = 200, showModel = Option(showModel), regWeight = regParameters, modelCoefficients = modelCoefficients)
      )
      val registered = model.instance(finalCoefficients)
      val avg = MeshMetrics.avgDistance(registered, targetMesh)
      val haus = MeshMetrics.hausdorffDistance(registered, targetMesh)
      print(s"registration eval, avg: ${avg}, haus: ${haus}")

      MeshIO.writeMesh(registered, new File(rawPath, s"computed/crop/premolar2/meshGradientBasedSmooth/${targetName}.ply"))
    }
  }
}
