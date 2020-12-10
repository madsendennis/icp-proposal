/*
 *  Copyright University of Basel, Graphics and Vision Research Group
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

package apps.molarSteps

import java.io.File

import apps.molar.Paths.alignedPath
import apps.molarSteps.Paths.rawPath
import apps.util.GeodesicHelper
import breeze.linalg.{DenseVector, max}
import breeze.numerics.tanh
import scalismo.common.DiscreteField.ScalarMeshField
import scalismo.common.interpolation.TriangleMeshInterpolator3D
import scalismo.common.{EuclideanSpace, Field, ScalarMeshField}
import scalismo.geometry.{EuclideanVector, Point, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.kernels.MatrixValuedPDKernel
import scalismo.mesh.TriangleMesh
import scalismo.statisticalmodel.{DiscreteLowRankGaussianProcess, GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

//case class ChangePointKernelSphere(pcaFull: StatisticalMeshModel, pcaPartial: StatisticalMeshModel, center: Point[_3D], radius: Double)
//  extends MatrixValuedPDKernel[_3D]() {
//
//  override def domain: EuclideanSpace[_3D] = EuclideanSpace[_3D]
//
//  val outputDim = 3
//
//  def s(p: Point[_3D]) = {
//    val distance = Math.max(0, Math.sqrt(Math.pow(p.x - center.x, 2) + Math.pow(p.y - center.y, 2) + Math.pow(p.z - center.z, 2)) - radius)
//    1.0 / (1 + math.exp(radius - distance))
//  }
//
//  def k(x: Point[_3D], y: Point[_3D]) = {
//    val xPartial = pcaPartial.referenceMesh.pointSet.findClosestPoint(x)
//    val yPartial = pcaPartial.referenceMesh.pointSet.findClosestPoint(y)
//    val sx = s(x)
//    val sy = s(y)
//    pcaPartial.gp.cov(xPartial.id, yPartial.id) * (1 - sx) * (1 - sy)
//  }
//}

case class ChangePointKernelGeoDistance(pcaFull: StatisticalMeshModel, pcaPartial: StatisticalMeshModel, distanceMesh: ScalarMeshField[Double], div: Double)
  extends MatrixValuedPDKernel[_3D]() {

  override def domain: EuclideanSpace[_3D] = EuclideanSpace[_3D]

  val outputDim = 3

  def s(id: Int) = {
    val dist = distanceMesh.data(id)
    max(0.0, tanh(dist / div))
  }

  def k(x: Point[_3D], y: Point[_3D]) = {
    val xPartial = pcaPartial.referenceMesh.pointSet.findClosestPoint(x)
    val yPartial = pcaPartial.referenceMesh.pointSet.findClosestPoint(y)
    val xFullId = pcaFull.referenceMesh.pointSet.findClosestPoint(x).id.id
    val yFullId = pcaFull.referenceMesh.pointSet.findClosestPoint(y).id.id
    val sx = s(xFullId)
    val sy = s(yFullId)
    pcaPartial.gp.cov(xPartial.id, yPartial.id) * sx * sy
  }
}

object Step8_CombinePCAmodels {

  //  def completePartial(partial: TriangleMesh3D, pcaFull: StatisticalMeshModel, pcaPartial: StatisticalMeshModel): TriangleMesh3D = {
  //    val fullRef = pcaFull.referenceMesh
  //    val crownRef = pcaPartial.referenceMesh
  //
  //    val idMap = crownRef.pointSet.pointIds.toIndexedSeq.map{cId =>
  //      val refP = crownRef.pointSet.point(cId)
  //      val meanP = pcaPartial.mean.pointSet.point(cId)
  //      val fId = fullRef.pointSet.findClosestPoint(refP).id
  //      (fId, cId, refP, meanP)
  //    }
  //    val posPoints = idMap.map{case (fId,cId, refP, meanP) => (fId, meanP)}
  //    pcaFull.posterior(posPoints, .1).mean
  //  }

  def myDist(p1: Point[_3D], p2: Point[_3D]): Double = {
    (p1 - p2).norm
  }

  def updateMeanTransform(combined: StatisticalMeshModel, full: StatisticalMeshModel, crown: StatisticalMeshModel, distanceMesh: ScalarMeshField[Double]): StatisticalMeshModel = {
    val gp = combined.gp

    def scale(dist: Double, div: Double): Double = {
      max(0.0, tanh(dist / div))
    }

    val maxDistToBoundary = distanceMesh.data.max
    println(s"Max dist to boundary: ${maxDistToBoundary}")

    //    val newMean = (full.referenceMesh.pointSet.points.toIndexedSeq zip full.mean.pointSet.points.toIndexedSeq).map{case(refP, meanP) => (meanP-refP)}
    val newMean = full.referenceMesh.pointSet.pointIds.toIndexedSeq.map { fId =>
      val pRef = full.referenceMesh.pointSet.point(fId)
      val pMean = full.mean.pointSet.point(fId)
      val cpCrown = crown.referenceMesh.pointSet.findClosestPoint(pRef)
      val fullMean = (pMean - pRef)
      if ((cpCrown.point - pRef).norm > 0.001)
        fullMean
      else {
        val crownMean = (crown.mean.pointSet.point(cpCrown.id) - pRef)
        //        val dist = boundaryPoints.map(p => geoHelp.distanceBetweenPoints(p, pRef)).min
        val dist = distanceMesh.data(fId.id)
        val s = scale(dist, maxDistToBoundary / 2)
        crownMean * s + fullMean * (1 - s)
      }
    }
    val newMeanVec: DenseVector[Double] = DenseVector(newMean.flatMap(_.toArray).toArray) //DenseVector(newMean.map(_.toArray).flatten.toArray)
    val newGp = DiscreteLowRankGaussianProcess[_3D, TriangleMesh, EuclideanVector[_3D]](combined.referenceMesh, newMeanVec, gp.variance, gp.basisMatrix)

    StatisticalMeshModel(combined.referenceMesh, newGp)
  }

  // ERROR IN GPA computation, remove above when it is fixed in scalismo (https://github.com/unibas-gravis/scalismo/issues/358)

  def computeDistanceMeshL2(full: TriangleMesh[_3D], partial: TriangleMesh[_3D]): ScalarMeshField[Double] = {
    val boundaryPointIds = partial.pointSet.pointIds.toIndexedSeq.filter(id => partial.operations.pointIsOnBoundary(id)).map(id => partial.pointSet.point(id)).map(p => full.pointSet.findClosestPoint(p).id)

    val data: Array[Double] = full.pointSet.pointIds.toArray.map { fId =>
      val fP = full.pointSet.point(fId)
      val pP = partial.pointSet.findClosestPoint(fP).point
      if (myDist(fP, pP) > 0.01)
        0.0
      else
        boundaryPointIds.map(id => myDist(full.pointSet.point(id), fP)).min
    }
    ScalarMeshField(full, data)
  }

  def computeDistanceMeshGeo(full: TriangleMesh[_3D], partial: TriangleMesh[_3D]): ScalarMeshField[Double] = {
    val boundaryPointIds = partial.pointSet.pointIds.toIndexedSeq.filter(id => partial.operations.pointIsOnBoundary(id)).map(id => partial.pointSet.point(id)).map(p => full.pointSet.findClosestPoint(p).id)
    println(s"Boundary points: ${boundaryPointIds.length}")
    val geoHelp = new GeodesicHelper(full, boundaryPointIds.toSet, 0)
    val t0 = System.currentTimeMillis()
    val data: Array[Double] = full.pointSet.pointIds.zipWithIndex.toArray.map { case (fId, i) =>
      if (i % 100 == 0) {
        val t1 = System.currentTimeMillis()
        println(s"${(t1 - t0) / 1.0e3} sec - Geo map, step: ${i}/${full.pointSet.numberOfPoints}")
      }
      val fP = full.pointSet.point(fId)
      val pP = partial.pointSet.findClosestPoint(fP).point
      if (myDist(fP, pP) > 0.01)
        0.0
      else
        geoHelp.distanceBetweenPoints(geoHelp.pathBetweenPointsTempDennis(fId))
    }
    ScalarMeshField(full, data)
  }


  def main(args: Array[String]) {
    scalismo.initialize()
    implicit val random: Random = Random(1024)

    val fullPCA = StatisticalModelIO.readStatisticalMeshModel(new File(alignedPath, "models/pca_basic.h5")).get
    val crownPCA = StatisticalModelIO.readStatisticalMeshModel(new File(alignedPath, "models/pca_basic_crown.h5")).get

    val distMeshL2File = new File(alignedPath, "models/ref_dist_l2.vtk")
    val distMeshL2: ScalarMeshField[Double] = MeshIO.readScalarMeshField[Double](distMeshL2File).getOrElse {
      val m = computeDistanceMeshL2(fullPCA.referenceMesh, crownPCA.referenceMesh)
      MeshIO.writeScalarMeshField[Double](m, distMeshL2File)
      m
    }
    val distMeshGeoFile = new File(alignedPath, "models/ref_dist_Geo.vtk")
    val distMeshGeo: ScalarMeshField[Double] = MeshIO.readScalarMeshField[Double](distMeshGeoFile).getOrElse {
      val m = computeDistanceMeshGeo(fullPCA.referenceMesh, crownPCA.referenceMesh)
      MeshIO.writeScalarMeshField[Double](m, distMeshGeoFile)
      m
    }

    val maxDistToBoundary = distMeshGeo.data.max


    //    val ui = ScalismoUI()
    //    val dumGroup = ui.createGroup("dum")
    //    ui.show(dumGroup, fullPCA.referenceMesh, "fullRef")
    ////    ui.show(dumGroup, fullPCA.mean, "fullMean")
    //    ui.show(dumGroup, crownPCA.referenceMesh, "crownRef")
    ////    ui.show(dumGroup, crownPCA.mean, "crownMean")
    //    ui.show(dumGroup, distMeshL2, "distL2")
    //    ui.show(dumGroup, distMeshGeo, "distGeo")


    //    //    val completed = completePartial(crownPCA.mean, fullPCA, crownPCA)
    //    //    ui.show(dumGroup, completed, "completed")


    println(s"Original FULL PCA, rank: ${fullPCA.rank}, vertices: ${fullPCA.referenceMesh.pointSet.numberOfPoints}")
    println(s"Original CROWN PCA, rank: ${crownPCA.rank}, vertices: ${crownPCA.referenceMesh.pointSet.numberOfPoints}")

    val modelLms = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_smooth_aligned.json")).get
    val centerLandmark = modelLms.find(_.id == "MiddlePit").get.point

    val fullRef = fullPCA.referenceMesh
    val crownRef = crownPCA.referenceMesh

    val ui = ScalismoUI()
    ui.show(fullRef, "full").opacity = 0.0
    ui.show(crownRef, "crown").opacity = 0.0

    //    val cpKernel = ChangePointKernelSphere(fullPCA, crownPCA, centerLandmark, 4.0)
    val cpKernel = ChangePointKernelGeoDistance(fullPCA, crownPCA, distMeshGeo, maxDistToBoundary / 2)
    val zeroMean = Field(EuclideanSpace[_3D], (_: Point[_3D]) => EuclideanVector.zeros[_3D])
    val gp = GaussianProcess[_3D, EuclideanVector[_3D]](zeroMean, cpKernel)

    val lowRankGP = LowRankGaussianProcess.approximateGPCholesky(fullPCA.referenceMesh, gp, relativeTolerance = 0.01, interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]())

    val combinedPCA = StatisticalMeshModel.augmentModel(fullPCA, lowRankGP)
    //    val combinedPCA = StatisticalMeshModel(fullPCA.referenceMesh, lowRankGP)

    val combinedMeanPCA = updateMeanTransform(combinedPCA, fullPCA, crownPCA, distMeshGeo)
    //    val combinedMeanPCA = updateMeanTransform(fullPCA, fullPCA, crownPCA, distMeshGeo)


    val rank = combinedPCA.rank
    println(s"Combined PCA model with rank: ${rank}")

    ui.show(ui.createGroup("full"), fullPCA, "full")
    ui.show(ui.createGroup("crown"), crownPCA, "crown")
    ui.show(ui.createGroup("combined"), combinedPCA, "combined")
    ui.show(ui.createGroup("comMean"), combinedMeanPCA, "comMean")

    StatisticalModelIO.writeStatisticalMeshModel(combinedMeanPCA, new File(alignedPath, "models/pca_combined.h5"))
  }
}
