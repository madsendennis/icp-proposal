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

import apps.molar.Paths.{alignedPath, rawPath}
import scalismo.common.DiscreteField
import scalismo.common.interpolation.NearestNeighborInterpolator
import scalismo.geometry._
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.mesh.{MeshMetrics, TriangleMesh}
import scalismo.registration.LandmarkRegistration
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.statisticalmodel.dataset.DataCollection.TriangleMeshDataCollection
import scalismo.transformations.Transformation
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

import scala.annotation.tailrec
import scala.collection.parallel.immutable.ParVector


object Step4_ComputePCAmodel {


  def meanTransformation(dc: TriangleMeshDataCollection[_3D]): Transformation[_3D] = {
    val fields = dc.fields(NearestNeighborInterpolator())

    Transformation { (pt: Point[_3D]) => {
      var meanPoint = EuclideanVector3D(0, 0, 0)

      for (field <- fields) {
        meanPoint = meanPoint + (pt + field(pt)).toVector
      }
      (meanPoint / fields.size).toPoint

    }
    }
  }

  def meanSurfaceFromDataCollection(dc: TriangleMeshDataCollection[_3D]): TriangleMesh[_3D] = {
    dc.reference.transform(meanTransformation(dc))
  }

  def gpa(dc: TriangleMeshDataCollection[_3D], maxIteration: Int = 5, haltDistance: Double = 1e-5)(
    implicit
    rng: Random
  ): TriangleMeshDataCollection[_3D] = {
    gpaComputation(dc, meanSurfaceFromDataCollection(dc), maxIteration, haltDistance)
  }

  @tailrec
  def gpaComputation(dc: TriangleMeshDataCollection[_3D],
                     meanShape: TriangleMesh[_3D],
                     maxIteration: Int,
                     haltDistance: Double)(implicit rng: Random): TriangleMeshDataCollection[_3D] = {

    if (maxIteration == 0) return dc

    val referencePoints = dc.reference.pointSet.points.toIndexedSeq
    val numberOfPoints = referencePoints.size
    val referenceCenterOfMass = referencePoints.foldLeft(Point3D(0, 0, 0))((acc, pt) => acc + (pt.toVector / numberOfPoints))

    val meanShapePoints = meanShape.pointSet.points.toIndexedSeq

    val fields = dc.fields(NearestNeighborInterpolator())

    // align all shape to it and create a transformations from the mean to the aligned shape
    val newDiscreteFields = new ParVector(fields.toVector).map { field =>
      val surface = dc.reference.transform(p => p + field(p))

      val transform = LandmarkRegistration.rigid3DLandmarkRegistration(
        surface.pointSet.points.toIndexedSeq.zip(meanShapePoints), referenceCenterOfMass
      )
      val newVecs = dc.reference.pointSet.points.toIndexedSeq.map(p => transform(p + field(p)) - p)
      new DiscreteField[_3D, TriangleMesh, EuclideanVector[_3D]](dc.reference, newVecs)
    }

    val newdc = DataCollection(newDiscreteFields.seq)
    val newMean = meanSurfaceFromDataCollection(newdc)

    if (MeshMetrics.procrustesDistance(meanShape, newMean) < haltDistance) {
      newdc
    } else {
      gpaComputation(newdc, newMean, maxIteration - 1, haltDistance)
    }
  }

  // ERROR IN GPA computation, remove above when it is fixed in scalismo (https://github.com/unibas-gravis/scalismo/issues/358)


  def main(args: Array[String]) {
    scalismo.initialize()
    implicit val random: Random = Random(1024)

    val alignedMeshesPath = new File(alignedPath, "mesh").listFiles(_.getName.endsWith(".ply"))
    val meshes = alignedMeshesPath.map { f => println(f); MeshIO.readMesh(f).get }
    val referenceMesh = MeshIO.readMesh(new File(rawPath, "reference/mesh/lowermolar_LowerJaw_full_mirrored_smooth_aligned.ply")).get

    println(s"Meshes: ${meshes.length}")

    val dc = DataCollection.fromTriangleMesh3DSequence(referenceMesh, meshes)
    println(s"Datacollection size: ${dc.size}")
    val alignedDC = gpa(dc)
    val pca = StatisticalMeshModel.createUsingPCA(alignedDC).get

    val ui = ScalismoUI()
    ui.show(pca, "model")

    val experimentPath = new File(alignedPath, "models")
    experimentPath.mkdir()

    StatisticalModelIO.writeStatisticalMeshModel(pca, new File(experimentPath, "pca_basic.h5"))
  }
}
