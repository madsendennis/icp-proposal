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

package apps.spine

import java.io.File

import apps.spine.Paths.generalPath
import breeze.linalg.DenseMatrix
import scalismo.common.interpolation.{NearestNeighborInterpolator, TriangleMeshInterpolator3D}
import scalismo.common.{EuclideanSpace, Field, RealSpace}
import scalismo.geometry._
import scalismo.io.{MeshIO, StatismoIO, StatisticalModelIO}
import scalismo.kernels._
import scalismo.mesh.TriangleMesh3D
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random


object CreateGPModel {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val ui = ScalismoUI()

    val referenceMesh = MeshIO.readMesh(new File(generalPath, "reference/mesh/01_L1.stl")).get

      println("Num of points in ref: " + referenceMesh.pointSet.numberOfPoints)

      val zeroMean = Field(EuclideanSpace[_3D], (_: Point[_3D]) => EuclideanVector.zeros[_3D])

      val k = DiagonalKernel(GaussianKernel[_3D](30) * 10.0, 3) +
        DiagonalKernel(GaussianKernel[_3D](10) * 2.0, 3)

      val gp = GaussianProcess[_3D, EuclideanVector[_3D]](zeroMean, k)

      val lowRankGP: LowRankGaussianProcess[_3D, EuclideanVector[_3D]] = LowRankGaussianProcess.approximateGPCholesky(referenceMesh , gp, relativeTolerance = 0.01, interpolator = TriangleMeshInterpolator3D[EuclideanVector[_3D]]())

      val rank = lowRankGP.rank

    val outputModelFile = new File(generalPath, s"/reference/gp_model_$rank-components.h5")

      println(lowRankGP.klBasis.map(_.eigenvalue))
      val mm = StatisticalMeshModel(referenceMesh, lowRankGP)
      val modelGroup = ui.createGroup(s"Model-$rank")
      ui.show(modelGroup, mm, "model")

      println(s"Writing GP model to: ${outputModelFile}")
      StatisticalModelIO.writeStatisticalMeshModel(mm, outputModelFile)
  }
}
