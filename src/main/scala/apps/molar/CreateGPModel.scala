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

package apps.molar

import java.io.File

import apps.molar.Paths.generalPath
import breeze.linalg.DenseMatrix
import breeze.linalg.svd.SVD
import scalismo.common.{Domain, NearestNeighborInterpolator, RealSpace, VectorField}
import scalismo.geometry._
import scalismo.io.{MeshIO, StatismoIO}
import scalismo.kernels._
import scalismo.mesh.TriangleMesh3D
import scalismo.numerics.UniformMeshSampler3D
import scalismo.statisticalmodel.{GaussianProcess, LowRankGaussianProcess, StatisticalMeshModel}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random


object CreateGPModel {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    val ui = ScalismoUI()

    val referenceMesh = MeshIO.readMesh(new File(generalPath, "reference/mesh/lowermolar_LowerJaw_full_mirrored_coarse.ply")).get

      println("Num of points in ref: " + referenceMesh.pointSet.numberOfPoints)

      val zeroMean = VectorField(RealSpace[_3D], (_: Point[_3D]) => EuclideanVector.zeros[_3D])

      val cov: MatrixValuedPDKernel[_3D] = new MatrixValuedPDKernel[_3D]() {
        // Adds more variance along the main direction of variation (the bone length)
        private val largeKernel = DiagonalKernel[_3D](GaussianKernel(6), 3) * 1.0
        private val midKernels = DiagonalKernel[_3D](GaussianKernel(3), 3) * 0.1
//        private val smallKernels = DiagonalKernel[_3D](GaussianKernel(10), 3) * 0.0

        override protected def k(x: Point[_3D], y: Point[_3D]): DenseMatrix[Double] = {
          largeKernel(x, y) + midKernels(x, y) //+ smallKernels(x, y)
        }

        override def outputDim = 3

        override def domain: Domain[_3D] = RealSpace[_3D]
      }

      val gp = GaussianProcess[_3D, EuclideanVector[_3D]](zeroMean, cov)

      val lowRankGP: LowRankGaussianProcess[_3D, EuclideanVector[_3D]] = LowRankGaussianProcess.approximateGPCholesky(referenceMesh.pointSet, gp, relativeTolerance = 0.01, interpolator = NearestNeighborInterpolator())

      val rank = lowRankGP.rank

    val outputModelFile = new File(generalPath, s"gp_model_$rank-components.h5")

      println(lowRankGP.klBasis.map(_.eigenvalue))
      val mm = StatisticalMeshModel(referenceMesh, lowRankGP)
      val modelGroup = ui.createGroup(s"Model-$rank")
      ui.show(modelGroup, mm, "model")

      println(s"Writing GP model to: ${outputModelFile}")
      StatismoIO.writeStatismoMeshModel(mm, outputModelFile)
  }
}
