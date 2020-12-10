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
import apps.molarSteps.Step4_ComputePCAmodel.gpa
import scalismo.io.{MeshIO, StatisticalModelIO}
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random


object Step7_ComputePCAmodel {

  // ERROR IN GPA computation, remove above when it is fixed in scalismo (https://github.com/unibas-gravis/scalismo/issues/358)

  def main(args: Array[String]) {
    scalismo.initialize()
    implicit val random: Random = Random(1024)

    val alignedMeshesPath = new File(rawPath, s"computed/crop/premolar2/meshGradientBasedSmooth_crown").listFiles(_.getName.endsWith(".ply"))
    val meshes = alignedMeshesPath.map { f => println(f); MeshIO.readMesh(f).get }
    val referenceMesh = MeshIO.readMesh(new File(rawPath, "reference/mesh/lowermolar_LowerJaw_full_mirrored_smooth_aligned_crown.ply")).get

    println(s"Meshes: ${meshes.length}")

    val dc = DataCollection.fromTriangleMesh3DSequence(referenceMesh, meshes)
    println(s"Datacollection size: ${dc.size}")
    val alignedDC = gpa(dc)
    val pca = StatisticalMeshModel.createUsingPCA(alignedDC).get

    val ui = ScalismoUI()
    ui.show(pca, "model")

    val experimentPath = new File(alignedPath, "models")
    experimentPath.mkdir()

    StatisticalModelIO.writeStatisticalMeshModel(pca, new File(experimentPath, "pca_basic_crown.h5"))
  }
}