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

import apps.molar.Paths.{alignedPath, rawPath}
import apps.util.FileUtils
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.registration.LandmarkRegistration
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.statisticalmodel.dataset.DataCollection
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

object ComputePCAmodel {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]) {
    scalismo.initialize()

    val alignedMeshesPath = new File(alignedPath, "detailed/mesh").listFiles(_.getName.endsWith(".ply"))
    val meshes = alignedMeshesPath.map(f => MeshIO.readMesh(f).get)
    val referenceMesh = MeshIO.readMesh(new File(rawPath, "reference/mesh/lowermolar_LowerJaw_full_mirrored.ply")).get

    val dc = DataCollection.fromTriangleMesh3DSequence(referenceMesh, meshes)
    val alignedDC = DataCollection.gpa(dc)
    val pca = StatisticalMeshModel.createUsingPCA(alignedDC).get

    val ui = ScalismoUI()
    ui.show(pca, "model")

    val experimentPath = new File(alignedPath, "detailed/pca")
    experimentPath.mkdir()

    StatisticalModelIO.writeStatisticalMeshModel(pca, new File(experimentPath, "pca_basic.h5"))
  }
}
