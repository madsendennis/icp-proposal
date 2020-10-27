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
import apps.util.{AlignmentTransforms, FileUtils}
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO}
import scalismo.registration.LandmarkRegistration
import scalismo.utils.Random

object AlignRegisteredShapes {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]) {
    scalismo.initialize()

    val computedPath = new File(rawPath, "computed/crop/premolar2/")
    val computedLMsPath = new File(computedPath, "landmarks")
    val computedMeshesPath = new File(computedPath, "mesh").listFiles(_.getName.endsWith(".ply"))

    alignedPath.mkdir()
    val alignedLMsPath = new File(alignedPath, "landmarks")
    alignedLMsPath.mkdir()
    val alignedMeshesPath = new File(alignedPath, "mesh")
    alignedMeshesPath.mkdir()

    val referenceMesh = MeshIO.readMesh(new File(rawPath, "reference/mesh/lowermolar_LowerJaw_full_mirrored_coarse.ply")).get
    val origin = Point3D(0, 0, 0)

    computedMeshesPath.foreach { f =>
      val basename = FileUtils.basename(f)
      val lmName = s"$basename.json"
      val mesh = MeshIO.readMesh(f).get
      val lms = LandmarkIO.readLandmarksJson[_3D](new File(computedLMsPath, lmName)).get
      val pointPairs = mesh.pointSet.points.toIndexedSeq zip referenceMesh.pointSet.points.toIndexedSeq
      val transform = LandmarkRegistration.rigid3DLandmarkRegistration(pointPairs, origin)
      MeshIO.writeMesh(mesh.transform(transform), new File(alignedMeshesPath, f.getName))
      LandmarkIO.writeLandmarksJson[_3D](lms.map(_.transform(transform)), new File(alignedLMsPath, lmName))
    }
  }
}
