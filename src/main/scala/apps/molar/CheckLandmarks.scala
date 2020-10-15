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

import apps.util.AlignmentTransforms
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO}
import scalismo.ui.api.ScalismoUI
import scalismo.ui.util.FileUtil

object CheckLandmarks extends App{
  scalismo.initialize()

  val refMesh = MeshIO.readMesh(new File(Paths.generalPath, "reference/mesh/lowermolar_LowerJaw_full_mirrored_coarse.ply")).get
  val refLM = LandmarkIO.readLandmarksJson[_3D](new File(Paths.generalPath, "reference/landmarks/ref.json")).get

  val targetMeshes = new File(Paths.generalPath, "targets/mesh/").listFiles(_.getName.endsWith(".stl"))
  val lmPath = new File(Paths.generalPath, "targets/landmarks/")

  println(s"refLM: ${refLM.length}")

  val ui = ScalismoUI()

  targetMeshes.foreach{meshFile =>
    val targetName = FileUtil.basename(meshFile)
    val mesh = MeshIO.readMesh(meshFile).get
    val lmFile = lmPath.listFiles().find(f => FileUtil.basename(f) == targetName).get
    val targetLMs = LandmarkIO.readLandmarksJson[_3D](lmFile).get

    val commonLmNames = refLM.map(_.id) intersect targetLMs.map(_.id)

    println(s"lmFile: ${lmFile}, ${targetLMs.length}, ${commonLmNames.length}")

    val trans = AlignmentTransforms.computeTransform(targetLMs, refLM, Point3D(0,0,0))

    val grp = ui.createGroup(targetName)
    ui.show(grp, mesh.transform(trans), targetName)
//    ui.show(grp, targetLMs, "landmarks")
  }
}