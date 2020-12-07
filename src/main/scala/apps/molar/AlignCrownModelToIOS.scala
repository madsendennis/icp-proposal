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
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.registration.LandmarkRegistration
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

// Align Crown GP model to intra oral laser scans (IOS)
object AlignCrownModelToIOS {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]) {
    scalismo.initialize()

    val gpCrown = StatisticalModelIO.readStatisticalMeshModel(new File(rawPath, "reference/gp_model_1503-components_smooth_crown.h5")).get
    val gpLM = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_crown_smooth.json")).get

    val targetPath = new File("/Volumes/storage/Dropbox/Workspace/uni-data/albert/surface/trainingsets_teeth/uk6er_extended/")

    val targetMesh = MeshIO.readMesh(new File(targetPath, "0012_36.stl")).get
    val targetLM = LandmarkIO.readLandmarksJson[_3D](new File(targetPath, "0012_36.json")).get

    val tf = AlignmentTransforms.computeTransform(gpLM, targetLM, Point3D(0,0,0))

    val model = gpCrown.transform(tf)

    val ui = ScalismoUI()
    val modelGroup = ui.createGroup("model")
    val targetGroup = ui.createGroup("target")
    ui.show(modelGroup, model, "model")
    ui.show(modelGroup, gpLM, "landmarks")
    ui.show(targetGroup, targetMesh, "target")

    StatisticalModelIO.writeStatisticalMeshModel(model, new File(rawPath, "reference/gp_model_1503-components_smooth_crown_aligned2IOS.h5"))
  }
}
