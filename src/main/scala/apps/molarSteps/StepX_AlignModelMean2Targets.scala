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

import apps.molarSteps.Paths.{alignedPath, rawPath}
import apps.util.AlignmentTransforms
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random

// Align Crown GP model to intra oral laser scans (IOS)
object StepX_AlignModelMean2Targets {
  implicit val random: Random = Random(1024)

  def main(args: Array[String]) {
    scalismo.initialize()

    val modelInit = StatisticalModelIO.readStatisticalMeshModel(new File(alignedPath, "models/pca_augmented-crown_611.h5")).get
    val modelLmInit = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_smooth_aligned.json")).get

    val modelLmMean = modelLmInit.map{lm =>
      val pId = modelInit.referenceMesh.pointSet.findClosestPoint(lm.point).id
      val p = modelInit.mean.pointSet.point(pId)
      lm.copy(point = p)
    }

    val targetPath = new File("/Volumes/storage/Dropbox/Workspace/uni-data/albert/surface/trainingsets_teeth/uk6er_extended/")

    val targetMesh = MeshIO.readMesh(new File(targetPath, "0012_36.stl")).get
    val targetLM = LandmarkIO.readLandmarksJson[_3D](new File(targetPath, "0012_36.json")).get

    val tf = AlignmentTransforms.computeTransform(modelLmMean, targetLM, Point3D(0,0,0))

    val modelLm = modelLmInit.map(_.transform(tf))


    val model = modelInit.transform(tf)

    val ui = ScalismoUI()
    val modelGroup = ui.createGroup("model")
    val targetGroup = ui.createGroup("target")
    ui.show(modelGroup, model, "model")
//    ui.show(modelGroup, modelLmInit, "landmarks")
    ui.show(targetGroup, targetMesh, "target")

    StatisticalModelIO.writeStatisticalMeshModel(model, new File(alignedPath, "models/pca_augmented-crown_611_aligned.h5"))
    LandmarkIO.writeLandmarksJson(modelLm, new File(alignedPath, "models/pca_augmented-crown_611_aligned.json"))
  }
}
