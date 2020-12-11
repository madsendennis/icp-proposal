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

import apps.molarSteps.Paths.{rawPath, alignedPath}
import scalismo.geometry._3D
import scalismo.io.{LandmarkIO, StatisticalModelIO}
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random


object Step9_PosteriorModel {

  def main(args: Array[String]) {
    scalismo.initialize()
    implicit val random: Random = Random(1024)

    val model = StatisticalModelIO.readStatisticalMeshModel(new File(alignedPath, "models/pca_combined.h5")).get
    val modelLms = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_smooth_aligned.json")).get

    //    val centerLandmark = modelLms.find(_.id == "MiddlePit").get.point
    val centerLandmark = modelLms.find(_.id == "MesialRoot").get.point
    //    val centerLandmark = modelLms.find(_.id == "DistalRoot").get.point


    val conditionPointsInit = model.referenceMesh.pointSet.pointIds.toIndexedSeq.filter { id =>
      val p = model.referenceMesh.pointSet.point(id)
      (p - centerLandmark).norm < 2.0
    }

    val conditionIds = scala.util.Random.shuffle(conditionPointsInit).take(math.max(conditionPointsInit.length, 100000))

    val sample = model.sample()
    val cp = conditionIds.map(id => (id, sample.pointSet.point(id)))

    val posteriorModel = model.posterior(cp, 0.01)

    val ui = ScalismoUI()
    val modelGroup = ui.createGroup("model")

    ui.show(modelGroup, posteriorModel, "model")
    ui.show(modelGroup, modelLms, "landmarks")


  }
}
