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

import apps.molar.Paths.alignedPath
import apps.molarSteps.Paths.rawPath
import apps.util.AlignmentTransforms
import scalismo.common.PointId
import scalismo.geometry.{Landmark, Point, Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.TriangleMesh3D
import scalismo.ui.api.ScalismoUI
import scalismo.utils.Random


object Step10_AnalyzeStuff {

  def main(args: Array[String]) {
    scalismo.initialize()
    implicit val random: Random = Random(1024)

    val crownReference = MeshIO.readMesh(new File(rawPath, "reference/mesh/lowermolar_LowerJaw_full_mirrored_smooth_aligned_crown.ply")).get

    val laserParisFiles = new File(rawPath, "supplied/cbctLaserPairs")
    val targetCBCTLmsInit: Seq[Landmark[_3D]] = LandmarkIO.readLandmarksJson3D(new File(laserParisFiles, "export/cbct_01.json")).get
    val targetIntraLms: Seq[Landmark[_3D]] = LandmarkIO.readLandmarksJson3D(new File(laserParisFiles, "DICOM_IOScanner/lower.json")).get
    val tCBCT = AlignmentTransforms.computeTransform(targetCBCTLmsInit, targetIntraLms, Point3D(0, 0, 0))
    val targetCBCTLms = targetCBCTLmsInit.map(_.transform(tCBCT))
    val targetCBCT: TriangleMesh3D = MeshIO.readMesh(new File(laserParisFiles, "export/cbct_01_reduced.stl")).get.transform(tCBCT)
    val targetMolar: TriangleMesh3D = MeshIO.readMesh(new File(laserParisFiles, "export/molar_1.stl")).get.transform(tCBCT)
    val targetIntra: TriangleMesh3D = MeshIO.readMesh(new File(laserParisFiles, "DICOM_IOScanner/8236_LowerJaw.stl")).get

    val modelLmsInit = LandmarkIO.readLandmarksJson[_3D](new File(rawPath, "reference/landmarks/ref_smooth_aligned.json")).get
    val tModel = AlignmentTransforms.computeTransform(modelLmsInit, targetIntraLms, Point3D(0, 0, 0))
    val modelInit = StatisticalModelIO.readStatisticalMeshModel(new File(alignedPath, "models/pca_combined.h5")).get
    val modelLms = modelLmsInit.map(_.transform(tModel))
    val model = modelInit.transform(tModel)

    val ref = model.referenceMesh

    val pIds: IndexedSeq[PointId] = ref.pointSet.pointIds.toIndexedSeq.filter { id =>
      val p = modelInit.referenceMesh.pointSet.point(id)
      val dist = (crownReference.pointSet.findClosestPoint(p).point - p).norm
      dist < 0.01
    }

    println(s"Points in crown area: ${pIds.length}")

    val commonLms = AlignmentTransforms.commonLandmarkPairs(modelLms, targetIntraLms)
    val posPoints: IndexedSeq[(PointId, Point[_3D])] = commonLms.map { case (lmM, lmT) => (model.referenceMesh.pointSet.findClosestPoint(lmM).id, lmT) }.toIndexedSeq
    val posterior = model.posterior(posPoints, 1.0)


    val ui = ScalismoUI()
    val modelGroup = ui.createGroup("model")
    val cbctGroup = ui.createGroup("CBCT")
    val intraGroup = ui.createGroup("Intra")
    val showModel = ui.show(modelGroup, posterior, "model")
    //    ui.show(modelGroup, modelLms, "landmarks")
    ui.show(cbctGroup, targetCBCT, "CBCT")
    ui.show(cbctGroup, targetMolar, "molar")
    //    ui.show(cbctGroup, targetCBCTLms, "Landmarks")
    ui.show(intraGroup, targetIntra, "Intra")
    //    ui.show(intraGroup, targetIntraLms, "Landmarks")


//    case class mySampler(mesh: TriangleMesh[_3D], numberOfPoints: Int)
    //      extends Sampler[_3D] {
    //      override val volumeOfSampleRegion = mesh.area
    //      private val p: Double = 1.0 / mesh.area
    //      val samplePoints = pIds.map { id => (mesh.pointSet.point(id), p) }
    //
    //      override def sample() = scala.util.Random.shuffle(samplePoints).take(numberOfPoints)
    //    }
    //
    //    val crownSampler = mySampler(ref, pIds.length)
    //
    //    val regWeights = Seq(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)
    //    val initialCoefficients = DenseVector.zeros[Double](model.rank)
    //    val finalCoefficients = regWeights.foldLeft(initialCoefficients) { (modelCoefficients, regParameters) =>
    //      println(s"Registration, current regularization weights: ${regParameters}")
    //      GradientBasedRegistration.fitting(posterior, targetIntra, Option(crownSampler), robustMetric = true, numOfIterations = 50, showModel = Option(showModel), regWeight = regParameters, modelCoefficients = modelCoefficients)
    //    }
    //    val fit = posterior.instance(finalCoefficients)
    //    ui.show(intraGroup, fit, "fit")
  }
}
