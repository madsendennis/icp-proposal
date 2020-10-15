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

import java.awt.Color
import java.io.File

import api.other.{LandmarkCorrespondence, ModelAndTargetSampling, ModelSampling, RegistrationComparison}
import api.sampling._
import api.sampling.evaluators.{ModelToTargetEvaluation, SymmetricEvaluation}
import apps.molar.Paths.generalPath
import apps.util.AlignmentTransforms
import scalismo.geometry.{Point3D, _3D}
import scalismo.io.{LandmarkIO, MeshIO, StatisticalModelIO}
import scalismo.mesh.{TriangleMesh, TriangleMesh3D}
import scalismo.sampling.DistributionEvaluator
import scalismo.sampling.proposals.MixtureProposal.ProposalGeneratorWithTransition
import scalismo.statisticalmodel.StatisticalMeshModel
import scalismo.ui.api.{ScalismoUI, StatisticalMeshModelViewControls}
import scalismo.ui.util.FileUtil
import scalismo.utils.Random.implicits._

object IcpProposalRegistration {

  def fitting(model: StatisticalMeshModel, targetMesh: TriangleMesh3D, evaluator: Map[String, DistributionEvaluator[ModelFittingParameters]], proposal: ProposalGeneratorWithTransition[ModelFittingParameters], numOfIterations: Int, showModel: Option[StatisticalMeshModelViewControls], log: File, initialParameters: Option[ModelFittingParameters] = None): ModelFittingParameters = {

    val samplingRegistration = new SamplingRegistration(model, targetMesh, showModel, modelUiUpdateInterval = 10, acceptInfoPrintInterval = 100)
    val t0 = System.currentTimeMillis()

    val best = samplingRegistration.runfitting(evaluator, proposal, numOfIterations, initialModelParameters = initialParameters, jsonName = log)

    val t1 = System.currentTimeMillis()
    println(s"ICP-Timing: ${(t1 - t0) / 1000.0} sec")
    best
  }

  def main(args: Array[String]): Unit = {
    scalismo.initialize()

    println(s"Starting Metropolis Hastings registrations with ICP-proposal!")

    val logPath = new File(generalPath, "log")

    val modelInit = StatisticalModelIO.readStatisticalMeshModel(new File(generalPath, "gp_model_200-components.h5")).get
    val modelLmsInit = LandmarkIO.readLandmarksJson[_3D](new File(generalPath, "reference/landmarks/ref.json")).get

    val targetMeshes = new File(Paths.generalPath, "targets/mesh/").listFiles(_.getName.endsWith("_1.stl"))
    val lmPath = new File(Paths.generalPath, "targets/landmarks/")

//    val meshFile = targetMeshes.head

    targetMeshes.foreach { meshFile =>

      val targetName = FileUtil.basename(meshFile)

      println(s"TARGET NAME: ${targetName}")

      val lmFile = lmPath.listFiles().find(f => FileUtil.basename(f) == targetName).get

      val targetMesh = MeshIO.readMesh(meshFile).get
      val targetLms = LandmarkIO.readLandmarksJson[_3D](lmFile).get

      val alignTransform = AlignmentTransforms.computeTransform(modelLmsInit, targetLms, Point3D(0, 0, 0))
      val model = modelInit.transform(alignTransform)
      val modelLms = modelLmsInit.map(_.transform(alignTransform))

      val numOfEvaluatorPoints = model.referenceMesh.pointSet.numberOfPoints / 2 // Used for the likelihood evaluator
      val numOfICPPointSamples = numOfEvaluatorPoints //model.rank*2 // Used for the ICP proposal
      val numOfSamples = 200 // Length of Markov Chain

      /** *** ***** ***** ***** ***** *****
        * Closest Point proposal configuration
        * projectionDirection:
        *  - TargetSampling (if registering partial meshes)
        *  - ModelSampling (if registering noisy meshes)
        *  - ModelAndTargetSampling (if registering clean complete meshes)
        *    **** ***** ***** ***** ***** **** */
      val proposalLM = MixedProposalDistributions.mixedProposalICP(model, targetMesh, modelLms, targetLms, 0, projectionDirection = LandmarkCorrespondence, tangentialNoise = 10.0, noiseAlongNormal = 2.0, stepLength = 0.1)
      val proposal = MixedProposalDistributions.mixedProposalICP(model, targetMesh, modelLms, targetLms, numOfICPPointSamples, projectionDirection = ModelAndTargetSampling, tangentialNoise = 3.0, noiseAlongNormal = 1.0, stepLength = 0.1)
      /* Uncomment below to use the standard "Random walk proposal" proposal */
      //    val proposal = MixedProposalDistributions.mixedProposalRandom(model)

      /** *** ***** ***** ***** ***** *****
        * Choosing the likelihood function
        *  - euclideanEvaluator (gaussian distribution): gives best L2 distance restults
        *  - hausdorffEvaluator (exponential distribution): gives best hausdorff result
        *    evaluationMode:
        *  - ModelToTargetEvaluation (if registering noisy meshes)
        *  - TargetToModelEvaluation (if registering partial meshes)
        *  - SymmetricEvaluation (if registering clean complete meshes)
        *    **** ***** ***** ***** ***** **** */
      val evaluator = ProductEvaluators.proximityAndIndependent(model, targetMesh, evaluationMode = ModelToTargetEvaluation, uncertainty = 2.0, numberOfEvaluationPoints = numOfEvaluatorPoints)
      /* Uncomment below to use the hausdorff likelihood function */
      //    val evaluator = ProductEvaluators.proximityAndHausdorff(model, targetMesh, uncertainty = 100.0)

      val ui = ScalismoUI(s"MH-ICP-proposal-registration: ${targetName}")
      val modelGroup = ui.createGroup("modelGroup")
      val targetGroup = ui.createGroup("targetGroup")
      val finalGroup = ui.createGroup("finalGroup")

      val showModel = ui.show(modelGroup, model, "model")
      ui.show(modelGroup, modelLms, "landmarks")
      val showTarget = ui.show(targetGroup, targetMesh, "target")
      ui.show(targetGroup, targetLms, "landmarks")
      showTarget.color = Color.YELLOW

      val lmReg = fitting(model, targetMesh, evaluator, proposalLM, 100, Option(showModel), new File(logPath, s"icpProposalRegistration.json"))
      val bestReg = fitting(model, targetMesh, evaluator, proposal, numOfSamples, Option(showModel), new File(logPath, s"icpProposalRegistration.json"), initialParameters = Option(lmReg))
      val bestRegistration = ModelFittingParameters.transformedMesh(model, bestReg)
      ui.show(finalGroup, bestRegistration, "best-fit")

      MeshIO.writeMesh(bestRegistration, new File(Paths.generalPath, s"targets/registered/${targetName}.vtk"))

      RegistrationComparison.evaluateReconstruction2GroundTruth("SAMPLE", bestRegistration, targetMesh)
    }
  }
}
