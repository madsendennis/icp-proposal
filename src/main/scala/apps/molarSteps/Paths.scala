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

object Paths {
  val userHome = System.getProperty("user.home")

  val generalPath = new File(userHome, "/export/skulls/projects/teeth/data")
//  val generalPath = new File(userHome, "Dropbox/Workspace/uni-data/tmp/molar")
  val rawPath = new File(generalPath, "raw")
  val alignedPath = new File(generalPath, "aligned")
}