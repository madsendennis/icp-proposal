package apps.util

import breeze.linalg.{CSCMatrix, DenseMatrix, DenseVector}
import breeze.numerics.sqrt
import scalismo.common.PointId
import scalismo.faces.mesh.DiscreteLaplaceBeltrami
import scalismo.faces.mesh.DiscreteLaplaceBeltrami.WeightFunction
import scalismo.faces.numerics.SparseCholesky
import scalismo.geometry
import scalismo.geometry.{Dim, Point, EuclideanVector, EuclideanVector3D, _3D}
import scalismo.mesh.{SurfacePointProperty, TriangleId, TriangleMesh}

/**
  * Geodesic distance from heat eqn
  * From Crane et al (2017): The heat method for distance computation
  * http://doi.acm.org/10.1145/3131280
  *
  * note: can be extended to tetrahedral meshes
  * note: voronoi area not yet implemented
  */
case class GeodesicDistanceCalculator(mesh: TriangleMesh[_3D],
                                      choleskyL: CSCMatrix[Double],
                                      choleskyLC: CSCMatrix[Double],
                                      pointAreas: Map[PointId, Double],
                                      gradientsCoefficients: Map[TriangleId, Map[PointId, EuclideanVector[_3D]]],
                                      divergenceCoefficients: Map[TriangleId, Map[PointId, EuclideanVector[_3D]]],
                                      initSourceTemp: Double = 1.0) {

  val nPoints: Int = mesh.pointSet.numberOfPoints

  def getTriangleGradient(tid: TriangleId, u: DenseVector[Double]): geometry.EuclideanVector[_3D] = {
    val pids = mesh.triangulation.triangle(tid).pointIds
    val vectors = pids.map { pid =>
      gradientsCoefficients(tid)(pid) * u(pid.id)
    }
    vectors.tail.foldLeft(vectors.head)(_ + _)
  }

  def getIntegratedDivergence(pid: PointId, X: Map[TriangleId, EuclideanVector[_3D]]): Double = {
    mesh.triangulation.adjacentTrianglesForPoint(pid).map { tid =>
      divergenceCoefficients(tid)(pid).dot(X(tid))
    }.sum
  }

  def distanceToPoint(pointId: PointId): SurfacePointProperty[Double] = {
    val gamma = DenseVector.zeros[Double](mesh.pointSet.numberOfPoints)
    gamma(pointId.id) = initSourceTemp * pointAreas(pointId)
    val u = SparseCholesky.substitutionSolver(choleskyL, gamma)

    val X = mesh.triangulation.triangleIds.map { tid =>
      val grad_u = getTriangleGradient(tid, u)
      val X_j = -1.0 *: grad_u.normalize
      tid -> X_j
    }.toMap

    val div = mesh.pointSet.pointIds.map {
      getIntegratedDivergence(_, X)
    }.toIndexedSeq
    val b = DenseVector[Double](div.toArray)
    val D = SparseCholesky.substitutionSolver(choleskyLC, b)

    SurfacePointProperty(mesh.triangulation, GeodesicDistanceCalculator.normalizeDistance(D.toArray, pointId))
  }


}

object GeodesicDistanceCalculator {

  def apply(mesh: TriangleMesh[_3D]): GeodesicDistanceCalculator = {

    val nPoints: Int = mesh.pointSet.numberOfPoints

    val t = calculateTimeStep(mesh)
    val triangleAreas: Map[TriangleId, Double] = calculateTriangleAreas(mesh)
    val anglesFromCosineRule = calculateAnglesFromCosineRule(mesh)

    //    val pointAreas: Map[PointId, Double] = calculatePointAreas(mesh, triangleAreas)
    val pointAreas: Map[PointId, Double] = calculateAreasVoronoi(mesh, triangleAreas, anglesFromCosineRule)
    val cotangentWeights: WeightFunction = DiscreteLaplaceBeltrami.cotangentWeight(mesh)

    val gradientCoefficients = calculateGradientCoefficients(mesh, triangleAreas)
    val divergenceCoefficients = calculateDivergenceCoefficients(mesh, anglesFromCosineRule)

    val sparseA = {
      val builder = new CSCMatrix.Builder[Double](nPoints, nPoints)
      builder.sizeHint(nPoints)
      for (a <- pointAreas) builder.add(a._1.id, a._1.id, a._2)
      builder.result
    }

    val sparseLC = getSparseCotanOperator(mesh, pointAreas, cotangentWeights)
    val choleskyL = SparseCholesky.sparseCholesky(sparseA - (sparseLC * t))
    val choleskyLC = SparseCholesky.sparseCholesky(-sparseLC)

    new GeodesicDistanceCalculator(mesh,
      choleskyL,
      choleskyLC,
      pointAreas,
      gradientCoefficients,
      divergenceCoefficients)
  }


  def calculatedTriangleAnglesFromCosineRule(mesh: TriangleMesh[_3D], tid: TriangleId): Map[PointId, Double] = {
    val triangle = mesh.triangulation.triangle(tid)

    val pta = mesh.pointSet.point(triangle.ptId1)
    val ptb = mesh.pointSet.point(triangle.ptId2)
    val ptc = mesh.pointSet.point(triangle.ptId3)

    val ab = (ptb - pta).normalize
    val bc = (ptc - ptb).normalize
    val ca = (pta - ptc).normalize

    val alpha = math.acos(ab.dot(-ca))
    val beta = math.acos(bc.dot(-ab))
    val gamma = math.acos(ca.dot(-bc))

    val angles = IndexedSeq(alpha, beta, gamma)
    triangle.pointIds.zip(angles).toMap
  }

  def calculateAnglesFromCosineRule(mesh: TriangleMesh[_3D]): Map[TriangleId, Map[PointId, Double]] = {
    mesh.triangulation.triangleIds.map { tid => tid -> calculatedTriangleAnglesFromCosineRule(mesh, tid) }.toMap
  }

  def calculateTriangleAreas(mesh: TriangleMesh[_3D]): Map[TriangleId, Double] = {
    mesh.triangulation.triangleIds.map(tid => tid -> getAreaHeron(mesh, tid)).toMap
  }

  def calculateDivergenceCoefficients(mesh: TriangleMesh[_3D],
                                      anglesFromCosineRule: Map[TriangleId, Map[PointId, Double]]): Map[TriangleId, Map[PointId, EuclideanVector[_3D]]] = {
    val permutations = Seq((0, 1, 2), (1, 2, 0), (2, 0, 1))
    mesh.triangulation.triangleIds.map { tid =>
      val triangle = mesh.triangulation.triangle(tid)
      tid -> permutations.map { case (a, b, c) => val pidA = triangle.pointIds(a)
        val pidB = triangle.pointIds(b)
        val pidC = triangle.pointIds(c)
        val pt = mesh.pointSet.point(pidA)

        def f(pid1: PointId, pid2: PointId): EuclideanVector[_3D] = {
          val p2 = mesh.pointSet.point(pid1)
          val angle = anglesFromCosineRule(tid)(pid2)
          (p2 - pt) / math.tan(angle)
        }

        pidA -> (f(pidB, pidC) + f(pidC, pidB))
      }.toMap
    }.toMap
  }

  def calculateGradientCoefficients(mesh: TriangleMesh[_3D],
                                    triangleAreas: Map[TriangleId, Double]): Map[TriangleId, Map[PointId, EuclideanVector3D]] = {
    val permutations = Seq((0, 1, 2), (1, 2, 0), (2, 0, 1))
    mesh.triangulation.triangleIds.map { tid =>
      val Af = triangleAreas(tid)
      val N = mesh.cellNormals(tid)
      val triangle = mesh.triangulation.triangle(tid)
      tid -> permutations.map { case (a, b, c) => val pidA = triangle.pointIds(a)
        val pidB = triangle.pointIds(b)
        val pidC = triangle.pointIds(c)
        val oppositeSide = mesh.pointSet.point(pidC) - mesh.pointSet.point(pidB)
        pidA -> (N.crossproduct(oppositeSide) * (1.0 / (2 * Af)))
      }.toMap
    }.toMap
  }

  /**
    * Get a list of half edges lengths.
    * This returns for each triangle each side length once. Inner edges occure therefore twice
    */
  def getHalfEdgesLengths(m: TriangleMesh[_3D]): IndexedSeq[Double] = {
    val circularEdgesInTriangle = Seq((0, 1), (1, 2), (2, 0))
    m.triangles.flatMap { t =>
      circularEdgesInTriangle.map { e =>
        (m.pointSet.point(t.pointIds(e._2)) - m.pointSet.point(t.pointIds(e._1))).norm
      }
    }
  }

  /** Calculates the time step used to evolve the heat equation.
    * t = m*h*h where m is a constant (1 in paper)
    * h is the mean edge length(if the mesh is nonuniform the maximum edge length should be used) */
  def calculateTimeStep(mesh: TriangleMesh[_3D], m: Double = 1.0, nonUniformMesh: Boolean = false): Double = {
    val spacings = getHalfEdgesLengths(mesh)
    val h = if (nonUniformMesh) spacings.max else spacings.sum / spacings.length
    m * h * h
  }

  def getAreaHeron[D <: Dim](point1: Point[D], point2: Point[D], point3: Point[D]): Double = {
    val a = (point1 - point2).norm
    val b = (point2 - point3).norm
    val c = (point3 - point1).norm
    val s = (a + b + c) / 2 //semi-perimeter
    sqrt(s * (s - a) * (s - b) * (s - c)) //Heron's formula
  }

  def getAreaHeron[D <: Dim](mesh: TriangleMesh[D], tri: TriangleId): Double = {
    val triangle = mesh.triangulation.triangle(tri)
    val point1 = mesh.pointSet.point(triangle.ptId1)
    val point2 = mesh.pointSet.point(triangle.ptId2)
    val point3 = mesh.pointSet.point(triangle.ptId3)
    getAreaHeron(point1, point2, point3)
  }

  def calculatePointAreas(mesh: TriangleMesh[_3D], triangleAreas: Map[TriangleId, Double]): Map[PointId, Double] = {
    mesh.pointSet.pointIds.map { pid =>
      val incidentTriangles = mesh.triangulation.adjacentTrianglesForPoint(pid)
      val areas = incidentTriangles.map { tri => triangleAreas(tri) }
      (pid, areas.sum / 3.0)
    }.toMap
  }

  def calculateAreasVoronoi(mesh: TriangleMesh[_3D],
                            triangleAreas: Map[TriangleId, Double],
                            anglesFromCosineRule: Map[TriangleId, Map[PointId, Double]]): Map[PointId, Double] = { //from Meyer et al. Discrete Differential Geometry Operators for Triangulated 2-Manifolds
    mesh.pointSet.pointIds.map { pid_i =>
      val neighborTriangles = mesh.triangulation.adjacentTrianglesForPoint(pid_i)
      val a = neighborTriangles.map { tri =>
        val anglesWithVertices = anglesFromCosineRule(tri)
        val nObtuse = anglesWithVertices.count(_._2 > Math.PI / 2) //should be 0 if acute, 1 if obtuse
        val t = if (nObtuse > 0) {
          if (anglesWithVertices(pid_i) > Math.PI / 2) {
            triangleAreas(tri) / 2.0
          } else {
            triangleAreas(tri) / 4.0
          }
        } else {
          getVoronoiRegionOfPinT(mesh, anglesFromCosineRule, tri, pid_i)
        }
        t
      }.sum
      (pid_i, a)
    }.toMap
  }

  def getVoronoiRegionOfPinT(mesh: TriangleMesh[_3D],
                             anglesFromCosineRule: Map[TriangleId, Map[PointId, Double]],
                             tri: TriangleId,
                             pid: PointId): Double = { //25.10.2019 from http://rodolphe-vaillant.fr/?e=20
    val ids = mesh.triangulation.triangle(tri).pointIds.filterNot(_.id == pid.id)
    val a = anglesFromCosineRule(tri)
    ((mesh.pointSet.point(pid) - mesh.pointSet.point(ids(0))).norm2 * 1 / math.tan(a(ids(1))) + (mesh.pointSet.point(pid) - mesh.pointSet.point(
      ids(1))).norm2 * 1 / math.tan(a(ids(0)))) / 8.0
  }


  def getCotanOperator(mesh: TriangleMesh[_3D],
                       pointAreas: Map[PointId, Double],
                       cotangentWeights: (PointId, PointId) => Double): DenseMatrix[Double] = {
    val nPoints = mesh.pointSet.numberOfPoints
    val Lc = DenseMatrix.fill[Double](nPoints, nPoints)(0.0)
    mesh.pointSet.pointIds.zipWithIndex.foreach { case (id_i, _) => val neighborVertices = mesh.triangulation.adjacentPointsForPoint(
      id_i)
      neighborVertices.foreach { id_j =>
        val cotValue = cotangentWeights(id_i, id_j)
        Lc(id_i.id, id_j.id) = cotValue
        Lc(id_i.id, id_i.id) = Lc(id_i.id, id_i.id) - cotValue
      }
    }
    Lc * 0.5 * 0.5
  }

  def getSparseCotanOperator(mesh: TriangleMesh[_3D],
                             pointAreas: Map[PointId, Double],
                             cotangentWeights: (PointId, PointId) => Double): CSCMatrix[Double] = {
    val n = mesh.pointSet.numberOfPoints
    val Lc = new CSCMatrix.Builder[Double](n, n)
    Lc.sizeHint(n * 10)
    mesh.pointSet.pointIds.zipWithIndex.foreach { case (id_i, _) => val neighborVertices = mesh.triangulation.adjacentPointsForPoint(
      id_i)
      neighborVertices.foreach { id_j =>
        val cotValue = cotangentWeights(id_i, id_j)
        Lc.add(id_i.id, id_j.id, cotValue)
        Lc.add(id_i.id, id_i.id, -cotValue)
      }
    }
    Lc.result * 0.5 * 0.5
  }

  def getTriangleGradient(mesh: TriangleMesh[_3D],
                          tid: TriangleId,
                          gradientCoefficients: Map[TriangleId, Map[PointId, EuclideanVector3D]],
                          u: SurfacePointProperty[Double]): geometry.EuclideanVector[_3D] = {
    val pids = mesh.triangulation.triangle(tid).pointIds
    val vectors = pids.map { pid =>
      gradientCoefficients(tid)(pid) * u(pid)
    }
    vectors.tail.foldLeft(vectors.head)(_ + _)
  }


  def getTriangleGradient(mesh: TriangleMesh[_3D],
                          tid: TriangleId,
                          gradientCoefficients: Map[TriangleId, Map[PointId, EuclideanVector3D]],
                          u: DenseVector[Double]): geometry.EuclideanVector[_3D] = {
    val pids = mesh.triangulation.triangle(tid).pointIds
    val vectors = pids.map { pid =>
      gradientCoefficients(tid)(pid) * u(pid.id)
    }
    vectors.tail.foldLeft(vectors.head)(_ + _)
  }


  def normalizeDistance(d: IndexedSeq[Double], pid: PointId): IndexedSeq[Double] = {
    val D = d.map(dist => math.abs(dist - d(pid.id))) //h)
    D
  }
}


