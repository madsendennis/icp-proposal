package apps.util

import breeze.linalg.DenseMatrix
import scalismo.geometry._

object MathHelp {
  def testOrtogonality(f1: IndexedSeq[EuclideanVector[_3D]], f2: IndexedSeq[EuclideanVector[_3D]]): Double = {
    f1.zip(f2).map(t => t._1.dot(t._2)).sum
  }

  /**
    * projects to the plane defined by origin and normal n. n is assumed to be normalized
    */
  def projectToPlane(p: Point[_3D], po: Point[_3D], n: EuclideanVector[_3D]): Point[_3D] = po + projectToPlane(p - po, n)

  def projectToPlane(v: EuclideanVector[_3D], n: EuclideanVector[_3D]): EuclideanVector[_3D] = v - n * n.dot(v) //.normalize)

  case class AngleInfo(angle: Double, axis: EuclideanVector[_3D])

  /**
    * returns in angle and rotation axis if 3 dimensional. axis is normalized
    */
  def getAngleInfo(v1: EuclideanVector[_3D], v2: EuclideanVector[_3D]): AngleInfo = {
    val a = v1.normalize
    val b = v2.normalize
    AngleInfo(math.acos(a.dot(b)), a.crossproduct(b))
  }

  def getInAngle[D](v1: EuclideanVector[D], v2: EuclideanVector[D]): Double = {
    val dot = v1.normalize.dot(v2.normalize)
    if (dot < -1.0) math.Pi else if (dot > 1.0) 0.0 else math.acos(dot)
  }

  /**
    * performs the edge check for all three edges.
    * also returns the int for which edge the point originates (0=12, 1=13, 2=23)
    */
  def intersectThreeEdges(t1: Point[_2D], t2: Point[_2D], t3: Point[_2D]): (Double, Int) = {
    Seq(checkEdge(t1, t2), checkEdge(t1, t3), checkEdge(t2, t3)).zipWithIndex.maxBy(_._1)
  }

  /**
    * performs the edge check for two edges 12,13 with the reason that 23 should not be considered to save time.
    */
  def intersectTwoEdges(t1: Point[_2D], t2: Point[_2D], t3: Point[_2D]): (Double, Int) = {
    Seq(checkEdge(t1, t2), checkEdge(t1, t3)).zipWithIndex.maxBy(_._1)
  }

  /**
    * returns the max x value of x-axis intersection of the given edge. returns -Inf if no intersections
    * exist.
    */
  def checkEdge(p1: Point[_2D], p2: Point[_2D]): Double = (p1.x, p1.y, p2.x, p2.y) match { //for speed changed to case class
    case (x1, 0.0, x2, 0.0) => math.max(x1, x2)
    case (_, _, x2, 0.0) => x2
    case (x1, 0.0, _, _) => x1
    case (x1, y1, x2, y2) => if (y1 * y2 < 0.0) {
      val k = -y1 / (y2 - y1)
      x1 + k * (x2 - x1)
    } else Double.NegativeInfinity
  }

  /**
    * t1 and t2 are on the xy plane and correspond to transformed triangle p1,p2,p3.
    * returns t3~p3 in xy plane. p3 is chosen to be on x+ side while the t points should be chosen such that the normal
    * under the same transformation is (0,0,1).
    */
  def projectLastPoint2D(t1: Point[_2D], t2: Point[_2D], p1: Point[_3D], p2: Point[_3D], p3: Point[_3D]): Point[_2D] = {
    val d12 = (p2 - p1).norm
    val d13 = (p3 - p1).norm
    val d23 = (p3 - p2).norm
    val x = (d12 * d12 + d13 * d13 - d23 * d23) / (2.0 * d12)
    val y = math.sqrt(d13 * d13 - x * x)
    val v = (t2 - t1).normalize
    //TODO handle v.y==0. probably handled in edge detection
    if (v.y == 0) println("problematic v.y value in projectLastPoint2D")
    val n = if (v.y > 0.0) EuclideanVector2D(v.y, -v.x) else EuclideanVector2D(-v.y, v.x)
    t1 + v * x + n * y
  }

  /**
    * shift t_i by -p then rotate around origin such that Rn = (0,0,1) and Rv = (1,0,0).
    * returns the resulting points expressed in two dims
    */
  def planeShift(p: Point[_3D], v: EuclideanVector[_3D], n: EuclideanVector[_3D], t1: Point[_3D], t2: Point[_3D], t3: Point[_3D]): (Point[_2D], Point[_2D], Point[_2D]) = {
    val fn = EuclideanVector3D(0.0, 0.0, 1.0)
    val fv = EuclideanVector3D(1.0, 0.0, 0.0)
    val r1 = rotMatrixFromAxisAngle(fn.crossproduct(n), -getInAngle(fn, n))
    val brv = r1 * v.toBreezeVector
    val av = EuclideanVector3D(brv(0), brv(1), brv(2))
    val r2 = rotMatrixFromAxisAngle(fv.crossproduct(av), -getInAngle(fv, av))
    val r = r2 * r1
    val p1 = r * (t1 - p).toBreezeVector
    val p2 = r * (t2 - p).toBreezeVector
    val p3 = r * (t3 - p).toBreezeVector
    (Point2D(p1(0), p1(1)), Point2D(p2(0), p2(1)), Point2D(p3(0), p3(1)))
  }

  /**
    * returns the resulting rotation matrix from turning around the axis 'axis' by the amount 'angle'
    */
  def rotMatrixFromAxisAngle(axis: EuclideanVector[_3D], angle: Double): DenseMatrix[Double] = {
    val a = axis.normalize
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (cos + a.x * a.x * (1 - cos), a.x * a.y * (1 - cos) - a.z * sin, a.x * a.z * (1 - cos) + a.y * sin),
      (a.y * a.x * (1 - cos) + a.z * sin, cos + a.y * a.y * (1 - cos), a.y * a.z * (1 - cos) - a.x * sin),
      (a.z * a.x * (1 - cos) - a.y * sin, a.z * a.y * (1 - cos) + a.x * sin, cos + a.z * a.z * (1 - cos))
    )
  }

  def rotX(angle: Double): DenseMatrix[Double] = {
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (1.0, 0.0, 0.0),
      (0.0, cos, -sin),
      (0.0, sin, cos)
    )
  }

  def rotY(angle: Double): DenseMatrix[Double] = {
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (cos, 0.0, sin),
      (0.0, 1.0, 0.0),
      (-sin, 0.0, cos)
    )
  }

  def rotZ(angle: Double): DenseMatrix[Double] = {
    val sin = math.sin(angle)
    val cos = math.cos(angle)
    DenseMatrix(
      (cos, -sin, 0.0),
      (sin, cos, 0.0),
      (0.0, 0.0, 1.0)
    )
  }
}
