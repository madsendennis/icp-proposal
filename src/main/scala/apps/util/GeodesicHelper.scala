package apps.util

import scalismo.common.{PointId, PointWithId, UnstructuredPointsDomain}
import scalismo.geometry.{EuclideanVector, Point, Point2D, _3D}
import scalismo.mesh.boundingSpheres.{ClosestPointInTriangle, ClosestPointIsVertex, ClosestPointOnLine, ClosestPointWithType}
import scalismo.mesh.{BarycentricCoordinates, TriangleCell, TriangleId, TriangleMesh}

import scala.collection.mutable

/**
  * variation of takashi kanai, hiromasa suzuki (2000)
  * alternatively 'fast exact and approximate geodesic paths on meshes' with interval propagation or the flat exact algo
  */
class GeodesicHelper(val original: TriangleMesh[_3D], borderids: Set[PointId], val edgePoints: Int = 1) {
  require(edgePoints >= 0)

  val border = UnstructuredPointsDomain(borderids.toIndexedSeq.map(_.id).map(domain))

  lazy val originalEdges = AdjacencyUtils.getEdgeList(original)
  lazy val originalAdjToEdge = originalEdges.map(t => (t,
    original.triangulation.adjacentPointsForPoint(t._1).intersect(
      original.triangulation.adjacentPointsForPoint(t._2))
  )).toMap
  lazy val originalAdjClockwise = original.pointSet.pointIds.map(pid => {
    //assumes that 'up' stays locally constant. returns all the points in clockwise order with an arbitrary start
    val p = original.pointSet.point(pid)
    val n = original.vertexNormals(pid)
    val adj = original.triangulation.adjacentPointsForPoint(pid)
    val ap = adj.map(original.pointSet.point)
    //finds the starting point where there are no adjacent overlap. there will always be one such point for normal shapes.
    //use this point to determine the clockwise order
    val aapids = adj.map(apid => if (pid.id < apid.id) originalAdjToEdge(pid, apid) else originalAdjToEdge(apid, pid))
    val (apid, next, _) = adj.zip(ap).zip(aapids).map { case ((apid, ap), aapids) => {
      val plane = (ap - p).crossproduct(n)
      val aasides = aapids.map(original.pointSet.point).map(aap => (aap - p).dot(plane))
      val s1 = math.signum(aasides.head)
      val s2 = math.signum(aasides.last)
      (apid, if (s1 > 0.0) aapids.head else aapids.last, s1 != s2)
    }
    }.find(_._3).get
    val aapidsMap = adj.zip(aapids).toMap
    val ids = Iterator.iterate((apid, next)) { case (apid, next) => {
      (next, aapidsMap(next).find(_.id != apid.id).get)
    }
    }.take(adj.length).map(_._1).toIndexedSeq
    (pid, ids)
  }).toMap
  lazy val pidsToTid = original.triangulation.triangleIds.map(tid => {
    val tr = original.triangulation.triangle(tid)
    val ordered = IndexedSeq(tr.ptId1.id, tr.ptId2.id, tr.ptId3.id).sorted //necessary to make pid-tr a unique map
    ((ordered(0), ordered(1), ordered(2)), tid)
  }).toMap

  def getTidToPids(p1: PointId, p2: PointId, p3: PointId): TriangleId = {
    val ordered = IndexedSeq(p1.id, p2.id, p3.id).sorted
    pidsToTid(ordered(0), ordered(1), ordered(2))
  }

  /**
    * this is a graph with edges. If viewed as mesh volume stays the same, faces are subdivided
    * to allow for a more accurate discrete shortest path calculation.
    */
  lazy val (domain, adjacency) = {
    if (edgePoints == 0) {
      (original.pointSet.points.toIndexedSeq,
        original.pointSet.pointIds.map(pid => (pid, original.triangulation.adjacentPointsForPoint(pid))).toMap)
    } else {
      val (newpoints, oldNewConnections) = originalEdges.zipWithIndex.map { case (edge, offset) => {
        val p1 = original.pointSet.point(edge._1)
        val p2 = original.pointSet.point(edge._2)
        val v = p2 - p1
        val newPoints = (1 to edgePoints).map(i => {
          val d = i / (edgePoints + 1.0)
          p1 + v * d
        })
        val others = original.triangulation.adjacentPointsForPoint(edge._1).intersect(original.triangulation.adjacentPointsForPoint(edge._2))

        //add new edges for every old point to each new point.
        val newInd = newPoints.indices.map(_ + (offset * edgePoints + original.pointSet.numberOfPoints)).map(PointId)
        val oldNewConnections = newInd.flatMap(np => /*oldPoints*/ others.map(op => (op, np)))
        val onEdgeOldNew = IndexedSeq((edge._1, newInd.head), (edge._2, newInd.last))
        val onEdgeNewNew = if (edgePoints > 1) newInd.sliding(2, 1).map(t => (t.head, t.last)).toIndexedSeq else IndexedSeq.empty[(PointId, PointId)]
        (newPoints, onEdgeOldNew ++ onEdgeNewNew ++ oldNewConnections)
      }
      }.unzip

      //add new edges among the new points
      val edgeIndex = originalEdges.zipWithIndex.toMap
      val newNewConnections = originalEdges.zip(newpoints).zipWithIndex.flatMap { case ((e, np), offset) => {
        val others = original.triangulation.adjacentPointsForPoint(e._1).intersect(original.triangulation.adjacentPointsForPoint(e._2))
        val (a, b) = if (others.head.id < others.last.id) (others.head, others.last) else (others.last, others.head) //edge order constraint could be relaxed here
        val otherEdges = IndexedSeq((e._1, a), (e._1, b), (e._2, a), (e._2, b)).map(t => if (t._1.id > t._2.id) t.swap else t).map(edgeIndex)
        val edgeIndexes = otherEdges.filter(_ > offset)
        np.indices.map(_ + offset * edgePoints + original.pointSet.numberOfPoints).map(PointId).flatMap(npid => edgeIndexes.flatMap(eind => {
          val bound = eind * edgePoints + original.pointSet.numberOfPoints
          (bound until bound + edgePoints).map(t => (npid, PointId(t)))
        }))
      }
      }

      val domain = original.pointSet.points.toIndexedSeq ++ newpoints.flatten
      //original edges are excluded at the moment to get a faster astar
      //generally also means that exploration along edges is the most expensive
      val connections = /*edges++*/ oldNewConnections.flatten ++ newNewConnections

      val out = connections.groupBy(_._1).mapValues(_.map(_._2))
      val in = connections.map(_.swap).groupBy(_._1).mapValues(_.map(_._2))
      val adjacency = (p: PointId) => {
        val i = in.get(p)
        val o = out.get(p)
        val e = IndexedSeq.empty[PointId]
        (if (i.isDefined) i.get else e) ++ (if (o.isDefined) o.get else e)
      }

      val adjMap = domain.indices.map(PointId).map(pid => (pid, adjacency(pid))).toMap

      (domain, adjMap)
    }
  }


  lazy val pointsOnEdge: Map[(PointId, PointId), IndexedSeq[PointId]] = {
    originalEdges.zipWithIndex.map { case (t, i) => {
      val offset = original.pointSet.numberOfPoints + i * edgePoints
      (t, (offset until offset + edgePoints).map(PointId))
    }
    }.toMap
  }

  def getPointsOnEdge(start: PointId, end: PointId): IndexedSeq[PointId] = {
    require(start.id < original.pointSet.numberOfPoints && end.id < original.pointSet.numberOfPoints)
    require(start.id != end.id)
    if (start.id < end.id) pointsOnEdge(start, end) else pointsOnEdge(end, start)
  }


  /**
    * simple path between known PointIds. Easiest use case. approximation of true geodesic
    */
  def pathBetweenPoints(start: PointId, end: PointId): IndexedSeq[PointId] = {
    if (start.id != end.id) {
      val ep = domain(end.id)

      def order(t: (Double, PointId, PointId)): Double = -(t._1 + (domain(t._2.id) - ep).norm)

      val queue = mutable.PriorityQueue[(Double, PointId, PointId)]()(Ordering.by(order))
      queue.+=((0.0, start, PointId(-1)))
      val closed = mutable.Map[PointId, PointId]()
      while (queue.nonEmpty) {
        val (dist, pid, origin) = queue.dequeue()
        val p = domain(pid.id)
        if (!closed.contains(pid)) {
          closed.put(pid, origin)
          if (pid == end) queue.clear() else {
            val apid = adjacency(pid)
            queue.enqueue(apid.map(a => (dist + (domain(a.id) - p).norm, a, pid)): _*)
          }
        }
      }
      Iterator.iterate(end)(closed).takeWhile(pid => pid.id != -1).toIndexedSeq.reverse
    } else IndexedSeq(start) //for a slight speedup by making the initializations unnecessary
  }

  def pathBetweenPointsTempDennis(start: PointId): IndexedSeq[PointId] = {
    //these two are provided by you outside the function
    if (! borderids.contains(start)){
      def order(t:(Double,PointId,PointId)):Double = {
        val state = domain(t._2.id)
        -(t._1+(state-border.pointSet.findClosestPoint(state).point).norm)
      }
      val queue = mutable.PriorityQueue[(Double, PointId, PointId)]()(Ordering.by(order))
      queue.+=((0.0,start,PointId(-1)))
      val closed = mutable.Map[PointId,PointId]()
      var foundEnd = PointId(-1)
      while(queue.nonEmpty){
        val (dist,pid,origin) = queue.dequeue()
        val p = domain(pid.id)
        if (!closed.contains(pid)){
          closed.put(pid,origin)
          if (borderids.contains(pid)) {
            foundEnd = pid
            queue.clear()
          } else {
            val apid = adjacency(pid)
            queue.enqueue(apid.map(a => (dist+(domain(a.id)-p).norm,a,pid)):_*)
          }
        }
      }
      Iterator.iterate(foundEnd)(closed).takeWhile(pid=>pid.id != -1).toIndexedSeq.reverse
    } else IndexedSeq(start) //for a slight speedup by making the initializations unnecessary
  }

  def distanceBetweenPoints(path: IndexedSeq[PointId]): Double = {
    if (path.length > 1) path.map(_.id).map(domain).sliding(2, 1).map(t => (t.head - t.last).norm).sum else 0.0
  }

  def distanceBetweenPoints(start: PointId, end: PointId): Double = distanceBetweenPoints(pathBetweenPoints(start, end))

  /**
    * simple path between two points. approximation of true geodesic. Apart from the closest points on the surface
    * of the start and end, also returns the pointids potentially encountered on the way. If direct path is on the
    * surface already then this will return an empty list for the intermediate pids. the returned pointids are referring
    * to the domain. The list excludes the start and end as they have no clearly associated pointid.
    */
  def pathBetweenPoints(start: Point[_3D], end: Point[_3D]): (ClosestPointWithType, IndexedSeq[PointId], ClosestPointWithType) = {
    val sclp = original.operations.closestPointOnSurface(start)
    val eclp = original.operations.closestPointOnSurface(end)
    (sclp, pathBetweenPoints(sclp, eclp), eclp)
  }

  def pathBetweenPoints(sclp: ClosestPointWithType, eclp: ClosestPointWithType): IndexedSeq[PointId] = {
    val sps = findAdjacentExistingPoints(sclp)
    val sfaces = findAdjacentFaces(sclp)
    val epsInd = findAdjacentExistingPoints(eclp)
    val eps = epsInd.toSet
    val efaces = findAdjacentFaces(eclp)
    if (sfaces.intersect(efaces).nonEmpty)
      IndexedSeq.empty
    else {
      def order(t: (Double, PointId, PointId)): Double = -(t._1 + (domain(t._2.id) - eclp.point).norm)

      val queue = mutable.PriorityQueue[(Double, PointId, PointId)]()(Ordering.by(order))
      val spid = PointId(-1)
      sps.foreach(pid => { //fill the queue by adding every adjacent starting point
        val p = domain(pid.id)
        queue.+=(((p - sclp.point).norm, pid, spid))
      })
      val closed = mutable.Map[PointId, PointId]()
      var endPid: PointId = spid //just a temp value that is overwritten once
      while (queue.nonEmpty) { //starting to expand nodes with the astar
        val (dist, pid, origin) = queue.dequeue()
        val p = domain(pid.id)
        if (!closed.contains(pid)) {
          closed.put(pid, origin)
          if (eps.contains(pid)) {
            endPid = pid
            queue.clear()
          } else {
            val apid = adjacency(pid)
            queue.enqueue(apid.map(a => (dist + (domain(a.id) - p).norm, a, pid)): _*)
          }
        }
      }
      Iterator.iterate(endPid)(closed).takeWhile(pid => pid.id != -1).toIndexedSeq.reverse
    }
  }

  def distanceBetweenPoints(start: ClosestPointWithType, path: IndexedSeq[PointId], end: ClosestPointWithType): Double = {
    (start.point +: path.map(_.id).map(domain) :+ end.point).sliding(2, 1).map(t => (t.head - t.last).norm).sum
  }

  def distanceBetweenPoints(start: Point[_3D], end: Point[_3D]): Double = {
    val (s, path, e) = pathBetweenPoints(start, end)
    distanceBetweenPoints(s, path, e)
  }

  /**
    * returns the existing 'adjacent' points to an arbitrary points. this is mainly used to help the p3d-p3d manifold search
    */
  def findAdjacentExistingPoints(p: Point[_3D]): IndexedSeq[PointId] = findAdjacentExistingPoints(original.operations.closestPointOnSurface(p))

  def findAdjacentExistingPoints(clp: ClosestPointWithType): IndexedSeq[PointId] = clp match {
    case ClosestPointIsVertex(_, _, pid) => IndexedSeq(pid)
    case ClosestPointOnLine(_, _, pids, _) => {
      val origAdj = IndexedSeq(pids._1, pids._2)
      val pOnEdge = getPointsOnEdge(pids._1, pids._2)
      origAdj ++ pOnEdge
    }
    case ClosestPointInTriangle(_, _, tid, _) => {
      val triangle = original.triangulation.triangle(tid)
      val origAdj = IndexedSeq(triangle.ptId1, triangle.ptId2, triangle.ptId3)
      val pOnEdges = IndexedSeq((triangle.ptId1, triangle.ptId2), (triangle.ptId2, triangle.ptId3), (triangle.ptId1, triangle.ptId3)).flatMap(t => getPointsOnEdge(t._1, t._2))
      origAdj ++ pOnEdges
    }
  }

  /**
    * the straight line in manifold begins on the start point and follows the direction vector for a certain length.
    * returns the end point of that line. first return is the start point projected onto the surface. if adjustLength
    * then the length will also be modified by the ratio of starting direction and original direction
    */
  def pathInDirection(start: Point[_3D], direction: EuclideanVector[_3D], length: Double, adjustLength: Boolean): (ClosestPointWithType, ClosestPointWithType) = {
    val clp = original.operations.closestPointOnSurface(start)
    (clp, pathInDirection(clp, direction, length, adjustLength))
  }

  def pathInDirection(clp: ClosestPointWithType, direction: EuclideanVector[_3D], lengthOrig: Double, adjustLength: Boolean): ClosestPointWithType = {
    //TODO maybe make a version that emits the edgepoints along the way
    val normal = findNormalToPoint(clp)
    val dir = MathHelp.projectToPlane(direction, normal)

    // returns the int of the first adjp used: 0->01, 5->56. builds ring by including connection of last to head
    def findStartingFaceWithPlanes(p: Point[_3D], adjpo: IndexedSeq[Point[_3D]]): Int = {
      val adjp = adjpo.map(ap => MathHelp.projectToPlane(ap, p, normal))
      //constructs planes. to be well defined assumes adjacent faces are not overlapping after projection.
      val planes = adjp.map(ap => ap - p).map(v => v.crossproduct(normal)) //the i+1 plane needs to be 'turned'
      val planeDots = planes.map(_.dot(dir))
      val resIt = (planeDots :+ planeDots.head).sliding(2, 1).map(tpl => {
        tpl.head >= 0.0 && tpl.last <= 0.0 //to handle some numerical cases where the direction is perfectly aligned -> 0.0 or -0.0
      }).toIndexedSeq
      resIt.zipWithIndex.find(_._1).get._2
    }

    val (startFace, tid) = clp match {
      case ClosestPointIsVertex(point, _, pid) => {
        val adjPidsOrdered = originalAdjClockwise(pid)
        val res = findStartingFaceWithPlanes(point, adjPidsOrdered.map(original.pointSet.point))
        val pid1 = adjPidsOrdered(res)
        val pid2 = adjPidsOrdered((res + 1) % adjPidsOrdered.length)
        (TriangleCell(pid, pid1, pid2), getTidToPids(pid, pid1, pid2))
      }
      case ClosestPointOnLine(point, _, pids, _) => {
        val adjpids = if (pids._1.id < pids._2.id) originalAdjToEdge(pids._1, pids._2) else originalAdjToEdge(pids._2, pids._1)
        val orderedPid1 = originalAdjClockwise(pids._1).filter(p => adjpids.exists(ap => p.id == ap.id))
        val ordered = IndexedSeq(pids._1, orderedPid1.head, pids._2, orderedPid1.last)
        val res = findStartingFaceWithPlanes(point, ordered.map(original.pointSet.point))
        if (res < 2)
          (TriangleCell(pids._1, pids._2, orderedPid1.head), getTidToPids(pids._1, pids._2, orderedPid1.head))
        else
          (TriangleCell(pids._1, pids._2, orderedPid1.last), getTidToPids(pids._1, pids._2, orderedPid1.last))
      }
      case ClosestPointInTriangle(_, _, tid, _) => (original.cells(tid.id), tid)
    }

    //a first complete test checking all three edges to start the general case
    val startNormal = original.cellNormals(tid)
    val startdirectionOrig = clp match {
      case ClosestPointInTriangle(_, _, _, _) => dir
      case _ => { //projects the dir onto the plane of the starting cell by shooting a ray along dir-direction
        val dd = dir - direction
        val intersect = dir.dot(startNormal) / dd.dot(startNormal)
        -dir + dd * intersect
      }
    }

    //the length is modified to reduce walked distance if direction is strongly orthogonal to manifold
    val length = if (adjustLength) lengthOrig * (startdirectionOrig.norm / direction.norm) else lengthOrig
    val startdirection = startdirectionOrig.normalize

    val p2d = MathHelp.planeShift(clp.point, MathHelp.projectToPlane(startdirection, startNormal), startNormal, domain(startFace.ptId1.id), domain(startFace.ptId2.id), domain(startFace.ptId3.id))
    var activeIntersect = MathHelp.intersectThreeEdges(p2d._1, p2d._2, p2d._3)
    var (activeEdge, lastPoint) = if (activeIntersect._2 == 0)
      (Seq(PointWithId(p2d._1, startFace.ptId1), PointWithId(p2d._2, startFace.ptId2)).sortBy(_.id.id),
        PointWithId(p2d._3, startFace.ptId3))
    else if (activeIntersect._2 == 1)
      (Seq(PointWithId(p2d._1, startFace.ptId1), PointWithId(p2d._3, startFace.ptId3)).sortBy(_.id.id),
        PointWithId(p2d._2, startFace.ptId2))
    else
      (Seq(PointWithId(p2d._2, startFace.ptId2), PointWithId(p2d._3, startFace.ptId3)).sortBy(_.id.id),
        PointWithId(p2d._1, startFace.ptId1))

    //continue x-axis while unfolding triangles
    while (activeIntersect._1 < length) {
      val ah = activeEdge.head
      val al = activeEdge.last
      val nextPoint = originalAdjToEdge((ah.id, al.id)).filter(p => p.id != lastPoint.id.id).head
      val next2d = MathHelp.projectLastPoint2D(ah.point, al.point, domain(ah.id.id), domain(al.id.id), domain(nextPoint.id))
      activeIntersect = MathHelp.intersectTwoEdges(next2d, ah.point, al.point)
      if (activeIntersect._2 == 0) {
        activeEdge = Seq(ah, PointWithId(next2d, nextPoint)).sortBy(_.id.id)
        lastPoint = al
      } else {
        activeEdge = Seq(al, PointWithId(next2d, nextPoint)).sortBy(_.id.id)
        lastPoint = ah
      }
    }

    //construct final point
    val ids = IndexedSeq(activeEdge.head.id.id, activeEdge.last.id.id, lastPoint.id.id)
    val idsSorted = ids.sorted
    val endtid = pidsToTid((idsSorted(0), idsSorted(1), idsSorted(2)))
    val bc = BarycentricCoordinates.pointInTriangle(Point2D(length, 0.0), activeEdge.head.point, activeEdge.last.point, lastPoint.point)
    val endPoint = bc.interpolateProperty(domain(activeEdge.head.id.id), domain(activeEdge.last.id.id), domain(lastPoint.id.id))
    //    if (bc.a<0.0||bc.b<0.0||bc.c<0.0 || math.abs(bc.a+bc.b+bc.c-1.0)>0.0001)
    //      println("problem detected in geodesic direction walk while constructing final point")
    //always create ClosestPointInTriangle. transformations to edge/vertex closest point need to be done manually.
    //the reason for that is the simple handling of ClosestPointInTriangle. with this chaining these operations is cheaper.
    ClosestPointInTriangle(endPoint, 0.0, endtid, bc)
  }

  def findInfoToPoint(p: Point[_3D]): (ClosestPointWithType, EuclideanVector[_3D]) = {
    val clp = original.operations.closestPointOnSurface(p)
    (clp, findNormalToPoint(clp))
  }

  def findNormalToPoint(p: ClosestPointWithType): EuclideanVector[_3D] = {
    p match {
      case ClosestPointIsVertex(point, distanceSquared, pid) => original.vertexNormals(pid)
      case ClosestPointOnLine(point, distanceSquared, pids, bc) => original.vertexNormals(pids._1) * bc + original.vertexNormals(pids._2) * (1.0 - bc)
      case ClosestPointInTriangle(point, distanceSquared, tid, bc) => original.cellNormals(tid)
    }
  }

  def findAdjacentFaces(p: ClosestPointWithType): IndexedSeq[TriangleId] = {
    p match {
      case ClosestPointIsVertex(point, distanceSquared, pid) => original.triangulation.adjacentTrianglesForPoint(pid)
      case ClosestPointOnLine(point, distanceSquared, pids, bc) => original.triangulation.adjacentTrianglesForPoint(pids._1).intersect(original.triangulation.adjacentTrianglesForPoint(pids._2))
      case ClosestPointInTriangle(point, distanceSquared, tid, bc) => IndexedSeq(tid)
    }
  }

}

object GeodesicHelper {

  /**
    * given starting points for mesh and precomputed adjacency, calculates the path of pointids. intrinsic l2 norm used
    * as heuristic to find shortest path on graph.
    * this can be used to limit the number of regarded triangles when finding the shortest path on the surface.
    * finds points that form with all their adjacent triangles a superset that contains the shortest path
    *
    * @deprecated use class methods instead
    */
  def useAstar(mesh: TriangleMesh[_3D], adj: TriangleMesh[_3D], start: Point[_3D], end: Point[_3D]): IndexedSeq[PointId] = {
    def order(t: (Double, PointId, PointId)): Double = -(t._1 + (mesh.pointSet.point(t._2) - end).norm)

    val queue = mutable.PriorityQueue[(Double, PointId, PointId)]()(Ordering.by(order))
    val closest = mesh.pointSet.findClosestPoint(start) //TODO find point in the right direction (currently points do not have to be in a triangle of this id -> superset incorrect)
    val endpid = mesh.pointSet.findClosestPoint(end).id //TODO find point in the right direction
    queue.+=(((closest.point - start).norm, closest.id, PointId(-1)))
    val closed = mutable.Map[PointId, PointId]()
    while (queue.nonEmpty) {
      val (dist, pid, origin) = queue.dequeue()
      val p = mesh.pointSet.point(pid)
      if (!closed.contains(pid)) {
        closed.put(pid, origin)
        if (pid == endpid) queue.clear() else {
          val apid = adj.triangulation.adjacentPointsForPoint(pid)
          queue.enqueue(apid.map(a => (dist + (mesh.pointSet.point(a) - p).norm, a, pid)): _*)
        }
      }
    }
    Iterator.iterate(endpid)(closed).takeWhile(pid => pid.id != -1).toIndexedSeq.reverse
  }

  /**
    * changes the correspondence collection to reflect the neighbor distances from the mesh. this is a simple first test
    * for comparison
    */
  def perpendicularNormalApprx(mesh: TriangleMesh[_3D], correspondence: IndexedSeq[Point[_3D]], target: TriangleMesh[_3D], lr: Double = 1.0): IndexedSeq[Point[_3D]] = {
    mesh.pointSet.pointIds.map(pid => { //uses flat vector interpolation everywhere -> should have no impact
      val mp = mesh.pointSet.point(pid)
      val adjpids = mesh.triangulation.adjacentPointsForPoint(pid)
      val adjp = adjpids.map(pid => mesh.pointSet.point(pid))
      val adjpMean = adjp.map(_.toVector).reduce(_ + _) * (1.0 / adjp.length)

      val cp = correspondence(pid.id)
      val cadjp = adjpids.map(_.id).map(correspondence)
      val cadjpMean = cadjp.map(_.toVector).reduce(_ + _) * (1.0 / cadjp.length)

      val n = mesh.vertexNormals(pid)
      val newp = mp + cadjpMean - adjpMean
      val rp = cp + (newp - cp) * lr

      val found = target.operations.getIntersectionPointsOnSurface(rp, n)
      found.map(t => {
        val tr = target.triangulation.triangle(t._1)
        val p = t._2.interpolateProperty(target.pointSet.point(tr.ptId1), target.pointSet.point(tr.ptId2), target.pointSet.point(tr.ptId3))
        p
      }).minBy(p => (cp - p).norm2)
    })
  }.toIndexedSeq
}
