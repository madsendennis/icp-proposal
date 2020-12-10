package apps.util

import scalismo.common.PointId
import scalismo.geometry._3D
import scalismo.mesh.TriangleMesh

object AdjacencyUtils {
  //  def apply[D](adj: IndexedSeq[IndexedSeq[Int]], res: IndexedSeq[D], f: (D,IndexedSeq[D]) => D): IndexedSeq[D] = {
  //    adj.indices.map(i => f(res(i),adj(i).map(res)))
  //  }

  /**
    * list of edges in a mesh. ordered by first then second pointid. no duplicates. directed to higher pointids
    */
  def getEdgeList[D](mesh: TriangleMesh[D]): IndexedSeq[(PointId, PointId)] = {
    mesh.pointSet.pointIds.flatMap(pid =>
      mesh.triangulation.adjacentPointsForPoint(pid).filter(apid => apid.id > pid.id).sortBy(_.id).map(apid => (pid, apid))
    ).toIndexedSeq
  }

  def getEdgeMap[D](mesh: TriangleMesh[D]): Map[PointId, IndexedSeq[PointId]] = getEdgeList(mesh).groupBy(_._1).map(t => (t._1, t._2.map(_._2)))

  def getLayersOfNeighborsLayered[D](mesh: TriangleMesh[D], degree: Int): IndexedSeq[IndexedSeq[(PointId, Int)]] = {
    require(degree >= 1)
    val all = (1 to degree).foldLeft(mesh.pointSet.pointIds.toIndexedSeq.map(IndexedSeq(_).map(pid => (pid, 0))))((a, layer) =>
      a.map(_.flatMap(t => mesh.triangulation.adjacentPointsForPoint(t._1)).distinct.map(pid => (pid, layer))))
    all.map(_.filter(t => t._2 > 0))
  }

  def getLayersOfNeighborsLayeredInt[D](mesh: TriangleMesh[D], degree: Int): IndexedSeq[IndexedSeq[(Int, Int)]] = {
    getLayersOfNeighborsLayered[D](mesh, degree).map(_.map(t => (t._1.id, t._2)))
  }

  def getLayersOfNeighbors[D](mesh: TriangleMesh[D], degree: Int): IndexedSeq[IndexedSeq[PointId]] = {
    require(degree >= 1)
    val all = (1 to degree).foldLeft(mesh.pointSet.pointIds.toIndexedSeq.map(IndexedSeq(_)))((a, _) =>
      a.map(_.flatMap(mesh.triangulation.adjacentPointsForPoint).distinct))
    all.zip(mesh.pointSet.pointIds.map(_.id).toIndexedSeq).map(t => t._1.filter(pid => pid.id != t._2))
  }

  def getLayersOfNeighborsInt[D](mesh: TriangleMesh[D], degree: Int): IndexedSeq[IndexedSeq[Int]] = {
    getLayersOfNeighbors[D](mesh, degree).map(_.map(_.id))
  }

  def getTotalArea(mesh: TriangleMesh[_3D]): Double = {
    mesh.triangulation.triangles.map(mesh.computeTriangleArea).sum
  }

  def averageSurfaceProperty(property: IndexedSeq[Double], adj: TriangleMesh[_3D], degree: Int): IndexedSeq[Double] = {
    AdjacencyUtils.getLayersOfNeighborsInt(adj, degree)
    adj.pointSet.pointIds.map(id => {

      val a = adj.triangulation.adjacentPointsForPoint(id)
      val sum = property(id.id) + a.map(id => property(id.id)).sum
      sum / (a.length + 1)
    }).toIndexedSeq
  }
}
