#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

/**
 * A class for dealin with the computation of Bezier curves and their properties
 **/
class Bezier
{

public:
  /**
 * An iterative implementation of the De Casteljau algorithm for computing a Bezier curve
 *
 * @param V  the vertices of the control polygon
 * @param t   an input parameter in [0, 1]
 *
 * @return    the point B(t), obtaining by evaluating the curve for the value 't'
 **/
  MatrixXd de_casteljau(const MatrixXd &V, double t)
  {
    // Number of rows
    int nV = V.rows();
    // Degree of the curve
    int degree = V.rows() - 1;

    // Condition of recursion function exit
    if (nV != 1){
      MatrixXd Res(nV-1,3);
      for (int i=0; i < nV-1; i++){
        Res(i,0) = (1 - t) * V(i,0) + t * V(i+1,0);
        Res(i,1) = (1 - t) * V(i,1) + t * V(i+1,1);
        Res(i,2) = (1 - t) * V(i,2) + t * V(i+1,2);
      }
      return de_casteljau(Res,t);
    }
    return V;
  }

  /**
	 * Plot the curve, for t=0, dt, 2*dt, ..., 1, with a given 'resolution' <br>
   * where dt=1/(resolution-1)
   *
   * @param resolution  number of points to be evaluated on the curve
	 */
  MatrixXd plot_curve(const MatrixXd &V, int resolution)
  {
    double dt = 1. / (resolution - 1);
    double t = 0.0d;
    MatrixXd Points(resolution-1,3);
    for (int i = 0; i < Points.rows(); i++){
      Points.row(i) = de_casteljau(V,t);
      t = t + dt;
    }
    return Points;
  }

  /**
	 * Perform the subdivision (once) of the Bezier curve (with parameter t) <br>
	 * Return two Bezier curves (with 'n' control points each)
	 */
  vector<MatrixXd> subdivide(const MatrixXd &V, double t)
  {
    vector<MatrixXd> curves{}; // the result: store the 2 curves obtained after subdivision
    MatrixXd division = de_casteljau(V,t);
    int nV = V.rows();
    vector<MatrixXd> result;
    MatrixXd res_cur(nV-1,3);
    for (int i = 0; i < nV-1; i++){
      res_cur(i,0) = (1 - t) * V(i,0) + t * V(i+1,0);
      res_cur(i,1) = (1 - t) * V(i,1) + t * V(i+1,1);
      res_cur(i,2) = (1 - t) * V(i,2) + t * V(i+1,2);
    }
    std::cout << "Result" << '\n';
    result.at(0) = res_cur;
    std::cout << result.at(0) << '\n';
    std::cout << "Try to get access to result(0)(0,0)" << '\n';
    for (int i = 1; i < nV-1; i++){
      MatrixXd res_cur(nV-1-i,3);
      for (int j = 0; j < nV-1-i; j++){
        res_cur(i,0) = (1 - t) * result.at(i-1)(i,0) + t * result.at(i-1)(i+1,0);
        res_cur(i,2) = (1 - t) * result.at(i-1)(i,2) + t * result.at(i-1)(i+1,2);
        res_cur(i,1) = (1 - t) * result.at(i-1)(i,1) + t * result.at(i-1)(i+1,1);
      }
      result.at(i) = res_cur;
    }

    return curves;
  }

  /**
	 * Plot the curve using recursive subdivision <br>
   *
   * @param levels  number of levels of subdivisions
   * @return  a polyline representing the curve to be rendered: this is obtained by concantenation of
   * the control points of all subdivided curves
	 */
  MatrixXd subdivision_plot(const MatrixXd &V, int levels)
  {
    std::cout << "computing recursive subdivision " << std::endl;
    list<MatrixXd> curves{};
    // to be completed

  }

  /**
 * Compute the tangent of a given curve c(t) for a given parameter t0
 *
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 *
 * @return    the tangent at c(t0)
 **/
  MatrixXd compute_tangent(const MatrixXd &V, double t0)
  {
    int n = V.rows(); // number of points in the control polygon
    // To be completed

  }

  /**
 * Compute the normal vector of a given curve c(t) for a given parameter t0
 *
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 *
 * @return    the normal at c(t0)
 **/
  MatrixXd compute_normal(const MatrixXd &V, double t0)
  {
    int n = V.rows(); // number of points in the control polygon
    // To be completed

  }

  /**
 * Compute a loop of points around a curve c(t) for a given parameter t0
 * The points belongs on a circle lying the hyperplane passing through c(t0) and orthogonal to tangent(t0)
 *
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 *
 * @return    a loop of vertices on the hyperplane passing through c(t0) and orthogonal to tangent(t0)
 **/
  MatrixXd compute_loop_of_vertices(const MatrixXd &V, double t0, int k, double radius)
  {
    int n = V.rows(); // number of points in the control polygon
    // To be completed
  }

};
