#include "Clothoid.hh"
#include "CubicRootsFlocke.hh"
#include <cmath>
#include <boost/python.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

Clothoid::d_vector
pointsOnClothoid(double x0, double y0, double theta0, double k, double dk, double L, int npts){
  Clothoid::d_vector C(npts);
  Clothoid::d_vector S(npts);
  Clothoid::d_vector tvec(npts);
  double tick = L / npts;
  for(int i=0; i<npts; ++i){
    double t = tick * i;
    tvec[i] = t;
    Clothoid::GeneralizedFresnelCS(dk*t*t, k*t, theta0, C[i], S[i]);
  }
  Clothoid::d_vector ret(3 * npts);
  for(int i=0; i<npts; ++i){
    double t = tick * i;
    ret[i] = x0 + t * C[i];
    ret[i+npts] = y0 + t * S[i];
  }
  return ret;
}

// Segment is the combination of clothoid and line
// return value: X[npts], Y[npts], dk, L
// format: [X[npts], Y[npts], dk, L]

BOOST_PYTHON_MODULE(clothoid2){
  using namespace boost::python;
  def("buildClothoid", &Clothoid::buildClothoidPython);
  def("GeneralizedFresnelCS", &Clothoid::GeneralizedFresnelCSPython);
  def("pointsOnClothoid", &pointsOnClothoid);

  class_<Clothoid::d_vector>("d_vector")
      .def(vector_indexing_suite<Clothoid::d_vector>());
}
