#include <Eigen/Core>
#include <chaiscript/chaiscript.hpp>

using namespace chaiscript;
using namespace std;
using namespace Eigen;

//Functions for Matrix
template<class Mat>
string Matrix_to_string(const Mat& m){
    stringstream ss;
    ss << m;
    return ss.str();
}

template<class Mat>
void Matrix_assign1(Mat& lhs, Mat& rhs){ lhs = rhs; }

template<class Mat>
void Matrix_assign2(Mat& lhs, const Mat& rhs){ lhs = rhs; }

template<class Mat>
float& Matrix_at(Mat& m, int i, int j){
    return m(i,j);
}

template<class Mat>
Mat Matrix_transpose(Mat& m){ Mat ret = m.transpose(); return ret; }


//Functions for Vectors
template<class Vec>
string Vector_to_string(const Vec& m){
    stringstream ss;
    ss << m.transpose();
    return ss.str();
}

template<class Vec>
float& Vector_index(Vec& m, int i){ return m(i); }

template<class Vec>
void Vector_assign1(Vec& lhs, Vec& rhs){ lhs = rhs; }

template<class Vec>
void Vector_assign2(Vec& lhs, const Vec& rhs){ lhs = rhs; }

template<class Vec>
Vector3f Vector_add1(const Vec& v1, const Vec& v2){ return v1+v2; }

template<class Vec>
Vector3f Vector_add2(const Vec& v1, double v2){ return v1.array()+v2; }

template<class Vec>
Vector3f Vector_add3(double v2, const Vec& v1){ return v1.array()+v2; }

template<class Vec>
Vector3f Vector_add4(const Vec& v1, int v2){ return v1.array()+v2; }

template<class Vec>
Vector3f Vector_add5(int v2, const Vec& v1){ return v1.array()+v2; }

template<class Vec>
Vector3f Vector_sub1(const Vec& v1, const Vec& v2){ return v1-v2; }

template<class Vec>
Vector3f Vector_sub2(const Vec& v1, double v2){ return v1.array()-v2; }

template<class Vec>
Vector3f Vector_sub3(double v2, const Vec& v1){ return v2-v1.array(); }

template<class Vec>
Vector3f Vector_sub4(const Vec& v1, int v2){ return v1.array()-v2; }

template<class Vec>
Vector3f Vector_sub5(int v2, const Vec& v1){ return v2-v1.array(); }

template<class Vec>
Vector3f Vector_mult1(const Vec& v1, const Vec& v2){ return v1.array()*v2.array(); }

template<class Vec>
Vector3f Vector_mult2(const Vec& v1, double v2){ return v1.array()*v2; }

template<class Vec>
Vector3f Vector_mult3(double v2, const Vec& v1){ return v1.array()*v2; }

template<class Vec>
Vector3f Vector_mult4(const Vec& v1, int v2){ return v1.array()*v2; }

template<class Vec>
Vector3f Vector_mult5(int v2, const Vec& v1){ return v1.array()*v2; }

template<class Vec>
void Vector_fill(Vec& vec, double v){ vec.fill(v); }

//-----------------------------------
template<class Mat>
void register_eigen_matrix_type(ChaiScript& chai, string name){
    chai.add(user_type<Mat>(), name);
    chai.add(constructor<Mat()>(), name);
    chai.add(constructor<Mat(const Mat &)>(), name);
    chai.add(fun<void(Mat&, Mat&)>(&Matrix_assign1<Mat>), "=");
    chai.add(fun<void(Mat&, const Mat&)>(&Matrix_assign2<Mat>), "=");
    chai.add(fun(&Matrix_to_string<Mat>), "to_string");
    chai.add(fun(&Matrix_at<Mat>), "at");
    chai.add(fun(&Matrix_transpose<Mat>), "transpose");
}

template<class Vec>
void register_eigen_vector_type(ChaiScript& chai, string name){
    chai.add(user_type<Vec>(), name);
    chai.add(constructor<Vec()>(), name);
    chai.add(constructor<Vec(const Vec &)>(), name);
    chai.add(fun<void(Vec&, Vec&)>(&Vector_assign1<Vec>), "=");
    chai.add(fun<void(Vec&, const Vec&)>(&Vector_assign2<Vec>), "=");
    chai.add(fun(&Vector_to_string<Vec>), "to_string");
    chai.add(fun(&Vector_index<Vec>), "[]");
    chai.add(fun(&Vector_index<Vec>), "at");
    chai.add(fun(&Vector_add1<Vec>), "+");
    chai.add(fun(&Vector_add2<Vec>), "+");
    chai.add(fun(&Vector_add3<Vec>), "+");
    chai.add(fun(&Vector_add4<Vec>), "+");
    chai.add(fun(&Vector_add5<Vec>), "+");
    chai.add(fun(&Vector_sub1<Vec>), "-");
    chai.add(fun(&Vector_sub2<Vec>), "-");
    chai.add(fun(&Vector_sub3<Vec>), "-");
    chai.add(fun(&Vector_sub4<Vec>), "-");
    chai.add(fun(&Vector_sub5<Vec>), "-");
    chai.add(fun(&Vector_mult1<Vec>), "*");
    chai.add(fun(&Vector_mult2<Vec>), "*");
    chai.add(fun(&Vector_mult3<Vec>), "*");
    chai.add(fun(&Vector_mult4<Vec>), "*");
    chai.add(fun(&Vector_mult5<Vec>), "*");
    chai.add(fun(&Vector_fill<Vec>), "fill");
}


void wrap_eigen(boost::shared_ptr<ChaiScript>& chai){
    register_eigen_matrix_type<Matrix3f>(*chai,"Matrix3f");
    register_eigen_matrix_type<MatrixXf>(*chai,"MatrixXf");
    chai->add(constructor<MatrixXf(int,int)>(), "MatrixXf");
    register_eigen_vector_type<Vector3f>(*chai,"Vector3f");
    chai->add(constructor<Vector3f(double,double,double)>(),"Vector3f");
    register_eigen_vector_type<Vector3f>(*chai,"VectorXf");
    chai->add(constructor<VectorXf(int)>(), "VectorXf");
    register_eigen_vector_type<Vector3f>(*chai,"Vector2i");
}
