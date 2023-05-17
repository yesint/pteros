#include "pteros/core/fitting.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"
#include "pteros/core/subset.h"

using namespace std;
using namespace Eigen;
using namespace pteros;

namespace pteros {

/// The code below is hacked from GROMACS 3.3.3
/// Used to compute the rotation matrix
/// It computes the rot matrix for two selections, which
/// are centerd at zero.
Matrix3f rot_transform_matrix(const Selection& sel1, const Selection& sel2){
    Matrix<float,6,6> omega,om;
    Matrix3f u,vh,vk;

    omega.fill(0.0);
    om.fill(0.0);
    int N = sel1.size();

    //Calculate the matrix U
    u.fill(0.0);
    for(int i=0;i<N;++i){ // Over atoms in selection
        u += sel1.xyz(i)*sel2.xyz(i).transpose()*sel1.mass(i);
    }

    //Construct omega
    for(int r=0; r<6; r++){
        for(int c=0; c<=r; c++){
            if (r>=3 && c<3) {
                omega(r,c)=u(r-3,c);
                omega(c,r)=u(r-3,c);
            } else {
                omega(r,c)=0;
                omega(c,r)=0;
            }
        }
    }

    //Finding eigenvalues of omega
    Eigen::SelfAdjointEigenSolver<Matrix<float,6,6>> solver(omega);
    om = solver.eigenvectors();

    /*  Copy only the first two eigenvectors
        The eigenvectors are already sorted ascending by their eigenvalues!
    */
    for(int j=0; j<2; j++){
        for(int i=0; i<3; i++) {
            vh(j,i)=sqrt(2.0)*om(i,5-j);
            vk(j,i)=sqrt(2.0)*om(i+3,5-j);
        }
    }

    // Calculate the last eigenvector as the cross-product of the first two.
    // This insures that the conformation is not mirrored and
    // prevents problems with completely flat reference structures.
    vh.row(2) = vh.row(0).cross(vh.row(1)) ;
    vk.row(2) = vk.row(0).cross(vk.row(1)) ;

    /* Determine rotational part */
    Matrix3f rot;
    for(int r=0; r<3; r++){
        for(int c=0; c<3; c++){
            rot(c,r) = vk(0,r)*vh(0,c) + vk(1,r)*vh(1,c) + vk(2,r)*vh(2,c);
        }
    }

    return rot;
}

/*
 *Code to validate multiplication:

#include <Eigen/Core>
#include <iostream>

using namespace Eigen;
using namespace std;

int main() {
    MatrixXf m1(3,4),m2(3,4);
    m1 << 1,4,7,1,
          2,5,8,2,
          3,6,9,3;
    m2 << 1,4,7,1,
          2,5,8,2,
          3,6,9,3;

    cout << m1*m2.transpose() << endl;

    Matrix3f res;
    res.fill(0);
    for(int i=0;i<4;++i){
        res += m1.col(i)*m2.col(i).transpose();
    }
    cout << res << endl;

    return 0;
}
 */

template <typename CrdT, typename MassT>
Matrix3f rot_transform_matrix_new(const MatrixBase<CrdT>& coord1,
                              const MatrixBase<CrdT>& coord2,
                              const MatrixBase<MassT>& masses){

    static_assert(masses.IsVectorAtCompileTime,"Vector is expected for masses");
    static_assert(coord1.RowsAtCompileTime==3,"Matrix<3,X> required!");
    static_assert(coord2.RowsAtCompileTime==3,"Matrix<3,X> required!");
    //fmt::print("{} {}\n",coord1.rows(),coord1.cols());

    Matrix<float,6,6> omega,om;
    Matrix3f u,vh,vk;

    omega.fill(0.0);
    om.fill(0.0);
    int N = coord1.size();

    //Calculate the matrix U
    /*
    u.fill(0.0);
    for(int i=0;i<N;++i){ // Over atoms in selection
        u += sel1.xyz(i)*sel2.xyz(i).transpose()*sel1.mass(i);
    }
    */
    u = (coord1.array().rowwise() * masses.array().transpose()).matrix() * coord2.transpose();

    //Construct omega
    /*
        1 4 7
        2 5 8
        3 6 9
    ==>
    0 0 0 1 2 3
    0 0 0 4 5 6
    0 0 0 7 8 9
    1 4 7 0 0 0
    2 5 8 0 0 0
    3 6 9 0 0 0

    for(int r=0; r<6; r++){
        for(int c=0; c<=r; c++){
            if (r>=3 && c<3) {
                omega(r,c)=u(r-3,c);
                omega(c,r)=u(r-3,c);
            } else {
                omega(r,c)=0;
                omega(c,r)=0;
            }
        }
    }
    */
    omega.block<3,3>(0,3) = u.transpose();
    omega.block<3,3>(3,0) = u;
    omega.block<3,3>(0,0) = Matrix3f::Zero();
    omega.block<3,3>(3,3) = Matrix3f::Zero();

    //Finding eigenvalues of omega
    Eigen::SelfAdjointEigenSolver<Matrix<float,6,6>> solver(omega);
    om = solver.eigenvectors();

    /*  Copy only the first two eigenvectors
        The eigenvectors are already sorted ascending by their eigenvalues!
    */
    for(int j=0; j<2; j++){
        for(int i=0; i<3; i++) {
            vh(j,i)=sqrt(2.0)*om(i,5-j);
            vk(j,i)=sqrt(2.0)*om(i+3,5-j);
        }
    }

    // Calculate the last eigenvector as the cross-product of the first two.
    // This insures that the conformation is not mirrored and
    // prevents problems with completely flat reference structures.
    vh.row(2) = vh.row(0).cross(vh.row(1)) ;
    vk.row(2) = vk.row(0).cross(vk.row(1)) ;

    /* Determine rotational part */
    Matrix3f rot;
    for(int r=0; r<3; r++){
        for(int c=0; c<3; c++){
            rot(c,r) = vk(0,r)*vh(0,c) + vk(1,r)*vh(1,c) + vk(2,r)*vh(2,c);
        }
    }

    return rot;
}


// Fitting transformation for two selections.
// If translate_to_zero=false it is presumed that selections are brough to zero already
Affine3f fit_transform(const Selection& sel1,
                       const Selection& sel2,
                       bool translate_to_zero){
    Affine3f tr; // result
    Vector3f cm1, cm2;
    int n1 = sel1.size();
    int n2 = sel2.size();

    if(n1!=n2) throw PterosError("Incompatible selections for fitting of sizes {} and {}", n1, n2);

    // Bring centers to zero
    if(translate_to_zero){
        cm1 = sel1.center(true);
        cm2 = sel2.center(true);
        const_cast<Selection&>(sel1).translate(-cm1);
        const_cast<Selection&>(sel2).translate(-cm2);
    }

    // Determine rotational part of transform
    tr.linear() = rot_transform_matrix(sel1,sel2);

    /*
    Subset s1(const_cast<Selection&>(sel1));
    Subset s2(const_cast<Selection&>(sel2));
    auto masses = sel1.get_mass();
    tr.linear() = rot_transform_matrix_new(s1.data(),s2.data(),
                    Map<VectorXf>(masses.data(),masses.size()));
    */

    //Clear translation part.
    tr.translation().fill(0.0);

    if(translate_to_zero){
        // Bring centers back
        const_cast<Selection&>(sel1).translate(cm1);
        const_cast<Selection&>(sel2).translate(cm2);
        // Add translation part to transform. Note reverse order of translations! This is important.
        tr = Translation3f(cm2) * tr * Translation3f(-cm1);
    }

    return tr;
}


float simple_rmsd(const Selection& sel1, const Selection& sel2) {
    float res = 0.0;
    size_t N = sel1.size();
    for(size_t i=0; i<N; ++i){
        res += (sel1.xyz(i)-sel2.xyz(i)).squaredNorm();
    }
    return sqrt(res/N);
}

/// Computes RMSD matrix over a trajectory
/// The trajectory have to be unwrapped properly before!
MatrixXf rmsd_matrix(Selection& sel)
{
    // Second selection
    Selection sel2 = sel;

    size_t N = sel.get_system()->num_frames();
    // Prepare matrix
    MatrixXf res(N,N);
    res.diagonal().fill(0.0);

    // Bring to zero for each frame
    // and save positions to restore later
    MatrixXf cm(3,N);
    for(size_t fr=0; fr<N; ++fr){
        sel.set_frame(fr);
        cm.col(fr) = sel.center(true);
        sel.translate(-cm.col(fr));
    }

    size_t Npairs = N*(N-1)/2;
    size_t cur=0;

    // Do fitting and rmsd
    for(size_t i=0; i<N-1; ++i){
        sel.set_frame(i);
        for(size_t j=i+1; j<N; ++j){
            sel2.set_frame(j);
            auto tr = fit_transform(sel2,sel,false); // Already brought to zero
            // Store coordinates
            auto crd = sel2.get_xyz();
            //sel2.apply_transform(tr);
            sel2.set_xyz( tr.linear()*crd ); // We dont have translation, so just apply rotation part of transform
            res(j,i) = res(i,j) = simple_rmsd(sel,sel2);
            // Restore coordinates
            sel2.set_xyz(crd);
            ++cur;
            if(cur%(Npairs/100)==0) LOG()->info("RMSD pair {} of {} ({}%)",cur,Npairs,100*cur/Npairs);
        }
    }

    // Restore centers
    for(size_t fr=0; fr<N; ++fr){
        sel.set_frame(fr);
        sel.translate(cm.col(fr));
    }

    return res;
}

} // namespace pteros
