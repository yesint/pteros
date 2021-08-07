/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#include <fstream>
#include <sstream>
#include <cmath>
#include "pteros/core/periodic_box.h"
#include "pteros/core/selection.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

using namespace std;
using namespace pteros;
using namespace Eigen;

PeriodicBox::PeriodicBox():_is_periodic(false),_is_triclinic(false),_box(Eigen::Matrix3f::Zero()){}

PeriodicBox::PeriodicBox(Matrix3f_const_ref m){
    set_matrix(m);
}

PeriodicBox::PeriodicBox(Vector3f_const_ref vectors, Vector3f_const_ref angles)
{
    from_vectors_angles(vectors,angles);
}

PeriodicBox::PeriodicBox(const PeriodicBox &other)
{
    if(&other == this) return;
    set_matrix(other._box);
}

PeriodicBox& PeriodicBox::operator=(const PeriodicBox& other){
    if(&other == this) return *this;
    set_matrix(other._box);
    return *this;
}

float PeriodicBox::get_element(int i, int j) const
{
    return _box(i,j);
}

void PeriodicBox::set_element(int i, int j, float val)
{
    _box(i,j) = val;
    recompute_internals();
}

void PeriodicBox::set_matrix(Matrix3f_const_ref matr)
{
    _box = matr;
    recompute_internals();
}

Vector3f PeriodicBox::get_vector(int i) const{
    return _box.col(i);
}

void PeriodicBox::set_vector(Vector3f_const_ref vec, int i)
{
    _box.col(i) = vec;
    recompute_internals();
}

Matrix3f PeriodicBox::get_matrix() const {
    return _box;
}

void PeriodicBox::scale_vectors(Vector3f_const_ref scale)
{
    if(!_is_periodic) throw PterosError("No periodicity! Can't scale!");
    for(int i=0;i<3;i++) _box.col(i) *= scale(i);
    _box_inv = _box.inverse();
}

Matrix3f PeriodicBox::get_inv_matrix() const {
    return _box_inv;
}

Vector3f PeriodicBox::lab_to_box(Vector3f_const_ref point) const {
    return _box_inv.colwise().normalized()*point;
}

Matrix3f PeriodicBox::lab_to_box_transform() const {
    return _box_inv.colwise().normalized();
}

Vector3f PeriodicBox::box_to_lab(Vector3f_const_ref point) const{
    return _box.colwise().normalized()*point;
}

Matrix3f PeriodicBox::box_to_lab_transform() const {
    return _box.colwise().normalized();
}

float PeriodicBox::extent(int i) const {
    return _box.col(i).norm();
}

Vector3f PeriodicBox::extents() const {
    return _box.colwise().norm();
}

float PeriodicBox::distance_squared(Vector3f_const_ref point1, Vector3f_const_ref point2, Array3i_const_ref pbc) const
{
    return shortest_vector(point1,point2,pbc).squaredNorm();
}

float PeriodicBox::distance(Vector3f_const_ref point1, Vector3f_const_ref point2, Array3i_const_ref pbc) const
{
    return shortest_vector(point1,point2,pbc).norm();
}

float PeriodicBox::volume(){
    return _box.col(1).cross( _box.col(2) ).dot( _box.col(0) );
}

void PeriodicBox::wrap_point(Vector3f_ref point, Array3i_const_ref pbc, Vector3f_const_ref origin) const
{
    point = _box_inv*(point-origin);
    for(int i=0;i<3;++i){
        if(pbc(i)!=0){
            point(i) -= round(point(i));
            if(point(i)<0) point(i) += 1.0;
        }
    }
    point = _box*point + origin;
}

bool PeriodicBox::in_box(Vector3f_const_ref point, Vector3f_const_ref origin) const
{
    Vector3f p = _box_inv*(point-origin);

    for(int i=0; i<3; ++i)
        if(p(i)<0 || p(i)>1.0) return false;

    return true;
}


Eigen::Vector3f PeriodicBox::closest_image(Vector3f_const_ref point, Vector3f_const_ref target, Array3i_const_ref pbc) const
{    
    return target + shortest_vector(target,point,pbc);
}

Vector3f PeriodicBox::shortest_vector(Vector3f_const_ref point1, Vector3f_const_ref point2, Array3i_const_ref pbc) const
{
    if(_is_periodic){
        Vector3f d = _box_inv*(point2-point1);
        for(int i=0;i<3;++i){
            if(pbc(i)!=0){
                d(i) -= round(d(i));
            }
        }
        return _box*d;
    } else {
        return point2-point1;
    }
}

// The code below is hacked from Gromacs 3.3.x and modified heavily

#define epbcSCREW 1111
#define XX 0
#define YY 1
#define ZZ 2
#define DEG2RAD (M_PI/180.0)
#define RAD2DEG (180.0/M_PI)

static inline float cos_angle_no_table(const Eigen::Vector3f a,const Eigen::Vector3f b)
{
  /* This version does not need the invsqrt lookup table */
  float   cos;
  int    m;
  double aa,bb,ip,ipa,ipb; /* For accuracy these must be double! */

  ip=ipa=ipb=0.0;
  for(m=0; (m<3); m++) {
    aa   = a[m];
    bb   = b[m];
    ip  += aa*bb;
    ipa += aa*aa;
    ipb += bb*bb;
  }
  cos=ip/sqrt(ipa*ipb);

  if (cos > 1.0)
    return  1.0;
  if (cos <-1.0)
    return -1.0;

  return cos;
}

void PeriodicBox::from_pdb_box(const char *line)
{
#define SG_SIZE 11
    char sa[12],sb[12],sc[12],sg[SG_SIZE+1],ident;
    double fa,fb,fc,alpha,beta,gamma,cosa,cosb,cosg,sing;
    int  syma,symb,symc;
    int  ePBC_file;
    Matrix3f box;

    sscanf(line,"%*s%s%s%s%lf%lf%lf",sa,sb,sc,&alpha,&beta,&gamma);

    ePBC_file = -1;
    if (strlen(line) >= 55) {
        strncpy(sg,line+55,SG_SIZE);
        sg[SG_SIZE] = '\0';
        ident = ' ';
        syma  = 0;
        symb  = 0;
        symc  = 0;
        sscanf(sg,"%c %d %d %d",&ident,&syma,&symb,&symc);
        if (ident == 'P' && syma == 21 && symb == 1 && symc == 1)
          ePBC_file = epbcSCREW;
    }

    fa = atof(sa);//*0.1;
    fb = atof(sb);//*0.1;
    fc = atof(sc);//*0.1;
    if (ePBC_file == epbcSCREW) {
      fa *= 0.5;
    }

    box.fill(0.0);
    box(XX,XX) = fa;

    if ((alpha!=90.0) || (beta!=90.0) || (gamma!=90.0)) {
      if (alpha != 90.0) {
    cosa = cos(alpha*DEG2RAD);
      } else {
    cosa = 0;
      }
      if (beta != 90.0) {
    cosb = cos(beta*DEG2RAD);
      } else {
    cosb = 0;
      }
      if (gamma != 90.0) {
    cosg = cos(gamma*DEG2RAD);
    sing = sin(gamma*DEG2RAD);
      } else {
    cosg = 0;
    sing = 1;
      }
      box(YY,XX) = fb*cosg;
      box(YY,YY) = fb*sing;
      box(ZZ,XX) = fc*cosb;
      box(ZZ,YY) = fc*(cosa - cosb*cosg)/sing;
      box(ZZ,ZZ) = sqrt(fc*fc
             - box(ZZ,XX)*box(ZZ,XX) - box(ZZ,YY)*box(ZZ,YY));
    } else {
      box(YY,YY) = fb;
      box(ZZ,ZZ) = fc;
    }

    // We obtained box as a set of row-vectors. Transform it to column vectors
    box.transposeInPlace();

    // Init object
    set_matrix(box);
}

// This function works with column-ordered box!
void PeriodicBox::to_vectors_angles(Vector3f_ref vectors, Vector3f_ref angles) const {
    //
    Eigen::Matrix3f boxT = _box.transpose();
    if (boxT.row(YY).squaredNorm()*boxT.row(ZZ).squaredNorm()!=0)
        angles(0) = RAD2DEG*acos(cos_angle_no_table(boxT.row(YY),boxT.row(ZZ)));
    else
        angles(0) = 90.0;
    if (boxT.row(XX).squaredNorm()*boxT.row(ZZ).squaredNorm()!=0)
        angles(1)  = RAD2DEG*acos(cos_angle_no_table(boxT.row(XX),boxT.row(ZZ)));
    else
        angles(1) = 90.0;
    if (boxT.row(XX).squaredNorm()*boxT.row(YY).squaredNorm()!=0)
        angles(2) = RAD2DEG*acos(cos_angle_no_table(boxT.row(XX),boxT.row(YY)));
    else
        angles(2) = 90.0;

    vectors(0) = boxT.row(XX).norm();
    vectors(1) = boxT.row(YY).norm();
    vectors(2) = boxT.row(ZZ).norm();
}

std::string PeriodicBox::to_pdb_box() const {
  float alpha,beta,gamma;
  char ch[80];

  Eigen::Vector3f vectors, angles;
  to_vectors_angles(vectors,angles);

  //sprintf(ch,"REMARK    THIS IS A SIMULATION BOX\n");

  //if (ePBC != epbcSCREW) {
  sprintf(ch,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
        vectors(0),vectors(1),vectors(2),angles(0),angles(1),angles(2),"P 1",1);
  /*
  } else {
    // Double the a-vector length and write the correct space group
    fprintf(out,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
        20*norm(box[XX]),10*norm(box[YY]),10*norm(box[ZZ]),
        alpha,beta,gamma,"P 21 1 1",1);

  }
  */
  return std::string(ch);
}

void PeriodicBox::from_vectors_angles(Vector3f_const_ref vectors, Vector3f_const_ref angles){
    double cosa,cosb,cosg,sing;

    Matrix3f box;

    box.fill(0.0);
    box(XX,XX) = vectors(0);

    if ((angles(0)!=90.0) || (angles(1)!=90.0) || (angles(2)!=90.0)) {
        if (angles(0) != 90.0) {
            cosa = cos(angles(0)*DEG2RAD);
        } else {
            cosa = 0;
        }
        if (angles(1) != 90.0) {
            cosb = cos(angles(1)*DEG2RAD);
        } else {
            cosb = 0;
        }
        if (angles(2) != 90.0) {
            cosg = cos(angles(2)*DEG2RAD);
            sing = sin(angles(2)*DEG2RAD);
        } else {
            cosg = 0;
            sing = 1;
        }
        box(YY,XX) = vectors(1)*cosg;
        box(YY,YY) = vectors(1)*sing;
        box(ZZ,XX) = vectors(2)*cosb;
        box(ZZ,YY) = vectors(2)*(cosa - cosb*cosg)/sing;
        box(ZZ,ZZ) = sqrt(vectors(2)*vectors(2)
                          - box(ZZ,XX)*box(ZZ,XX) - box(ZZ,YY)*box(ZZ,YY));
    } else {
        box(YY,YY) = vectors(1);
        box(ZZ,ZZ) = vectors(2);
    }

    // We obtained box as a set of row-vectors. Transform it to column vectors
    box.transposeInPlace();

    // Recompute internals
    set_matrix(box);
}

void PeriodicBox::recompute_internals()
{
    _is_periodic = (_box.array()!=0).any() ? true : false;
    if(!_is_periodic) return;
    _box_inv = _box.inverse();
    _is_triclinic = (_box(0,1)||_box(0,2)||_box(1,0)||_box(1,2)||_box(2,0)||_box(2,1));
}




