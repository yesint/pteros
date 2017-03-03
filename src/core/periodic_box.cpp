/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#include <iostream>
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

Periodic_box::Periodic_box():_is_periodic(false),_is_triclinic(false),_box(Eigen::Matrix3f::Zero()){}

Periodic_box::Periodic_box(Matrix3f_const_ref box){
    modify(box);
}

Periodic_box::Periodic_box(Vector3f_const_ref vectors, Vector3f_const_ref angles)
{
    from_vectors_angles(vectors,angles);
}

Periodic_box& Periodic_box::operator=(Periodic_box other){
    modify(other._box);
    return *this;
}

void Periodic_box::modify(Matrix3f_const_ref box)
{
    _box = box;
    _is_periodic = (_box.array().abs().sum()>0) ? true : false;
    if(!_is_periodic) return;
    _box_inv = _box.inverse();
    _is_triclinic = (_box(0,1)||_box(0,2)||_box(1,0)||_box(1,2)||_box(2,0)||_box(2,1));
}

Vector3f Periodic_box::get_vector(int i){
    return _box.col(i);
}

Matrix3f Periodic_box::get_matrix() const {
    return _box;
}

void Periodic_box::scale_vectors(Vector3f_const_ref scale)
{
    if(!_is_periodic) throw Pteros_error("No periodicity! Can't scale!");
    for(int i=0;i<3;i++) _box.col(i) *= scale(i);
    _box_inv = _box.inverse();
}

Matrix3f Periodic_box::get_inv_matrix() const {
    return _box_inv;
}

Vector3f Periodic_box::lab_to_box(Vector3f_const_ref point) const {
    return _box_inv.colwise().normalized()*point;
}

Matrix3f Periodic_box::lab_to_box_transform() const {
    return _box_inv.colwise().normalized();
}

Vector3f Periodic_box::box_to_lab(Vector3f_const_ref point) const{
    return _box.colwise().normalized()*point;
}

Matrix3f Periodic_box::box_to_lab_transform() const {
    return _box.colwise().normalized();
}

float Periodic_box::extent(int i) const {
    return _box.col(i).norm();
}

Vector3f Periodic_box::extents() const {
    return _box.colwise().norm();
}

float Periodic_box::distance_squared(Vector3f_const_ref point1, Vector3f_const_ref point2, Vector3i_const_ref dims) const
{
    return shortest_vector(point1,point2,dims).squaredNorm();
}

float Periodic_box::distance(Vector3f_const_ref point1, Vector3f_const_ref point2, Vector3i_const_ref dims) const
{
    return shortest_vector(point1,point2,dims).norm();
}

float Periodic_box::volume(){
    return _box.col(1).cross( _box.col(2) ).dot( _box.col(0) );
}

void Periodic_box::wrap_point(Vector3f_ref point, Vector3i_const_ref dims) const
{
    point = _box_inv*point;
    for(int i=0;i<3;++i){
        if(dims(i)!=0){
            point(i) -= round(point(i));
            if(point(i)<0) point(i) += 1.0;
        }
    }
    point = _box*point;
}

bool Periodic_box::in_box(Vector3f_const_ref point, Vector3f_const_ref origin) const
{
    Vector3f p = _box_inv*(point-origin);

    for(int i=0; i<3; ++i)
        if(p(i)<0 || p(i)>1.0) return false;

    return true;
}


Eigen::Vector3f Periodic_box::get_closest_image(Vector3f_const_ref point, Vector3f_const_ref target, Vector3i_const_ref dims) const
{    
    return target + shortest_vector(target,point,dims);
}

Vector3f Periodic_box::shortest_vector(Vector3f_const_ref point1, Vector3f_const_ref point2, Vector3i_const_ref dims) const
{
    if(_is_periodic){
        Vector3f d = _box_inv*(point2-point1);
        for(int i=0;i<3;++i){
            if(dims(i)!=0){
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

void Periodic_box::read_pdb_box(const char *line)
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
    modify(box);
}

// This function works with column-ordered box!
void Periodic_box::to_vectors_angles(Vector3f_ref vectors, Vector3f_ref angles) const {
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

std::string Periodic_box::write_pdb_box() const {
  float alpha,beta,gamma;
  char ch[80];

  Eigen::Vector3f vectors, angles;
  to_vectors_angles(vectors,angles);

  sprintf(ch,"REMARK    THIS IS A SIMULATION BOX\n");

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

void Periodic_box::from_vectors_angles(Vector3f_const_ref vectors, Vector3f_const_ref angles){
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
    modify(box);
}
