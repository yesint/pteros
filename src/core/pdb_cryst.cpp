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

#include "pteros/core/pdb_cryst.h"
#include <ctype.h>
#include <cstdio>
#include <Eigen/Geometry> // For M_PI constant

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

void pteros::read_pdb_box(const char *line, Eigen::Matrix3f& box)
{
#define SG_SIZE 11
    char sa[12],sb[12],sc[12],sg[SG_SIZE+1],ident;
    double fa,fb,fc,alpha,beta,gamma,cosa,cosb,cosg,sing;
    int  syma,symb,symc;
    int  ePBC_file;

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
}

// This function works with column-ordered box!
void pteros::box_to_vectors_angles(const Eigen::Matrix3f& box,
                               Eigen::Vector3f& vectors, Eigen::Vector3f& angles){
    //
    Eigen::Matrix3f boxT = box.transpose();
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

std::string pteros::write_pdb_box(const Eigen::Matrix3f& box)
{
  float alpha,beta,gamma;
  char ch[80];

  Eigen::Vector3f vectors, angles;
  box_to_vectors_angles(box,vectors,angles);

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
