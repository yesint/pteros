/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
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

#include "pteros/core/utilities.h"
#include "pteros/core/pteros_error.h"
#include <fstream>

using namespace std;
using namespace pteros;
using namespace Eigen;

float pteros::angle_between_vectors(Vector3f_const_ref vec1, Vector3f_const_ref vec2)
{
    float ang = vec1.dot(vec2)/vec1.norm()/vec2.norm();
    if(ang>1.0) ang = 1.0;
    if(ang<-1.0) ang = -1.0;
    return acos(ang);
}


Vector3f pteros::project_vector(Vector3f_const_ref vec1, Vector3f_const_ref vec2)
{
    return (vec1.dot(vec2)/vec2.dot(vec2))*vec2;
}


float pteros::rad_to_deg(float ang)
{
    return ang*3.141592/180.0;
}

float deg_to_rad(float ang)
{
    return ang*180.0/3.141592;
}

Histogram::Histogram(float minval, float maxval, int n): minv(minval), maxv(maxval), nbins(n), normalized(false)
{
    val.resize(nbins);
    val.fill(0.0);
    pos.resize(nbins);
    d = (maxv-minv)/float(nbins);
    for(int i=0;i<nbins;++i) pos(i) = minv+0.5*d+d*i;
}

void Histogram::add(float v)
{
    if(normalized) throw Pteros_error("Can't add value to normalized histogram!");
    int b = floor((v-minv)/d);
    if(b>=0 && b<nbins) val(b) += 1.0;
}

void Histogram::normalize()
{
    val.normalize();
    normalized = true;
}

float Histogram::value(int i) const
{
    return val[i];
}

float Histogram::position(int i) const
{
    return pos[i];
}

const VectorXd &Histogram::values() const
{
    return val;
}

const VectorXd &Histogram::positions() const
{
    return pos;
}

int Histogram::num_bins() const
{
    return nbins;
}

void Histogram::save_to_file(const string &fname)
{
    ofstream f(fname);
    for(int i=0;i<nbins;++i) f << pos(i) << " " << val(i) << endl;
    f.close();
}


