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


#include "pteros/core/utilities.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/version.h"
#include <fstream>
#include <iostream>
#include <string>
// Periodic table from VMD molfile plugins
#include "periodic_table.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

namespace pteros {

void guess_element(const string& name, int& anum, float& mass){
    anum = get_pte_idx(name.c_str());
    mass = get_pte_mass(anum);
}

float angle_between_vectors(Vector3f_const_ref vec1, Vector3f_const_ref vec2)
{
    float ang = vec1.dot(vec2)/vec1.norm()/vec2.norm();
    if(ang>1.0) ang = 1.0;
    if(ang<-1.0) ang = -1.0;
    return acos(ang);
}


Vector3f project_vector(Vector3f_const_ref vec1, Vector3f_const_ref vec2)
{
    return (vec1.dot(vec2)/vec2.dot(vec2))*vec2;
}


float deg_to_rad(float ang)
{
    return ang*3.141592/180.0;
}

float rad_to_deg(float ang)
{
    return ang*180.0/3.141592;
}

Affine3f rotation_transform(Vector3f_const_ref pivot, Vector3f_const_ref axis, float angle)
{
    Affine3f m;
    m = AngleAxisf(angle,axis.normalized());
    m.translation().fill(0.0);
    return Translation3f(pivot)*m*Translation3f(-pivot);
}

string get_element_name(int elnum){
    return (elnum<nr_pte_entries && elnum>=0) ? string(pte_label[elnum]) : "X";
}

float get_vdw_radius(int elnum, const string &name) {
    if(elnum<=0){
        switch(name[0]){
        case 'H': return  0.12;
        case 'C': return  0.17;
        case 'N': return  0.155;
        case 'O': return  0.152;
        case 'S': return  0.18;
        case 'P': return  0.18;
        case 'F': return  0.147;
        default:  return  0.15;
        }
    } else {
        // Use periodic table from VMD plugins
        return (elnum<nr_pte_entries) ? 0.1*pte_vdw_radius[elnum] : 0.15;
    }
}

int get_element_number(const string &name)
{
    return get_pte_idx(name.c_str());
}

static const map<string,char> code_3to1 = {
    {"ALA",	'A'},
    {"ARG",	'R'},
    {"ASN",	'N'},
    {"ASP",	'D'},
    {"ASX",	'B'},
    {"CYS",	'C'},
    {"GLU",	'E'},
    {"GLN",	'Q'},
    {"GLX",	'Z'},
    {"GLY", 'G'},
    {"HIS",	'H'},
    {"ILE",	'I'},
    {"LEU",	'L'},
    {"LYS",	'K'},
    {"MET",	'M'},
    {"PHE",	'F'},
    {"PRO",	'P'},
    {"SER",	'S'},
    {"THR",	'T'},
    {"TRP",	'W'},
    {"TYR",	'Y'},
    {"VAL", 'V'}
};

static const map<char,string> code_1to3 = {
    {'A', "ALA"},
    {'R', "ARG"},
    {'N', "ASN"},
    {'D', "ASP"},
    {'B', "ASX"},
    {'C', "CYS"},
    {'E', "GLU"},
    {'Q', "GLN"},
    {'Z', "GLX"},
    {'G', "GLY"},
    {'H', "HIS"},
    {'I', "ILE"},
    {'L', "LEU"},
    {'K', "LYS"},
    {'M', "MET"},
    {'F', "PHE"},
    {'P', "PRO"},
    {'S', "SER"},
    {'T', "THR"},
    {'W', "TRP"},
    {'Y', "TYR"},
    {'V', "VAL"}
};

char resname_1char(const string &code)
{
    if(code_3to1.count(code))
        return code_3to1.at(code);
    else
        return 'X';
}

string resname_3char(char code)
{
    if(code_1to3.count(code))
        return code_1to3.at(code);
    else
        return "XXX";
}

void greet_string(const string& s, int sz){
    cout << "| " << s;  for(int i=0;i<sz-3-s.size();++i) cout << " "; cout << "|" << endl;
}

void greeting(string tool_name)
{
    string s1 = "|***    Pteros molecular modeling library    ***|";
    string ver_str = fmt::format("Version: {} ({})", _version_tag, _git_revision);
    string timestamp = fmt::format("Built at: {}", _build_time);
    string author = "(C) Semen Yesylevskyy, 2021";
    string cite = "Please cite: 10.1002/jcc.23943";
    string web = "Web: github.com/yesint/pteros";
    // Line
    string line="+";  for(int i=0;i<s1.size()-2;++i) line+="-";  line+="+";

    cout << line << endl;
    cout << s1 << endl;
    cout << line << endl;
    greet_string(ver_str, s1.size());
    greet_string(timestamp, s1.size());
    cout << line << endl;
    greet_string(author ,s1.size());
    greet_string(web ,s1.size());
    greet_string(cite ,s1.size());
    cout << line << endl;

    if(tool_name.size()>0){
        greet_string("Tool: "+tool_name, s1.size());
        cout << line << endl;
    }
}

void get_element_from_atom_name(const string& name, int &anum, float &mass){
    // Find first character, which is not digit to account for cases like 21C2
    int i = name.find_first_not_of("1234567890");

         if(name[i]=='C') { mass = 12.0; anum = 6; }
    else if(name[i]=='O') { mass = 16.0; anum = 8; }
    else if(name[i]=='N') { mass = 14.0; anum = 7; }
    else if(name[i]=='S') { mass = 32.0; anum = 16; }
    else if(name[i]=='H') { mass = 1.0;  anum = 1; }
    else if(name[i]=='P') { mass = 31.0; anum = 15; }
    else if(name[i]=='F') { mass = 19.0; anum = 9; }
    else if(name[i]=='B') { mass = 11.0; anum = 5; }
         else { mass = 1.0; anum = 0; } //default
}

} // namespace pteros



Histogram::Histogram(float minval, float maxval, int n): normalized(false)
{
    create(minval,maxval,n);
}

void Histogram::create(float minval, float maxval, int n)
{
    nbins = n;
    minv = minval;
    maxv = maxval;

    val.resize(nbins);
    val.fill(0.0);
    pos.resize(nbins);
    d = (maxv-minv)/float(nbins);
    for(int i=0;i<nbins;++i) pos(i) = minv+0.5*d+d*i;

    normalized = false;
}

int Histogram::get_bin(float v)
{
    return floor((v-minv)/d);
}

void Histogram::add(float v,float  weight)
{
    if(normalized) throw PterosError("Can't add value to normalized histogram!");
    int b = floor((v-minv)/d);
    if(b>=0 && b<nbins) val(b) += weight;
}

void Histogram::add(const std::vector<float>& v)
{
    for(auto& val: v) add(val);
}

void Histogram::add_cylindrical(float r, float w, float sector, float cyl_h)
{
    if(normalized) throw PterosError("Can't add value to normalized histogram!");
    int b = floor((r-minv)/d);
    if(b>=0 && b<nbins){
        float r1 = pos(b)-0.5*d;
        float r2 = pos(b)+0.5*d;
        val(b) += w/(cyl_h*sector*(r2*r2-r1*r1));
    }
}

void Histogram::normalize(float norm)
{
    if(norm){
        val /= norm;
    } else {
        val /= val.sum()*d;
    }
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

VectorXd &Histogram::values()
{
    return val;
}

Eigen::VectorXd &Histogram::positions()
{
    return pos;
}

int Histogram::num_bins() const
{
    return nbins;
}

void Histogram::save_to_file(const string &fname, float x_shift)
{
    ofstream f(fname);
    for(int i=0;i<nbins;++i) f << pos(i)+x_shift << " " << val(i) << endl;
    f.close();
}

//----------------------------------

Histogram2D::Histogram2D(float minval1, float maxval1, int n1, float minval2, float maxval2, int n2): normalized(false)
{
    create(minval1, maxval1, n1, minval2, maxval2, n2);
}

void Histogram2D::create(float minval1, float maxval1, int n1, float minval2, float maxval2, int n2)
{
    minv = Vector2f(minval1,minval2);
    maxv = Vector2f(maxval1,maxval2);
    nbins = Vector2i(n1,n2);

    val.resize(nbins(0),nbins(1));
    val.fill(0.0);

    d(0) = (maxv(0)-minv(0))/float(nbins(0));
    d(1) = (maxv(1)-minv(1))/float(nbins(1));

    normalized = false;
}

void Histogram2D::add(float v1, float v2, float weight)
{
    if(normalized) throw PterosError("Can't add value to normalized histogram!");
    int b1 = floor((v1-minv(0))/d(0));
    int b2 = floor((v2-minv(1))/d(1));
    if(b1>=0 && b1<nbins(0) && b2>=0 && b2<nbins(1)) val(b1,b2) += weight;
}

void Histogram2D::normalize(float norm)
{
    if(norm){
        val /= norm;
    } else {
        val /= val.sum()*d(0)*d(1);
    }
    normalized = true;
}

float Histogram2D::value(int i, int j) const
{
    return val(i,j);
}

void Histogram2D::save_to_file(const string &fname)
{
    ofstream f(fname);
    for(int i=0;i<nbins(0);++i){
        for(int j=0;j<nbins(1);++j)
            f << val(i,j) << " ";
        f << endl;
    }
    f.close();
}

namespace pteros {

int str_to_int(const string& str){
    /*
    int val;
    const auto res = std::from_chars(str.data(), str.data()+str.size(),val);
    if (res.ec == std::errc::invalid_argument) throw PterosError("Value '{}' is not INT!",str);
    return val;
    */
    //from_chars is not implemented for floats in gcc...
    stringstream ss(str);
    int val;
    ss >> val;
    return val;
}

float str_to_float(const string& str){
    /*
    float val;
    const auto res = std::from_chars(str.data(), str.data()+str.size(),val);
    if (res.ec == std::errc::invalid_argument) throw PterosError("Value '{}' is not FLOAT!",str);
    return val;
    */
    //from_chars is not implemented for floats in gcc...
    stringstream ss(str);
    float val;
    ss >> val;
    return val;
}

void str_to_lower_in_place(string& str){
    std::locale loc;
    for(auto elem: str) tolower(elem,loc);
}

void str_to_upper_in_place(string& str){
    std::locale loc;
    for(auto elem: str) toupper(elem,loc);
}


string str_to_lower_copy(const string &str)
{
    string s(str);
    str_to_lower_in_place(s);
    return s;
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
void str_trim_in_place(std::string &s) {
    ltrim(s);
    rtrim(s);
}

} // namespace pteros
