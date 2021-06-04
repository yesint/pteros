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


#include "pteros/python/compiled_plugin.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/distance_search.h"
#include <boost/lexical_cast.hpp>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

using namespace std;
using namespace pteros;
using namespace Eigen;

class rdf: public pteros::Compiled_plugin_base {
public:
    rdf(pteros::Trajectory_processor* pr, const pteros::Options& opt): Compiled_plugin_base(pr,opt) {}

    string help(){
        return  "Purpose:\n"
                "\tComputes radial distribution function of two selections.\n"
                "Output:\n"
                "\tFile <label>.dat containing the following columns:\n"
                "\ttime RMSD\n"
                "\tAlso reports mean RMSD in the file header.\n"
                "Options:\n"
                "\t-selection <string>\n"
                "\t\tSelection text\n"
                "\t-nojump <true|false>. Default: true\n"
                "\t\tRemove jumps of atoms over periodic box boundary.\n"
                "\t\tSetting this to false is only meaningful is very special cases."
                ;
    }

protected:

    void pre_process(){                        
        sel1.modify(system,options("sel1").as_string() );
        sel2.modify(system,options("sel2").as_string() );
        do_cm = options("cm","false").as_bool();
        n_bins = options("bins","50").as_int();
        min_dist = options("min","0.0").as_float();
        // Max dist will be set later when we get the box dimension
        mass_weighted = options("mass_weighted","true").as_bool();

        data.resize(n_bins);
        data.fill(0.0);

        max_dist = 0.5*options("max",
                           boost::lexical_cast<string>(system.Box(0).extents().minCoeff())
                          ).as_float();
        bin_sz = (max_dist-min_dist)/float(n_bins);
    }     

    void process_frame(const pteros::Frame_info &info){        
        float d;

        sel1.apply();
        sel2.apply();

        density = (sel1.size()*sel2.size())/system.Box(0).volume();

        vector<Vector2i> bon;
        vector<float> dist_vec;
        search_contacts(max_dist,sel1,sel2,bon,true,true,&dist_vec);


        // Atom-atom rdf
        int bin;
        float r;
        for(int i=0;i<dist_vec.size();++i){
            bin = int(floor(dist_vec[i]/bin_sz));
            if(bin<n_bins){
                r = (0.5+bin)*bin_sz;
                data[bin] += 1.0/( 4.0*M_PI*r*r*density*bin_sz );
            }
        }

    }

    void post_process(const pteros::Frame_info &info){        
        // Output
        string fname = label+".dat";                      

        ofstream f(fname.c_str());
        f << "# RDF of selections [" << sel1.get_text() << "] and ["
          << sel2.get_text() << "]"
          << endl;        
        for(int i=0; i<data.size(); ++i){
            float r = 0.5*bin_sz+i*bin_sz;
            f << r << " " << data[i]/(float)info.valid_frame << endl;
        }
        f.close();
    }

private:
    VectorXf data;
    pteros::Selection sel1, sel2;
    bool do_cm;
    bool mass_weighted;
    float bin_sz;
    int n_bins;
    float min_dist, max_dist;
    float density;
};


CREATE_COMPILED_PLUGIN(rdf)




