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

#include "pteros/python/compiled_plugin.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/grid_search.h"
#include <fstream>

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
                "\t\tSelection text"
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
    }     

    void process_frame(const pteros::Frame_info &info){
        if(info.valid_frame==0){
            max_dist = 0.5*options("max",
                               boost::lexical_cast<string>(system.Box(0).extents().minCoeff())
                              ).as_float();
            bin_sz = (max_dist-min_dist)/float(n_bins);
        }

        float d;

        sel1.apply();
        sel2.apply();

        density = (sel1.size()*sel2.size())/system.Box(0).volume();

        vector<Vector2i> bon;
        vector<float> dist_vec;
        Grid_searcher(max_dist,sel1,sel2,bon,true,true,&dist_vec);


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
