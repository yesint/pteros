#include "pteros/extras/membrane/lipid_tail_descr.h"
#include "pteros/core/pteros_error.h"
#include "pteros/core/logging.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

void LipidTailDescr::init(const Selection& lipid_sel, const string &tail_descr_str)
{
    // Parse tail description string
    bond_orders.clear();
    c_names.clear();
    int token_start = 0;
    for(int cur=0; cur<tail_descr_str.size(); ++cur){
        if(tail_descr_str[cur]=='-'){
            bond_orders.push_back(1);
            c_names.push_back(tail_descr_str.substr(token_start,cur-token_start));
            token_start = cur+1;
        } else if(tail_descr_str[cur]=='=') {
            bond_orders.push_back(2);
            c_names.push_back(tail_descr_str.substr(token_start,cur-token_start));
            token_start = cur+1;
        }

        if(cur==tail_descr_str.size()-1){
            // Last token
            c_names.push_back(tail_descr_str.substr(token_start,cur-token_start+1));
            // Includes last name symbol!
        }
    }


    //LOG()->debug("{}",c_names.join());

    // Set offsets of carbon atoms to use as local selection indexes
    c_offsets.resize(c_names.size());
    for(int i=0;i<c_offsets.size();++i){
        auto c_sel = lipid_sel("name "+c_names[i]);
        if(c_sel.size()!=1) throw PterosError("Incorrect tail carbon atom selection: {}!",c_names[i]);
        c_offsets[i] = c_sel.index(0) - lipid_sel.index(0);
    }
}
