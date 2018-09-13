/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * (C) 2009-2018, Semen Yesylevskyy
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


#include "pteros/extras/substructure_search.h"
#include "pteros/core/pteros_error.h"
#include <set>

#include <openbabel/mol.h>
#include <openbabel/isomorphism.h>
#include "openbabel/query.h"
#include "openbabel/obconversion.h"


using namespace std;
using namespace pteros;
using namespace Eigen;

namespace pteros {

std::vector<std::vector<int>> find_equivalent_atoms(const Selection& sel)
{
    std::vector<std::vector<int>> res;

    OpenBabel::OBMol mol;
    selection_to_obmol(sel,mol);

    OpenBabel::OBConversion conv;
    conv.WriteFile(&mol,"a.smi");

    std::vector<OpenBabel::OBIsomorphismMapper::Mapping> aut;
    OpenBabel::FindAutomorphisms(&mol,aut);

    vector<set<int>> sym(sel.size());

    for(int i=0;i<aut.size();++i){
        for(int j=0;j<aut[i].size();++j){
            sym[aut[i][j].first].insert(aut[i][j].second);
        }
    }

    for(int i=0;i<sym.size();++i){
        for(int a: sym[i]){
            if(i!=a) sym[a].clear();
        }
    }

    for(int i=0;i<sym.size();++i){
        if(!sym[i].empty()){
            res.emplace_back();
            copy(sym[i].begin(),sym[i].end(),back_inserter(res.back()));
        }
    }

    return res;
}


vector<vector<int>> find_substructures(const Selection& source, const Selection& query, bool find_all)
{
    vector<vector<int>> res;

    OpenBabel::OBMol src,sample;
    selection_to_obmol(source,src);
    selection_to_obmol(query,sample);

    OpenBabel::OBQuery* obquery = OpenBabel::CompileMoleculeQuery(&sample);
    OpenBabel::OBIsomorphismMapper *mapper = OpenBabel::OBIsomorphismMapper::GetInstance(obquery);
    OpenBabel::OBIsomorphismMapper::Mappings maps;
    if(find_all){
        mapper->MapUnique(&src,maps);
    } else {
        maps.resize(1);
        mapper->MapFirst(&src,maps[0]);
    }

    delete mapper;
    delete obquery;

    res.resize(maps.size());
    for(int i=0;i<maps.size();++i){
        res[i].reserve(maps[i].size());
        for(auto& el: maps[i]){
            res[i].push_back(el.second);
        }
    }

    return res;
}

} // namespace
