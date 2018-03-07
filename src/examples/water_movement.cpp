#include <fstream>
#include "pteros/core/grid.h"
#include "pteros/core/pteros_error.h"
#include "pteros/analysis/trajectory_reader.h"
#include "pteros/analysis/task_plugin.h"

using namespace std;
using namespace pteros;
using namespace Eigen;

TASK_SERIAL(Water_processor)
protected:
    void pre_process() override {
        // Create selection for water        
        water.modify(system,std::string("resname SOL and name OW"));

        // Allocate a grid
        // Get size of the cell
        float cell_size = options("cell_size").as_float();
        // Get extents of the box
        Vector3f box_dim = system.box(0).extents();


        NgridX = floor(box_dim(0)/cell_size);
        NgridY = floor(box_dim(1)/cell_size);
        NgridZ = floor(box_dim(2)/cell_size);

        if(NgridX==0) NgridX = 1;
        if(NgridY==0) NgridY = 1;
        if(NgridZ==0) NgridZ = 1;

        //NgridX = NgridY = NgridZ = 3;

        searcher.resize(NgridX, NgridY, NgridZ);

        // Create grid for result
        grid.resize(boost::extents[NgridX][NgridY][NgridZ]);
        // Clear it
        int i,j,k;
        for(i=0;i<NgridX;++i)
            for(j=0;j<NgridY;++j)
                for(k=0;k<NgridZ;++k)
                    grid[i][j][k] = 0.0;


    }


    void process_frame(const Frame_info& info) override {
        int i,j,k,n,w,ind;

        if(info.valid_frame==0){
            water.get_xyz(last_pos);
            return;
        }

        // Assign waters to grid
        searcher.populate_periodic(water,false);
        // Cycle over grid cells
        for(i=0;i<NgridX;++i)
            for(j=0;j<NgridY;++j)
                for(k=0;k<NgridZ;++k){
                    // Get waters in this cell
                    n = searcher.cell(i,j,k).size();
                    for(w=0;w<n;++w){
                        ind = searcher.cell(i,j,k)[w].index;
                        // For current water get delta of prev and cur coorfinates
                        // And add it to result grid

                        grid[i][j][k] += pow( system.box(0).distance(water.xyz(ind),last_pos.col(ind)) ,2);
                        //grid[i][j][k] += 1;
                    }
                }

        water.get_xyz(last_pos);        
    }


    void post_process(const Frame_info& info) override {
        int i,j,k;
        for(i=0;i<NgridX;++i)
            for(j=0;j<NgridY;++j)
                for(k=0;k<NgridZ;++k){
                    grid[i][j][k] /= (float)info.valid_frame;
                }

        // Output
        ofstream f("map.dx");
        f << "# Made by Pteros" << endl;
        f << "object 1 class gridpositions counts " << NgridX << " " << NgridY << " " << NgridZ << endl;
        f << "origin 0 0 0" << endl;
        f << "delta 1 0 0" << endl;
        f << "delta 0 1 0" << endl;
        f << "delta 0 0 1" << endl;
        f << "object 2 class gridconnections counts " << NgridX << " " << NgridY << " " << NgridZ << endl;
        f << "object 3 class array type double rank 0 items " << NgridX*NgridY*NgridZ << " data follows" << endl;

        int num = 0;
        for(i=0;i<NgridX;++i)
            for(j=0;j<NgridY;++j)
                for(k=0;k<NgridZ;++k){
                    f << grid[i][j][k] << " ";
                    ++num;
                    if(num==3){
                        f << endl;
                        num = 0;
                    }
                }
        f<< "\nobject \"WaterMap\" class field" << endl;
        f.close();
    }

    Options options;
    Grid searcher;
    boost::multi_array<float,3> grid;
    int NgridX, NgridY, NgridZ;
    Selection water;
    MatrixXf last_pos; // Previous positions of water
};

int main(int argc, char** argv){
    try {
        Options opt;
        parse_command_line(argc,argv,opt);

        Trajectory_reader proc(opt);
        proc.add_task( new Water_processor(opt) );

        proc.run();


    } catch(const Pteros_error& e){
        cout << e.what() << endl;
    }
}

