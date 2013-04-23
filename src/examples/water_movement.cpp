#include "pteros/analysis/trajectory_processor.h"
#include "pteros/analysis/consumer.h"
#include <fstream>
#include "pteros/core/grid_search.h"


using namespace std;
using namespace pteros;
using namespace Eigen;

class Water_processor: public Consumer {
public:
    Water_processor(Trajectory_processor* pr, Options_tree* opt): Consumer(pr){
        options = opt;
    }
protected:
    virtual void pre_process(){       
        // Create selection for water
        water.modify(system,"resname SOL and name OW");

        // Allocate a grid
        // Get size of the cell
        float cell_size = options->get_value<double>("cell_size");
        // Get extents of the box
        Vector3f box_dim = system.Box(0).colwise().norm();


        NgridX = floor(box_dim(0)/cell_size);
        NgridY = floor(box_dim(1)/cell_size);
        NgridZ = floor(box_dim(2)/cell_size);

        if(NgridX==0) NgridX = 1;
        if(NgridY==0) NgridY = 1;
        if(NgridZ==0) NgridZ = 1;

        //NgridX = NgridY = NgridZ = 3;

        searcher.create_custom_grid(NgridX, NgridY, NgridZ);

        // Create grid for result
        grid.resize(boost::extents[NgridX][NgridY][NgridZ]);
        // Clear it
        int i,j,k;
        for(i=0;i<NgridX;++i)
            for(j=0;j<NgridY;++j)
                for(k=0;k<NgridZ;++k)
                    grid[i][j][k] = 0.0;


    }


    virtual void process_frame(const Frame_info& info){
        int i,j,k,n,w,ind;

        if(info.valid_frame==0){
            water.get_xyz(last_pos);
            return;
        }

        // Assign waters to grid
        searcher.fill_custom_grid(water,false);
        // Cycle over grid cells
        for(i=0;i<NgridX;++i)
            for(j=0;j<NgridY;++j)
                for(k=0;k<NgridZ;++k){
                    // Get waters in this cell
                    n = searcher.cell_of_custom_grid(i,j,k).size();
                    for(w=0;w<n;++w){
                        ind = searcher.cell_of_custom_grid(i,j,k)[w];
                        // For current water get delta of prev and cur coorfinates
                        // And add it to result grid

                        grid[i][j][k] += pow( system.distance(water.XYZ(ind),last_pos.col(ind),0,true) ,2);
                        //grid[i][j][k] += 1;
                    }
                }

        water.get_xyz(last_pos);        
    }


    virtual void post_process(const Frame_info& info){
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

    Options_tree* options;
    Grid_searcher searcher;
    boost::multi_array<float,3> grid;
    int NgridX, NgridY, NgridZ;
    Selection water;
    MatrixXf last_pos; // Previous positions of water
};

int main(int argc, char** argv){
    try {
        Options_tree opt;
        opt.from_command_line(argc,argv);

        Trajectory_processor proc(opt);
        Water_processor wp(&proc,&opt);

        proc.run();


    } catch(Pteros_error e){
        e.print();
    }
}

