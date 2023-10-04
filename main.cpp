#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>


//1-11 17- end
using namespace std;
//External libraries
//GSL
//Google Test

//Twenty lines per function maximum
//Pointer arithmetic forbidden
//output.write(reninterpret_cast<char*>(&x),sizeof(x)); //NOLINT(..)
//Vectores iniciales - fuerzas resultantes -
//Unique pointer or shared point

//Vector of blocks where each block is a list or vector of particles



struct Particle{
    double px;
    double py;
    double pz;
    double hvx;
    double hvy;
    double hvz;
    double vx;
    double vy;
    double vz;
    double a = 9.81;
    double p = 0;
};


struct Block{
    int x;
    int y;
    int z;
    vector<Particle> particles;
};



double distance_squared(Particle p1, Particle p2){
    double dx = p1.px - p2.px;
    double dy = p1.py - p2.py;
    double dz = p1.pz - p2.pz;

    return dx * dx + dy * dy + dz * dz;
};

/*Introducir en cada bloque las partículas que le corresponden(en proceso)
 Esta función va a ser para, cuando queramos acceder a un bloque al principio, para meter las partículas
 que le corresponden a cada bloque según su posición, y cuando querramos acceder a ese bloque desde el grid,
     cuyas posiciones están organizadas de esta manera*/

int find_block(Particle particle,vector<int> blockNumber,vector<double> block_size,vector<double> min_coord){
    int block_x = floor((particle.px - min_coord[0])/block_size[0]);
    int block_y = floor((particle.py - min_coord[1])/block_size[1]);
    int block_z = floor((particle.pz - min_coord[2])/block_size[2]);
    if (block_x < 0){
        block_x = 0;
    } else if (block_x > blockNumber[0]) {
        block_x = blockNumber[0]-1;
    }if (block_y < 0){
        block_y = 0;
    } else if (block_y > blockNumber[1]) {
        block_y = blockNumber[1]-1;
    }if (block_z < 0){
        block_z = 0;
    } else if (block_z > blockNumber[2]) {
        block_z = blockNumber[2]-1;
    }


    cout << "This is the x block " << block_x << ", y block " << block_y << ", z block " << block_z;
    int num_block = block_z + block_y*blockNumber[2] + block_x*blockNumber[2]*blockNumber[1];
    return num_block;
}

int main(int argc, char** argv) {

    float ppm_float;
    int np;

    float px;
    float py;
    float pz;
    float hvx;
    float hvy;
    float hvz;
    float vx;
    float vy;
    float vz;
    //Contsants intialization
    double r = 1.695;
    double p = 1000;
    double ps = 3.0;
    double sc = 30000;
    double dv = 128.0;
    double nu = 0.4;
    double dp = 0.0002;
    double time_step = 0.001;
    int particle_num = 2;




    //Constants for grid creation
    //Remove normalization
    vector<double> bmax = {0.08, 0.1, 0.08};
    vector<double> bmin = {-0.08, -0.09, -0.08};



    //Parte de lectura de Juan crear todas las particulas
    ifstream file("/Users/raulpineda/iCloud Drive (Archive)/Document/Arquitectura/Lab1/small.fld", ios::binary);

    if (!file.is_open()) {
        cerr << "Error opening file." << std::endl;
        return 1;
    }

    // Read ppm and np

    //cap 10-11
    file.read(reinterpret_cast<char*>(&ppm_float), sizeof(float));
    file.read(reinterpret_cast<char*>(&np), sizeof(int));
    double ppm = static_cast<float>(ppm_float);
    double m = p*pow(ppm,3);
    double h = r/ppm;

    std::cout << "ppm: " << ppm_float << ", np: " << np << std::endl;

    // Create a vector to store particles
    std::vector<Particle> particles;

    for (int i = 0; i < np; ++i) {
        Particle particle;
        file.read(reinterpret_cast<char*>(&px), sizeof(float));
        file.read(reinterpret_cast<char*>(&py), sizeof(float));
        file.read(reinterpret_cast<char*>(&pz), sizeof(float));
        file.read(reinterpret_cast<char*>(&hvx), sizeof(float));
        file.read(reinterpret_cast<char*>(&hvy), sizeof(float));
        file.read(reinterpret_cast<char*>(&hvz), sizeof(float));
        file.read(reinterpret_cast<char*>(&vx), sizeof(float));
        file.read(reinterpret_cast<char*>(&vy), sizeof(float));
        file.read(reinterpret_cast<char*>(&vz), sizeof(float));


        particle.px = static_cast<double>(trunc(px));
        particle.py = static_cast<double>(py);
        particle.pz = static_cast<double>(pz);
        particle.hvx = static_cast<double>(hvx);
        particle.hvy = static_cast<double>(hvy);
        particle.hvz = static_cast<double>(hvz);
        particle.vx = static_cast<double>(vx);
        particle.vy = static_cast<double>(vy);
        particle.vz = static_cast<double>(vz);
        particles.push_back(particle);
    }

    file.close();
    int counter;
    //Access particles
    /*for (const Particle& particle : particles) {
        std::cout << "Particle Data:" << counter  << std::endl;
        std::cout << "px: " << particle.px << ", py: " << particle.py << ", pz: " << particle.pz << std::endl;
        std::cout << "hvx: " << particle.hvx << ", hvy: " << particle.hvy << ", hvz: " << particle.hvz << std::endl;
        std::cout << "vx: " << particle.vx << ", vy: " << particle.vy << ", vz: " << particle.vz << std::endl;
        counter +=1;
    }*/

    double boxx = bmax[0] - bmin[0];
    double boxy = bmax[1] - bmin[1];
    double boxz = bmax[2] - bmin[2];

    m = p/pow(ppm,3);
    h = r/ppm;
    cout << boxx;


    int nx = floor(boxx/h);
    int ny = floor(boxy/h);
    int nz = floor(boxz/h);
    cout << "\n r: " << r << " ppm: " << ppm;

    double sx = boxx/nx;
    double sy = boxy/ny;
    double sz = boxz/nz;
    cout << "\n sx: " << sx << " sy " << sy << " sz " << sz;
    int NumberofBlocks = (nx-1)*(ny-1)*(nz-1);




    //Placing particles in blocks
    //Create a grid which is made of blocks
    std::vector<Block> grid;
    grid.reserve(NumberofBlocks);
    int counterBlock= 0;

    cout << "\nnx " << nx << " ny " << ny << " nz " << nz << " Number of blocks " << NumberofBlocks;

    for (int x = 0; x < nx-1; x++){
        for (int y = 0; y < ny-1; y++){
            for (int z = 0; z < nz-1; z++){
                counterBlock+=1;
                Block block;
                block.x = x;
                block.y = y;
                block.z = z;
                grid.push_back(block);

            }
        }
    }


    //Grid organization
    cout << "\n Grid organization";
    int blockPositionx;
    int blockPositiony;
    int blockPositionz;

    vector<int> blockAmmount = {nx,ny,nz};
    vector<double> blockSize = {sx,sy,sz};
    int blockNumber;
    counter = 0;
    for (auto particle = particles.begin(); particle != particles.end(); particle++){
        blockNumber =  find_block(*particle,blockAmmount,blockSize,bmin);
        cout << "\nThis is particle number: " << counter <<  "Block number" << blockNumber;
        counter +=1;
        //For every y there are nz number of z blocks and for every x there are ny * nz number of blocks

        grid[blockNumber].particles.push_back(*particle);
    }




    //Density computation
    double d;
    double i_density;
    double j_density;
    double density_inc;
    double acceleration_inc;
    double density_constant_transformation_part = 64*M_PI*pow(h,9);
    double densityConstantTransformation =  315 * m / density_constant_transformation_part;

    for (auto i = particles.begin();i !=particles.end(); i++){
        for (auto j = i +1;j != particles.end(); j++){
            d = distance_squared(*i,*j);

            if (d < h){
                //Density
                density_inc =pow((h*h-d),3);
                i_density = i->p + density_inc;
                j_density = i->p +density_inc;
                //Linear transformation
                i->p =(i_density + pow(h,6))*densityConstantTransformation;
                //Check in optimization
                j->p =(j_density + pow(h,6))*densityConstantTransformation;
                //Acceleration
                //Local to each iteration
                //Check formula later
                i->hvx = i->hvx + ((i->px-j->px)*15*m*(h-d)*(h-d)*(i->p+j->p -p)+45*(i->vx-j->vx)*nu*m)/(M_PI*pow(h,6)*i->p*j->p);
                i->hvy = i->hvy + ((i->py-j->py)*15*m*(h-d)*(h-d)*(i->p+j->p -p)+45*(i->vy-j->vy)*nu*m)/(M_PI*pow(h,6)*i->p*j->p);
                i->hvz = i->hvz + ((i->pz-j->pz)*15*m*(h-d)*(h-d)*(i->p+j->p -p)+45*(i->vz-j->vz)*nu*m)/(M_PI*pow(h,6)*i->p*j->p);
            }

        }
    }

//Collisions


//Write floats use standard conversion static cast or narrow cast


}
