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
//Structure for initially read variables ppm,np,m,h
struct Initial_values{
    float ppm;
    int np;
    double m;
    double h;
};

struct GridSize{ /// Struct to store grid size, used later to ease arguments passing in fucntions
    int nx;
    int ny;
    int nz;
    double sx;
    double sy;
    double sz;
};

//Global Variables //porque está los componentes de una partícula como variables globales??!
float px;
float py;
float pz;
float hvx;
float hvy;
float hvz;
float vx;
float vy;
float vz;
//Constants intialization
double r = 1.695;
double p = 1000;
double ps = 3.0;
double sc = 30000;
double dv = 128.0;
double nu = 0.4;
double dp = 0.0002;
double time_step = 0.001;
int particle_num = 2;
std::string address = "./small.fld";
vector<double> bmax = {0.08, 0.1, 0.08};
vector<double> bmin = {-0.08, -0.09, -0.08};


int particles_statistics(vector<vector <Particle>> &grid){
    int maximum = 0;
    int minimum = 0;
    int particles_num;
    int average= 0 ;
    for (auto & i : grid){
        particles_num = i.size();
        if (particles_num > maximum){
            maximum = particles_num;
        }else if(particles_num < minimum){
            minimum = particles_num;
        }
        average += particles_num;
    }
    int block_number = grid.size();
    average = average/block_number;
    cout << "\nThe average number of particles is: " << average << "\nThe maximum is: " << maximum << "\n The minimum is: " << minimum;
}

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

int find_block(Particle particle,GridSize gridSize){
    int block_x = floor((particle.px - bmin[0])/gridSize.sx);
    int block_y = floor((particle.py - bmin[1])/gridSize.sy);
    int block_z = floor((particle.pz - bmin[2])/gridSize.sz);
    if (block_x < 0){
        block_x = 0;
    } else if (block_x >= gridSize.nx) {
        block_x = gridSize.nx-1;
    }if (block_y < 0){
        block_y = 0;
    } else if (block_y >= gridSize.ny) {
        block_y = gridSize.ny-1;
    }if (block_z < 0){
        block_z = 0;
    } else if (block_z >= gridSize.nz) {
        block_z = gridSize.nz-1;
    }
    cout << "This is the x block " << block_x << ", y block " << block_y << ", z block " << block_z;
    int num_block = block_z + block_y*gridSize.nz + block_x*gridSize.nz*gridSize.ny;
    return num_block;
}


///FUNCTIONS FOR SIMULATION


void reposition_particles(std::vector<Particle> &particles, vector<vector<Particle>> &grid){
    for (int i = 0; i < particles.size(); i++){
        
    }
}


void simulate(int nsteps, std::vector<Particle> &particles, vector<vector<Particle>> &grid){

}


///------------------------


Initial_values read_general_info(ifstream &file){
    // Read ppm and np
    //cap 10-11
    Initial_values initialValues;
    file.read(reinterpret_cast<char*>(&initialValues.ppm), sizeof(float));
    file.read(reinterpret_cast<char*>(&initialValues.np), sizeof(int));
    initialValues.ppm = static_cast<float>(initialValues.ppm);
    initialValues.m = p*pow(initialValues.ppm,3);
    initialValues.h = r/initialValues.ppm;
    std::cout << "ppm: " << initialValues.ppm << ", np: " << initialValues.np << std::endl;
    return initialValues;
}
int read_particle_info(ifstream &file,vector<Particle> &particles,Initial_values initialValues){
    int counter;
    for (int i = 0; i < initialValues.np; ++i) {
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
        counter = i;
    }
    return counter;
}

vector<Particle> initial_read(std::string file_address,Initial_values &initialValues){
    //Parte de lectura de Juan crear todas las particulas
    ifstream file(file_address, ios::binary);
    if (!file.is_open()) {
        cerr << "Error opening file." << std::endl;
        exit(1);
    }
    std::vector<Particle> particles;
    initialValues = read_general_info(file);
    int counter = read_particle_info(file,particles,initialValues);
    if(counter == 0){
        cout<< "Error : Invalid number of particles: " << counter <<".";
    }
    else if(counter != initialValues.np){
        cout<<"Error : Number of particles mismatch. Header " << initialValues.np << " Found " << counter <<".";
    }
    file.close();
    return particles;
}

int main(int argc, char** argv) {
    Initial_values initialValues;
    std::vector<Particle> particles = initial_read(address,initialValues);
    //Constants for grid creation
    //Remove normalization

    // Create a vector to store particles

    GridSize gridSize;

    double boxx = bmax[0] - bmin[0];
    double boxy = bmax[1] - bmin[1];
    double boxz = bmax[2] - bmin[2];

    //cout << boxx;

    gridSize.nx = floor(boxx/initialValues.h);
    gridSize.ny = floor(boxy/initialValues.h);
    gridSize.nz = floor(boxz/initialValues.h);

    //cout << "\n r: " << r << " ppm: " << ppm;

    gridSize.sx = boxx/gridSize.nx;
    gridSize.sy = boxy/gridSize.ny;
    gridSize.sz = boxz/gridSize.nz;


    //cout << "\n sx: " << sx << " sy " << sy << " sz " << sz;

    //Placing particles in blocks
    //Create a grid which is made of blocks
    std::vector<vector <int>> grid;

    //cout << "\nnx " << nx << " ny " << ny << " nz " << nz << " Number of blocks " << NumberofBlocks;

    for (int x = 0; x < (gridSize.nx-1)*(gridSize.ny-1)*(gridSize.nz-1); x++){
        vector <int> new_vector;
        grid.push_back(new_vector);
    }


    ///Comento esta sección porque esto ya se va a hacer por la funcion reposition_particles
    /*
    //Grid organization
    cout << "\n Grid organization";

    vector<int> blockAmmount = {nx,ny,nz};
    vector<double> blockSize = {sx,sy,sz};
    int blockNumber;
    int count = 0;
    for (auto particle = particles.begin(); particle != particles.end(); particle++){
        blockNumber =  find_block(*particle,blockAmmount,blockSize,bmin);
        cout << "\nThis is particle number: " << count <<  "Block number" << blockNumber;
        count +=1;
        //For every y there are nz number of z blocks and for every x there are ny * nz number of blocks

        grid[blockNumber].push_back(*particle); ///ESTO SE PUDE USAR???
    }
    */


    ///Raul explicanos esta parte cuando lo leas, te queremos
    //Density computation
    
    double m = initialValues.m;
    double h = initialValues.h;
    
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
