
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
//output.write(reinterpret_cast<char*>(&x),sizeof(x)); //NOLINT(..)
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

struct Grid{
    GridSize size;
    std::vector<std::vector<int>> blocks;
};

struct Acceleration{
    double ax;
    double ay;
    double az;
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
//Constants initialization
double r = 1.695;
double p = 1000;
double ps = 3.0;
double sc = 30000;
double dv = 128.0;
double nu = 0.4;
double dp = 0.0002; //Particle Size
double time_step = 0.001;
int particle_num = 2;
double g = 9.8;
//std::string address = "../new.fld";
vector<double> bmax = {0.08, 0.1, 0.08};
vector<double> bmin = {-0.08, -0.09, -0.08};
double h;
double m;

///ESTA FUNCION PARA QUE SIRVE?
/*
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
}*/

double distance_squared(Particle p1, Particle p2){
    double dx = p1.px - p2.px;
    double dy = p1.py - p2.py;
    double dz = p1.pz - p2.pz;

    return dx * dx + dy * dy + dz * dz;
}
///FUNCTIONS FOR SIMULATION



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
     //cout << "This is the x block " << block_x << ", y block " << block_y << ", z block " << block_z;
    int num_block = block_z + block_y*gridSize.nz + block_x*gridSize.nz*gridSize.ny;
    return num_block;
}

std::vector<int> get_contiguous_blocks(int current_block, GridSize gsize){
    std::vector<int> contiguous_blocks;
    int bx = current_block / (gsize.nz*gsize.ny);
    int by = (current_block - bx*(gsize.nz*gsize.ny)) / gsize.nz;
    int bz = current_block - (bx*gsize.nz*gsize.ny) - (by*gsize.nz);

    for (int x = -1; x < 2; x++){
        if ( (bx + x>=0) && (bx + x <= gsize.nx-1)) {
            for (int y = -1; y < 2; y++) {
                if ( (by + y>=0) && (by + y <= gsize.ny-1)) {
                    for (int z = -1; z < 2; z++) {
                        if ( (bz + z>=0) && (bz + z <= gsize.nz-1)){
                            int computed_block = (bz + z) + (by + y)*gsize.nz + (bx+x)*gsize.nz*gsize.ny;
                            contiguous_blocks.push_back(computed_block);
                        }
                    }
                }
            }
        }
    }

    return contiguous_blocks;
}

void AllBlocksArrayCreation(vector<Particle> &particles,GridSize gridSize){
    vector<vector<int>> all_blocks;
    for (int x = 0; x < (gridSize.nx-1)*(gridSize.ny-1)*(gridSize.nz-1); x++){
        vector <int> new_vector;
        all_blocks.push_back(new_vector);

    }
    /*
    int insertion_block =0;
    for (int i = 0; i < particles.size(); i++){
        insertion_block = find_block(particles[i],xxxxxx);
        all_blocks[insertion_block].push_back(particles.at(i));
    }
     */
}

void densities_increase(std::vector<Particle> &particles, Grid &grid, vector<double> &densities){ /// Cambiar p por part, porque ya hay una varibale gloabl p
    for (int i = 0; i < grid.blocks.size();i++){ ///Go through all blocks
        vector<int> contiguous_blocks = get_contiguous_blocks(i,grid.size); ///Get contiguous blocks to current block
        for (int p = 0; p < grid.blocks[i].size(); p++){ ///Go through each particle of the current block
            int particle_i_index = grid.blocks[i][p];
            Particle pi = particles[particle_i_index];
            for (int b = 0; b < contiguous_blocks.size(); b++){ ///Traverse the contiguous blocks
                int c_block_index = contiguous_blocks[b]; /// Get the index of the contiguous block to traverse
                for (int j = 0; j < grid.blocks[c_block_index].size();j++){ /// Go through each particle in the contiguous block
                    int particle_j_index = grid.blocks[c_block_index][j];
                    Particle pj = particles[particle_j_index];
                    if (particle_i_index != particle_j_index) { /// Check pi != pj
                        if (distance_squared(pi, pj) < (pow(h, 2))) {
                            densities[particle_i_index] += pow((pow(h, 2) - distance_squared(pi, pj)), 3);
                        }
                    }
                }
            }
        }
    }
}

void densities_transform(vector<double> &densities){
    for (int i = 0; i < densities.size(); i++){
        densities[i] = (densities[i] + pow(h,6))* (315*m)/(64*M_PI* pow(h,9));
    }
}

void acceleration_transfer(std::vector<Particle> &particles, Grid &grid, vector<double> &densities, vector<Acceleration> &accelerations){
    for (int i = 0; i < grid.blocks.size();i++){ ///Go through all blocks
        vector<int> contiguous_blocks = get_contiguous_blocks(i,grid.size); ///Get contiguous blocks to current block
        for (int part = 0; part < grid.blocks[i].size(); part++){ ///Go through each particle of the current block
            int particle_i_index = grid.blocks[i][part];
            Particle pi = particles[particle_i_index];
            for (int b = 0; b < contiguous_blocks.size(); b++){ ///Traverse the contiguous blocks
                int c_block_index = contiguous_blocks[b]; /// Get the index of the contiguous block to traverse
                for (int j = 0; j < grid.blocks[c_block_index].size();j++){ /// Go through each particle in the contiguous block
                    int particle_j_index = grid.blocks[c_block_index][j];
                    Particle pj = particles[particle_j_index];
                    if (particle_i_index != particle_j_index) { /// Check pi != pj
                        double dist_squared = distance_squared(pi, pj);
                        if (dist_squared < (pow(h, 2))) {
                            double distij = sqrt(max(dist_squared, pow(10,-12))); /// In these 4 lines calculate distij as stated in project and update accelerations
                            accelerations[particle_i_index].ax += ((pi.px - pj.px) * (15 / (M_PI*pow(h,6))) * (3 * m * ps/2) * pow(h-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*p) + (pj.vx - pi.vx) * (45/(M_PI*pow(h,6)) )* nu * m)/(densities[particle_i_index] * densities[particle_j_index]);
                            accelerations[particle_i_index].ay += ((pi.py - pj.py) * (15 / (M_PI*pow(h,6))) * (3 * m * ps/2) * pow(h-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*p) + (pj.vy - pi.vy) * (45/(M_PI*pow(h,6)) )* nu * m)/(densities[particle_i_index] * densities[particle_j_index]);
                            accelerations[particle_i_index].az += ((pi.pz - pj.pz) * (15 / (M_PI*pow(h,6))) * (3 * m * ps/2) * pow(h-distij,2)/distij * (densities[particle_i_index] + densities[particle_j_index] - 2*p) + (pj.vz - pi.vz) * (45/(M_PI*pow(h,6)) )* nu * m)/(densities[particle_i_index] * densities[particle_j_index]);
                        }
                    }
                }
            }
        }
    }
}

void update_particle_block(Particle &particle, Grid &grid, vector<int> offset){

}

///FUNCTIONS FOR SIMULATION


void reposition_particles(std::vector<Particle> &particles, Grid &grid){
    for (int i = 0; i < particles.size(); i++){
        int index = find_block(particles[i],grid.size);
        grid.blocks[index].push_back(i);
    }
}

void accelerations_computation(std::vector<Particle> &particles, Grid &grid,std::vector <double> &densities, std::vector <Acceleration> &accelerations){
    //initialization of densities and accelerations
    densities.clear();
    accelerations.clear();
    for (int i = 0; i < particles.size(); i++){
        densities.push_back(0);
        struct Acceleration a;
        a.ax = 0;
        a.ay = -g;
        a.az = 0;
        accelerations.push_back(a);
    }

    densities_increase(particles,grid,densities);
    densities_transform(densities);
    acceleration_transfer(particles,grid,densities,accelerations);

}
void particle_collision_with_X_axis(std::vector<Particle> &particles,  Grid &grid, std::vector <Acceleration> &accelerations) {

    double X_increment =0;
    //pared x=0
    for (int i = 0; i < grid.size.nz*grid.size.ny; i++){     //the number of particles in the x axix is ny * nz twice x min and xmax
            grid.blocks[i];
            for (int j = 0; j < grid.blocks[i].size();j++) {
                X_increment = dp-(particles[j].px-bmin[1]);  //dp − (x − xmin)
                if(X_increment<= pow(10,-10))
                    accelerations[j].ax = accelerations[j].ax + sc * X_increment - dv * particles[j].vx;       //ax + sc · ∆x − dv · vx

            }
    }
}

void particle_collision_with_Y_axis(std::vector<Particle> &particles , Grid &grid) {
    for (int i = 0; i < grid.size.nx*grid.size.nz; i++){
        int index = find_block(particles[i],);
        grid.blocks[index].push_back(i);
    }
}
void particle_collision_with_Z_axis(std::vector<Particle> &particles, Grid &grid) {
    for (int i = 0; i < grid.size.nx*grid.size.ny; i++){
        int index = find_block(particles[i],);
        grid.blocks[index].push_back(i);
    }
}
void particle_collision(std::vector<Particle> &particles, GridSize &gridSize){


}

void particles_motion(std::vector<Particle> &particles, std::vector <Acceleration> &accelerations){
    for (int i = 0; i < particles.size(); i++){

        double move_x = particles[i].hvx*time_step + accelerations[i].ax*pow(time_step,2);
        double move_y = particles[i].hvy*time_step + accelerations[i].ay*pow(time_step,2);
        double move_z = particles[i].hvz*time_step + accelerations[i].az*pow(time_step,2);

        particles[i].px += move_x;
        particles[i].py += move_y;
        particles[i].pz += move_z;
        particles[i].vx = particles[i].hvx + (accelerations[i].ax*time_step)/2;
        particles[i].vy = particles[i].hvy + (accelerations[i].ay*time_step)/2;
        particles[i].vz = particles[i].hvz + (accelerations[i].az*time_step)/2;
        particles[i].hvx = particles[i].hvx + accelerations[i].ax*time_step;
        particles[i].hvy = particles[i].hvy + accelerations[i].ay*time_step;
        particles[i].hvz = particles[i].hvz + accelerations[i].az*time_step;
    }
}
void simulate(int nsteps, std::vector<Particle> &particles, Grid &grid){
    vector <double> densities;
    vector <Acceleration> accelerations;
    //Stages of Simulation
    //Stage 1: Reposition Particles
    reposition_particles(particles,grid);
    //Stage 2: Accelerations computation
    accelerations_computation(particles,grid,densities,accelerations);
    //stage 4: Particles motion
    particles_motion(particles,accelerations);
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
int read_particle_info(ifstream &file,vector<Particle> &particles){
    int counter = 0;
    while(!file.eof()){
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
        counter ++;
    }
    return counter;}

vector<Particle> initial_read(std::string file_address,Initial_values &initialValues){
    //Read file
    ifstream file(file_address, ios::binary);
    if (!file.is_open()) { //Check error opening
        cout<<"Error: Cannot open " << file_address <<" for reading";
        exit (-3);
    }
    std::vector<Particle> particles; //Iniatilize a vector of particles to be stored
    initialValues = read_general_info(file);//call to a function to read parameters
    //
    int counter = read_particle_info(file,particles);
    if(counter == 0){
        cout<< "Error : Invalid number of particles: " << counter <<".";
    }
    else if(counter != initialValues.np){
        cout<<"Error : Number of particles mismatch. Header " << initialValues.np << " Found " << counter <<".";
    }
    file.close();
    return particles;
}

void check_command_errors(int argc,char** argv) {
    if (argc != 3) {
        cout << "Error: Invalid number of arguments: " << argc << ".";
        exit(-1);
    }
    int i=0;
    while (argv[1] != "\0") {
        if (!isdigit((argv[1])[i])) {
            exit(-1);
        }
        i++;
    }
    if((argv[1])[0] == '-'){
        exit (-2);
    }

}

void write_to_file(std::string output_file_address,std::vector<Particle> particles,Initial_values initialValues){


}

int main(int argc, char** argv) {
    check_command_errors(argc,argv);
    Initial_values initialValues;
    std::vector<Particle> particles = initial_read(argv[2],initialValues);
    h = initialValues.h;
    m = initialValues.m;
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
    std::vector<vector <int>> all_blocks;

    //cout << "\nnx " << nx << " ny " << ny << " nz " << nz << " Number of blocks " << NumberofBlocks;

    /*
    GridSize g;
    g.nx = 3;
    g.ny = 5;
    g.nz = 5;

    std::vector<int> conti;
    conti = get_contiguous_blocks(31,g);
    for(int i = 0; i<conti.size();i++)
        cout << conti[i]<< " ";

    cout<<"\n";
    */


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
    /*
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
*/
//Collisions


//Write floats use standard conversion static cast or narrow cast

///Looooooooooooooooooooooooooooool

}
