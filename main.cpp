
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

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
double sc = 30000;   //stiffness collision
double dv = 128.0;
double nu = 0.4;
double dp = 0.0002; //Particle Size
double time_step = 0.001;
int particle_num = 2;
double g = 9.8;
//std::string address = "../new.fld";
vector<double> bmax = {0.065, 0.1, 0.065};
vector<double> bmin = {-0.065, -0.08, -0.065};
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



///Functions for trace debug ----- $$$$$$$$$$$$$$$$$$$$$$$$$$

void compare_particle(Particle &p1, Particle &p2,long id){
    if (p1.px != p2.px){
        cout<<"id = "<<id<<" "<<"Particles x pos differ, p1.px = "<<p1.px<<" p2.px = "<<p2.px<<'\n';
        exit(-1);
    }
    if (p1.py != p2.py){
        cout<<"id = "<<id<<" "<<"Particles y pos differ, p1.py = "<<p1.py<<" p2.py = "<<p2.py<<'\n';
        exit(-1);
    }
    if (p1.pz != p2.pz){
        cout<<"id = "<<id<<" "<<"Particles z pos differ, p1.pz = "<<p1.pz<<" p2.pz = "<<p2.pz<<'\n';
        exit(-1);
    }
    if (p1.hvx != p2.hvx){
        cout<<"id = "<<id<<" "<<"Particles hvx pos differ, p1.hvx = "<<p1.hvx<<" p2.hvx = "<<p2.hvx<<'\n';
        exit(-1);
    }
    if (p1.hvy != p2.hvy){
        cout<<"id = "<<id<<" "<<"Particles hvy pos differ, p1.hvy = "<<p1.hvy<<" p2.hvy = "<<p2.hvy<<'\n';
        exit(-1);
    }
    if (p1.hvz != p2.hvz){
        cout<<"id = "<<id<<" "<<"Particles hvz pos differ, p1.hvz = "<<p1.hvz<<" p2.hvz = "<<p2.hvz<<'\n';
        exit(-1);
    }
    if (p1.vx != p2.vx){
        cout<<"id = "<<id<<" "<<"Particles vx pos differ, p1.vx = "<<p1.vx<<" p2.vx = "<<p2.vx<<'\n';
        exit(-1);
    }
    if (p1.vy != p2.vy){
        cout<<"id = "<<id<<" "<<"Particles vy pos differ, p1.vy = "<<p1.vy<<" p2.vy = "<<p2.vy<<'\n';
        exit(-1);
    }
    if (p1.vz != p2.vz){
        cout<<"id = "<<id<<" "<<"Particles vz pos differ, p1.vz = "<<p1.vz<<" p2.vz = "<<p2.vz<<'\n';
        exit(-1);
    }

}

void compare_accelerations(Acceleration &a1, Acceleration &a2, long id){
    if (a1.ax!=a2.ax){
        cout<<"id = "<<id<<" "<<"Accelerations ax differ, a1.ax = "<<a1.ax<<" a2.ax = "<<a2.ax<<'\n';
        exit(-1);
    }
    if (a1.ay!=a2.ay){
        cout<<"id = "<<id<<" "<<"Accelerations ay differ, a1.ay = "<<a1.ay<<" a2.ay = "<<a2.ay<<'\n';
        exit(-1);
    }
    if (a1.az!=a2.az){
        cout<<"id = "<<id<<" "<<"Accelerations az differ, a1.az = "<<a1.az<<" a2.az = "<<a2.az<<'\n';
        exit(-1);
    }
}

void find_elem(int e, std::vector<int> &v){
    ///Only works if both vectors have same size
    int found = 0;
    for (int i = 0; i<v.size();i++){
        if(e == v[i]){
            found = 1;
        }
    }
    if (found==0){
        cout<<"Id "<<e<<" not found in grid block\n";
        exit(-1);
    }
}

void check_trace(string trz, Grid &grid, vector<Particle> &particles, vector<double> &densities, vector<Acceleration> &accelerations){
    ifstream file(trz, ios::binary);
    if (!file.is_open()) { //Check error opening
        cout<<"Error: Cannot open trace file: " << trz <<" for reading";
        exit (-1);
    }

    int num_blocks;
    file.read(reinterpret_cast<char*>(&num_blocks), sizeof(int));

    if (num_blocks != grid.blocks.size()){
        cout<<"Number of blocks differ from trace:\n"<<"trz_blocks = "<<num_blocks<<'\t'<<"grid_blocks = "<<grid.blocks.size()<<'\n';
        //cout<<"Pruebita"<<(grid.size.nx-1)*(grid.size.ny-1)*(grid.size.nz-1)<<'\n';
        exit(-1);
    }

    long particles_in_block;
    long id;
    Particle part;
    double d;
    Acceleration a;

    for (int i = 0; i < num_blocks;i++){
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));
        cout<<"Entering block: "<<i<<" particles in block = "<<particles_in_block<<'\n';

        if(grid.blocks[i].size()!=particles_in_block){
            cout<<"Number of particles for block "<<i<<" mismatch: "<<"grid.blocks["<<i<<"].size() = "<<grid.blocks[i].size()<<" particles in block = "<<particles_in_block<<'\n';
            exit(-1);
        }

        for (int p = 0; p<particles_in_block;p++){
            file.read(reinterpret_cast<char*>(&id), sizeof(long));
            find_elem(id,grid.blocks[i]);
            file.read(reinterpret_cast<char*>(&part.px), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.py), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.pz), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.hvx), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.hvy), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.hvz), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.vx), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.vy), sizeof(double));
            file.read(reinterpret_cast<char*>(&part.vz), sizeof(double));
            file.read(reinterpret_cast<char*>(&d), sizeof(double));
            file.read(reinterpret_cast<char*>(&a.ax), sizeof(double));
            file.read(reinterpret_cast<char*>(&a.ay), sizeof(double));
            file.read(reinterpret_cast<char*>(&a.az), sizeof(double));


            compare_particle(particles[id],part,id);

            if (densities[id]!=d){
                cout<<"Densities for particle "<<id<<" differ, d = "<<d<<" densities["<<id<<"] = "<<densities[id]<<'\n';
                exit(-1);
            }

            compare_accelerations(accelerations[id],a,id);

        }
    }


    cout<<"\nTrace is equal to current state of the simulation\n";

    file.close();
}

void load_trace(string trz, Grid &grid, vector<Particle> &particles, vector<double> &densities, vector<Acceleration> &accelerations, Initial_values &i_v){
    ifstream file(trz, ios::binary);
    if (!file.is_open()) { //Check error opening
        cout<<"Error: Cannot open trace file: " << trz <<" for reading";
        exit (-1);
    }

    int num_blocks;
    file.read(reinterpret_cast<char*>(&num_blocks), sizeof(int));

    std::vector<std::vector<int>> blocks(num_blocks);
    grid.blocks = blocks;

    particles = std::vector<Particle> (i_v.np);
    densities = std::vector<double> (i_v.np);
    accelerations = std::vector<Acceleration> (i_v.np);

    long particles_in_block;
    long id;
    Particle part;
    double d;
    Acceleration a;

    cout<<"Grid.size = "<<grid.blocks.size()<<" total particles = "<<particles.size()<<'\n';

    for (int i = 0; i<num_blocks;i++){
        file.read(reinterpret_cast<char*>(&particles_in_block), sizeof(long));
        //cout<<"Loading block: "<<i<<" particles in block = "<<particles_in_block<<'\n';

        for (int p = 0; p<particles_in_block;p++){
            file.read(reinterpret_cast<char*>(&id), sizeof(long));
            grid.blocks[i].push_back(id);
            file.read(reinterpret_cast<char*>(&particles[id].px), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].py), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].pz), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].hvx), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].hvy), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].hvz), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].vx), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].vy), sizeof(double));
            file.read(reinterpret_cast<char*>(&particles[id].vz), sizeof(double));
            file.read(reinterpret_cast<char*>(&densities[id]), sizeof(double));
            file.read(reinterpret_cast<char*>(&accelerations[id].ax), sizeof(double));
            file.read(reinterpret_cast<char*>(&accelerations[id].ay), sizeof(double));
            file.read(reinterpret_cast<char*>(&accelerations[id].az), sizeof(double));
        }
    }

    cout<<"\nTrace loaded\n";
}


///------------------------
///-$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




double distance_squared(Particle p1, Particle p2){
    double dx = p1.px - p2.px;
    double dy = p1.py - p2.py;
    double dz = p1.pz - p2.pz;

    return dx * dx + dy * dy + dz * dz;
}

//FILE FUNCTIONS
Initial_values read_general_info(ifstream &file){
    // Read ppm and np
    //cap 10-11
    Initial_values initialValues;
    file.read(reinterpret_cast<char*>(&initialValues.ppm), sizeof(float));
    file.read(reinterpret_cast<char*>(&initialValues.np), sizeof(int));
    initialValues.ppm = static_cast<float>(initialValues.ppm); ///Esto no debería ser double en vez de float??????
    initialValues.m = p*pow(initialValues.ppm,3);
    initialValues.h = r/initialValues.ppm;
    std::cout << "ppm: " << initialValues.ppm << ", np: " << initialValues.np << std::endl;
    return initialValues;
}

int read_particle_info(ifstream &file,vector<Particle> &particles){
    int counter = 0;
    while(file.peek()!=EOF){
        counter ++;
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
        particle.px = static_cast<double>(px);
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
    cout<<"Lee bien ppm np";
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



void write_to_file(std::string output_file_address,std::vector<Particle> particles, Initial_values &initialValues){
    //Write to the file all the new values
    ofstream output_file(output_file_address, ios::binary);
    if (!output_file.is_open()) { //Check error opening
        cout<<"Error: Cannot open " << output_file_address <<" for writing";
        exit (-3);
    }

    output_file.write(reinterpret_cast<char*>(&initialValues.ppm), sizeof(float));
    output_file.write(reinterpret_cast<char*>(&initialValues.np), sizeof(int));

    for (int i = 0; i < particles.size(); i++){
        px = static_cast<float>(particles[i].px);
        py = static_cast<float>(particles[i].py);
        pz = static_cast<float>(particles[i].pz);
        hvx = static_cast<float>(particles[i].hvx);
        hvy = static_cast<float>(particles[i].hvy);
        hvz = static_cast<float>(particles[i].hvz);
        vx = static_cast<float>(particles[i].vx);
        vy = static_cast<float>(particles[i].vy);
        vz = static_cast<float>(particles[i].vz);
        output_file.write(reinterpret_cast<const char*>(&px), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&py), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&pz), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&hvx), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&hvy), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&hvz), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&vx), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&vy), sizeof(float));
        output_file.write(reinterpret_cast<const char*>(&vz), sizeof(float));
        //output_file.write("\n", sizeof(char ));
    }
    output_file.close();
}
vector <vector<int>> gridCreation(vector<Particle> &particles,GridSize gridSize){
    vector<vector<int>> blocks;
    for (int x = 0; x < (gridSize.nx)*(gridSize.ny)*(gridSize.nz); x++){
        vector <int> new_vector;
        blocks.push_back(new_vector);
    }
    cout<<"Before : "<< blocks.size()<<'\n';
    return blocks;
}



void check_command_errors(int argc,char** argv) {
    if (argc != 4) {
        cout << "Error: Invalid number of arguments: " << argc << ".";
        exit(-1);
    }
    int i=0;
    while (argv[1][0] != '\0') {
        if (!isdigit((argv[1])[i])) {
            exit(-1);
        }
        i++;
    }
    if((argv[1])[0] == '-'){
        exit (-2);
    }

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

/*
void reposition_particles(std::vector<Particle> &particles, Grid &grid){
    for (int i = 0; i < particles.size(); i++){
        int index = find_block(particles[i],grid.size);
        grid.blocks[index].push_back(i);
    }
}
*/


void densities_increase(std::vector<Particle> &particles, Grid &grid, vector<double> &densities){ /// Cambiar p por part, porque ya hay una varibale gloabl p
    vector<int> contiguous_blocks;int particle_i_index;
    Particle pi;Particle pj;
    int c_block_index;
    int particle_j_index;
    for (int i = 0; i < grid.blocks.size();i++){ ///Go through all blocks
        contiguous_blocks = get_contiguous_blocks(i,grid.size); ///Get contiguous blocks to current block
        for (int p = 0; p < grid.blocks[i].size(); p++){ ///Go through each particle of the current block
            particle_i_index = grid.blocks[i][p];
            pi = particles[particle_i_index];
            for (int b = 0; b < contiguous_blocks.size(); b++){ ///Traverse the contiguous blocks
                c_block_index = contiguous_blocks[b]; /// Get the index of the contiguous block to traverse
                for (int j = 0; j < grid.blocks[c_block_index].size();j++){ /// Go through each particle in the contiguous block
                    particle_j_index = grid.blocks[c_block_index][j];
                    pj = particles[particle_j_index];
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





void accelerations_computation(std::vector<Particle> &particles, Grid &grid,std::vector <double> &densities, std::vector <Acceleration> &accelerations){
    //initialization of densities and accelerations
    densities.clear();
    accelerations.clear();
    struct Acceleration a;
    for (int i = 0; i < particles.size(); i++){
        densities.push_back(0);
        a.ax = 0;
        a.ay = -g;
        a.az = 0;
        accelerations.push_back(a);
    }

    //check_trace("../trz/small/densinc-base-1.trz",grid,particles,densities,accelerations);
    densities_increase(particles,grid,densities);
    check_trace("../trz/small/densinc-base-1.trz",grid,particles,densities,accelerations);
    densities_transform(densities);
    acceleration_transfer(particles,grid,densities,accelerations);

}

void particle_collision_with_X_axis(std::vector<Particle> &particles,  Grid &grid, std::vector <Acceleration> &accelerations) {

    double x_param;
    double increment;
    for (int i = 0; i < grid.size.nz*grid.size.ny; i++){   //pared x_0 //the number of particles in the x axis is ny * nz twice x min and xmax
        for (int j : grid.blocks[i]) {
            x_param = particles[j].px + particles[j].hvx * time_step;        //  x = px + hvx · ∆t
            increment = dp-(x_param-bmin[0]);                                //  dp − (x − xmin)
            if(increment> pow(10,-10))                                 //  ax + (cs · ∆x − dv · vx)
                accelerations[j].ax = accelerations[j].ax + (sc * increment - dv * particles[j].vx);
        }
    }
    for (int i = (grid.size.nz*grid.size.ny*grid.size.nx - grid.size.nz*grid.size.ny)-1; i < grid.size.nz*grid.size.ny*grid.size.nx; i++){    //pared x_max
        for (int j : grid.blocks[i]) {
            x_param = particles[j].px + particles[j].hvx * time_step;        //  x = px + hvx · ∆t
            increment = dp-(bmax[0]-x_param);                                //  dp − (xmax − x)
            if (increment > pow(10, -10))                              //  ax − (cs · ∆x + dv · vx)
                accelerations[j].ax = accelerations[j].ax - (sc * increment + dv * particles[j].vx);
        }
    }
}

void particle_collision_with_Y_axis(std::vector<Particle> &particles , Grid &grid, std::vector <Acceleration> &accelerations) {

    double y_param;
    double increment;
    for (int i = 0; i < grid.size.nx; i++){     //the number of particles in the y axis is nx * nz, twice y min and xmax
        for (int j = 0, k=grid.size.nz*(grid.size.ny-1); j < grid.size.nz; j++,k++) {//pared Y_0
            for (int l : grid.blocks[ j + i * grid.size.nz * grid.size.ny]) {//  block_index = j + i * grid.size.nz * grid.size.ny;
                y_param = particles[l].py  + particles[l].hvy * time_step;   //  y = py + hvy · ∆t
                increment = dp - (y_param - bmin[1]);                        //  dp − (y − ymin)grid.blocks[block_index]
                if (increment > pow(10, -10))                          //  ay + (sc · ∆y − dv · vy)
                    accelerations[l].ay = accelerations[l].ay + (sc * increment - dv * particles[l].vy);
            }                                                                //  pared Y_max
            for (int l : grid.blocks[k + i * grid.size.nz * grid.size.ny]) { //  block_index = k + i * grid.size.nz * grid.size.ny;
                y_param = particles[l].py + particles[l].hvy * time_step;    //  y = py + hvy · ∆t
                increment = dp - (bmax[1] - y_param);                        //  dp − (ymay − y) i
                if (increment > pow(10, -10))                          //  ay − (sc · ∆y + dv · vy)
                    accelerations[l].ay = accelerations[l].ay - (sc * increment + dv * particles[l].vy);
            }
        }
    }
}

void particle_collision_with_Z_axis(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &accelerations) {

    double z_param ;
    double increment;
    for (int i = 0,j =grid.size.nz-1; i < grid.size.nz*grid.size.ny*grid.size.nx; i+=grid.size.nz, j+=grid.size.nz){
        for (int l : grid.blocks[i]) {                                       //pared Z_0
            z_param = particles[l].pz + particles[l].hvz * time_step;        //  z = pz + hvz · ∆t
            increment =dp - (z_param - bmin[2]);                             //  dp − (z − zmin)
            if (increment > pow(10, -10))                              //az + (sc · ∆z − dv · vz)
                accelerations[l].az =accelerations[l].az + (sc * increment - dv * particles[l].vz);
        }
        for (int l : grid.blocks[j]) {                                       //pared Z_max
            z_param = particles[l].pz + particles[l].hvz * time_step;        //z = pz + hvz · ∆t
            increment =dp - (bmax[2] - z_param);                             //dp − (zmax− z) i
            if (increment > pow(10, -10))                              //az − (sc · ∆z + dv · vz)
                accelerations[l].az =accelerations[l].az - (sc * increment + dv * particles[l].vz);
        }
    }
}


void X_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_x;
    for (int i = 0; i < grid.size.nz*grid.size.ny; i++){//pared x_0
        for (int j : grid.blocks[i]) {
            distance_x = particles[j].px-bmin[0];
            if (distance_x < 0){
                particles[j].px = bmin[0] - distance_x;
                particles[j].vx = -particles[j].vx;
                particles[j].hvx = -particles[j].hvx;
            }
        }
    }                                                   //pared x_max
    for (int i = (grid.size.nz*grid.size.ny*grid.size.nx - grid.size.nz*grid.size.ny); i < grid.size.nz*grid.size.ny*grid.size.nx; i++){    
        for (int j : grid.blocks[i]) {
            distance_x = bmax[0]-particles[i].px;
            if (distance_x < 0){
                particles[j].px = bmax[0] + distance_x;
                particles[j].vx = -particles[j].vx;
                particles[j].hvx = -particles[j].hvx;
            }
        }
    }
}

void Y_boundary_interaction(std::vector<Particle> &particles,  Grid &grid) {

    double distance_y;
    for (int i = 0; i < grid.size.nx; i++){
        for (int j = 0, k=grid.size.nz*(grid.size.ny-1); j < grid.size.nz; j++,k++) {  //pared Y_0
            for (int l : grid.blocks[ j + i * grid.size.nz * grid.size.ny]) {//  block_index = j + i * grid.size.nz * grid.size.ny;
                distance_y = particles[l].py - bmin[1];
                if (distance_y < 0) {
                    particles[l].py = bmin[1] - distance_y;
                    particles[l].vy = -particles[l].vy;
                    particles[l].hvy = -particles[l].hvy;
                }
            }                                                                          //pared Y_max
            for (int l : grid.blocks[k + i * grid.size.nz * grid.size.ny]) { //  block_index = k + i * grid.size.nz * grid.size.ny;
                distance_y = bmax[1] - particles[l].py;
                if (distance_y < 0) {
                    particles[l].py = bmax[1] + distance_y;
                    particles[l].vy = -particles[l].vy;
                    particles[l].hvy = -particles[l].hvy;
                }
            }
        }
    }
}

void Z_boundary_interaction(std::vector<Particle> &particles, Grid &grid) {

    double distance_z;
    for (int i = 0,j =grid.size.nz-1; i < grid.size.nz*grid.size.ny*grid.size.nx; i+=grid.size.nz, j+=grid.size.nz){
        for (int l : grid.blocks[i]) {                                       //pared Z_0
            distance_z = particles[l].py - bmin[2];
            if (distance_z < 0) {
                particles[l].py = bmin[2] - distance_z;
                particles[l].vz = -particles[l].vz;
                particles[l].hvz = -particles[l].hvz;
            }
        }
        for (int l : grid.blocks[j]) {                                       //pared Z_max
            distance_z = bmax[2] - particles[l].py;
            if (distance_z < 0) {
                particles[l].py = bmax[2] + distance_z;
                particles[l].vz = -particles[l].vz;
                particles[l].hvz = -particles[l].hvz;
            }
        }
    }
}


void particle_collision(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &accelerations){
    particle_collision_with_X_axis(particles,grid,accelerations);
    particle_collision_with_Y_axis(particles,grid,accelerations);
    particle_collision_with_Z_axis(particles,grid,accelerations);
}
void boundary_collision(std::vector<Particle> &particles, Grid &grid){
    X_boundary_interaction(particles,grid);
    Y_boundary_interaction(particles,grid);
    Z_boundary_interaction(particles,grid);
}


void particles_motion(std::vector<Particle> &particles, Grid &grid, std::vector <Acceleration> &accelerations){
    for (int i = 0; i < particles.size(); i++){

        double move_x = particles[i].hvx*time_step + accelerations[i].ax*pow(time_step,2);
        double move_y = particles[i].hvy*time_step + accelerations[i].ay*pow(time_step,2);
        double move_z = particles[i].hvz*time_step + accelerations[i].az*pow(time_step,2);

        int old_block = find_block(particles[i],grid.size);

        particles[i].px += move_x;
        particles[i].py += move_y;
        particles[i].pz += move_z;
        particles[i].vx = particles[i].hvx + (accelerations[i].ax*time_step)/2;
        particles[i].vy = particles[i].hvy + (accelerations[i].ay*time_step)/2;
        particles[i].vz = particles[i].hvz + (accelerations[i].az*time_step)/2;
        particles[i].hvx = particles[i].hvx + accelerations[i].ax*time_step;
        particles[i].hvy = particles[i].hvy + accelerations[i].ay*time_step;
        particles[i].hvz = particles[i].hvz + accelerations[i].az*time_step;

        int new_block = find_block(particles[i],grid.size);

        if (old_block!=new_block){
            grid.blocks[new_block].push_back(i);

            auto x = grid.blocks[old_block].begin();
            while ( *x != i){x++;}

            grid.blocks[old_block].erase(x);

        }

    }
}

Grid grid_initialization(Initial_values &initialValues,vector<Particle> &particles){
    double boxx = bmax[0] - bmin[0];
    double boxy = bmax[1] - bmin[1];
    double boxz = bmax[2] - bmin[2];
    GridSize gridSize;
    gridSize.nx = floor(boxx/initialValues.h);
    gridSize.ny = floor(boxy/initialValues.h);
    gridSize.nz = floor(boxz/initialValues.h);

    //cout << "\n r: " << r << " ppm: " << ppm;

    gridSize.sx = boxx/gridSize.nx;
    gridSize.sy = boxy/gridSize.ny;
    gridSize.sz = boxz/gridSize.nz;
    std::vector<vector <int>> blocks = gridCreation(particles,gridSize);
    Grid grid;
    grid.size = gridSize;
    grid.blocks = blocks;
    cout<<"After : "<<grid.blocks.size()<<'\n';

    for (int i = 0; i < particles.size(); i++){
        int index = find_block(particles[i],grid.size);
        //cout<<i<<' ';
        grid.blocks[index].push_back(i);

    }

    return grid;
}


///------------------------



void simulate(int nsteps, std::vector<Particle> &particles, Grid &grid){
    for(int i=0;i < nsteps; i++) {
        vector<double> densities;
        vector<Acceleration> accelerations;

        //Stages of Simulation
        //Stage 2: Accelerations computation
        accelerations_computation(particles, grid, densities, accelerations);
        //Stage 3: Particle Collisions
        particle_collision(particles,grid,accelerations);
        //stage 4: Particles motion
        particles_motion(particles, grid, accelerations);
        //Stage 5: Boundary collisions
        boundary_collision(particles,grid);
    }
}


int main(int argc, char** argv) {
    /*
    cout<<"h1\n";
    //check_command_errors(argc,argv);
    cout<<"h2\n";
    Initial_values initialValues;
    cout<<"h3\n";
    std::vector<Particle> particles = initial_read(argv[1],initialValues);
    cout<<"h4\n";
    h = initialValues.h;
    m = initialValues.m;
    Grid grid = grid_initialization(initialValues,particles);
    cout<<"h5\n";
    simulate(stoi(argv[2]),particles,grid);
    cout<<"h6\n";
    write_to_file(argv[3],particles);
    cout<<"h7\n";
     */

    Initial_values initialValues;
    std::vector<Particle> myparticles = initial_read(argv[2],initialValues);
    cout<<"\nNum particles: "<<myparticles.size()<<'\n';
    //write_to_file(argv[3],myparticles, initialValues);

    Grid grid = grid_initialization(initialValues,myparticles);


    cout<<"grid.blocks.size() = "<<grid.blocks.size()<<'\n';
    int saved_particles = 0;
    for (int i = 0; i< grid.blocks.size();i++) {

        if (grid.blocks[i].size()>0) {
            cout << '\n';
            cout<<"grid.blocks["<<i<<"].size() = "<<grid.blocks[i].size()<<'\n';
            saved_particles += grid.blocks[i].size();
        }
        for (int x = 0; x < grid.blocks[i].size(); x++) {
            cout << grid.blocks[i][x] << ' ';
        }
    }

    cout<<'\n'<<"saved_particles = "<<saved_particles;
    cout<<"nx  = "<< grid.size.nx<<"    ny  = "<< grid.size.ny<<"     nz  = "<< grid.size.nz <<'\n';

    cout<<"\n\n\n Nuevooooo.... \n\n\n";
    std::vector<double> densities;
    std::vector<Acceleration> accelerations;

    //check_trace("../trz/small/densinc-base-1.trz",grid,myparticles,densities,accelerations);

    /*/ TRAZE  PARTICLE COLLITION
    load_trace("../trz/small/acctransf-base-1.trz",grid,myparticles,densities,accelerations,initialValues);
    int num_block = 0;
    Particle mypart101;
    mypart101= myparticles[101];
    num_block = find_block( mypart101, grid.size); // ver que la particula 101 en el traze pertence al bloque 1 pero para nosotros al 315
    cout<<"num block is "<<num_block;
    particle_collision(myparticles,grid,accelerations);
    cout<<"nx  = "<< grid.size.nx<<"    ny  = "<< grid.size.ny<<"     nz  = "<< grid.size.nz <<'\n';
    check_trace("../trz/small/partcol-base-1.trz",grid,myparticles,densities,accelerations);
    */


    //simulate(1,myparticles,grid);
    write_to_file(argv[3],myparticles,initialValues);

    return 0;
}
