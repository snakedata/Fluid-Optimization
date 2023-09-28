#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
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
double distance_squared(Particle p1, Particle p2){
    double dx = p1.px - p2.px;
    double dy = p1.py - p2.py;
    double dz = p1.pz - p2.pz;

    return dx * dx + dy * dy + dz * dz;
};

int main(int argc, char** argv) {

    double ppm = 2;
    int np;

    //Parte de lectura de Juan crear todas las particulas
    ifstream file("/Users/raulpineda/iCloud Drive (Archive)/Document/Arquitectura/Lab1/small.fld", ios::binary);

    if (!file.is_open()) {
        cerr << "Error opening file." << std::endl;
        return 1;
    }

    // Read ppm and np


    file.read(reinterpret_cast<char*>(&ppm), sizeof(float));
    file.read(reinterpret_cast<char*>(&np), sizeof(int));

    std::cout << "ppm: " << ppm << ", np: " << np << std::endl;

    // Create a vector to store particles
    std::vector<Particle> particles;

    for (int i = 0; i < np; ++i) {
        Particle particle;
        file.read(reinterpret_cast<char*>(&particle), sizeof(Particle));

        particles.push_back(particle);
    }

    file.close();
    int counter;
    // Access particles
    for (const Particle& particle : particles) {
        std::cout << "Particle Data:" << counter  << std::endl;
        std::cout << "px: " << particle.px << ", py: " << particle.py << ", pz: " << particle.pz << std::endl;
        std::cout << "hvx: " << particle.hvx << ", hvy: " << particle.hvy << ", hvz: " << particle.hvz << std::endl;
        std::cout << "vx: " << particle.vx << ", vy: " << particle.vy << ", vz: " << particle.vz << std::endl;
        counter +=1;
    }

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

    double mass = p*pow(ppm,3);
    double h = r/ppm;


    //Parte del grid creation

    int b_min[3] = {1,2,3};
    int b_max[3] = {4,5,6};
    int boxx = b_max[0] - b_min[0];
    int boxy = b_max[1] - b_min[1];
    int boxz = b_max[2] - b_min[2];




    //Density computation
    double d;
    double i_density;
    double j_density;
    double density_inc;
    double acceleration_inc;
    double density_constant_transformation_part = 64*M_PI*pow(h,9);
    double densityConstantTransformation =  315 * mass / density_constant_transformation_part;

    for (auto i = particles.begin();i !=particles.end(); i++){
        for (auto j = i++;j != particles.end(); j++){
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
                i->hvx = i->hvx + ((i->px-j->px)*15*mass*(h-d)*(h-d)*(i->p+j->p -p)+45*(i->vx-j->vx)*nu*mass)/(M_PI*pow(h,6)*i->p*j->p);
                i->hvy = i->hvy + ((i->py-j->py)*15*mass*(h-d)*(h-d)*(i->p+j->p -p)+45*(i->vy-j->vy)*nu*mass)/(M_PI*pow(h,6)*i->p*j->p);
                i->hvz = i->hvz + ((i->pz-j->pz)*15*mass*(h-d)*(h-d)*(i->p+j->p -p)+45*(i->vz-j->vz)*nu*mass)/(M_PI*pow(h,6)*i->p*j->p);
            }
        }
    }

//Collisions



//Write floats use standard conversion static cast or narrow cast






}
