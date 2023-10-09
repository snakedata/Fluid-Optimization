#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>


int main() {
    std::ofstream outfile ("new.fld",std::ofstream::binary);

    int numParticles = 10.0;
    float ppm = 2.3;

    outfile.write(reinterpret_cast<const char *>(&ppm), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&numParticles), sizeof(int));

    for (int i = 0; i < numParticles; i++){
        for (int x = 0; i<10; i++){
            float f = static_cast <float>(i);
            outfile.write(reinterpret_cast<const char*>(&f), sizeof(float));

        }
    }
    outfile.close();
    return 0;
}
