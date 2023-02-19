#include <iostream>
#include <fstream>
#include <random>
#include <vector>

int main() {
    using namespace std;
    const int N = 1000;
    double points_per_length = 100.0;
    vector<double> potential;
    vector<double> to;

    random_device gen;
    uniform_real_distribution<double> dist;

    for (double i = 0; i < 2*N; ++i) {
        potential.push_back(pow((i - N) / points_per_length, 2));
    }
    
    // ‘‚¢‚Ä
        std::fstream data_file_write("C:/Users/harun/experiment/discrete.csv", ios_base::out | ios_base::binary);
        data_file_write.seekp(0);
        for (int i = 0; i < 2*N; ++i) {
            double val = potential[i];
            data_file_write.write(reinterpret_cast<char*>(&val), sizeof(val));
        }
    
    data_file_write.close();
    
    // “Ç‚Þ
        std::fstream data_file_read("C:/Users/harun/experiment/discrete.csv", ios_base::in | ios_base::binary);
        data_file_read.seekg(0);
        for (int i = 0; i < 2*N; ++i) {
            double val;
            data_file_read.read(reinterpret_cast<char*>(&val), sizeof(val));
            to.push_back(val);
        }
    
    
    
    for (int i = 0; i < 2*N; ++i) {
        cout << potential[i] << '\t' << to[i];
        cout << (potential[i] == to[i] ? " ok" : " ??") << endl;
    }

    data_file_read.close();
    
}

