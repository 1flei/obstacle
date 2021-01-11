#pragma once

#include <vector>
#include <string>
#include <fstream>

inline void write_dataset(const std::string& filename, const std::vector<std::vector<double>>& datasets)
{
    int n = datasets.size();
    int dim = datasets[0].size();

    std::ofstream ofs(filename, std::ofstream::binary);

    ofs.write((const char*)&n, sizeof(int));
    ofs.write((const char*)&dim, sizeof(int));
    for(int i=0;i<n;i++) {
        ofs.write((const char*)&datasets[i][0], sizeof(double)*dim);
    }
}

inline std::vector<std::vector<double>> read_dataset(const std::string& filename, int nrow=-1)
{
    std::vector<std::vector<double>> ret;

    std::ifstream ifs(filename, std::ifstream::binary);

    int n, dim;

    ifs.read((char*)&n, sizeof(int));
    ifs.read((char*)&dim, sizeof(int));

    if(nrow > 0 && nrow < n){
        n = nrow;
    }

    ret.reserve(n);

    for(int i=0;i<n;i++) {
        ret.emplace_back(dim);
        ifs.read(( char*)&ret.back()[0], sizeof(double)*dim);
    }
    return ret;
}

inline void write_ground_truth(const std::string& filename, const std::vector<double>& ground_truth)
{
    int n = ground_truth.size();
    std::ofstream ofs(filename);
    ofs << n << "\n";
    for(int i=0;i<n;i++) {
        ofs << ground_truth[i] << " ";
    }
}

inline std::vector<double> read_ground_truth(const std::string& filename, int nrow=-1)
{
    int n = 0;
    std::ifstream ifs(filename, std::ofstream::binary);
    ifs >> n;

    if(nrow>0 && nrow < n){
        n = nrow;
    }
    std::vector<double> ret;
    ret.resize(n);
    for(int i=0;i<n;i++){
        ifs >> ret[i];
    }
    return ret;
}