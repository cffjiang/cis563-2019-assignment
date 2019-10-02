#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "SimulationDriver.h"

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T,dim,1>;

    SimulationDriver<T,dim> driver;

    // set up mass spring system
    T youngs_modulus = 2;
    T damping_coeff = 2; // 0

    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    std::vector<bool> node_is_fixed;

    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> rest_length;

    int N = 64;
    int N_points = N*N;
    T dx = (T)1/(N-1);
    m.resize(N_points);
    x.resize(N_points);
    v.resize(N_points);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            int id = i*N+j;
            m[id] = (T)1/N_points;
            x[id](0) = (i-1)*dx;
            x[id](2) = (j-1)*dx;
            x[id](1) = 1;
            v[id] = TV::Zero();
        }
    }
    node_is_fixed.resize(N_points,false);
    // for(int i=0; i<N; i++)
    //     node_is_fixed[i] = true;

    
//    node_is_fixed[0] = true;
//    node_is_fixed[(N-1)*N] = true;
//    node_is_fixed[N_points-1] = true;
//    node_is_fixed[(N-1)] = true;
//    for(int i=(N-1)*N; i<N_points; i++)
//        node_is_fixed[i] = true;


    // structure
    for(int i=0; i<N-1; i++){
        for(int j=0; j<N; j++){
            Eigen::Matrix<int,2,1> seg;
            int p = i*N+j, q=(i+1)*N+j;
            seg << p,q;
            segments.push_back(seg);
            rest_length.push_back((x[p]-x[q]).norm());
        }
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N-1; j++){
            Eigen::Matrix<int,2,1> seg;
            int p = i*N+j, q=i*N+j+1;
            seg << p,q;
            segments.push_back(seg);
            rest_length.push_back((x[p]-x[q]).norm());
        }
    }

    // shear
    for(int i=0; i<N-1; i++){
        for(int j=0; j<N-1; j++){
            Eigen::Matrix<int,2,1> seg;
            int p = i*N+j, q=(i+1)*N+j+1;
            seg << p,q;
            segments.push_back(seg);
            rest_length.push_back((x[p]-x[q]).norm());
        }
    }
    for(int i=0; i<N-1; i++){
        for(int j=0; j<N-1; j++){
            Eigen::Matrix<int,2,1> seg;
            int p = (i+1)*N+j, q=i*N+j+1;
            seg << p,q;
            segments.push_back(seg);
            rest_length.push_back((x[p]-x[q]).norm());
        }
    }

    //bending
    for(int i=0; i<N-2; i++){
        for(int j=0; j<N; j++){
            Eigen::Matrix<int,2,1> seg;
            int p = i*N+j, q=(i+2)*N+j;
            seg << p,q;
            segments.push_back(seg);
            rest_length.push_back((x[p]-x[q]).norm());
        }
    }
    for(int i=0; i<N; i++){
        for(int j=0; j<N-2; j++){
            Eigen::Matrix<int,2,1> seg;
            int p = i*N+j, q=i*N+j+2;
            seg << p,q;
            segments.push_back(seg);
            rest_length.push_back((x[p]-x[q]).norm());
        }
    }


    // simulate
    driver.ms.segments = segments;
    driver.ms.m = m;
    driver.ms.v = v;
    driver.ms.x = x;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.damping_coeff = damping_coeff;
    driver.ms.node_is_fixed = node_is_fixed;
    driver.ms.rest_length = rest_length;

    driver.run(240);

    return 0;
}
