#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MassSpringSystem.h"


template<class T, int dim>
class SimulationDriver{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MassSpringSystem<T,dim> ms;
    T dt;
    TV gravity;

    TV sphere_center;
    T sphere_radius;
    T ground;
    T collision_stiffness;

    SimulationDriver()
    : dt((T)0.00001) // 0.0015 for implicit
    {
        gravity.setZero();
        gravity(1) = -9.8;

        sphere_center = TV::Ones()*0.5;
        sphere_radius = 0.2;
        ground = 0.1;
        collision_stiffness = 1.1*0;
    }

    void run(const int max_frame)
    {
        for(int frame=1; frame<max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                std::cout << "Step " << step << std::endl;
                advanceOneStepExplicitIntegration();
            }
            mkdir("output/", 0777);
            std::string filename = "output/" + std::to_string(frame) + ".poly";
            ms.dumpPoly(filename);
            std::cout << std::endl;
        }
    }

    void advanceOneStepExplicitIntegration()
    {
        int N_points = ms.x.size();
        int N_dof = dim*N_points;
	std::vector<TV> f_spring;
        ms.evaluateSpringForces(f_spring);
	std::vector<TV> f_damping;
	ms.evaluateDampingForces(f_damping);
	
	for(int p=0; p<N_points; p++){
            if(ms.node_is_fixed[p]){
	      ms.v[p] = TV::Zero();
            }
	    else{
	      ms.v[p] += ((f_spring[p]+f_damping[p])/ms.m[p]+gravity)*dt;
	      ms.x[p] += ms.v[p]*dt;
	    }
        }
    }


};
