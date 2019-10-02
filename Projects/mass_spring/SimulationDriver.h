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

    SimulationDriver()
      // : dt((T)0.00001) 
      : dt((T)0.0015)  // 150 times bigger dt than explicit. We can't go arbitrarily large because we are still doing approximations to the non-linear problem using taylor expansion.
    {
        gravity.setZero();
        gravity(1) = -9.8;

        sphere_center = TV::Ones()*0.5;
        sphere_radius = 0.2;
        ground = 0.1;
    }

    void run(const int max_frame)
    {
        for(int frame=1; frame<max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int)(((T)1/24)/dt);
            for (int step = 1; step <= N_substeps; step++) {
                std::cout << "Step " << step << std::endl;
                // advanceOneStepExplicitIntegration();
		advanceOneStepImplicitIntegration();
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
    
    void advanceOneStepImplicitIntegration()
    {
        int N_points = ms.x.size();
        int N_dof = dim*N_points;
        SpMat A(N_dof,N_dof);
        A.reserve(Eigen::VectorXi::Constant(N_dof,dim*20)); // estimate non-zero entries per column

        // build right hand side
        std::vector<TV> fn;
        ms.evaluateSpringForces(fn);

	Vec rhs(N_dof);
        rhs.setZero();
	///////////////////////////////////////////////
	//  ASSIGNMENT ////////////////////////////////
	//  Add you code to build f /////////////////// 
	///////////////////////////////////////////////
  
        // build the matrix
        // Mass matrix contribution (assembly to global)
        for(int p=0; p<N_points; p++) {
            for (int d = 0; d < dim; d++) {
                int i = p * dim + d; // global dof index
                A.coeffRef(i, i) += ms.m[p] / (dt * dt);
            }
        }

        for(size_t e=0; e<ms.segments.size(); e++)
        {
            int particle[2]; // global particle index
            particle[0] = ms.segments[e](0);
            particle[1] = ms.segments[e](1);


            T l0 = ms.rest_length[e];
            T E = ms.youngs_modulus;
            T l = (ms.x[particle[0]]-ms.x[particle[1]]).norm();
            TV n = (ms.x[particle[0]]-ms.x[particle[1]])/l;
            T b = ms.damping_coeff;

            // Damping matrix contribution
            Eigen::Matrix<T,dim,dim> bnn = b * n * n.transpose();
            Eigen::Matrix<T,dim*2,dim*2> G_local;
            G_local.template block<dim,dim>(0,0) = -bnn;
            G_local.template block<dim,dim>(dim,0) = bnn;
            G_local.template block<dim,dim>(0,dim) = bnn;
            G_local.template block<dim,dim>(dim,dim) = -bnn;

	    ////////////////////////////////////////////////////////////////// 
	    //  ASSIGNMENT /////////////////////////////////////////////////// 
	    //  Add you code to construct local elasticity matrix K_local
	    /////////////////////////////////////////////// //////////////////

	    ////////////////////////////////////////////////////////////////// 
	    //  ASSIGNMENT /////////////////////////////////////////////////// 
	    //  Add you code to add contributions of elasticity and damping to A
	    //  Note that you need to take care of dirichlet-0 nodes in the
	    //     corresponding row and columns (by keeping those entries 0)
	    /////////////////////////////////////////////// //////////////////
        }

        // process dirichlet-0 nodes at the diagonal of A
        for (size_t p=0; p<ms.node_is_fixed.size(); p++){
            if(ms.node_is_fixed[p]){
                for(int d=0; d<dim; d++) {
                    A.coeffRef(dim * p + d, dim * p + d) = 1;
                    rhs[p*dim+d] = 0;
                }
            }
        }

        // Eigen::ConjugateGradient<SpMat , Eigen::Lower|Eigen::Upper> krylov;
        Eigen::MINRES<SpMat, Eigen::Lower|Eigen::Upper> krylov;

        krylov.setTolerance((T)1e-7);
        krylov.compute(A);
        Vec dx = krylov.solve(rhs);
        std::cout << "#iterations:     " << krylov.iterations() << std::endl;
        std::cout << "estimated error: " << krylov.error()      << std::endl;

        for(int p=0; p<N_points; p++){

            TV tentative_dx;
            for(int d=0; d<dim; d++)
                tentative_dx(d) = dx(dim*p+d);

            TV x0 = ms.x[p];
            TV v0 = tentative_dx/dt;
            TV xt = x0 + tentative_dx;

            // ground impulse
            if(xt(1) < ground){
                T v0n = (ground - x0(1))/dt;

		// This is some hack friction. Real friction force should depend on normal force magnitude.
		T friction_coefficient = 0.1;
		v0(0)*=friction_coefficient;
		v0(2)*=friction_coefficient;
		
                v0(1) = v0n;
                xt = x0 + v0*dt;
                tentative_dx = xt-x0;
            }

            // sphere impulse
            T distance_to_sphere_center = (xt-sphere_center).norm();
            if (distance_to_sphere_center < sphere_radius)
            {
                TV inward_normal = (sphere_center-xt).normalized();
                TV v0normal = v0.dot(inward_normal) * inward_normal;
                TV v0tangent = v0 - v0normal;
                v0normal -= (sphere_radius-distance_to_sphere_center)/dt*inward_normal;
                v0 = v0normal + v0tangent;
                xt = x0 + v0*dt;
                tentative_dx = xt-x0;
            }

            ms.x[p] += tentative_dx;
            ms.v[p] = tentative_dx/dt;
        }

    }



};
