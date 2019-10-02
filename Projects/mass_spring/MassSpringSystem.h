#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class MassSpringSystem{
public:
    using TV = Eigen::Matrix<T,dim,1>;
    
    std::vector<Eigen::Matrix<int,2,1> > segments;
    std::vector<T> m;
    std::vector<TV> x;
    std::vector<TV> v;
    T youngs_modulus;
    T damping_coeff;
    std::vector<bool> node_is_fixed;
    std::vector<T> rest_length;

    MassSpringSystem()
    {}

    void evaluateSpringForces(std::vector<TV >& f)
    {
        f.clear();
        f.resize(x.size(), TV::Zero());
        for (size_t i=0; i<segments.size(); ++i)
        {
            int A = segments[i](0);
            int B = segments[i](1);
            T l = (x[A]-x[B]).norm();
            TV n = (x[A]-x[B]).normalized();
            TV force = -youngs_modulus*(l/rest_length[i]-(T)1)*n;
            f[A] += force;
            f[B] -= force;
        }
    }

    void evaluateDampingForces(std::vector<TV >& f)
    {
        f.clear();
        f.resize(x.size(), TV::Zero());
        for (size_t i=0; i<segments.size(); ++i)
        {
            int A = segments[i](0);
            int B = segments[i](1);
            TV n = (x[A]-x[B]).normalized();
            Eigen::Matrix<T,dim,dim> nn = n * n.transpose();
            TV force = -damping_coeff * nn * (v[A]-v[B]);
            f[A] += force;
            f[B] -= force;
        }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : x) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        for (const Eigen::Matrix<int, 2, 1>& seg : segments)
            fs << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
        fs << "END\n";
        fs.close();
    }
    
};
