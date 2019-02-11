#ifndef Mesh3D_H
#define Mesh3D_H
#include "../DEC/face3d.h"
#include "../DEC/edge3d.h"
#include "../Graphics/Vertex/Vertex.h"
#include <map>
#include <set>
#include <Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <assimp/scene.h>

class Mesh3D
{
public:
    Mesh3D();
    Mesh3D(float resolution,aiMesh* data);
    ~Mesh3D();
    unsigned int getResolution();
    unsigned int getNumVertices();
    unsigned int getNumEdges();
    unsigned int getNumFaces();
    unsigned int getNumVoxels();
    Eigen::VectorXf getVelocityField();
    void buildLaplace();

    friend Eigen::SparseMatrix<float> hodge2(Mesh3D& mesh,bool inverse);
    friend Eigen::SparseMatrix<float> hodge1(Mesh3D& mesh,bool inverse);
    friend Eigen::SparseMatrix<float> hodge0(Mesh3D& mesh,bool inverse);

    friend Eigen::SparseMatrix<float> derivative0(Mesh3D& mesh,bool dual);
    friend Eigen::SparseMatrix<float> derivative1(Mesh3D& mesh,bool dual);

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& getBasisField();
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& getBasisCoefficients();

    void integrate();
private:
    float resolution;

    void voxelize();
    bool checkVoxel(unsigned int x,unsigned int y,unsigned int z);

    aiMesh *data;

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenValues;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenVectors;

    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenVortValues;
    Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> eigenVortVectors;

public:
    Eigen::VectorXf velocityField;
    std::vector<glm::vec2> velocityField2;
    std::vector<glm::vec2> vorticityField2;
    Eigen::VectorXf vorticityField;
    std::vector<Eigen::MatrixXd> advection;
    Eigen::SparseMatrix<double> curl;
    Eigen::SparseMatrix<double> b2;
    Eigen::SparseMatrix<double> b1;

    std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int,
                        unsigned int,unsigned int,unsigned int,unsigned int>> voxels;
    std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>> faces;
    std::set<std::tuple<unsigned int,unsigned int>> edges;
    std::set<unsigned int> points;

    std::map<unsigned int,Vertex> vertex;
};

#endif // Mesh3D_H
