#ifndef DEC_H
#define DEC_H
#include <Eigen/Eigen>
#include "decmesh3d.h"
#include "../Geometry/mesh3d.h"

Eigen::SparseMatrix<double> hodge2(DECMesh3D& mesh,bool dual=false);
Eigen::SparseMatrix<double> hodge1(DECMesh3D& mesh,bool dual=false);
Eigen::SparseMatrix<double> hodge0(DECMesh3D& mesh,bool dual=false);

Eigen::SparseMatrix<double> derivative2(DECMesh3D& mesh,bool dual=false);
Eigen::SparseMatrix<double> derivative1(DECMesh3D& mesh,bool dual=false);
Eigen::SparseMatrix<double> derivative0(DECMesh3D& mesh,bool dual=false);



#endif // DEC_H
