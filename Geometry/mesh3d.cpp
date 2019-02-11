#include "mesh3d.h"
#include <iostream>
#include <eigen3/unsupported/Eigen/ArpackSupport>
#include "../DEC/dec.h"
typedef Eigen::SparseMatrix<float> SparseMat;
typedef Eigen::SimplicialLDLT<SparseMat> SparseChol;
typedef Eigen::ArpackGeneralizedSelfAdjointEigenSolver <SparseMat, SparseChol> Arpack;

Mesh3D::Mesh3D()
{
}

Mesh3D::Mesh3D(float resolution,aiMesh* data) : Mesh3D()
{
    this->data = data;
    this->resolution = resolution;
    voxelize();
    buildLaplace();
}

Mesh3D::~Mesh3D()
{
}

unsigned int Mesh3D::getNumVertices()
{
    return vertex.size();
}

unsigned int Mesh3D::getNumEdges()
{
    return edges.size();
}

unsigned int Mesh3D::getNumFaces()
{
    return faces.size();
}

unsigned int Mesh3D::getNumVoxels()
{
    return voxels.size();
}

unsigned int Mesh3D::getResolution()
{
    return resolution;
}

Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& Mesh3D::getBasisField()
{
    return eigenVectors;
}

Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>& Mesh3D::getBasisCoefficients()
{
    return eigenValues;
}

Eigen::VectorXf Mesh3D::getVelocityField()
{
    return velocityField;
}

void Mesh3D::voxelize()
{
/*    vertex.resize((height/resY+1)*(width/resX+1));
    edges.resize(2*(height/resY+2)*(width/resX+1)-(height/resY+2)-(width/resX+1));
    faces.resize(height/resY*width/resX);*/

    //faces.resize(width/w*height/h,width/w*height/h);
    unsigned int iv1,iv2,iv3,iv4;
    int e1,e2,e3,e4;
    int f;
    int y;

    unsigned int depth=0;
    unsigned int height=0;
    unsigned int width=0;

    for(int z=0;z<depth;z++)
    {
        for(int y=0;y<height;y++)
        {
            for(int x=0;x<width;x++)
            {
                iv1 = y*(width/resolution)+x+y; //UP LEFT Vert
                iv2 = y*(width/resolution)+x+1+y; //UP RIGHT Vert
                iv3 = (y+1)*(width/resolution)+x+(y+1); //DOWN LEFT vert
                iv4 = (y+1)*(width/resolution)+x+1+(y+1); //DOWN RIGHT vert*/

                f = y*(width/resolution)+x; //FACE
                if(checkVoxel(x,y,z))
                {
                    Vertex v1,v2,v3,v4;
                    v1.pos.x = x*resolution;
                    v1.pos.y = y*resolution;
                    v2.pos.x = (x+1)*resolution;
                    v2.pos.y = y*resolution;
                    v3.pos.x = x*resolution;
                    v3.pos.y = (y+1)*resolution;
                    v4.pos.x = (x+1)*resolution;
                    v4.pos.y = (y+1)*resolution;
                    /*vertex[iv1] = v1;
                    vertex[iv2] = v2;
                    vertex[iv3] = v3;
                    vertex[iv4] = v4;*/
                    std::set<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>>::const_iterator itf=faces.emplace(std::make_tuple(iv1,iv3,iv4,iv2)).first;
                    std::set<std::tuple<unsigned int,unsigned int>>::const_iterator it;
                    std::set<unsigned int>::const_iterator pit1,pit2;

                    //First Triangle
                    if((it=edges.find(std::make_tuple(iv1,iv3)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        pit1=points.find(iv1);
                        pit2=points.find(iv3);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else if((it=edges.find(std::make_tuple(iv3,iv1)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                        pit1=points.find(iv3);
                        pit2=points.find(iv1);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;

                    }
                    else
                    {
                        it = edges.emplace(std::make_tuple(iv1,iv3)).first;
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        if((pit1=points.find(iv1))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        }
                        else
                        {
                            pit1 = points.emplace(iv1).first;
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                            vertex[iv1] = v1;
                        }
                        if((pit2=points.find(iv3))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        }
                        else
                        {
                            pit2 = points.emplace(iv3).first;
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                            vertex[iv3] = v3;
                        }
                    }

                    if((it=edges.find(std::make_tuple(iv3,iv4)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        pit1=points.find(iv3);
                        pit2=points.find(iv4);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else if((it=edges.find(std::make_tuple(iv4,iv3)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                        pit1=points.find(iv4);
                        pit2=points.find(iv3);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        it = edges.emplace(std::make_tuple(iv3,iv4)).first;
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        if((pit1=points.find(iv3))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        }
                        else
                        {
                            pit1 = points.emplace(iv3).first;
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                            vertex[iv3] = v3;
                        }
                        if((pit2=points.find(iv4))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        }
                        else
                        {
                            pit2 = points.emplace(iv4).first;
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                            vertex[iv4] = v4;
                        }
                    }

                    if((it=edges.find(std::make_tuple(iv4,iv2)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        pit1=points.find(iv4);
                        pit2=points.find(iv2);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else if((it=edges.find(std::make_tuple(iv2,iv4)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                        pit1=points.find(iv2);
                        pit2=points.find(iv4);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        it = edges.emplace(std::make_tuple(iv4,iv2)).first;
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        if((pit1=points.find(iv4))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        }
                        else
                        {
                            pit1 = points.emplace(iv4).first;
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                            vertex[iv4] = v4;
                        }
                        if((pit2=points.find(iv2))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        }
                        else
                        {
                            pit2 = points.emplace(iv2).first;
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                            vertex[iv2] = v2;
                        }
                    }

                    if((it=edges.find(std::make_tuple(iv2,iv1)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        pit1=points.find(iv2);
                        pit2=points.find(iv1);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else if((it=edges.find(std::make_tuple(iv1,iv2)))!=edges.end())
                    {
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = -1;
                        pit1=points.find(iv1);
                        pit2=points.find(iv2);
                        //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                    }
                    else
                    {
                        it = edges.emplace(std::make_tuple(iv2,iv1)).first;
                        //b2.insert(std::distance(edges.begin(),it),std::distance(faces.begin(),itf)) = 1;
                        if((pit1=points.find(iv2))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                        }
                        else
                        {
                            pit1 = points.emplace(iv2).first;
                            //b1.insert(std::distance(points.begin(),pit1),std::distance(edges.begin(),it)) = 1;
                            vertex[iv2] = v2;
                        }
                        if((pit2=points.find(iv1))!=points.end())
                        {
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                        }
                        else
                        {
                            pit2 = points.emplace(iv1).first;
                            //b1.insert(std::distance(points.begin(),pit2),std::distance(edges.begin(),it)) = -1;
                            vertex[iv1] = v1;
                        }
                    }
                }
            }
        }
    }
}

bool Mesh3D::checkVoxel(unsigned int x,unsigned int y,unsigned int z)
{
    return false;
}

void Mesh3D::buildLaplace()
{
    //Build Face Matrix
    Eigen::SparseMatrix<double> mat,bound;
    mat.resize(edges.size(),edges.size());
    //curl.resize(edges.size(),edges.size());
    //mat = derivative1(*this).transpose();
    //derivative1(*this)*
    //Calculate Laplace Operator
    //mat = -1.0f*derivative1(*this)*hodge1(*this,true)*derivative1(*this).transpose()*hodge2(*this,false);
    //std::cout<<"HODGE DIMS"<<hodge1(*this,false).cols()<<" "<<hodge1(*this,false).rows()<<std::endl;
    //std::cout<<"DERIVATIVE DIMS"<<derivative1(*this,false).cols()<<" "<<derivative0(*this,false).rows()<<std::endl;

/*    curl = -1.0f*hodge2(*this,true)*derivative1(*this,true)*hodge1(*this,false);
    mat = -1.0f*derivative0(*this,false)*hodge2(*this,true)*derivative1(*this,true)*hodge1(*this,false);
    bound = derivative1(*this,false);
    for(int k=0;k<bound.outerSize();k++)
    {
        unsigned int nFaces=0;
        for(Eigen::SparseMatrix<double>::InnerIterator it(bound,k);it;++it)
        {
            nFaces++;
        }
        if(nFaces!=2)
        {
            mat.prune([k](int i,int j,float v){return i!=k;});
        }
        std::cout<<"FACES: "<<nFaces<<std::endl;
    }

    //mat = -1.0f*derivative0(*this)*hodge2(*this,true)*derivative1(*this,true)*hodge1(*this,false);

    std::cout<<"MAT DIMS: "<<mat.cols()<<" "<<mat.rows()<<std::endl;
    //Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<float>> eigenSolver;
    //eigenSolver.compute(A);
    //mat = hodge1(*this,true);//derivative1(*this).transpose()*hodge2(*this,false);
    //std::cout<<eigenSolver.eigenvalues()<<std::endl;

    std::cout<<mat<<std::endl;

    Arpack arpack;
    arpack.compute(mat,16,"LM");
    eigenValues = arpack.eigenvalues();
    eigenVectors = arpack.eigenvectors();

    std::cout<<eigenValues.transpose()<<std::endl;
    std::cout<<eigenVectors<<std::endl;

    std::cout<<"Velocity"<<std::endl;
    std::cout.precision(5);

    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    velocityField = Eigen::VectorXf::Zero(getNumEdges());
    vorticityField = Eigen::VectorXf::Zero(getNumEdges());

    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        Eigen::VectorXf e = eigenVectors.col(k);
        velocityField += (eigenValues(k,0)*e).transpose();
    }

    Eigen::VectorXd vort = curl*velocityField;
    for(std::set<std::tuple<unsigned int,unsigned int>>::iterator eit=edges.begin();eit!=edges.end();eit++)
    {
        int idx = std::distance(edges.begin(),eit);
        std::set<unsigned int>::iterator pit1,pit2;
        pit1 = points.find(std::get<0>(*eit));
        pit2 = points.find(std::get<1>(*eit));
        int pidx1 = std::distance(points.begin(),pit1);
        int pidx2 = std::distance(points.begin(),pit2);

        vorticityField(idx) = vort(pidx2) - vort(pidx1);
    }

    std::cout<<"CURL: "<<curl<<std::endl;
    std::cout.flush();

    std::cout<<vorticityField<<std::endl;
    std::cout.flush();

    std::cout.flush();*/

/*    for(unsigned int i=0;i<eigenValues.rows();i++)
    {
        for(unsigned int h=0;h<vortictyBasisField.rows();h++)
        {
            for(unsigned int g=0;g<vortictyBasisField.rows();g++)
            {
                advection(g,h,i) = (curl*())*vortictyBasisField.row(g);
            }
        }
    }*/

    for(unsigned int k=0;k<velocityField.rows();k++)
    {
        std::cout.width(10);
        std::cout<<velocityField(k)<<" ";
        if((k+1)%16==0)
        {
            std::cout<<std::endl;
        }
    }
}

void Mesh3D::integrate()
{
    std::cout<<"Advect"<<std::endl;
    /*for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        //std::cout<<(eigenValues(k,0)*eigenVectors.col(k))<<" ";
        if(k%15==0)
        {
            std::cout<<std::endl;
        }
    }*/
    float e1 = 0.0f;
    float e2 = 0.0f;
    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        e1 += eigenValues(k,0)*eigenValues(k,0);
    }
    Eigen::VectorXf vel;
    vel.resize(16);
    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        //vel(k)=(eigenValues.transpose()*advection(k)*eigenValues)(0);
    }
    eigenValues += vel*(1.0);
    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        e2 += eigenValues(k,0)*eigenValues(k,0);
    }
    eigenValues *= std::sqrt(e1/e2);
    /*velocityField.setZero();
    for(unsigned int k=0;k<eigenValues.rows();k++)
    {
        Eigen::VectorXf e = eigenVectors.col(k);
        velocityField += (eigenValues(k,0)*e).transpose();
    }*/
    std::cout<<e1<<std::endl;
}
