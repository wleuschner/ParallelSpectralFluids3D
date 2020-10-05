#include"Model.h"
#include <glm/gtc/matrix_transform.hpp>
#include <GL/glew.h>
#include <QtOpenGL/QGLFramebufferObject>
#include <assimp/Importer.hpp>
#include <assimp/postprocess.h>
#include <assimp/mesh.h>
#include <assimp/scene.h>
#include <cmath>
#include <iostream>
#include "../FrameBufferObject/FrameBufferObject.h"
ShaderProgram* Model::voxelProgram = NULL;

Model::Model()
{
    if(voxelProgram==NULL)
    {
        voxelProgram = new ShaderProgram();
        Shader vert(GL_VERTEX_SHADER,"Res/Voxel/voxel.vert");
        if(!vert.compile())
        {
            std::cout<<vert.compileLog()<<std::endl;
        }
        Shader frag(GL_FRAGMENT_SHADER,"Res/Voxel/voxel.frag");
        if(!frag.compile())
        {
            std::cout<<frag.compileLog()<<std::endl;
        }
        Shader geom(GL_GEOMETRY_SHADER,"Res/Voxel/voxel.geom");
        if(!geom.compile())
        {
            std::cout<<geom.compileLog()<<std::endl;
        }
        voxelProgram->attachShader(vert);
        voxelProgram->attachShader(geom);
        voxelProgram->attachShader(frag);
        if(!voxelProgram->link())
        {
            std::cout<<voxelProgram->linkLog()<<std::endl;
        }
        voxelProgram->bind();
    }
}

Model::~Model()
{
}

bool Model::load(std::string path)
{
    Assimp::Importer importer;
    modelMat = glm::mat4(1.0f);
    const aiScene* scene = importer.ReadFile(path,aiProcessPreset_TargetRealtime_Quality);
    for(unsigned int j=0;j</*scene->mNumMeshes*/1;j++)
    {
        const aiMesh* mesh = scene->mMeshes[j];
        aiColor3D mat_ambient;
        aiColor3D mat_diffuse;
        aiColor3D mat_specular;

        int mat_index = mesh->mMaterialIndex;
        const aiMaterial* mat = scene->mMaterials[mat_index];

        float mat_shininess;

        mat->Get(AI_MATKEY_COLOR_AMBIENT,mat_ambient);
        mat->Get(AI_MATKEY_COLOR_DIFFUSE,mat_diffuse);
        mat->Get(AI_MATKEY_COLOR_SPECULAR,mat_specular);
        mat->Get(AI_MATKEY_SHININESS,mat_shininess);
        position.reserve(mesh->mNumVertices);
        for(unsigned int i=0;i<mesh->mNumVertices;i++)
        {
            aiVector3D v = mesh->mVertices[i];
            glm::vec3 pos = glm::vec3(v.x,v.y,v.z);
            aabb.max.x = glm::max(pos.x,aabb.max.x);
            aabb.max.y = glm::max(pos.y,aabb.max.y);
            aabb.max.z = glm::max(pos.z,aabb.max.z);
            aabb.min.x = glm::min(pos.x,aabb.min.x);
            aabb.min.y = glm::min(pos.y,aabb.min.y);
            aabb.min.z = glm::min(pos.z,aabb.min.z);
            position.push_back(pos);
            if(mesh->HasNormals())
            {
                aiVector3D n = mesh->mNormals[i];
                normal.push_back(glm::vec3(n.x,n.y,n.z));
            }
            //if(mesh->HasTextureCoords())
            //{
            //    aiVector2D uv = mesh->mTextureCoords[i];
            //    uv_coords.push_back(glm::vec2(uv[0],uv[1]));
            //}
            Vertex vert;
            vert.pos = position[i];
            vert.normal = normal[i];
            vert.uv = glm::vec2(0,0);
            vertices.push_back(vert);
        }
        for(int i=0;i<mesh->mNumFaces;i++)
        {
            indices.push_back(mesh->mFaces[i].mIndices[0]);
            indices.push_back(mesh->mFaces[i].mIndices[1]);
            indices.push_back(mesh->mFaces[i].mIndices[2]);
        }
    }
    createVBO();
    createIndex();
    return true;
}

bool Model::release()
{
    delete vbo;
    delete index;
    indices.clear();
    vertices.clear();
    uv_coords.clear();
    normal.clear();
    return true;
}

void Model::bind()
{
    vbo->bind();
    index->bind();
}

void Model::update()
{
    vbo->bind();
    vbo->upload(vertices);
    bind();
}

std::vector<Vertex>& Model::getVertices()
{
    return vertices;
}

std::vector<unsigned int>& Model::getIndices()
{
    return indices;
}

void Model::setModelMat(const glm::mat4& mat)
{
    this->modelMat = mat;
}

glm::mat4 Model::getModelMat()
{
    return glm::translate(this->modelMat,aabb.getCenter());
}

AABB Model::getAABB()
{
    return aabb;
}

bool Model::createVBO()
{
    vbo = new VertexBuffer();
    vbo->bind();
    vbo->upload(vertices);
    return true;
}

bool Model::createIndex()
{
    index = new IndexBuffer();
    index->bind();
    index->upload(indices);
    return true;
}

Material Model::getMaterial() const
{
    return material;
}

void Model::setMaterial(const Material &value)
{
    material = value;
}

void Model::draw()
{
    /*bind();

    glm::vec3 ambient = material.getAmbient();
    glm::vec3 diffuse = material.getDiffuse();
    glm::vec3 specular = material.getSpecular();
    float shininess = material.getShininess();*/

    //shader->uploadVec3("material.amb",ambient);
    //shader->uploadVec3("material.dif",diffuse);
    //shader->uploadVec3("material.spec",specular);
    //shader->uploadScalar("material.shininess",shininess);


    glDrawElements(GL_TRIANGLES,indices.size(),GL_UNSIGNED_INT,(void*)0);
}

DECMesh3D Model::voxelize(float resolution)
{
    glm::vec3 ext = getAABB().getExtent();

    unsigned int layers = ceil(ext.z/(resolution*32*4));
    unsigned int depth = ceil(ext.z/(resolution));
    unsigned int width = ceil(ext.x/(resolution));
    unsigned int height = ceil(ext.y/(resolution));

    glm::mat4 orthProj = glm::ortho(-ext.x/2,ext.x/2,-ext.y/2,ext.y/2);
    glm::mat4 pv = orthProj*glm::translate(glm::mat4(1.0f),glm::vec3(-getAABB().getCenter().x,-getAABB().getCenter().y,-getAABB().getCenter().z));


    DECMesh3D decMesh = DECMesh3D(resolution,glm::uvec3(width,height,depth),resolution,aabb.min);
    FrameBufferObject fbo;
    TextureArray texArray;
    texArray.bind(0);
    /*ext.x/particleSize,ext.y/particleSize,((ext.z/32)/particleSize)*/
    texArray.createRenderArray(width,height,layers);
    texArray.unbind(0);
    //texArray.createRenderArray(10,10,512);

    fbo.bind();
    fbo.attachColorArray(texArray,0);
    //fbo.attachDepthArray(depthArray);
    fbo.setRenderBuffer({GL_COLOR_ATTACHMENT0});
    if(!fbo.isComplete())
    {
        std::cout<<"FBO incomplete"<<std::endl;
    }
    fbo.bind();
    bind();
    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();
    voxelProgram->bind();
    voxelProgram->uploadMat4("projection",pv);
    voxelProgram->uploadScalar("zOffset",getAABB().getExtent().z/2);
    voxelProgram->uploadScalar("resolution",resolution);
    texArray.bind(0);

    glViewport(0,0,width,height);
    glEnable(GL_DEPTH_CLAMP);
    glDepthRange(-ext.z/2,ext.z/2);
    glEnable(GL_COLOR_LOGIC_OP);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_BLEND);
    glLogicOp(GL_XOR);
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
    draw();
    glEnable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_LOGIC_OP);
    glDepthRange(-1.0,1.0);
    glDisable(GL_DEPTH_CLAMP);
    glFlush();
    fbo.unbind();
    vertices.clear();
    vertices.resize(decMesh.getNumPoints());


    //Read Voxels from Texture Array
    unsigned int* buffer = new unsigned int[4*width*height*layers];
    //glReadBuffer(GL_COLOR_ATTACHMENT0);
    //glReadPixels(0,0,width,height,GL_RGBA,GL_UNSIGNED_INT,buffer);
    texArray.bind(0);
    glGetTexImage(GL_TEXTURE_2D_ARRAY,0,GL_RGBA_INTEGER,GL_UNSIGNED_INT,buffer);

    unsigned int v1,v2,v3,v4,v5,v6,v7,v8;
    for(unsigned int l=0;l<layers;l++)
    {
        //glGetTextureSubImage(GL_TEXTURE_2D_ARRAY,0,0,0,l,width,height,0,GL_RGBA,GL_UNSIGNED_INT,4*width*height*layers,buffer);
        for(unsigned int y=0;y<height;y++)
        {
            for(unsigned int x=0;x<width;x++)
            {

                //Red Channel
                for(unsigned int d=0;d<32;d++)
                {
                    if(buffer[(x*4)+width*4*y+(width*4*height)*l]&1)
                    {
                        unsigned int z = (32*4*l+d);
                        v1 = decMesh.getPointIndex(x,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1);
                        v2 = decMesh.getPointIndex(x+1,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1+1);
                        v3 = decMesh.getPointIndex(x+1,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1+1);
                        v4 = decMesh.getPointIndex(x,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1);

                        v5 = decMesh.getPointIndex(x,y,z);//((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1);
                        v6 = decMesh.getPointIndex(x+1,y,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1+1);
                        v7 = decMesh.getPointIndex(x+1,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1+1);
                        v8 = decMesh.getPointIndex(x,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1);
                        decMesh.setVoxelInside(Voxel3D(decMesh.getVoxelIndex(x,y,z),GridState::INSIDE));
                        float xp = getAABB().min.x+x*resolution;
                        float yp = getAABB().min.y+y*resolution;
                        float zp = getAABB().min.z+z*resolution;
                        glm::vec3 p1(xp,yp+resolution,zp);
                        glm::vec3 p2(xp+resolution,yp+resolution,zp);
                        glm::vec3 p3(xp+resolution,yp+resolution,zp+resolution);
                        glm::vec3 p4(xp,yp+resolution,zp+resolution);

                        glm::vec3 p5(xp,yp,zp);
                        glm::vec3 p6(xp+resolution,yp,zp);
                        glm::vec3 p7(xp+resolution,yp,zp+resolution);
                        glm::vec3 p8(xp,yp,zp+resolution);

                        vertices[v1] = p1;
                        vertices[v2] = p2;
                        vertices[v3] = p3;
                        vertices[v4] = p4;
                        vertices[v5] = p5;
                        vertices[v6] = p6;
                        vertices[v7] = p7;
                        vertices[v8] = p8;
                    }
                    buffer[(x*4)+width*4*y+(width*4*height)*l]>>=1;
                }
                //Green Channel
                for(unsigned int d=0;d<32;d++)
                {
                    if(buffer[(x*4)+width*4*y+(width*4*height)*l+1]&1)
                    {
                        unsigned int z = (32*4*l+d+32);
                        v1 = decMesh.getPointIndex(x,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1);
                        v2 = decMesh.getPointIndex(x+1,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1+1);
                        v3 = decMesh.getPointIndex(x+1,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1+1);
                        v4 = decMesh.getPointIndex(x,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1);

                        v5 = decMesh.getPointIndex(x,y,z);//((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1);
                        v6 = decMesh.getPointIndex(x+1,y,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1+1);
                        v7 = decMesh.getPointIndex(x+1,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1+1);
                        v8 = decMesh.getPointIndex(x,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1);
                        decMesh.setVoxelInside(Voxel3D(decMesh.getVoxelIndex(x,y,z),GridState::INSIDE));
                        float xp = getAABB().min.x+x*resolution;
                        float yp = getAABB().min.y+y*resolution;
                        float zp = getAABB().min.z+z*resolution;
                        glm::vec3 p1(xp,yp+resolution,zp);
                        glm::vec3 p2(xp+resolution,yp+resolution,zp);
                        glm::vec3 p3(xp+resolution,yp+resolution,zp+resolution);
                        glm::vec3 p4(xp,yp+resolution,zp+resolution);

                        glm::vec3 p5(xp,yp,zp);
                        glm::vec3 p6(xp+resolution,yp,zp);
                        glm::vec3 p7(xp+resolution,yp,zp+resolution);
                        glm::vec3 p8(xp,yp,zp+resolution);

                        vertices[v1] = p1;
                        vertices[v2] = p2;
                        vertices[v3] = p3;
                        vertices[v4] = p4;
                        vertices[v5] = p5;
                        vertices[v6] = p6;
                        vertices[v7] = p7;
                        vertices[v8] = p8;
                    }
                    buffer[(x*4)+width*4*y+(width*4*height)*l+1]>>=1;
                }
                //Blue Channel
                for(unsigned int d=0;d<32;d++)
                {
                    if(buffer[(x*4)+width*4*y+(width*4*height)*l+2]&1)
                    {
                        unsigned int z = (32*4*l+d+64);
                        v1 = decMesh.getPointIndex(x,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1);
                        v2 = decMesh.getPointIndex(x+1,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1+1);
                        v3 = decMesh.getPointIndex(x+1,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1+1);
                        v4 = decMesh.getPointIndex(x,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1);

                        v5 = decMesh.getPointIndex(x,y,z);//((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1);
                        v6 = decMesh.getPointIndex(x+1,y,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1+1);
                        v7 = decMesh.getPointIndex(x+1,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1+1);
                        v8 = decMesh.getPointIndex(x,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1);
                        decMesh.setVoxelInside(Voxel3D(decMesh.getVoxelIndex(x,y,z),GridState::INSIDE));
                        float xp = getAABB().min.x+x*resolution;
                        float yp = getAABB().min.y+y*resolution;
                        float zp = getAABB().min.z+z*resolution;
                        glm::vec3 p1(xp,yp+resolution,zp);
                        glm::vec3 p2(xp+resolution,yp+resolution,zp);
                        glm::vec3 p3(xp+resolution,yp+resolution,zp+resolution);
                        glm::vec3 p4(xp,yp+resolution,zp+resolution);

                        glm::vec3 p5(xp,yp,zp);
                        glm::vec3 p6(xp+resolution,yp,zp);
                        glm::vec3 p7(xp+resolution,yp,zp+resolution);
                        glm::vec3 p8(xp,yp,zp+resolution);

                        vertices[v1] = p1;
                        vertices[v2] = p2;
                        vertices[v3] = p3;
                        vertices[v4] = p4;
                        vertices[v5] = p5;
                        vertices[v6] = p6;
                        vertices[v7] = p7;
                        vertices[v8] = p8;
                    }
                    buffer[(x*4)+width*4*y+(width*4*height)*l+2]>>=1;
                }
                //Alpha Channel
                for(unsigned int d=0;d<32;d++)
                {
                    if(buffer[(x*4)+width*4*y+(width*4*height)*l+3]&1)
                    {
                        unsigned int z = (32*4*l+d+96);
                        v1 = decMesh.getPointIndex(x,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1);
                        v2 = decMesh.getPointIndex(x+1,y+1,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1+1))+(x+1+1);
                        v3 = decMesh.getPointIndex(x+1,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1+1);
                        v4 = decMesh.getPointIndex(x,y+1,z+1);//((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1+1))+(x+1);

                        v5 = decMesh.getPointIndex(x,y,z);//((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1);
                        v6 = decMesh.getPointIndex(x+1,y,z);//(((width+1+2)*(height+1+2))*(32*4*l+d+1))+((width+1+2)*(y+1))+(x+1+1);
                        v7 = decMesh.getPointIndex(x+1,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1+1);
                        v8 = decMesh.getPointIndex(x,y,z+1);//(((width+1+2)*(height+1+2))*((32*4*l+d+1+1)))+((width+1+2)*(y+1))+(x+1);
                        decMesh.setVoxelInside(Voxel3D(decMesh.getVoxelIndex(x,y,z),GridState::INSIDE));
                        float xp = getAABB().min.x+x*resolution;
                        float yp = getAABB().min.y+y*resolution;
                        float zp = getAABB().min.z+z*resolution;
                        glm::vec3 p1(xp,yp+resolution,zp);
                        glm::vec3 p2(xp+resolution,yp+resolution,zp);
                        glm::vec3 p3(xp+resolution,yp+resolution,zp+resolution);
                        glm::vec3 p4(xp,yp+resolution,zp+resolution);

                        glm::vec3 p5(xp,yp,zp);
                        glm::vec3 p6(xp+resolution,yp,zp);
                        glm::vec3 p7(xp+resolution,yp,zp+resolution);
                        glm::vec3 p8(xp,yp,zp+resolution);

                        vertices[v1] = p1;
                        vertices[v2] = p2;
                        vertices[v3] = p3;
                        vertices[v4] = p4;
                        vertices[v5] = p5;
                        vertices[v6] = p6;
                        vertices[v7] = p7;
                        vertices[v8] = p8;
                    }
                    buffer[(x*4)+width*4*y+(width*4*height)*l+3]>>=1;
                }
            }
        }
    }
    texArray.unbind(0);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);

    delete[] buffer;
    return decMesh;
}
