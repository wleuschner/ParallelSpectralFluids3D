#ifndef MODEL_H
#define MODEL_H
#include<string>
#include<vector>
#include<map>
#include<set>
#include<glm/glm.hpp>
#include"../../Graphics/Shader/Shader.h"
#include"../../Graphics/ShaderProgram/ShaderProgram.h"
#include"../../Graphics/VertexBuffer/VertexBuffer.h"
#include"../../Graphics/IndexBuffer/IndexBuffer.h"
#include"../../Graphics/Material/Material.h"
#include"../../Graphics/AABB/AABB.h"
#include"../../DEC/decmesh3d.h"

class Model
{
public:
    Model();
    ~Model();
    bool load(std::string path);
    bool release();
    void bind();
    void draw();

    DECMesh3D voxelize(float resolution);

    Material getMaterial() const;
    void setMaterial(const Material &value);
    std::vector<Vertex>& getVertices();
    std::vector<unsigned int>& getIndices();

    void setModelMat(const glm::mat4& mat);
    glm::mat4 getModelMat();

    AABB getAABB();

    void update();

private:
    bool createVBO();
    bool createIndex();

    std::string name;

    static ShaderProgram* voxelProgram;

    AABB aabb;
    glm::mat4 modelMat;

    Material material;
    std::vector<unsigned int> indices;

    std::vector<glm::vec3> position;
    std::vector<glm::vec3> normal;
    std::vector<glm::vec2> uv_coords;

    std::vector<Vertex> vertices;

    VertexBuffer* vbo;
    IndexBuffer* index;
};
#endif
