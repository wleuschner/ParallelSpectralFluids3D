#ifndef __SHADER_H_
#define __SHADER_H_
#include<string>

class Shader
{
    friend class ShaderProgram;
public:
    Shader(unsigned int shaderType);
    Shader(unsigned int shaderType,const std::string& fileName);
    ~Shader();
    bool compile();
    std::string compileLog();
private:
    unsigned int id;
    unsigned int shaderType;
};

#endif
