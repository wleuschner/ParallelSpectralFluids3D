#include"Shader.h"
#include<GL/glew.h>
#include<vector>
#include<string>
#include<fstream>
#include<streambuf>

Shader::Shader(unsigned int shaderType)
{
    this->shaderType = shaderType;
    id = glCreateShader(shaderType);
}

Shader::Shader(unsigned int shaderType,const std::string& fileName)
{
    this->shaderType = shaderType;
    id = glCreateShader(shaderType);
    std::ifstream t(fileName.c_str());
    std::string source((std::istreambuf_iterator<char>(t)),
                     std::istreambuf_iterator<char>());
    const char* c_str = source.c_str();
    glShaderSource(id,1,&c_str,NULL);
}

Shader::~Shader()
{
    glDeleteShader(id);
}

bool Shader::compile()
{
    glCompileShader(id);

    int compiled;
    glGetShaderiv(id, GL_COMPILE_STATUS, &compiled);
    return compiled==GL_TRUE;
}

std::string Shader::compileLog()
{
    GLint maxLength = 0;
    glGetShaderiv(id, GL_INFO_LOG_LENGTH, &maxLength);
    std::vector<GLchar> infoLog(maxLength);
    glGetShaderInfoLog(id, maxLength, &maxLength, &infoLog[0]);
    std::string ret(infoLog.begin(),infoLog.end());
    return ret;
}
