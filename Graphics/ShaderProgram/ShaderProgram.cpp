#include"ShaderProgram.h"
#include<GL/glew.h>
#include<glm/gtc/type_ptr.hpp>

ShaderProgram::ShaderProgram()
{
    id = glCreateProgram();
}

void ShaderProgram::attachShader(const Shader &shader)
{
    glAttachShader(id,shader.id);
}

bool ShaderProgram::link()
{
    glLinkProgram(id);
    GLint linked = 0;
    glGetProgramiv(id, GL_LINK_STATUS, &linked);
    return linked==GL_TRUE;
}

void ShaderProgram::bind()
{
    glUseProgram(id);
}

std::string ShaderProgram::linkLog()
{
    GLint maxLength = 0;
    glGetProgramiv(id, GL_INFO_LOG_LENGTH, &maxLength);

    //The maxLength includes the NULL character
    std::vector<GLchar> infoLog(maxLength);
    glGetProgramInfoLog(id, maxLength, &maxLength, &infoLog[0]);
    std::string ret(infoLog.begin(),infoLog.end());
    return ret;
}

unsigned int ShaderProgram::getAttribLocation(const std::string& name)
{
    return glGetAttribLocation(id,name.c_str());
}

void ShaderProgram::setAttribute(unsigned int loc,unsigned int type,unsigned int offset,unsigned int n,unsigned int size)
{
    glVertexAttribPointer(loc,3,GL_FLOAT,GL_FALSE,24,(void*)0);
}

void ShaderProgram::setAttribute(const std::string& name,unsigned int type,unsigned int offset,unsigned int n,unsigned int size)
{
    glVertexAttribPointer(glGetAttribLocation(id,name.c_str()),3,GL_FLOAT,GL_FALSE,24,(void*)0);
}

void ShaderProgram::enableAttribute(unsigned int loc)
{
    glEnableVertexAttribArray(loc);
}

void ShaderProgram::enableAttribute(const std::string& name)
{
    glEnableVertexAttribArray(glGetAttribLocation(id,name.c_str()));
}

void ShaderProgram::uploadInt(const std::string &var, unsigned int val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform1i(loc,val);
}

void ShaderProgram::uploadUnsignedInt(const std::string &var, unsigned int val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform1ui(loc,val);
}

void ShaderProgram::uploadScalar(const std::string& var,float val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform1f(loc,val);
}

void ShaderProgram::uploadVec2(const std::string& var,glm::vec2 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform2fv(loc,1,glm::value_ptr(val));
}

void ShaderProgram::uploadVec3(const std::string& var,glm::vec3 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform3fv(loc,1,glm::value_ptr(val));
}

void ShaderProgram::uploadVec4(const std::string& var,glm::vec4 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform4fv(loc,1,glm::value_ptr(val));
}

void ShaderProgram::uploadIVec2(const std::string& var,glm::ivec2 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform2iv(loc,1,glm::value_ptr(val));
}

void ShaderProgram::uploadIVec3(const std::string& var,glm::ivec3 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform3iv(loc,1,glm::value_ptr(val));
}

void ShaderProgram::uploadIvec4(const std::string& var,glm::ivec4 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniform4iv(loc,1,glm::value_ptr(val));
}

void ShaderProgram::uploadMat2(const std::string& var,glm::mat2 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniformMatrix2fv(loc,1,GL_FALSE,glm::value_ptr(val));
}

void ShaderProgram::uploadMat3(const std::string& var,glm::mat3 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniformMatrix3fv(loc,1,GL_FALSE,glm::value_ptr(val));
}

void ShaderProgram::uploadMat4(const std::string& var,glm::mat4 val)
{
    unsigned int loc = glGetUniformLocation(id,var.c_str());
    glUniformMatrix4fv(loc,1,GL_FALSE,glm::value_ptr(val));
}

void ShaderProgram::uploadLight(const std::string& var,const Light& val,const glm::mat4& view)
{
    glm::vec3 pos = val.pos;
    glm::vec3 ldir = glm::vec3(view*glm::vec4(val.pos,1.0));
    glm::vec3 ambient = val.amb;
    glm::vec3 diffuse = val.diff;
    glm::vec3 specular = val.spec;

    std::string pos_string = std::string("light.pos");
    std::string ldir_string = std::string("light.ldir");
    std::string ambient_string = std::string("light.amb");
    std::string diffuse_string = std::string("light.dif");
    std::string specular_string = std::string("light.spec");

    uploadVec3(pos_string.c_str(),pos);
    uploadVec3(ldir_string.c_str(),ldir);
    uploadVec3(ambient_string.c_str(),ambient);
    uploadVec3(diffuse_string.c_str(),diffuse);
    uploadVec3(specular_string.c_str(),specular);
}

glm::ivec3 ShaderProgram::getMaxWorkGroupSize()
{
    glm::ivec3 size;
    glGetIntegeri_v( GL_MAX_COMPUTE_WORK_GROUP_SIZE, 0, &size.x );
    glGetIntegeri_v( GL_MAX_COMPUTE_WORK_GROUP_SIZE, 1, &size.y );
    glGetIntegeri_v( GL_MAX_COMPUTE_WORK_GROUP_SIZE, 2, &size.z );
    return size;
}

glm::ivec3 ShaderProgram::getMaxWorkGroups()
{
    glm::ivec3 size;
    glGetIntegeri_v( GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &size.x );
    glGetIntegeri_v( GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &size.y );
    glGetIntegeri_v( GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &size.z );
    return size;
}

void ShaderProgram::dispatch(unsigned int gwx,unsigned int gwy,unsigned int gwz,unsigned int lwx,unsigned int lwy,unsigned int lwz)
{
    glDispatchComputeGroupSizeARB(gwx,gwy,gwz,lwx,lwy,lwz);
}
