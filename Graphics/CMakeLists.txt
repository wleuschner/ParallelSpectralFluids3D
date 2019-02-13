cmake_minimum_required(VERSION 2.8.8)

add_subdirectory("Shader")
add_subdirectory("ShaderProgram")
add_subdirectory("VertexBuffer")
add_subdirectory("IndexBuffer")
add_subdirectory("FrameBufferObject")
add_subdirectory("Texture")
add_subdirectory("TextureArray")
add_subdirectory("Light")
add_subdirectory("Vertex")
add_subdirectory("Particle")
#add_subdirectory("ParticleBuffer")
add_subdirectory("Camera")
add_subdirectory("Model")
add_subdirectory("Material")
add_subdirectory("AABB")
add_subdirectory("CubeMap")

set(HEADERS
   ${HEADERS}
   PARENT_SCOPE
)

set(SOURCE
   ${SOURCE}
   PARENT_SCOPE
)
