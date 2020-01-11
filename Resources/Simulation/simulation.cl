#pragma extension cl_khr_gl_sharing : enable
#define WARP_SHIFT 4
#define GRP_SHIFT 4
#define BANK_OFFSET(n)     (((n) >> WARP_SHIFT) + ((n) >> GRP_SHIFT))

uint getVoxelIndex(uint x,uint y,uint z,uint4 dims)
{
    return z*(dims.x*dims.y)+y*dims.x+x;
}

uint getZFaceIndex(uint x,uint y,uint z,uint4 dims)
{
    return (z*(dims.x*dims.y)+y*dims.x+x);
}

uint getYFaceIndex(uint x,uint y,uint z,uint4 dims)
{
    return (z*((dims.y+1)*(dims.x))+y*(dims.x)+x)+(dims.z+1)*(dims.x*dims.y);
}

uint getXFaceIndex(uint x,uint y,uint z,uint4 dims)
{
    return (z*((dims.y)*(dims.x+1))+y*(dims.x+1)+x+(dims.z+1)*(dims.x*dims.y)+(dims.y+1)*(dims.x*dims.z));
}

int getFaceSignum(uint vid,uint fidx,__global const uint* signBitString)
{
    int idx=vid/5;
    int offset=vid%5;
    return signBitString[idx]>>(offset*6+(fidx))&1?1:-1;
}

__attribute__((reqd_work_group_size(1024, 1, 1)))
__kernel void advection(
                        __global float4* inParticles,
                        __global const double* velField,
                        __global const uint* signBitString,
                        const uint4 dims,
                        const float4 aabb_min,
                        const float4 aabb_max,
                        const float4 aabb_extent,
                        const float resolution,
                        const float timestep,
                        const uint3 randoms,
                        const float lifeTime)
{

    uint idx = get_global_id(0);

    uint seed = randoms.x * idx;
    uint t = seed ^ (seed << 11);
    uint rand1 = seed >> 16;
    seed = (seed * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
    uint rand2 = seed >> 16;
    seed = (seed * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
    uint rand3 = seed >> 16;
    seed = (seed * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
    uint rand4 = seed >> 16;

    int yOfs = floor((aabb_min.y+aabb_extent.y+inParticles[idx].y)/resolution);
    int xOfs = floor((aabb_min.x+aabb_extent.x+inParticles[idx].x)/resolution);
    int zOfs = floor((aabb_min.z+aabb_extent.z+inParticles[idx].z)/resolution);

    int zOfsMinus = floor((aabb_min.z+aabb_extent.z+inParticles[idx].z-resolution/2)/resolution);
    int zOfsPlus = floor((aabb_min.z+aabb_extent.z+inParticles[idx].z+resolution/2)/resolution);
    int yOfsMinus = floor((aabb_min.y+aabb_extent.y+inParticles[idx].y-resolution/2)/resolution);
    int yOfsPlus = floor((aabb_min.y+aabb_extent.y+inParticles[idx].y+resolution/2)/resolution);
    int xOfsMinus = floor((aabb_min.x+aabb_extent.x+inParticles[idx].x-resolution/2)/resolution);
    int xOfsPlus = floor((aabb_min.x+aabb_extent.x+inParticles[idx].x+resolution/2)/resolution);

    xOfs = clamp(xOfs,0,convert_int(dims.x-1));
    yOfs = clamp(yOfs,0,convert_int(dims.y-1));
    zOfs = clamp(zOfs,0,convert_int(dims.z-1));

    xOfsMinus = clamp(xOfsMinus,0,convert_int(dims.x-1));
    xOfsPlus = clamp(xOfsPlus,0,convert_int(dims.x-1));
    yOfsMinus = clamp(yOfsMinus,0,convert_int(dims.y-1));
    yOfsPlus = clamp(yOfsPlus,0,convert_int(dims.y-1));
    zOfsMinus = clamp(zOfsMinus,0,convert_int(dims.z-1));
    zOfsPlus = clamp(zOfsPlus,0,convert_int(dims.z-1));

    uint v1xIdx = getVoxelIndex(xOfs,yOfsMinus,zOfsMinus,dims);
    uint v2xIdx = getVoxelIndex(xOfs,yOfsPlus,zOfsMinus,dims);
    uint v3xIdx = getVoxelIndex(xOfs,yOfsMinus,zOfsPlus,dims);
    uint v4xIdx = getVoxelIndex(xOfs,yOfsPlus,zOfsPlus,dims);

    uint v1yIdx = getVoxelIndex(xOfsMinus,yOfs,zOfsMinus,dims);
    uint v2yIdx = getVoxelIndex(xOfsPlus,yOfs,zOfsMinus,dims);
    uint v3yIdx = getVoxelIndex(xOfsMinus,yOfs,zOfsPlus,dims);
    uint v4yIdx = getVoxelIndex(xOfsPlus,yOfs,zOfsPlus,dims);

    uint v1zIdx = getVoxelIndex(xOfsMinus,yOfsMinus,zOfs,dims);
    uint v2zIdx = getVoxelIndex(xOfsPlus,yOfsMinus,zOfs,dims);
    uint v3zIdx = getVoxelIndex(xOfsMinus,yOfsPlus,zOfs,dims);
    uint v4zIdx = getVoxelIndex(xOfsPlus,yOfsPlus,zOfs,dims);

    uint v1xf5 = getXFaceIndex(xOfs,yOfsMinus,zOfsMinus,dims);
    uint v1xf6 = getXFaceIndex(xOfs+1,yOfsMinus,zOfsMinus,dims);
    uint v2xf5 = getXFaceIndex(xOfs,yOfsPlus,zOfsMinus,dims);
    uint v2xf6 = getXFaceIndex(xOfs+1,yOfsPlus,zOfsMinus,dims);
    uint v3xf5 = getXFaceIndex(xOfs,yOfsMinus,zOfsPlus,dims);
    uint v3xf6 = getXFaceIndex(xOfs+1,yOfsMinus,zOfsPlus,dims);
    uint v4xf5 = getXFaceIndex(xOfs,yOfsPlus,zOfsPlus,dims);
    uint v4xf6 = getXFaceIndex(xOfs+1,yOfsPlus,zOfsPlus,dims);

    uint v1yf3 = getYFaceIndex(xOfsMinus,yOfs,zOfsMinus,dims);
    uint v1yf4 = getYFaceIndex(xOfsMinus,yOfs+1,zOfsMinus,dims);
    uint v2yf3 = getYFaceIndex(xOfsPlus,yOfs,zOfsMinus,dims);
    uint v2yf4 = getYFaceIndex(xOfsPlus,yOfs+1,zOfsMinus,dims);
    uint v3yf3 = getYFaceIndex(xOfsMinus,yOfs,zOfsPlus,dims);
    uint v3yf4 = getYFaceIndex(xOfsMinus,yOfs+1,zOfsPlus,dims);
    uint v4yf3 = getYFaceIndex(xOfsPlus,yOfs,zOfsPlus,dims);
    uint v4yf4 = getYFaceIndex(xOfsPlus,yOfs+1,zOfsPlus,dims);

    uint v1zf1 = getZFaceIndex(xOfsMinus,yOfsMinus,zOfs,dims);
    uint v1zf2 = getZFaceIndex(xOfsMinus,yOfsMinus,zOfs+1,dims);
    uint v2zf1 = getZFaceIndex(xOfsPlus,yOfsMinus,zOfs,dims);
    uint v2zf2 = getZFaceIndex(xOfsPlus,yOfsMinus,zOfs+1,dims);
    uint v3zf1 = getZFaceIndex(xOfsMinus,yOfsPlus,zOfs,dims);
    uint v3zf2 = getZFaceIndex(xOfsMinus,yOfsPlus,zOfs+1,dims);
    uint v4zf1 = getZFaceIndex(xOfsPlus,yOfsPlus,zOfs,dims);
    uint v4zf2 = getZFaceIndex(xOfsPlus,yOfsPlus,zOfs+1,dims);

    float3 cf1x = aabb_min.xyz+resolution*(float3)(xOfs,yOfsMinus,zOfsMinus)+0.5f*(float3)(0.0f,resolution,resolution);
    float3 cf1y = aabb_min.xyz+resolution*(float3)(xOfsMinus,yOfs,zOfsMinus)+0.5f*(float3)(resolution,0.0f,resolution);
    float3 cf1z = aabb_min.xyz+resolution*(float3)(xOfsMinus,yOfsMinus,zOfs)+0.5f*(float3)(resolution,resolution,0.0f);


    float vel1x,vel2x,vel3x,vel4x,vel5x,vel6x,vel7x,vel8x;
    vel1x = -getFaceSignum(v1xIdx,4,signBitString)*convert_float(velField[v1xf5]);
    vel3x = getFaceSignum(v1xIdx,5,signBitString)*convert_float(velField[v1xf6]);
    vel2x = -getFaceSignum(v2xIdx,4,signBitString)*convert_float(velField[v2xf5]);
    vel4x = getFaceSignum(v2xIdx,5,signBitString)*convert_float(velField[v2xf6]);
    vel5x = -getFaceSignum(v3xIdx,4,signBitString)*convert_float(velField[v3xf5]);
    vel7x = getFaceSignum(v3xIdx,5,signBitString)*convert_float(velField[v3xf6]);
    vel6x = -getFaceSignum(v4xIdx,4,signBitString)*convert_float(velField[v4xf5]);
    vel8x = getFaceSignum(v4xIdx,5,signBitString)*convert_float(velField[v4xf6]);

    float vel1y,vel2y,vel3y,vel4y,vel5y,vel6y,vel7y,vel8y;
    vel1y = -getFaceSignum(v1yIdx,2,signBitString)*convert_float(velField[v1yf3]);
    vel2y = getFaceSignum(v1yIdx,3,signBitString)*convert_float(velField[v1yf4]);
    vel3y = -getFaceSignum(v2yIdx,2,signBitString)*convert_float(velField[v2yf3]);
    vel4y = getFaceSignum(v2yIdx,3,signBitString)*convert_float(velField[v2yf4]);
    vel5y = -getFaceSignum(v3yIdx,2,signBitString)*convert_float(velField[v3yf3]);
    vel6y = getFaceSignum(v3yIdx,3,signBitString)*convert_float(velField[v3yf4]);
    vel7y = -getFaceSignum(v4yIdx,2,signBitString)*convert_float(velField[v4yf3]);
    vel8y = getFaceSignum(v4yIdx,3,signBitString)*convert_float(velField[v4yf4]);

    float vel1z,vel2z,vel3z,vel4z,vel5z,vel6z,vel7z,vel8z;
    vel1z = getFaceSignum(v1zIdx,0,signBitString)*convert_float(velField[v1zf1]);
    vel2z = getFaceSignum(v3zIdx,0,signBitString)*convert_float(velField[v3zf1]);
    vel3z = getFaceSignum(v2zIdx,0,signBitString)*convert_float(velField[v2zf1]);
    vel4z = getFaceSignum(v4zIdx,0,signBitString)*convert_float(velField[v4zf1]);
    vel5z = -getFaceSignum(v1zIdx,1,signBitString)*convert_float(velField[v1zf2]);
    vel6z = -getFaceSignum(v3zIdx,1,signBitString)*convert_float(velField[v3zf2]);
    vel7z = -getFaceSignum(v2zIdx,1,signBitString)*convert_float(velField[v2zf2]);
    vel8z = -getFaceSignum(v4zIdx,1,signBitString)*convert_float(velField[v4zf2]);

    float3 vel=(float3)(0.0,0.0,0.0);

    float3 particleNormalizedX = (1.0f/resolution)*(inParticles[idx].xyz-cf1x);

    //if(particleNormalizedX.x<0.0f) printf("ERROR X 1 %f %f %f %f %f %f %f %f %f",particleNormalizedX.x,particleNormalizedX.y,particleNormalizedX.z,inParticles[idx].x,inParticles[idx].y,inParticles[idx].z,cf1x.x,cf1x.y,cf1x.z);
    //if(particleNormalizedX.y<0.0f) printf("ERROR X 2");
    //if(particleNormalizedX.z<0.0f) printf("ERROR X 3");

    float4 xVelXInterp = mix((float4)(vel1x,vel3x,vel5x,vel7x),(float4)(vel2x,vel4x,vel6x,vel8x),clamp(particleNormalizedX.y,0.0f,1.0f));
    float2 xVelYInterp = mix((float2)(xVelXInterp.x,xVelXInterp.z),(float2)(xVelXInterp.y,xVelXInterp.w),clamp(particleNormalizedX.x,0.0f,1.0f));
    vel.x = mix(xVelYInterp.x,xVelYInterp.y,clamp(particleNormalizedX.z,0.0f,1.0f));

    float3 particleNormalizedY = (1.0f/resolution)*(inParticles[idx].xyz-cf1y);

    //if(particleNormalizedY.x<0.0f) printf("ERROR Y 1");
    //if(particleNormalizedY.y<0.0f) printf("ERROR Y 2");
    //if(particleNormalizedY.z<0.0f) printf("ERROR Y 3");

    float4 yVelXInterp = mix((float4)(vel1y,vel3y,vel5y,vel7y),(float4)(vel2y,vel4y,vel6y,vel8y),clamp(particleNormalizedY.y,0.0f,1.0f));
    float2 yVelYInterp = mix((float2)(yVelXInterp.x,yVelXInterp.z),(float2)(yVelXInterp.y,yVelXInterp.w),clamp(particleNormalizedY.x,0.0f,1.0f));
    vel.y = mix(yVelYInterp.x,yVelYInterp.y,clamp(particleNormalizedY.z,0.0f,1.0f));

    float3 particleNormalizedZ = (1.0f/resolution)*(inParticles[idx].xyz-cf1z);

    //if(particleNormalizedZ.x<0.0f) printf("ERROR Z 1");
    //if(particleNormalizedZ.y<0.0f) printf("ERROR Z 2");
    //if(particleNormalizedZ.z<0.0f) printf("ERROR Z 3");

    float4 zVelXInterp = mix((float4)(vel1z,vel3z,vel5z,vel7z),(float4)(vel2z,vel4z,vel6z,vel8z),clamp(particleNormalizedZ.y,0.0f,1.0f));
    float2 zVelYInterp = mix((float2)(zVelXInterp.x,zVelXInterp.z),(float2)(zVelXInterp.y,zVelXInterp.w),clamp(particleNormalizedZ.x,0.0f,1.0f));
    vel.z = mix(zVelYInterp.x,zVelYInterp.y,clamp(particleNormalizedZ.z,0.0f,1.0f));

    float4 newPart;
    float4 aabb_center = aabb_min+0.5*(aabb_max-aabb_min);
    newPart.x = ((rand1%1024)/1024.0-0.5)*0.25;
    newPart.y = aabb_min.y+0.2f;
    //newPart.y =((rand4%1024)/1024.0-0.5)*0.5;
    newPart.z = -aabb_center.z+((rand3%1024)/1024.0-0.5)*0.25;
    //newPart.z = ((rand3%1024)/1024.0-0.5)*0.5;
    newPart.w = lifeTime*((rand2%1024)/1024.0);
    float4 oldPart = inParticles[idx];

    inParticles[idx] = mix(oldPart,newPart,convert_float(inParticles[idx].w<0.0));

    inParticles[idx].xyz = (clamp(inParticles[idx].xyz+timestep*vel,aabb_min.xyz,aabb_max.xyz));
    inParticles[idx].w -= 1.0;
}

__kernel void normalization_viscocity_gravity(__global double* e1,__global double* e2,__global double4* basisCoeff,__global double4* eigenValues,__global double4* gravity,double visc,double timestep,double gravityOn)
{
    uint idx = get_global_id(0);
    basisCoeff[idx] *= sqrt(e1[0]/e2[0]);
    basisCoeff[idx] *= exp(-visc*eigenValues[idx]*timestep);
    basisCoeff[idx] += timestep*gravity[idx]*gravityOn;
}

__kernel void update_vel(__global double4* basisCoeff,__global double4* vel,double timestep)
{
    uint idx = get_global_id(0);
    basisCoeff[idx] += timestep*vel[idx];
}

__kernel void advection_reduce_x(__global double4* advection_xyz,
                                 __global double* advection_yz,
                                 __global double4* baseCoeff,
                                 uint3 dims,
                                 __local double* temp)
{
    uint x_offset = get_global_id(0);
    uint y_offset = get_global_id(1);
    uint z_offset = get_global_id(2);

    uint advection_offset = z_offset*(dims.y*dims.x) + y_offset*(dims.x);

    uint lid = get_local_id(0);
    uint n = get_local_size(0) * 2;

    int ai = lid;
    int bi = lid + n/2;
    int bankOffsetA = BANK_OFFSET(ai);
    int bankOffsetB = BANK_OFFSET(bi);

    temp[ai + bankOffsetA] = dot(advection_xyz[advection_offset/4 + ai],baseCoeff[ai]);
    temp[bi + bankOffsetB] = dot(advection_xyz[advection_offset/4 + bi],baseCoeff[bi]);


    int offset = 1;
    for(int d = n/2;d>0;d/=2)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if(lid<d)
        {
            int ai = offset * (2*lid+1)-1;
            int bi = offset * (2*lid+2)-1;
            ai += BANK_OFFSET(ai);
            bi += BANK_OFFSET(bi);
            temp[bi] += temp[ai];
        }
        offset*=2;
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    if(lid==0)
    {
        advection_yz[z_offset*(dims.y)+y_offset] = temp[n-1 + BANK_OFFSET(n-1)];
    }
}

__kernel void advection_reduce_y(__global double4* advection_yz,
                                 __global double* advection_z,
                                 __global double4* baseCoeff,
                                 uint3 dims,
                                 __local double* temp)
{
    uint x_offset = get_global_id(0);
    uint y_offset = get_global_id(1);

    uint advection_offset = y_offset*(dims.y);

    uint lid = get_local_id(0);
    uint n = get_local_size(0) * 2;

    int ai = lid;
    int bi = lid + n/2;
    int bankOffsetA = BANK_OFFSET(ai);
    int bankOffsetB = BANK_OFFSET(bi);


    double4 a = advection_yz[advection_offset/4 + ai];
    double4 b = baseCoeff[ai];

    temp[ai + bankOffsetA] = dot(advection_yz[advection_offset/4 + ai],baseCoeff[ai]);
    temp[bi + bankOffsetB] = dot(advection_yz[advection_offset/4 + bi],baseCoeff[bi]);

    int offset = 1;
    for(int d = n/2;d>0;d/=2)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if(lid<d)
        {
            int ai = offset * (2*lid+1)-1;
            int bi = offset * (2*lid+2)-1;
            ai += BANK_OFFSET(ai);
            bi += BANK_OFFSET(bi);

            temp[bi] += temp[ai];
        }
        offset*=2;
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    if(lid==0)
    {
        advection_z[y_offset] = temp[n-1 + BANK_OFFSET(n-1)];
    }
}
