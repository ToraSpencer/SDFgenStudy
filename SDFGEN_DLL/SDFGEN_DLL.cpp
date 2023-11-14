#include "SDFGEN_DLL.h"
#include "makelevelset2.h"
#include "makelevelset3.h"

template<typename TV, typename TI>
bool genTriMeshSDF3D(SDF_RESULT& result, const TRIANGLE_MESH::triMesh<TV, TI>& mesh, \
        const float step, const int interCounts)
{
    result.interCounts = 0;
    result.step = 0;
    result.origin = verF{0, 0, 0};
    result.stepsCount.clear();
    result.SDFvalues.clear();

    // 0.
    int interCounts0 = interCounts;
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    if (interCounts < 1)
        interCounts0 = 1;

    // 1. 读取输入网格，生成网格包围盒:  
    std::vector<Vec3f> vers;                  // 输入网格顶点
    std::vector<Vec3ui> tris;                 // 输入网格三角片；
    for(const auto& ver: mesh.vertices)
    { 
        Vec3f point{static_cast<float>(ver.x), static_cast<float>(ver.y), static_cast<float>(ver.z)};
        vers.push_back(point);
        update_minmax(point, min_box, max_box); 
    }

    for(const auto& tri: mesh.triangles) 
        tris.push_back(Vec3ui{static_cast<unsigned int>(tri.x), static_cast<unsigned int>(tri.y), \
            static_cast<unsigned int>(tri.z)});
  
    // 2. Add interCounts around the box.
    Vec3f unit(1, 1, 1);
    min_box -= interCounts0 * step * unit;
    max_box += interCounts0 * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step); 

    // 3. 计算距离场： 
    SDF_GEN::Array3f SDFvalues;
    make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], SDFvalues);

    // 4. 输出：
    const Array1<float>& SDFdata = SDFvalues.a;
    result.step = step;
    result.interCounts = interCounts0;
    result.origin = verF{ min_box[0], min_box[1], min_box[2] };
    result.stepsCount = std::vector<unsigned>{ sizes[0], sizes[1], sizes[2] };
    result.SDFvalues.resize(SDFdata.n);
    std::memcpy(&result.SDFvalues[0], SDFdata.data, sizeof(float) * SDFdata.n);

	return true;
}


 // 输入三角网格，生成3D情形下的符号距离场； 
 SDFGEN_DLL_API bool genSDF3D(SDF_RESULT& result, \
     const triMeshF& mesh, const float step, const int interCounts)
 {
     return genTriMeshSDF3D<float, int>(result, mesh, step, interCounts);
 }

 SDFGEN_DLL_API bool genSDF3D(SDF_RESULT& result, \
     const triMeshD& mesh, const float step, const int interCounts)
 {
     return genTriMeshSDF3D<double, int>(result, mesh, step, interCounts);
 }
