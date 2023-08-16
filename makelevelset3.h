#ifndef MAKELEVELSET3_H
#define MAKELEVELSET3_H

#include "array3.h"
#include "vec.h"
#include <windows.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <thread>
#include <mutex>

#include "myEigen.h"

#define USE_MULTITHREADS

#ifdef USE_MULTITHREADS
    static std::mutex g_mutex;
#endif
 


/*
     tri is a list of triangles in the mesh, 

     and x is the positions of the vertices absolute distances will be nearly correct for triangle soup, 
     but a closed mesh is needed for accurate signs. 
     Distances for all grid cells within exact_band cells of a triangle should be exact; 
     further away a distance is calculated but it might not be to the closest triangle - just one nearby.
*/
void make_level_set3(const std::vector<Vec3ui> &tri, const std::vector<Vec3f> &x,
                     const Vec3f &origin, float dx, int nx, int ny, int nz,
                     SDF_GEN::Array3f &phi, const int exact_band=1);

 
float point_segment_distance(const Vec3f& x0, const Vec3f& x1, const Vec3f& x2);
float point_triangle_distance(const Vec3f& x0, const Vec3f& x1, const Vec3f& x2, const Vec3f& x3);
void check_neighbour(const std::vector<Vec3ui>& tris, const std::vector<Vec3f>& vers,
    SDF_GEN::Array3f& DFvalues, SDF_GEN::Array3i& closest_tris,
    const Vec3f& ver0, int i0, int j0, int k0, int i1, int j1, int k1);
void sweep(const std::vector<Vec3ui>& tris, const std::vector<Vec3f>& vers,
    SDF_GEN::Array3f& DFvalues, SDF_GEN::Array3i& closest_tris, const Vec3f& origin, float step,
    int di, int dj, int dk);
int orientation(double x1, double y1, double x2, double y2, double& twice_signed_area);
bool point_in_triangle_2d(double x0, double y0,
    double x1, double y1, double x2, double y2, double x3, double y3,
    double& a, double& b, double& c);



///////////////////////////////////////////////////////////////////// 一些debug接口：
  
// 打印三维矩阵的某一层切面：
template <typename T, class ArrayT = std::vector<T>>
void dispArrSlice(const SDF_GEN::Array3<T, ArrayT>& arr, const int k)
{
    for (unsigned i = 0; i<arr.ni; ++i) 
    {
        for (unsigned j = 0; j < arr.nj; ++j)
            std::cout << arr(i, j, k) << ", ";
        std::cout << std::endl;
    }
}




#endif
