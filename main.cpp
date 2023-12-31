#include "makelevelset2.h"
#include "makelevelset3.h"
#include "config.h"
#include <list> 


// 静态库：
#include "myEigenIO/myEigenIO.h"
#pragma comment(lib,"myEigenIO.lib")	

#include "myEigenBasicMath/myEigenBasicMath.h"
#pragma comment(lib, "myEigenBasicMath.lib")

#include "myEigenModeling/myEigenModeling.h"
#pragma comment(lib, "myEigenModeling.lib")

#include "myEigenPMP/myEigenPMP.h"
#pragma comment(lib, "myEigenPMP.lib")

 
// 在指定空间内生成三角网格的符号距离场；
/*
    SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
    Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
    ...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
    This code is public domain. Feel free to mess with it, let me know if you like it.
*/


// 栅格坐标系： 用顶点到坐标栅格起点的步数来表示顶点的位置；
/*
    约定栅格坐标系一律用大写字母表示，X0 = x0 - ox / step;   
            X0为顶点ver0(x0, y0, z0)在栅格坐标系下的坐标的X分量；
            其中ox为栅格起点的x坐标，step为栅格步长；
 */


// 序列化的距离场数据：
/*
    距离场数据按照x优先y其次z最后的顺序存储，
        即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), .

*/


////////////////////////////////////////////////////////////////////////////////////////////// 表象转换、DEBUG 接口、辅助接口
namespace MY_DEBUG
{
    static std::string g_debugPath = "E:/";


    // 笛卡尔坐标系表示的向量转换为栅格坐标系：
    template <typename T, typename Tc, typename To>
    void CC2GC(Vec<3, T>& vecGC, const Vec<3, Tc>& vec, const Vec<3, To>& origin, const float step)
    {        
        const T ox = static_cast<T>(origin[0]);
        const T oy = static_cast<T>(origin[1]);
        const T oz = static_cast<T>(origin[2]);
        const T cx = static_cast<T>(vec[0]);
        const T cy = static_cast<T>(vec[1]);
        const T cz = static_cast<T>(vec[2]);
        const T stepT = static_cast<T>(step);
        vecGC[0] = (cx - ox) / stepT;
        vecGC[1] = (cy - oy) / stepT;
        vecGC[2] = (cz - oz) / stepT;
    }


    template <typename T>
    Eigen::Matrix<T, 1, 3> vec2Vec(const Vec<3, T>& vec)
    {
        return Eigen::Matrix<T, 1, 3>{vec[0], vec[1], vec[2]};
    }


    template <typename T>
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> vers2mat(const std::vector<Vec<3, T>>& vers)
    {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> versMat;
        versMat.resize(vers.size(), 3);
        for (unsigned i = 0; i < vers.size(); ++i)
        {
            versMat(i, 0) = vers[i][0];
            versMat(i, 1) = vers[i][1];
            versMat(i, 2) = vers[i][2];
        }
        return versMat;
    }

    template <unsigned int N, typename Derived, typename T>
    void mat2vers(std::vector<Vec<N, T>>& vers, const Eigen::PlainObjectBase<Derived>& versMat)
    {
        assert(N == versMat.cols(), "Assert!!! versMat and vers should have the same dimensions.");
        unsigned versCount = versMat.rows();
        vers.resize(versCount);
        for (unsigned i = 0; i < versCount; ++i)
            for (unsigned k = 0; k < N; ++k)
                vers[i][k] = static_cast<T>(versMat(i, k)); 
    }

    template <typename T>
    Eigen::MatrixXi tris2mat(const std::vector<Vec<3, T>>& tris)
    {
        Eigen::MatrixXi trisMat;
        trisMat.resize(tris.size(), 3);
        for (unsigned i = 0; i < tris.size(); ++i)
        {
            trisMat(i, 0) = static_cast<int>(tris[i][0]);
            trisMat(i, 1) = static_cast<int>(tris[i][1]);
            trisMat(i, 2) = static_cast<int>(tris[i][2]);
        }
        return trisMat;
    }


    template <unsigned int N, typename T, typename Derived>
    void eigenVec2sdfVec(Vec<N, T>& vec, const Eigen::PlainObjectBase<Derived>& eigenVec)
    {
        assert(1 == eigenVec.rows() || 1 == eigenVec.cols(), "assert!!! eigenVec should be a vector.");
        assert(N == eigenVec.rows() * eigenVec.cols(), "assert!!! vec and eigenVec should have the same size.");

        for (unsigned i = 0; i < N; ++i)
            vec[i] = static_cast<T>(eigenVec(i));
    }


    template <unsigned int N, typename T, typename Derived>
    void sdfVec2eigenVec(Eigen::PlainObjectBase<Derived>& eigenVec, const Vec<N, T>& vec)
    {
        using Scalar = typename Derived::Scalar;
        assert(1 == eigenVec.rows() || 1 == eigenVec.cols(), "assert!!! eigenVec should be a vector.");
        assert(N == eigenVec.rows() * eigenVec.cols(), "assert!!! vec and eigenVec should have the same size.");
        for (unsigned i = 0; i < N; ++i)
            eigenVec(i) = static_cast<Scalar>(vec[i]);
    }


    template <typename T, typename Derived>
    void eigenVec2sdfArray3(SDF_GEN::Array3<T, SDF_GEN::Array1<T>>& arry, const Eigen::PlainObjectBase<Derived>& eigenVec)
    {
        assert(1 == eigenVec.rows() || 1 == eigenVec.cols(), "assert!!! eigenVec should be a vector.");

        unsigned elemCount = static_cast<unsigned>(eigenVec.rows() * eigenVec.cols());
        assert(static_cast<unsigned>(arry.ni * arry.nj * arry.nk) == elemCount, "assert!!! vec and eigenVec should have the same size.");

        for (unsigned i = 0; i < elemCount; ++i)
            arry.a[i] = static_cast<T>(eigenVec(i));
    }


    template <typename T, typename Derived>
    void sdfArray32eigenVec(Eigen::PlainObjectBase<Derived>& eigenVec, const SDF_GEN::Array3<T, SDF_GEN::Array1<T>>& arry)
    {
        assert(1 == eigenVec.rows() || 1 == eigenVec.cols(), "assert!!! eigenVec should be a vector.");
        assert(arry.ni * arry.nj * arry.nk == eigenVec.rows() * eigenVec.cols(), "assert!!! vec and eigenVec should have the same size.");

        std::vector<T>& values = arry.a;
        eigenVec2Vec(values, eigenVec);
    }

    template <typename T, typename Derived>
    void sdfArray22eigenVec(Eigen::PlainObjectBase<Derived>& eigenVec, const SDF_GEN::Array2<T, SDF_GEN::Array1<T>>& arry)
    {
        assert(1 == eigenVec.rows() || 1 == eigenVec.cols(), "assert!!! eigenVec should be a vector.");
        assert(arry.ni * arry.nj  == eigenVec.rows() * eigenVec.cols(), "assert!!! vec and eigenVec should have the same size.");

        std::vector<T>& values = arry.a;
        eigenVec2Vec(values, eigenVec);
    }


    // 读取OBJ网格，生成计算距离场程序所用的网格数据向量：
    void readOBJ(std::vector<Vec3f>& vers, std::vector<Vec3ui>& tris, const char* path)
    {
        std::ifstream infile(path);
        if (!infile)
        {
            std::cerr << "Failed to open. Terminating.\n";
            exit(-1);
        }

        int ignored_lines = 0;
        std::string line;

        while (!infile.eof())
        {
            std::getline(infile, line);
            if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn")) //.obj files sometimes contain vertex normals indicated by "vn"
            {
                std::stringstream data(line);
                char c;
                Vec3f point;
                data >> c >> point[0] >> point[1] >> point[2];
                vers.push_back(point);
            }
            else if (line.substr(0, 1) == std::string("f"))
            {
                std::stringstream data(line);
                char c;
                int v0, v1, v2;
                data >> c >> v0 >> v1 >> v2;
                tris.push_back(Vec3ui(v0 - 1, v1 - 1, v2 - 1));
            }
            else if (line.substr(0, 2) == std::string("vn"))
            {
                std::cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
                exit(-2);
            }
            else
                ++ignored_lines;
        }
        infile.close();
        if (ignored_lines > 0)
            std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";
        std::cout << "Read in " << vers.size() << " vertices and " << tris.size() << " faces." << std::endl;
    }


    static void debugDisp()			// 递归终止
    {						//		递归终止设为无参或者一个参数的情形都可以。
        std::cout << std::endl;
        return;
    }


    template <typename T, typename... Types>
    static void debugDisp(const T& firstArg, const Types&... args)
    {
        std::cout << firstArg << " ";
        debugDisp(args...);
    }


    template <typename T, int M, int N>
    static void dispData(const Eigen::Matrix<T, M, N>& m)
    {
        auto dataPtr = m.data();
        unsigned elemsCount = m.size();

        for (unsigned i = 0; i < elemsCount; ++i)
            std::cout << dataPtr[i] << ", ";

        std::cout << std::endl;
    }


    template <typename Derived>
    static void dispData(const Eigen::PlainObjectBase<Derived>& m)
    {
        int m0 = m.RowsAtCompileTime;
        int n0 = m.ColsAtCompileTime;

        auto dataPtr = m.data();
        unsigned elemsCount = m.size();

        for (unsigned i = 0; i < elemsCount; ++i)
            std::cout << dataPtr[i] << ", ";

        std::cout << std::endl;
    }


    template <typename Derived>
    static void dispElem(const Eigen::MatrixBase<Derived>& m)
    {
        const Derived& mm = m.derived();
        std::cout << mm(1, 1) << std::endl;
    }


    template<typename DerivedV>
    static void debugWriteVers(const char* name, const Eigen::PlainObjectBase<DerivedV>& vers)
    {
        char path[512] = { 0 };
        sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
        objWriteVerticesMat(path, vers);
    }


    template<typename T>
    static void debugWriteMesh(const char* name, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& vers, const Eigen::MatrixXi& tris)
    {
        char path[512] = { 0 };
        sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
        objWriteMeshMat(path, vers, tris);
    }


    template<typename DerivedV>
    static void debugWriteEdges(const char* name, const Eigen::MatrixXi& edges, const Eigen::PlainObjectBase<DerivedV>& vers)
    {
        char path[512] = { 0 };
        sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), name);
        objWriteEdgesMat(path, edges, vers);
    }


    // 自定义计时器，使用WINDOWS计时API
    class tiktok
    {
    private:
        tiktok() = default;
        tiktok(const tiktok&) {}
        ~tiktok() = default;

    public:
        DWORD startTik;
        DWORD endTik;
        unsigned recordCount;
        std::vector<DWORD> records;

        static tiktok& getInstance()
        {
            static tiktok tt_instance;
            return tt_instance;
        }

        void start()
        {
            this->startTik = GetTickCount();
            this->recordCount = 0;
            this->records.clear();
        }

        void endCout(const char* str)
        {
            this->endTik = GetTickCount();
            std::cout << str << endTik - startTik << std::endl;
        }

        DWORD endGetCount()
        {
            this->endTik = GetTickCount();
            return endTik - startTik;
        }

        bool endWrite(const char* fileName, const char* str)
        {
            this->endTik = GetTickCount();
            std::ofstream file(fileName, std::ios_base::out | std::ios_base::app);
            if (!file)
                return false;

            file << str << endTik - startTik << std::endl;
            file.close();
            return true;
        }

        void takeArecord()
        {
            this->records.push_back(GetTickCount());
            recordCount++;
        }
    };
}
using namespace MY_DEBUG;


// 生成loopVersCount个点组成的环路的边数据
Eigen::MatrixXi getLoopEdges(const int loopVersCount)
{
    Eigen::MatrixXi edges(loopVersCount, 2);
    edges.col(0) = Eigen::VectorXi::LinSpaced(loopVersCount, 0, loopVersCount - 1);
    edges.block(0, 1, loopVersCount - 1, 1) = Eigen::VectorXi::LinSpaced(loopVersCount - 1, 1, loopVersCount - 1);
    edges(loopVersCount - 1, 1) = 0;

    return edges;
}


// 原make_level_set3()函数；
void calcSDF(const std::vector<Vec3ui>& tris, const std::vector<Vec3f>& vers,
        const Vec3f& startPos, const float step, const int ni, const int nj, const int nk,
        SDF_GEN::Array3f& SDFvalues, const int exact_band = 1)
{
    /*
    void make_level_set3(
            const std::vector<Vec3ui> &tris,                网格三角片
            const std::vector<Vec3f>    &vers,             网格点云
            const Vec3f &origin,                                坐标栅格原点
            float step,                                                 坐标栅格步长
            int ni, int nj, int nk,                                   坐标栅格三个维度的步数；
            SDF_GEN::Array3f &SDFvalues,                输出的距离场三维矩阵；
            const int exact_band                                默认值为1，求三角片包围盒时外扩的栅格数；
            )
    */

    //  起点坐标：
    const float ox = startPos[0];
    const float oy = startPos[1];
    const float oz = startPos[2]; 

    SDFvalues.resize(ni, nj, nk);
    SDFvalues.assign((ni + nj + nk) * step);                           // upper bound on distance
    SDF_GEN::Array3i closest_tris(ni, nj, nk, -1);                     // 每个采样点最近的三角片的索引，初始化为-1；

    // 判断距离符号时有用；intersection_count(i,j,k) is # of tris intersections in (i-1,i]vers{j}vers{k} 
    SDF_GEN::Array3i intersection_count(ni, nj, nk, 0);         //  每个栅格与网格三角片交点的数目；
                                                               
    // we begin by initializing distances near the mesh, and figuring out intersection counts

    // 1. 对三角片的遍历 
#ifdef USE_MULTITHREADS

    PARALLEL_FOR(0, tris.size(), [&](unsigned triIdx)
        {
            unsigned int vaIdx, vbIdx, vcIdx;                  // 当前三角片的三个顶点索引；
            assign(tris[triIdx], vaIdx, vbIdx, vcIdx);

            // 1f.1 用顶点到坐标栅格起点的步数来表示顶点的位置； 
            Vec3f verA = vers[vaIdx];
            Vec3f verB = vers[vbIdx];
            Vec3f verC = vers[vcIdx];
            Vec3d verA_gc, verB_gc, verC_gc;
            CC2GC(verA_gc, verA,  startPos, step);
            CC2GC(verB_gc, verB, startPos, step);
            CC2GC(verC_gc, verC, startPos, step);
            double Xa =verA_gc[0], Ya = verA_gc[1], Za = verA_gc[2]; 
            double Xb = verB_gc[0], Yb = verB_gc[1], Zb = verB_gc[2];
            double Xc = verC_gc[0], Yc= verC_gc[1], Zc = verC_gc[2];
    
            // 1f.2 求三角片的包围盒——栅格坐标系下包围盒起始点取整数，且扩张exact_band的尺寸;
            int Xmin = clamp(int(min(Xa, Xb, Xc)) - exact_band, 0, ni - 1), Xmax = clamp(int(max(Xa, Xb, Xc)) + exact_band + 1, 0, ni - 1);
            int Ymin = clamp(int(min(Ya, Yb, Yc)) - exact_band, 0, nj - 1), Ymax = clamp(int(max(Ya, Yb, Yc)) + exact_band + 1, 0, nj - 1);
            int Zmin = clamp(int(min(Za, Zb, Zc)) - exact_band, 0, nk - 1), Zmax = clamp(int(max(Za, Zb, Zc)) + exact_band + 1, 0, nk - 1);

            // 1f.3 计算三角片包围盒内的所有栅格采样点，到当前三角片的距离，找出最小距离，作为该采样点暂时的SDF值；
            for (int k = Zmin; k <= Zmax; ++k)
                for (int j = Ymin; j <= Ymax; ++j)
                    for (int i = Xmin; i <= Xmax; ++i)
                    {
                        Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);                                         // 当前栅格采样点
                        float dis = point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);  // 当前栅格采样点到当前三角片的距离
                        if (dis < SDFvalues(i, j, k))       // 若当前计算的距离比之前计算的小，更新距离数据；需要加互斥锁的写操作；
                        {
                            std::lock_guard<std::mutex> guard(g_mutex);
                            SDFvalues(i, j, k) = dis;
                            closest_tris(i, j, k) = triIdx;
                        }
                    }

            // 1f.4 相交检测：  
            Xmin = clamp((int)std::ceil(min(Xa, Xb, Xc)), 0, ni - 1);                // 当前三角片y最小值对应的y栅格下标；若超出栅格区域则取栅格边界；
            Xmax = clamp((int)std::floor(max(Xa, Xb, Xc)), 0, ni - 1);
            Ymin = clamp((int)std::ceil(min(Ya, Yb, Yc)), 0, nj - 1);
            Ymax = clamp((int)std::floor(max(Ya, Yb, Yc)), 0, nj - 1);

            //      遍历分析当前三角片在每个XOY平面上的投影；相当于做射线检测（Z方向的射线），确定射线和当前三角片的交点在哪个栅格；
            for (int Xp = Xmin; Xp <= Xmax; ++Xp)
                for (int Yp = Ymin; Yp <= Ymax; ++Yp)
                {
                    // 相当于做射线检测，当前射线为经过点(Xp, Yp, 0), 方向平行于Z方向的直线；
                    double alpha, beta, gamma;                       // 当前栅格点在投影三角片中的重心坐标；
                    if (point_in_triangle_2d(Xp, Yp, Xa, Ya, Xb, Yb, Xc, Yc, alpha, beta, gamma))
                    {
                        std::lock_guard<std::mutex> guard(g_mutex);
                        double Z_isct = alpha * Za + beta * Zb + gamma * Zc;      // 射线交点的z坐标（栅格坐标系）；
                        int Zp = int(std::ceil(Z_isct));                                                // intersection is in (Zp-1,Zp]

                        if (Zp < 0)
                            ++intersection_count(Xp, Yp, 0);    // we enlarge the first interval to include everything to the -vers direction
                        else if (Zp < nk)
                            ++intersection_count(Xp, Yp, Zp);  // we ignore intersections that are beyond the +vers side of the grid
                    }
                }
        });

#else

    for (unsigned int triIdx = 0; triIdx < tris.size(); ++triIdx)
    {
        unsigned int vaIdx, vbIdx, vcIdx;                  // 当前三角片的三个顶点索引；
        assign(tris[triIdx], vaIdx, vbIdx, vcIdx);

        // 1.1 用顶点到坐标栅格起点的步数来表示顶点的位置；
        double Xa = ((double)vers[vaIdx][0] - ox) / step, Ya = ((double)vers[vaIdx][1] - oy) / step, Za = ((double)vers[vaIdx][2] - oz) / step;
        double Xb = ((double)vers[vbIdx][0] - ox) / step, Yb = ((double)vers[vbIdx][1] - oy) / step, Zb = ((double)vers[vbIdx][2] - oz) / step;
        double Xc = ((double)vers[vcIdx][0] - ox) / step, Yc = ((double)vers[vcIdx][1] - oy) / step, Zc = ((double)vers[vcIdx][2] - oz) / step;

        // 求三角片的包围盒，用坐标栅格中的索引表示——(i0, j0, k0), (i1, j1, k1)确定的长方体空间；
        int i0 = clamp(int(min(Xa, Xb, Xc)) - exact_band, 0, ni - 1), i1 = clamp(int(max(Xa, Xb, Xc)) + exact_band + 1, 0, ni - 1);
        int j0 = clamp(int(min(Ya, Yb, Yc)) - exact_band, 0, nj - 1), j1 = clamp(int(max(Ya, Yb, Yc)) + exact_band + 1, 0, nj - 1);
        int k0 = clamp(int(min(Za, Zb, Zc)) - exact_band, 0, nk - 1), k1 = clamp(int(max(Za, Zb, Zc)) + exact_band + 1, 0, nk - 1);

        // 计算三角片包围盒内的所有坐标栅格点，到当前三角片的距离，找出最小距离；
        for (int k = k0; k <= k1; ++k) for (int j = j0; j <= j1; ++j) for (int i = i0; i <= i1; ++i)
        {
            Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);           // 当前栅格坐标点；
            float d = point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);
            if (d < SDFvalues(i, j, k))
            {
                SDFvalues(i, j, k) = d;                // 若当前计算的距离比之前计算的小，更新距离数据
                closest_tris(i, j, k) = triIdx;        // 
            }
        }

        // and do intersection counts
        // int i0 = clamp(int(min(Xa, Xb, Xc)) - exact_band, 0, ni - 1)
        j0 = clamp((int)std::ceil(min(Ya, Yb, Yc)), 0, nj - 1);                // 当前三角片y最小值对应的y栅格下标；
        j1 = clamp((int)std::floor(max(Ya, Yb, Yc)), 0, nj - 1);
        k0 = clamp((int)std::ceil(min(Za, Zb, Zc)), 0, nk - 1);            // 当前三角片z最小值对应的z栅格下标；
        k1 = clamp((int)std::floor(max(Za, Zb, Zc)), 0, nk - 1);

        for (int k = k0; k <= k1; ++k) for (int j = j0; j <= j1; ++j)
        {
            double x_bc, y_bc, z_bc;       // 
            if (point_in_triangle_2d(j, k, Ya, Za, Yb, Zb, Yc, Zc, x_bc, y_bc, z_bc))
            {
                double fi = x_bc * Xa + y_bc * Xb + z_bc * Xc;      // intersection i coordinate
                int i_interval = int(std::ceil(fi));               // intersection is in (i_interval-1,i_interval]

                if (i_interval < 0)
                    ++intersection_count(0, j, k);          // we enlarge the first interval to include everything to the -vers direction
                else if (i_interval < ni)
                    ++intersection_count(i_interval, j, k);  // we ignore intersections that are beyond the +vers side of the grid
            }
        }
    }
#endif

    // 2. and now we fill in the rest of the distances with fast sweeping；貌似是计算那些三角片包围盒外的采样点的距离场值；
    for (unsigned int i = 0; i < 2; ++i)
    {
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, +1, +1);
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, -1, -1);
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, +1, -1);
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, -1, +1);
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, -1, +1);
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, +1, -1);
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, +1, -1, -1);
        sweep(tris, vers, SDFvalues, closest_tris, startPos, step, -1, +1, +1);
    }

    // 3. 符号判断；then figure out signs (inside/outside) from intersection counts 
#ifdef USE_MULTITHREADS 
    PARALLEL_FOR(0, ni, [&](int i)
        {
            for (int j = 0; j < nj; ++j)
            {
                int total_count = 0;
                // 平行z轴方向上遍历一条直线上的栅格，通过交点数来确定距离场值的符号；
                for (int k = 0; k < nk; ++k)
                {
                    // 在起点k == 0上，采样点必然在网格外，随着采样点沿着z轴上升，累计交点+1则表示进入网格，再+1表示出网格，以此类推
                    total_count += intersection_count(i, j, k);
                    if (total_count % 2 == 1)                                    // 若该栅格内交点数为奇数，则该栅格在网格内；
                    {
                        std::lock_guard<std::mutex> guard(g_mutex);
                        SDFvalues(i, j, k) = -SDFvalues(i, j, k);            // we are inside the mesh
                    }
                }
            }
        });
#else
    for (int k = 0; k < nk; ++k) for (int j = 0; j < nj; ++j)
    {
        int total_count = 0;
        // 平行x轴方向上遍历一条直线上的栅格，通过交点数来确定距离场值的符号；
        for (int i = 0; i < ni; ++i)
        {
            total_count += intersection_count(i, j, k);
            if (total_count % 2 == 1)                            // if parity of intersections so far is odd,
                SDFvalues(i, j, k) = -SDFvalues(i, j, k);            // we are inside the mesh
        }
    }
#endif

}


void calcSDF2Dconst(const std::vector<Vec2f>& vers,
    const Vec2f& startPos, const float step, const int ni, const int nj,
    SDF_GEN::Array3f& SDFvalues, const int exact_band = 1) {}


// 写.sdf文本文件：
bool writeSDF(const char* fileName, const SDF_GEN::Array3f& SDFvalues, const Vec3f& origin, const float step)
{
    // 输出数据写入到SDF文件中：
    std::ofstream outfile(fileName);
    outfile << SDFvalues.ni << " " << SDFvalues.nj << " " << SDFvalues.nk << std::endl;        // 第一行：
    outfile << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;   // 第二行：
    outfile << step << std::endl;                                                 // 第三行：grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < SDFvalues.a.size(); ++i)                // 第三行之后：
        outfile << SDFvalues.a[i] << std::endl;
    outfile.close();

    return true;
}


template<typename T>
bool writeSDF2D(const char* fileName, const SDF_GEN::Array2<T, Array1<T> >& SDFvalues, \
    const Vec<2, T>& origin, const T step)
{
    // 输出数据写入到SDF文件中：
    std::ofstream outfile(fileName);
    outfile << SDFvalues.ni << " " << SDFvalues.nj << std::endl;        // 第一行：
    outfile << origin[0] << " " << origin[1] << std::endl;   // 第二行：
    outfile << step << std::endl;                                                 // 第三行：grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < SDFvalues.a.size(); ++i)                // 第三行之后：
        outfile << SDFvalues.a[i] << std::endl;
    outfile.close();

    return true;
}


// 解析.sdf文本文件：
double parseSDF(std::vector<int>& stepCounts, Eigen::RowVector3d& gridsOri, Eigen::VectorXd& SDF, const char* filePath)
{
    /*
        double parseSDF(												返回距离场的采样步长；
                std::vector<int>& stepCounts,					xyz三个方向上的步数
                Eigen::RowVector3d& gridsOri,					距离场空间的原点
                Eigen::VectorXd& SDF,								距离场数据
                const char* filePath										SDF文件目录
                )

    */
    double SDFstep = -1;
    stepCounts.resize(3);
    std::string readStr(1024, '\0');
    std::ifstream sdfFile(filePath);
    if (!sdfFile)
        return SDFstep;

    // 第一行：步数
    {
        std::string tmpStr;
        sdfFile.getline(&readStr[0], 1024);

        unsigned index = 0;
        for (const auto& ch : readStr)
        {
            if (ch >= '0' && ch <= '9' || ch == '.' || ch == '+' || ch == '-')
                tmpStr.push_back(ch);
            else
            {
                if (tmpStr.size() > 0)
                {
                    stepCounts[index] = std::stoi(tmpStr);
                    index++;
                    tmpStr.clear();
                }
            }
        }
    }

    // 第二行：栅格原点
    {
        std::string tmpStr;
        sdfFile.getline(&readStr[0], 1024);

        unsigned index = 0;
        for (const auto& ch : readStr)
        {
            if (ch >= '0' && ch <= '9' || ch == '.' || ch == '+' || ch == '-')
                tmpStr.push_back(ch);
            else
            {
                if (tmpStr.size() > 0)
                {
                    gridsOri(index) = std::stod(tmpStr);
                    index++;
                    tmpStr.clear();
                }
            }
        }
    }

    // 第三行：距离场空间步长：
    sdfFile.getline(&readStr[0], 1024);
    SDFstep = std::stod(readStr);

    // 第三行之后：距离场数据：
    unsigned dataCount = stepCounts[0] * stepCounts[1] * stepCounts[2];
    SDF.resize(dataCount);
    for (unsigned i = 0; i < dataCount; ++i)
    {
        sdfFile.getline(&readStr[0], 1024);
        SDF(i) = std::stod(readStr);
    }
    sdfFile.close();

    return SDFstep;
}



////////////////////////////////////////////////////////////////////////////////////////////// 测试函数：

// 原main函数：
int test0(int argc, char* argv[])
{
    if (argc != 4)
    {
        std::cout << "SDFGen - A utility for converting closed oriented triangle meshes into grid-based signed distance fields.\n";
        std::cout << "\nThe output file format is:";
        std::cout << "<ni> <nj> <nk>\n";
        std::cout << "<origin_x> <origin_y> <origin_z>\n";
        std::cout << "<dx>\n";
        std::cout << "<value_1> <value_2> <value_3> [...]\n\n";

        std::cout << "(ni,nj,nk) are the integer dimensions of the resulting distance field.\n";
        std::cout << "(origin_x,origin_y,origin_z) is the 3D position of the grid origin.\n";
        std::cout << "<dx> is the grid spacing.\n\n";
        std::cout << "<value_n> are the signed distance data values, in ascending order of i, then j, then k.\n";

        std::cout << "The output filename will match that of the input, with the OBJ suffix replaced with SDF.\n\n";

        std::cout << "Usage: SDFGen <filename> <dx> <interCounts>\n\n";
        std::cout << "Where:\n";
        std::cout << "\t<filename> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix \".obj\".\n";
        std::cout << "\t<dx> specifies the length of grid cell in the resulting distance field.\n";
        std::cout << "\t<interCounts> specifies the number of cells worth of interCounts between the object bound box and the boundary of the distance field grid. Minimum is 1.\n\n";

        exit(-1);
    }

    std::string filename(argv[1]);            // 参数1——网格文件地址
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    std::stringstream arg2(argv[2]);      // 参数2——栅格中单个方块(grid cell)的尺寸；
    float dx;
    arg2 >> dx;

    std::stringstream arg3(argv[3]);      // 参数3——设定的栅格距离场边界到网格包围盒的最小距离，用方块数量表示，最小为1；
    int interCounts;
    arg3 >> interCounts;

    if (interCounts < 1)
        interCounts = 1;

    // start with a massive inside out bound box.
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    std::cout << "Reading data.\n";

    // 1. 读取输入网格，生成网格包围盒:
    std::ifstream infile(argv[1]);
    if (!infile)
    {
        std::cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }

    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f> vers;                  // 输入网格顶点
    std::vector<Vec3ui> tris;                 // 输入网格三角片；

    while (!infile.eof())
    {
        std::getline(infile, line);
        if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn")) //.obj files sometimes contain vertex normals indicated by "vn"
        {
            std::stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            vers.push_back(point);
            update_minmax(point, min_box, max_box);
        }
        else if (line.substr(0, 1) == std::string("f"))
        {
            std::stringstream data(line);
            char c;
            int v0, v1, v2;
            data >> c >> v0 >> v1 >> v2;
            tris.push_back(Vec3ui(v0 - 1, v1 - 1, v2 - 1));
        }
        else if (line.substr(0, 2) == std::string("vn"))
        {
            std::cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
            exit(-2);
        }
        else
            ++ignored_lines;
    }
    infile.close();

    if (ignored_lines > 0)
        std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";
    std::cout << "Read in " << vers.size() << " vertices and " << tris.size() << " faces." << std::endl;

    // 2. Add interCounts around the box.
    Vec3f unit(1, 1, 1);
    min_box -= interCounts * dx * unit;
    max_box += interCounts * dx * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / dx);
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 3. 计算距离场：
    std::cout << "Computing signed distance field.\n";
    SDF_GEN::Array3f SDFvalues;
    make_level_set3(tris, vers, min_box, dx, sizes[0], sizes[1], sizes[2], SDFvalues);

    std::string outname;

#ifdef HAVE_VTK
    // If compiled with VTK, we can directly output a volumetric image format instead
    //Very hackily strip off file suffix.
    outname = filename.substr(0, filename.size() - 4) + std::string(".vti");
    std::cout << "Writing results to: " << outname << "\n";
    vtkSmartPointer<vtkImageData> output_volume = vtkSmartPointer<vtkImageData>::New();

    output_volume->SetDimensions(SDFvalues.ni, SDFvalues.nj, SDFvalues.nk);
    output_volume->SetOrigin(SDFvalues.ni * dx / 2, SDFvalues.nj * dx / 2, SDFvalues.nk * dx / 2);
    output_volume->SetSpacing(dx, dx, dx);

    vtkSmartPointer<vtkFloatArray> distance = vtkSmartPointer<vtkFloatArray>::New();

    distance->SetNumberOfTuples(SDFvalues.a.size());

    output_volume->GetPointData()->AddArray(distance);
    distance->SetName("Distance");

    for (unsigned int i = 0; i < SDFvalues.a.size(); ++i) {
        distance->SetValue(i, SDFvalues.a[i]);
    }

    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(outname.c_str());

#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(output_volume);
#else
    writer->SetInputData(output_volume);
#endif
    writer->Write();

#else

    // if VTK support is missing, default back to the original ascii file-dump. Very hackily strip off file suffix.
    outname = filename.substr(0, filename.size() - 4) + std::string(".sdf");
    std::cout << "Writing results to: " << outname << "\n";

    // 输出数据写入到SDF文件中：
    std::ofstream outfile(outname.c_str());
    outfile << SDFvalues.ni << " " << SDFvalues.nj << " " << SDFvalues.nk << std::endl;        // 第一行：
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // 第二行：
    outfile << dx << std::endl;                 // 第三行：dx——grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < SDFvalues.a.size(); ++i)         // 第三行之后：
        outfile << SDFvalues.a[i] << std::endl;
    outfile.close();
#endif

    std::cout << "Processing complete.\n";
    return 0;
}


// 手动读取文件计算距离场：
int test1(int argc, char* argv[])
{
    tiktok& tt = tiktok::getInstance();
    tt.start();

    // 参数1——网格文件地址
    std::string filePath("E:/材料/");
    std::string filename("twoTeeth.obj");
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    // 参数2—— 距离场采样步长(grid cell) 
    float step = 0.2;

    // 参数3——设定的采样空间边界到网格包围盒的最小距离，用步数表示，最小为1；
    int interCounts = 3;
    if (interCounts < 1)
        interCounts = 1;

    //      网格轴向包围盒，初始尺寸为无穷大；
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    // 1. 读取输入网格，生成网格包围盒:
    std::ifstream infile((filePath + filename).c_str());
    if (!infile)
    {
        std::cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }

    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f> vers;                  // 输入网格顶点
    std::vector<Vec3ui> tris;                 // 输入网格三角片；

    while (!infile.eof())
    {
        std::getline(infile, line);
        if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn")) //.obj files sometimes contain vertex normals indicated by "vn"
        {
            std::stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            vers.push_back(point);
            update_minmax(point, min_box, max_box);         // 更新网格包围盒数据；
        }
        else if (line.substr(0, 1) == std::string("f"))
        {
            std::stringstream data(line);
            char c;
            int v0, v1, v2;
            data >> c >> v0 >> v1 >> v2;
            tris.push_back(Vec3ui(v0 - 1, v1 - 1, v2 - 1));
        }
        else if (line.substr(0, 2) == std::string("vn"))
        {
            std::cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
            exit(-2);
        }
        else
            ++ignored_lines;
    }
    infile.close();

    if (ignored_lines > 0)
        std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";
    std::cout << "Read in " << vers.size() << " vertices and " << tris.size() << " faces." << std::endl;

    // 2. 生成采样点；
    Vec3f unit(1, 1, 1);
    //      包围盒增加空隙的尺寸，生成坐标栅格；
    min_box -= interCounts * step * unit;                       // 坐标栅格的起点坐标；
    max_box += interCounts * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step);        // 坐标栅格三个维度上的步数；
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 3. 计算距离场：
    SDF_GEN::Array3f SDFvalues;
    // make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], SDFvalues);
    calcSDF(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], SDFvalues);
    tt.endCout("读网格、计算距离场耗时：");

    // 写输出数据——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::string outname;
    outname = std::string{ "E:/" } + filename.substr(0, filename.size() - 4) + std::string(".sdf");
    std::ofstream outfile(outname.c_str());
    outfile << SDFvalues.ni << " " << SDFvalues.nj << " " << SDFvalues.nk << std::endl;        // 第一行：
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // 第二行：
    outfile << step << std::endl;                                                    // 第三行：step——grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < SDFvalues.a.size(); ++i)         // 第三行之后：
        outfile << SDFvalues.a[i] << std::endl;
    outfile.close();

    // 4. 生成距离场小于0的点云——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::list<Vec3f> versList;
    for (unsigned i = 0; i < SDFvalues.ni; ++i) 
        for (unsigned j = 0; j < SDFvalues.nj; ++j)
            for (unsigned k = 0; k < SDFvalues.nk; ++k)
                if (SDFvalues(i, j, k) < 0)
                    versList.push_back(Vec3f{ min_box[0] + i * step, min_box[1] + j * step, min_box[2] + k * step });

    // 5. 写输出数据
    Eigen::MatrixXf versOut;
    int versOutCount = versList.size();
    auto iter = versList.begin();
    versOut.resize(versOutCount, 3);
    for (int i = 0; i < versOutCount; ++i)
    {
        versOut(i, 0) = iter->v[0];
        versOut(i, 1) = iter->v[1];
        versOut(i, 2) = iter->v[2];
        iter++;
    }

    debugWriteVers("versOut", versOut);
     
    std::cout << "finished.\n";

    return 0;
}


// 尝试设定栅格原点和采样空间尺寸，生成网格SDF数据：
int test2(int argc, char* argv[])
{
    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f> vers;                  // 输入网格顶点
    std::vector<Vec3ui> tris;                 // 输入网格三角片；
    Eigen::MatrixXf versMat;
    Eigen::MatrixXi trisMat;

    // 参数1——网格文件地址
    std::string filePath("E:/");
    std::string filename("meshOuter.obj");
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    // 参数2—— 距离场采样步长(grid cell) 
    float step = 0.5;

    // 参数3——设定的采样空间边界到网格包围盒的最小距离，用步数表示，最小为1；
    int interCounts = 3;
    if (interCounts < 1)
        interCounts = 1;

    //      网格轴向包围盒，初始尺寸为无穷大；
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    // 1. 读取输入网格，生成网格包围盒:
    std::ifstream infile((filePath + filename).c_str());
    if (!infile)
    {
        std::cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }
    while (!infile.eof())
    {
        std::getline(infile, line);
        if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn")) //.obj files sometimes contain vertex normals indicated by "vn"
        {
            std::stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            vers.push_back(point);
            update_minmax(point, min_box, max_box);         // 更新网格包围盒数据；
        }
        else if (line.substr(0, 1) == std::string("f"))
        {
            std::stringstream data(line);
            char c;
            int v0, v1, v2;
            data >> c >> v0 >> v1 >> v2;
            tris.push_back(Vec3ui(v0 - 1, v1 - 1, v2 - 1));
        }
        else if (line.substr(0, 2) == std::string("vn"))
        {
            std::cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
            exit(-2);
        }
        else
            ++ignored_lines;
    }
    infile.close();
    if (ignored_lines > 0)
        std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";
    std::cout << "Read in " << vers.size() << " vertices and " << tris.size() << " faces." << std::endl;

    versMat = vers2mat(vers);
    trisMat = tris2mat(tris);
    objWriteMeshMat("E:/meshInput.obj", versMat, trisMat);


    // 2. 生成采样点；
    Vec3f unit(1, 1, 1);
    //      包围盒增加空隙的尺寸，生成坐标栅格；
    min_box -= interCounts * step * unit;                           // 坐标栅格的起点坐标；
    max_box += interCounts * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step);        // 坐标栅格三个维度上的步数；
    unsigned xCount = sizes[0];
    unsigned yCount = sizes[1];
    unsigned zCount = sizes[2];
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 2+: 添加约束——要求指定点位置必须有采样点，否则平移网格；
    Vec3f pos0(0.3, 0.4, 0.5);
    Vec3f arrow = pos0 - min_box;
    Vec3f posN = min_box + Vec3f{step * std::ceil(arrow[0]/step), step * std::ceil(arrow[1] / step) , step * std::ceil(arrow[2] / step) };

    Vec3f offsetArrow = pos0 - posN;                // 所有栅格点需要平移的向量；
    min_box = min_box + offsetArrow;            // 只需要平移min_box即可；

    // 3. 计算距离场：
    SDF_GEN::Array3f SDFvalues;
    make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], SDFvalues);

    std::string outname;
    outname = std::string{ "E:/" } + filename.substr(0, filename.size() - 4) + std::string(".sdf");

    // for debug——打印所有栅格点：
    Eigen::MatrixXf gridCenters;
    Eigen::RowVector3f origin = vec2Vec(min_box);

    genGrids(gridCenters, origin, step, std::vector<unsigned>{xCount, yCount, zCount});
    objWriteVerticesMat("E:/gridCenters.obj", gridCenters);

#if 1

    // 写输出数据——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::ofstream outfile(outname.c_str());
    outfile << SDFvalues.ni << " " << SDFvalues.nj << " " << SDFvalues.nk << std::endl;        // 第一行：
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // 第二行：
    outfile << step << std::endl;                                                    // 第三行：step——grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < SDFvalues.a.size(); ++i)            // 第三行之后：
        outfile << SDFvalues.a[i] << std::endl;
    outfile.close();

#endif
    std::cout << "finished.\n";

    return 0;
} 


// 读取网格生成距离场数据，进而生成网格内部点云；
void test3()
{
    // 1. 读取网格数据：
    tiktok& tt = tiktok::getInstance();
    tt.start();
    float step = 0.5;                    // 参数2—— 距离场采样步长(grid cell) 
    int interCounts = 3;              // 参数3——设定的采样空间边界到网格包围盒的最小距离，用步数表示，最小为1；
    std::vector<Vec3f> vers;
    std::vector<Vec3ui> tris;
    Eigen::MatrixXf versIn;
    Eigen::MatrixXi trisIn;
    DWORD tikCount = 0;
    readOBJ(vers, tris, "E:/材料/jawMesh.obj");
    versIn = vers2mat(vers);
    trisIn = tris2mat(tris);
    tikCount += tt.endGetCount();
    debugWriteMesh("meshInput", versIn, trisIn);

    // 2. 生成包围盒数据：
    tt.start();
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), \
        std::numeric_limits<float>::max()), \
        max_box(-std::numeric_limits<float>::max(), \
            - std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());            //      网格轴向包围盒，初始尺寸为无穷大；
    for (const auto& point : vers)
        update_minmax(point, min_box, max_box);         // 更新网格包围盒数据；

    // 2. 包围盒增加空隙的尺寸，生成坐标栅格；
    Vec3f unit(1, 1, 1);
    if (interCounts < 1)
        interCounts = 1;
    min_box -= interCounts * step * unit;                               // 坐标栅格的起点坐标；
    max_box += interCounts * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step);        // 坐标栅格三个维度上的步数；
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 3. 计算距离场：
    SDF_GEN::Array3f SDFvalues;
    calcSDF(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], SDFvalues);
    tikCount += tt.endGetCount();
    debugDisp("SDFGen：读网格+计算SDF总耗时：", tikCount);

    // 4. 生成距离场小于0的点云——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::list<Vec3f> versList;
    for (unsigned i = 0; i < SDFvalues.nk; ++i)
        for (unsigned j = 0; j < SDFvalues.nj; ++j)
            for (unsigned k = 0; k < SDFvalues.nk; ++k)
                if (SDFvalues(i, j, k) >= 0 && SDFvalues(i, j, k) <= 0.2)
                    versList.push_back(Vec3f{ min_box[0] + i * step, min_box[1] + j * step, min_box[2] + k * step });

    // 5. 写输出数据
    Eigen::MatrixXf versOut;
    int versOutCount = versList.size();
    auto iter = versList.begin();
    versOut.resize(versOutCount, 3);
    for (int i = 0; i < versOutCount; ++i)
    {
        versOut(i, 0) = iter->v[0];
        versOut(i, 1) = iter->v[1];
        versOut(i, 2) = iter->v[2];
        iter++;
    }

    debugWriteVers("versOut", versOut);
    debugDisp("finished.");
}


// testIO:
void test4() 
{
    std::vector<int> stepCounts;
    Eigen::RowVector3d gridsOri;
    Eigen::VectorXd SDF;
    double step = parseSDF(stepCounts, gridsOri, SDF, "E:/材料/tooth.sdf");
    const int xCount = stepCounts[0];
    const int yCount = stepCounts[1];
    const int zCount = stepCounts[2];

    Vec3f origin;
    SDF_GEN::Array3f SDFmat(xCount, yCount, zCount);
    eigenVec2sdfVec(origin, gridsOri);
    eigenVec2sdfArray3(SDFmat, SDF);

    writeSDF("E:/tooth.sdf", SDFmat, origin, static_cast<float>(step));

    debugDisp("finished.");
}


// test make_level_set2()
void test5()
{
    Eigen::MatrixXd versMat, versMat2D;
    std::vector<Vec2d> vers2D;
    objReadVerticesMat(versMat, "E:/材料/bdryUpper2D.obj");
    versMat2D = versMat.leftCols(2);
    mat2vers(vers2D, versMat2D);

    // 生成边数据：
    const unsigned versCount = vers2D.size();
    const unsigned edgesCount = versCount;
    std::vector<Vec2ui> edges(versCount);
    for (unsigned i = 0; i < edgesCount; ++i)
    {
        edges[i][0] = i;
        edges[i][1] = i + 1;
    }
    edges.rbegin()->operator[](1) = 0;

    // 生成采样点；
    double step = 0.5; 
    int interCounts = 3;  
    Vec2d min_box(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()),\
        max_box(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    Vec2d unit(1, 1); 
    for (const auto& ver : vers2D)
        update_minmax(ver, min_box, max_box);                    // 更新网格包围盒数据；

    min_box -= interCounts * step * unit;                               // 坐标栅格的起点坐标；
    max_box += interCounts * step * unit;
    Vec2ui sizes = Vec2ui((max_box - min_box) / step);        // 坐标栅格三个维度上的步数；
    unsigned xCount = sizes[0];
    unsigned yCount = sizes[1];   

    // 生成SDF:
    Array2d SDFvalues;
    make_level_set2(edges, vers2D, min_box, step, xCount, yCount, SDFvalues);

    // 输出：
    writeSDF2D("E:/tmpSDF2D.sdf", SDFvalues, min_box, step);

    debugDisp("finished.");
}


////////////////////////////////////////////////////////////////////////////////////////////// 生成控制台程序：
int cmdGenSDF3D(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cout << "SDFGen - A utility for converting closed oriented triangle meshes into grid-based signed distance fields.\n";
        std::cout << "\nThe output file format is:";
        std::cout << "<ni> <nj> <nk>\n";
        std::cout << "<origin_x> <origin_y> <origin_z>\n";
        std::cout << "<step>\n";
        std::cout << "<value_1> <value_2> <value_3> [...]\n\n";

        std::cout << "(ni,nj,nk) are the integer dimensions of the resulting distance field.\n";
        std::cout << "(origin_x,origin_y,origin_z) is the 3D position of the grid origin.\n";
        std::cout << "<step> is the grid spacing.\n\n";
        std::cout << "<value_n> are the signed distance data values, in ascending order of i, then j, then k.\n";

        std::cout << "Usage: SDFGen <sdfPath><objPath> <step><interCounts>\n\n";
        std::cout << "Where:\n";
        std::cout << "\t<objPath> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix \".obj\".\n";
        std::cout << "\t<step> specifies the length of grid cell in the resulting distance field.\n";
        std::cout << "\t<interCounts> specifies the number of cells worth of interCounts between the object bound box and the boundary of the distance field grid. Minimum is 1.\n\n";

        exit(-1);
    }

    std::string sdfPath(argv[1]);            // 参数2——网格文件地址
    if (sdfPath.size() < 5 || sdfPath.substr(sdfPath.size() - 4) != std::string(".sdf"))
    {
        std::cerr << "Error: Expected sdfPath of the form <name>.sdf.\n";
        exit(-1);
    }

    std::string objPath(argv[2]);            // 参数2——网格文件地址
    if (objPath.size() < 5 || objPath.substr(objPath.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with objPath of the form <name>.obj.\n";
        exit(-1);
    }

    std::stringstream arg2(argv[3]);      // 参数3——栅格中单个方块(grid cell)的尺寸；
    float step;
    arg2 >> step;

    std::stringstream arg3(argv[4]);      // 参数4——设定的栅格距离场边界到网格包围盒的最小距离，用方块数量表示，最小为1；
    int interCounts;
    arg3 >> interCounts;

    if (interCounts < 1)
        interCounts = 1;

    // start with a massive inside out bound box.
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    std::cout << "Reading data.\n";

    // 1. 读取输入网格，生成网格包围盒:
    std::ifstream infile(objPath.c_str());
    if (!infile)
    {
        std::cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }

    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f> vers;                  // 输入网格顶点
    std::vector<Vec3ui> tris;                 // 输入网格三角片；

    while (!infile.eof())
    {
        std::getline(infile, line);
        if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn")) //.obj files sometimes contain vertex normals indicated by "vn"
        {
            std::stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            vers.push_back(point);
            update_minmax(point, min_box, max_box);
        }
        else if (line.substr(0, 1) == std::string("f"))
        {
            std::stringstream data(line);
            char c;
            int v0, v1, v2;
            data >> c >> v0 >> v1 >> v2;
            tris.push_back(Vec3ui(v0 - 1, v1 - 1, v2 - 1));
        }
        else if (line.substr(0, 2) == std::string("vn"))
        {
            std::cerr << "Obj-loader is not able to parse vertex normals, please strip them from the input file. \n";
            exit(-2);
        }
        else
            ++ignored_lines;
    }
    infile.close();

    if (ignored_lines > 0)
        std::cout << "Warning: " << ignored_lines << " lines were ignored since they did not contain faces or vertices.\n";
    std::cout << "Read in " << vers.size() << " vertices and " << tris.size() << " faces." << std::endl;

    // 2. Add interCounts around the box.
    Vec3f unit(1, 1, 1);
    min_box -= interCounts * step * unit;
    max_box += interCounts * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step);
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 3. 计算距离场：
    std::cout << "Computing signed distance field.\n";
    SDF_GEN::Array3f SDFvalues;
    make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], SDFvalues);
     
    std::cout << "Writing results to: " << sdfPath << "\n";

    // 输出数据写入到SDF文件中：
    writeSDF(sdfPath.c_str(), SDFvalues, min_box, step);
 
    std::cout << "Processing complete.\n";
    return 0;
}


int cmdGenSDF2D(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cout << "SDFGen - A utility for converting closed oriented triangle meshes into grid-based signed distance fields.\n";
        std::cout << "Notice that the loop must be in the XOY plane.\n";
        std::cout << "\nThe output file format is:";
        std::cout << "<ni> <nj> <nk>\n";
        std::cout << "<origin_x> <origin_y> <origin_z>\n";
        std::cout << "<step>\n";
        std::cout << "<value_1> <value_2> <value_3> [...]\n\n";

        std::cout << "(ni,nj,nk) are the integer dimensions of the resulting distance field.\n";
        std::cout << "(origin_x,origin_y,origin_z) is the 3D position of the grid origin.\n";
        std::cout << "<step> is the grid spacing.\n\n";
        std::cout << "<value_n> are the signed distance data values, in ascending order of i, then j, then k.\n";

        std::cout << "Usage: SDFGen <sdfPath><objPath> <step><interCounts>\n\n";
        std::cout << "Where:\n";
        std::cout << "\t<objPath> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix \".obj\".\n";
        std::cout << "\t<step> specifies the length of grid cell in the resulting distance field.\n";
        std::cout << "\t<interCounts> specifies the number of cells worth of interCounts between the object bound box and the boundary of the distance field grid. Minimum is 1.\n\n";

        exit(-1);
    }

    std::string sdfPath(argv[1]);            // 参数1——输出SDF文件路径
    if (sdfPath.size() < 5 || sdfPath.substr(sdfPath.size() - 4) != std::string(".sdf"))
    {
        std::cerr << "Error: Expected sdfPath of the form <name>.sdf.\n";
        exit(-1);
    }

    std::string objPath(argv[2]);            // 参数2——网格文件路径
    if (objPath.size() < 5 || objPath.substr(objPath.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with objPath of the form <name>.obj.\n";
        exit(-1);
    }

    std::stringstream arg2(argv[3]);      // 参数3——栅格中单个方块(grid cell)的尺寸；
    double step;
    arg2 >> step;

    std::stringstream arg3(argv[4]);      // 参数4——设定的栅格距离场边界到网格包围盒的最小距离，用方块数量表示，最小为1；
    int interCounts;
    arg3 >> interCounts;

    if (interCounts < 1)
        interCounts = 1;
     
    // 1. 读取输入顶点，生成输入数据
    std::ifstream infile(objPath.c_str());
    if (!infile)
    {
        std::cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }
    int ignored_lines = 0;
    std::string line;
    std::vector<Vec2d> vers2D;                  // 输入网格顶点 
    while (!infile.eof())
    {
        std::getline(infile, line);
        if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn")) //.obj files sometimes contain vertex normals indicated by "vn"
        {
            std::stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            vers2D.push_back(Vec2d{point[0], point[1]});
        } 
        else
            ++ignored_lines;
    }
    infile.close(); 
     
    //      生成边数据：
    const unsigned versCount = vers2D.size();
    const unsigned edgesCount = versCount;
    std::vector<Vec2ui> edges(versCount);
    for (unsigned i = 0; i < edgesCount; ++i)
    {
        edges[i][0] = i;
        edges[i][1] = i + 1;
    }
    edges.rbegin()->operator[](1) = 0;

    //      生成采样点； 
    Vec2d min_box(std::numeric_limits<double>::max(), std::numeric_limits<double>::max()), \
        max_box(-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
    Vec2d unit(1, 1);
    for (const auto& ver : vers2D)
        update_minmax(ver, min_box, max_box);                    // 更新网格包围盒数据；

    min_box -= interCounts * step * unit;                               // 坐标栅格的起点坐标；
    max_box += interCounts * step * unit;
    Vec2ui sizes = Vec2ui((max_box - min_box) / step);        // 坐标栅格三个维度上的步数；
    unsigned xCount = sizes[0];
    unsigned yCount = sizes[1];

    // 3. 计算距离场：
    std::cout << "Computing signed distance field.\n";
    SDF_GEN::Array2d SDFvalues;
    make_level_set2(edges, vers2D, min_box, step, sizes[0], sizes[1], SDFvalues);
    std::cout << "Writing results to: " << sdfPath << "\n";

    // 输出数据写入到SDF文件中：
    writeSDF2D(sdfPath.c_str(), SDFvalues, min_box, step);

    std::cout << "Processing complete.\n";
    return 0;
}


int main(int argc, char* argv[]) 
{
    // return cmdGenSDF2D(argc, argv);

    test1(argc, argv);


    debugDisp("main() finished.");
 
    return 0;
}
 
 