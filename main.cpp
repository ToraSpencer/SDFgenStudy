// 在指定空间内生成三角网格的符号距离场；
/*
    SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
    Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
    ...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
    This code is public domain. Feel free to mess with it, let me know if you like it.
*/
#include "makelevelset3.h"
#include "config.h"
#include "myEigen.h"

#include <list>
#ifdef HAVE_VTK
  #include <vtkImageData.h>
  #include <vtkFloatArray.h>
  #include <vtkXMLImageDataWriter.h>
  #include <vtkPointData.h>
  #include <vtkSmartPointer.h>
#endif


// 栅格坐标系： 用顶点到坐标栅格起点的步数来表示顶点的位置；
/*
    约定栅格坐标系一律用大写字母表示，X0 = x0 - ox / step;   
            X0为顶点ver0(x0, y0, z0)在栅格坐标系下的坐标的X分量；
            其中ox为栅格起点的x坐标，step为栅格步长；
 */


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
    for (unsigned i = 0; i<vers.size(); ++i) 
    {
        versMat(i, 0) = vers[i][0];
        versMat(i, 1) = vers[i][1];
        versMat(i, 2) = vers[i][2];
    }
    return versMat;
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


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG 接口
namespace MY_DEBUG
{
    static std::string g_debugPath = "E:/";


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
}
using namespace MY_DEBUG;


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
    SDF_GEN::Array3f DFvalues_grid;
    make_level_set3(tris, vers, min_box, dx, sizes[0], sizes[1], sizes[2], DFvalues_grid);

    std::string outname;

#ifdef HAVE_VTK
    // If compiled with VTK, we can directly output a volumetric image format instead
    //Very hackily strip off file suffix.
    outname = filename.substr(0, filename.size() - 4) + std::string(".vti");
    std::cout << "Writing results to: " << outname << "\n";
    vtkSmartPointer<vtkImageData> output_volume = vtkSmartPointer<vtkImageData>::New();

    output_volume->SetDimensions(DFvalues_grid.ni, DFvalues_grid.nj, DFvalues_grid.nk);
    output_volume->SetOrigin(DFvalues_grid.ni * dx / 2, DFvalues_grid.nj * dx / 2, DFvalues_grid.nk * dx / 2);
    output_volume->SetSpacing(dx, dx, dx);

    vtkSmartPointer<vtkFloatArray> distance = vtkSmartPointer<vtkFloatArray>::New();

    distance->SetNumberOfTuples(DFvalues_grid.a.size());

    output_volume->GetPointData()->AddArray(distance);
    distance->SetName("Distance");

    for (unsigned int i = 0; i < DFvalues_grid.a.size(); ++i) {
        distance->SetValue(i, DFvalues_grid.a[i]);
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
    outfile << DFvalues_grid.ni << " " << DFvalues_grid.nj << " " << DFvalues_grid.nk << std::endl;        // 第一行：
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // 第二行：
    outfile << dx << std::endl;                 // 第三行：dx——grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < DFvalues_grid.a.size(); ++i)         // 第三行之后：
        outfile << DFvalues_grid.a[i] << std::endl;
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
    std::string filename("jawMeshDense.obj");
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    // 参数2—— 距离场采样步长(grid cell) 
    float step = 0.25;

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
    SDF_GEN::Array3f DFvalues_grid;
    make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], DFvalues_grid);
    tt.endCout("读网格、计算距离场耗时：");

    // 写输出数据——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::string outname;
    outname = std::string{ "E:/" } + filename.substr(0, filename.size() - 4) + std::string(".sdf");
    std::ofstream outfile(outname.c_str());
    outfile << DFvalues_grid.ni << " " << DFvalues_grid.nj << " " << DFvalues_grid.nk << std::endl;        // 第一行：
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // 第二行：
    outfile << step << std::endl;                                                    // 第三行：step——grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < DFvalues_grid.a.size(); ++i)         // 第三行之后：
        outfile << DFvalues_grid.a[i] << std::endl;
    outfile.close();

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
    objWriteVerticesMat("E:/min_box.obj", vec2Vec(min_box));
    objWriteVerticesMat("E:/pos0.obj", vec2Vec(pos0));
    objWriteVerticesMat("E:/posN.obj", vec2Vec(posN));

    Vec3f offsetArrow = pos0 - posN;                // 所有栅格点需要平移的向量；
    min_box = min_box + offsetArrow;            // 只需要平移min_box即可；

    // 3. 计算距离场：
    SDF_GEN::Array3f DFvalues_grid;
    make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], DFvalues_grid);

    std::string outname;
    outname = std::string{ "E:/" } + filename.substr(0, filename.size() - 4) + std::string(".sdf");

    // for debug——打印所有栅格点：
    Eigen::MatrixXf gridCenters;
    Eigen::RowVector3f origin = vec2Vec(min_box);
    genGrids(gridCenters, origin, step, std::vector<unsigned>{xCount, yCount, zCount});
    objWriteVerticesMat("E:/gridCenters.obj", gridCenters);

    // 写输出数据——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::ofstream outfile(outname.c_str());
    outfile << DFvalues_grid.ni << " " << DFvalues_grid.nj << " " << DFvalues_grid.nk << std::endl;        // 第一行：
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // 第二行：
    outfile << step << std::endl;                                                    // 第三行：step——grad spacing，即栅格的尺寸；
    for (unsigned int i = 0; i < DFvalues_grid.a.size(); ++i)            // 第三行之后：
        outfile << DFvalues_grid.a[i] << std::endl;
    outfile.close();

    std::cout << "finished.\n";

    return 0;
}
 

// 原make_level_set3()函数；
void calcSDF(const std::vector<Vec3ui>& tris, const std::vector<Vec3f>& vers,
    const Vec3f& startPos, float step, int ni, int nj, int nk,
    SDF_GEN::Array3f& SDFvalues, const int exact_band = 1)
{
    /*
    void make_level_set3(
            const std::vector<Vec3ui> &tris,            网格三角片
            const std::vector<Vec3f> &vers,             网格点集
            const Vec3f &origin,                                坐标栅格原点
            float step,                                                 坐标栅格步长
            int ni, int nj, int nk,                                   坐标栅格三个维度的步数；
            SDF_GEN::Array3f &SDFvalues,                                  输出的距离场三维矩阵；
            const int exact_band                                默认值为1，求三角片包围盒时外扩的栅格数；
            )


    */
    //  起点坐标：
    const float ox = startPos[0];
    const float oy = startPos[1];
    const float oz = startPos[2];

    SDFvalues.resize(ni, nj, nk);
    SDFvalues.assign((ni + nj + nk) * step);                           // upper bound on distance
    SDF_GEN::Array3i closest_tris(ni, nj, nk, -1);                     // 每个坐标栅格点最近的三角片的索引，初始化为-1；
    SDF_GEN::Array3i intersection_count(ni, nj, nk, 0);         //  每个栅格与网格三角片交点的数目；
                                                                              // 判断距离符号时有用；intersection_count(i,j,k) is # of tris intersections in (i-1,i]vers{j}vers{k} 

    // we begin by initializing distances near the mesh, and figuring out intersection counts

    // 1. 对三角片的遍历
    tiktok& tt = tiktok::getInstance(); 
    PARALLEL_FOR(0, tris.size(), [&](unsigned triIdx)
        {
            unsigned int vaIdx, vbIdx, vcIdx;                  // 当前三角片的三个顶点索引；
            assign(tris[triIdx], vaIdx, vbIdx, vcIdx);

            // 1.1 用顶点到坐标栅格起点的步数来表示顶点的位置；
            double Xa = ((double)vers[vaIdx][0] - ox) / step, Ya = ((double)vers[vaIdx][1] - oy) / step, Za = ((double)vers[vaIdx][2] - oz) / step;
            double Xb = ((double)vers[vbIdx][0] - ox) / step, Yb = ((double)vers[vbIdx][1] - oy) / step, Zb = ((double)vers[vbIdx][2] - oz) / step;
            double Xc = ((double)vers[vcIdx][0] - ox) / step, Yc = ((double)vers[vcIdx][1] - oy) / step, Zc = ((double)vers[vcIdx][2] - oz) / step;

            // 1.2 求三角片的包围盒——(Xmin, Ymin, Zmin), (Xmax, Ymax, Zmax)确定的长方体空间；
            int Xmin = clamp(int(min(Xa, Xb, Xc)) - exact_band, 0, ni - 1), Xmax = clamp(int(max(Xa, Xb, Xc)) + exact_band + 1, 0, ni - 1);
            int Ymin = clamp(int(min(Ya, Yb, Yc)) - exact_band, 0, nj - 1), Ymax = clamp(int(max(Ya, Yb, Yc)) + exact_band + 1, 0, nj - 1);
            int Zmin = clamp(int(min(Za, Zb, Zc)) - exact_band, 0, nk - 1), Zmax = clamp(int(max(Za, Zb, Zc)) + exact_band + 1, 0, nk - 1);

            // 1.3 计算三角片包围盒内的所有栅格采样点，到当前三角片的距离，找出最小距离，作为该采样点暂时的SDF值；
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

            // 1.4 相交检测：  
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

    // 2. 看不懂，好像没什么屌用；and now we fill in the rest of the distances with fast sweeping
#if 0
    tt.start();
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
    tt.endCout("elapsed time of step2: ");
#endif

    // 3. 符号判断；then figure out signs (inside/outside) from intersection counts 
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
    SDF_GEN::Array3f DFvalues_grid;
    calcSDF(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], DFvalues_grid);
    tikCount += tt.endGetCount();
    debugDisp("SDFGen：读网格+计算SDF总耗时：", tikCount);

    // 4. 生成距离场小于0的点云——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::list<Vec3f> versList;
    for (unsigned i = 0; i < DFvalues_grid.nk; ++i)
        for (unsigned j = 0; j < DFvalues_grid.nj; ++j)
            for (unsigned k = 0; k < DFvalues_grid.nk; ++k)
                if (DFvalues_grid(i, j, k) >= 0 && DFvalues_grid(i, j, k) <= 0.2)
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


Eigen::MatrixXi getLoopEdges(const int loopVersCount)
{
    Eigen::MatrixXi edges(loopVersCount, 2);
    edges.col(0) = Eigen::VectorXi::LinSpaced(loopVersCount, 0, loopVersCount - 1);
    edges.block(0, 1, loopVersCount - 1, 1) = Eigen::VectorXi::LinSpaced(loopVersCount - 1, 1, loopVersCount - 1);
    edges(loopVersCount - 1, 1) = 0;

    return edges;
}


int main(int argc, char* argv[]) 
{
    // test1(argc, argv);
    test3();

    debugDisp("main() finished.");
 
    return 0;
}

 