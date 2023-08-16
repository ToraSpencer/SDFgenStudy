// ��ָ���ռ���������������ķ��ž��볡��
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
 

template <typename T>
Eigen::Matrix<T, 1, 3> vec2Vec(const Vec<3, T>& vec)
{
    return Eigen::Matrix<T, 1, 3>{vec[0], vec[1], vec[2]};
}
 
template <typename T>
Eigen::MatrixXd vers2mat(const std::vector<Vec<3, T>>& vers)
{
    Eigen::MatrixXd versMat;
    versMat.resize(vers.size(), 3);
    for (unsigned i = 0; i<vers.size(); ++i) 
    {
        versMat(i, 0) = static_cast<double>(vers[i][0]);
        versMat(i, 1) = static_cast<double>(vers[i][1]);
        versMat(i, 2) = static_cast<double>(vers[i][2]);
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


// ��ȡOBJ�������ɼ�����볡�������õ���������������
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


////////////////////////////////////////////////////////////////////////////////////////////// DEBUG �ӿ�
namespace MY_DEBUG
{
    static std::string g_debugPath = "E:/";


    static void debugDisp()			// �ݹ���ֹ
    {						//		�ݹ���ֹ��Ϊ�޲λ���һ�����������ζ����ԡ�
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


////////////////////////////////////////////////////////////////////////////////////////////// ���Ժ�����

// ԭmain������
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

    std::string filename(argv[1]);            // ����1���������ļ���ַ
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    std::stringstream arg2(argv[2]);      // ����2����դ���е�������(grid cell)�ĳߴ磻
    float dx;
    arg2 >> dx;

    std::stringstream arg3(argv[3]);      // ����3�����趨��դ����볡�߽絽�����Χ�е���С���룬�÷���������ʾ����СΪ1��
    int interCounts;
    arg3 >> interCounts;

    if (interCounts < 1)
        interCounts = 1;

    // start with a massive inside out bound box.
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    std::cout << "Reading data.\n";

    // 1. ��ȡ�����������������Χ��:
    std::ifstream infile(argv[1]);
    if (!infile)
    {
        std::cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }

    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f> vers;                  // �������񶥵�
    std::vector<Vec3ui> tris;                 // ������������Ƭ��

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

    // 3. ������볡��
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

    // �������д�뵽SDF�ļ��У�
    std::ofstream outfile(outname.c_str());
    outfile << DFvalues_grid.ni << " " << DFvalues_grid.nj << " " << DFvalues_grid.nk << std::endl;        // ��һ�У�
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // �ڶ��У�
    outfile << dx << std::endl;                 // �����У�dx����grad spacing����դ��ĳߴ磻
    for (unsigned int i = 0; i < DFvalues_grid.a.size(); ++i)         // ������֮��
        outfile << DFvalues_grid.a[i] << std::endl;
    outfile.close();
#endif

    std::cout << "Processing complete.\n";
    return 0;
}


// �ֶ���ȡ�ļ�������볡��
int test1(int argc, char* argv[])
{
    // ����1���������ļ���ַ
    std::string filePath("E:/����/");
    std::string filename("tooth.obj");
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    // ����2���� ���볡��������(grid cell) 
    float step = 0.5;

    // ����3�����趨�Ĳ����ռ�߽絽�����Χ�е���С���룬�ò�����ʾ����СΪ1��
    int interCounts = 3;
    if (interCounts < 1)
        interCounts = 1;

    //      ���������Χ�У���ʼ�ߴ�Ϊ�����
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    // 1. ��ȡ�����������������Χ��:
    std::ifstream infile((filePath + filename).c_str());
    if (!infile)
    {
        std::cerr << "Failed to open. Terminating.\n";
        exit(-1);
    }

    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f> vers;                  // �������񶥵�
    std::vector<Vec3ui> tris;                 // ������������Ƭ��

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
            update_minmax(point, min_box, max_box);         // ���������Χ�����ݣ�
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


    // 2. ���ɲ����㣻

    Vec3f unit(1, 1, 1);
    //      ��Χ�����ӿ�϶�ĳߴ磬��������դ��
    min_box -= interCounts * step * unit;                       // ����դ���������ꣻ
    max_box += interCounts * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step);        // ����դ������ά���ϵĲ�����
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 3. ������볡��
    SDF_GEN::Array3f DFvalues_grid;
    make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], DFvalues_grid);

    std::string outname;
    outname = std::string{"E:/"} + filename.substr(0, filename.size() - 4) + std::string(".sdf");

    // д������ݡ������볡���ݰ���x����y���z����˳��洢��������������������洢����f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::ofstream outfile(outname.c_str());
    outfile << DFvalues_grid.ni << " " << DFvalues_grid.nj << " " << DFvalues_grid.nk << std::endl;        // ��һ�У�
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // �ڶ��У�
    outfile << step << std::endl;                                                    // �����У�step����grad spacing����դ��ĳߴ磻
    for (unsigned int i = 0; i < DFvalues_grid.a.size(); ++i)         // ������֮��
        outfile << DFvalues_grid.a[i] << std::endl;
    outfile.close();

    std::cout << "finished.\n";

    return 0;
}


// �����趨դ��ԭ��Ͳ����ռ�ߴ磬��������SDF���ݣ�
int test2(int argc, char* argv[])
{
    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f> vers;                  // �������񶥵�
    std::vector<Vec3ui> tris;                 // ������������Ƭ��
    Eigen::MatrixXd versMat;
    Eigen::MatrixXi trisMat;

    // ����1���������ļ���ַ
    std::string filePath("E:/");
    std::string filename("meshOuter.obj");
    if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj"))
    {
        std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
        exit(-1);
    }

    // ����2���� ���볡��������(grid cell) 
    float step = 0.5;

    // ����3�����趨�Ĳ����ռ�߽絽�����Χ�е���С���룬�ò�����ʾ����СΪ1��
    int interCounts = 3;
    if (interCounts < 1)
        interCounts = 1;

    //      ���������Χ�У���ʼ�ߴ�Ϊ�����
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()),
        max_box(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    // 1. ��ȡ�����������������Χ��:
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
            update_minmax(point, min_box, max_box);         // ���������Χ�����ݣ�
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


    // 2. ���ɲ����㣻
    Vec3f unit(1, 1, 1);
    //      ��Χ�����ӿ�϶�ĳߴ磬��������դ��
    min_box -= interCounts * step * unit;                           // ����դ���������ꣻ
    max_box += interCounts * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step);        // ����դ������ά���ϵĲ�����
    unsigned xCount = sizes[0];
    unsigned yCount = sizes[1];
    unsigned zCount = sizes[2];
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 2+: ���Լ������Ҫ��ָ����λ�ñ����в����㣬����ƽ������
    Vec3f pos0(0.3, 0.4, 0.5);
    Vec3f arrow = pos0 - min_box;
    Vec3f posN = min_box + Vec3f{step * std::ceil(arrow[0]/step), step * std::ceil(arrow[1] / step) , step * std::ceil(arrow[2] / step) };
    objWriteVerticesMat("E:/min_box.obj", vec2Vec(min_box));
    objWriteVerticesMat("E:/pos0.obj", vec2Vec(pos0));
    objWriteVerticesMat("E:/posN.obj", vec2Vec(posN));

    Vec3f offsetArrow = pos0 - posN;                // ����դ�����Ҫƽ�Ƶ�������
    min_box = min_box + offsetArrow;            // ֻ��Ҫƽ��min_box���ɣ�

    // 3. ������볡��
    SDF_GEN::Array3f DFvalues_grid;
    make_level_set3(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], DFvalues_grid);

    std::string outname;
    outname = std::string{ "E:/" } + filename.substr(0, filename.size() - 4) + std::string(".sdf");

    // for debug������ӡ����դ��㣺
    Eigen::MatrixXf gridCenters;
    Eigen::RowVector3f origin = vec2Vec(min_box);
    genGrids(gridCenters, origin, step, std::vector<unsigned>{xCount, yCount, zCount});
    objWriteVerticesMat("E:/gridCenters.obj", gridCenters);

    // д������ݡ������볡���ݰ���x����y���z����˳��洢��������������������洢����f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::ofstream outfile(outname.c_str());
    outfile << DFvalues_grid.ni << " " << DFvalues_grid.nj << " " << DFvalues_grid.nk << std::endl;        // ��һ�У�
    outfile << min_box[0] << " " << min_box[1] << " " << min_box[2] << std::endl;   // �ڶ��У�
    outfile << step << std::endl;                                                    // �����У�step����grad spacing����դ��ĳߴ磻
    for (unsigned int i = 0; i < DFvalues_grid.a.size(); ++i)            // ������֮��
        outfile << DFvalues_grid.a[i] << std::endl;
    outfile.close();

    std::cout << "finished.\n";

    return 0;
}
 

void foo(const std::vector<Vec3ui>& tris, const std::vector<Vec3f>& vers,
    const Vec3f& startPos, float step, int ni, int nj, int nk,
    SDF_GEN::Array3f& DFvalues, const int exact_band = 1)
{
    /*
    void make_level_set3(
            const std::vector<Vec3ui> &tris,            ��������Ƭ
            const std::vector<Vec3f> &vers,             ����㼯
            const Vec3f &origin,                                ����դ��ԭ��
            float step,                                                 ����դ�񲽳�
            int ni, int nj, int nk,                                   ����դ������ά�ȵĲ�����
            SDF_GEN::Array3f &DFvalues,                                  ����ľ��볡��ά����
            const int exact_band                                Ĭ��ֵΪ1��������Ƭ��Χ��ʱ������դ������
            )


    */
    //  ������꣺
    const float ox = startPos[0];
    const float oy = startPos[1];
    const float oz = startPos[2];

    DFvalues.resize(ni, nj, nk);
    DFvalues.assign((ni + nj + nk) * step);                                    // upper bound on distance
    SDF_GEN::Array3i closest_tris(ni, nj, nk, -1);                     // ÿ������դ������������Ƭ����������ʼ��Ϊ-1��
    SDF_GEN::Array3i intersection_count(ni, nj, nk, 0);       // �жϾ������ʱ���ã�intersection_count(i,j,k) is # of tris intersections in (i-1,i]vers{j}vers{k}

    // we begin by initializing distances near the mesh, and figuring out intersection counts

    // 1. ������Ƭ�ı���
    tiktok& tt = tiktok::getInstance();
    tt.start();

    PARALLEL_FOR(0, tris.size(), [&](unsigned triIdx)
        {
            unsigned int vaIdx, vbIdx, vcIdx;                  // ��ǰ����Ƭ����������������
            assign(tris[triIdx], vaIdx, vbIdx, vcIdx);

            // 1.1 �ö��㵽����դ�����Ĳ�������ʾ�����λ�ã�
            double vax_sc = ((double)vers[vaIdx][0] - ox) / step, vay_sc = ((double)vers[vaIdx][1] - oy) / step, vaz_sc = ((double)vers[vaIdx][2] - oz) / step;
            double vbx_sc = ((double)vers[vbIdx][0] - ox) / step, vby_sc = ((double)vers[vbIdx][1] - oy) / step, vbz_sc = ((double)vers[vbIdx][2] - oz) / step;
            double vcx_sc = ((double)vers[vcIdx][0] - ox) / step, vcy_sc = ((double)vers[vcIdx][1] - oy) / step, vcz_sc = ((double)vers[vcIdx][2] - oz) / step;

            // ������Ƭ�İ�Χ�У�������դ���е�������ʾ����(i0, j0, k0), (i1, j1, k1)ȷ���ĳ�����ռ䣻
            int i0 = clamp(int(min(vax_sc, vbx_sc, vcx_sc)) - exact_band, 0, ni - 1), i1 = clamp(int(max(vax_sc, vbx_sc, vcx_sc)) + exact_band + 1, 0, ni - 1);
            int j0 = clamp(int(min(vay_sc, vby_sc, vcy_sc)) - exact_band, 0, nj - 1), j1 = clamp(int(max(vay_sc, vby_sc, vcy_sc)) + exact_band + 1, 0, nj - 1);
            int k0 = clamp(int(min(vaz_sc, vbz_sc, vcz_sc)) - exact_band, 0, nk - 1), k1 = clamp(int(max(vaz_sc, vbz_sc, vcz_sc)) + exact_band + 1, 0, nk - 1);

            // ��������Ƭ��Χ���ڵ���������դ��㣬����ǰ����Ƭ�ľ��룬�ҳ���С���룻
            for (int k = k0; k <= k1; ++k)
                for (int j = j0; j <= j1; ++j)
                    for (int i = i0; i <= i1; ++i)
                    {
                        Vec3f ver0(i * step + ox, j * step + oy, k * step + oz);           // ��ǰդ������㣻
                        float d = point_triangle_distance(ver0, vers[vaIdx], vers[vbIdx], vers[vcIdx]);
                        if (d < DFvalues(i, j, k))
                        {
                            std::lock_guard<std::mutex> guard(g_mutex);
                            DFvalues(i, j, k) = d;                // ����ǰ����ľ����֮ǰ�����С�����¾�������
                            closest_tris(i, j, k) = triIdx;        // 
                        }
                    }

            // and do intersection counts 
            j0 = clamp((int)std::ceil(min(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);                // ��ǰ����Ƭy��Сֵ��Ӧ��yդ���±ꣻ
            j1 = clamp((int)std::floor(max(vay_sc, vby_sc, vcy_sc)), 0, nj - 1);
            k0 = clamp((int)std::ceil(min(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);            // ��ǰ����Ƭz��Сֵ��Ӧ��zդ���±ꣻ
            k1 = clamp((int)std::floor(max(vaz_sc, vbz_sc, vcz_sc)), 0, nk - 1);

            for (int k = k0; k <= k1; ++k)
                for (int j = j0; j <= j1; ++j)
                {
                    double x_bc, y_bc, z_bc;       // 
                    if (point_in_triangle_2d(j, k, vay_sc, vaz_sc, vby_sc, vbz_sc, vcy_sc, vcz_sc, x_bc, y_bc, z_bc))
                    {
                        std::lock_guard<std::mutex> guard(g_mutex);
                        double fi = x_bc * vax_sc + y_bc * vbx_sc + z_bc * vcx_sc;      // intersection i coordinate
                        int i_interval = int(std::ceil(fi));               // intersection is in (i_interval-1,i_interval]

                        if (i_interval < 0)
                            ++intersection_count(0, j, k);          // we enlarge the first interval to include everything to the -vers direction
                        else if (i_interval < ni)
                            ++intersection_count(i_interval, j, k);  // we ignore intersections that are beyond the +vers side of the grid
                    }
                }
        });

    tt.endCout("elapsed time of step1: ");

    // 2. and now we fill in the rest of the distances with fast sweeping
    tt.start();
    for (unsigned int pass = 0; pass < 2; ++pass)
    {
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, +1, +1);
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, -1, -1);
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, +1, -1);
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, -1, +1);
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, -1, +1);
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, +1, -1);
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, +1, -1, -1);
        sweep(tris, vers, DFvalues, closest_tris, startPos, step, -1, +1, +1);
    }
    tt.endCout("elapsed time of step2: ");

    // 3. �����жϣ�then figure out signs (inside/outside) from intersection counts
    tt.start();
    PARALLEL_FOR(0, nk, [&](int k)
        {
            for (int j = 0; j < nj; ++j)
            {
                int total_count = 0;
                // ƽ��x�᷽���ϱ���һ��ֱ���ϵ�դ��ͨ����������ȷ�����볡ֵ�ķ��ţ�
                for (int i = 0; i < ni; ++i)
                {
                    total_count += intersection_count(i, j, k);
                    if (total_count % 2 == 1)                            // if parity of intersections so far is odd,
                    {
                        std::lock_guard<std::mutex> guard(g_mutex);
                        DFvalues(i, j, k) = -DFvalues(i, j, k);            // we are inside the mesh
                    }
                }
            }
        });
    tt.endCout("elapsed time of step3: ");
}


// ��ȡ�������ɾ��볡���ݣ��������������ڲ����ƣ�
void test3()
{
    // 1. ��ȡ�������ݣ�
    float step = 0.5;                    // ����2���� ���볡��������(grid cell) 
    int interCounts = 3;              // ����3�����趨�Ĳ����ռ�߽絽�����Χ�е���С���룬�ò�����ʾ����СΪ1��
    std::vector<Vec3f> vers;
    std::vector<Vec3ui> tris;
    readOBJ(vers, tris, "E:/����/tooth.obj");

    // 2. ���ɰ�Χ�����ݣ�
    Vec3f min_box(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), \
        std::numeric_limits<float>::max()), \
        max_box(-std::numeric_limits<float>::max(), \
            - std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());            //      ���������Χ�У���ʼ�ߴ�Ϊ�����
    for (const auto& point : vers)
        update_minmax(point, min_box, max_box);         // ���������Χ�����ݣ�

    // 2. ��Χ�����ӿ�϶�ĳߴ磬��������դ��
    Vec3f unit(1, 1, 1);
    if (interCounts < 1)
        interCounts = 1;
    min_box -= interCounts * step * unit;                               // ����դ���������ꣻ
    max_box += interCounts * step * unit;
    Vec3ui sizes = Vec3ui((max_box - min_box) / step);        // ����դ������ά���ϵĲ�����
    std::cout << "Bound box size: (" << min_box << ") to (" << max_box << ") with dimensions " << sizes << "." << std::endl;

    // 3. ������볡��
    SDF_GEN::Array3f DFvalues_grid;
    foo(tris, vers, min_box, step, sizes[0], sizes[1], sizes[2], DFvalues_grid);

    // 4. ���ɾ��볡С��0�ĵ��ơ������볡���ݰ���x����y���z����˳��洢��������������������洢����f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
    std::list<Vec3f> versList;
    for (unsigned i = 0; i < DFvalues_grid.nk; ++i)
        for (unsigned j = 0; j < DFvalues_grid.nj; ++j)
            for (unsigned k = 0; k < DFvalues_grid.nk; ++k)
                if (DFvalues_grid(i, j, k) <= 0)
                    versList.push_back(Vec3f{ min_box[0] + i * step, min_box[1] + j * step, min_box[2] + k * step });

    // 5. д�������
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
    
    //Eigen::MatrixXd vers;
    //objReadVerticesMat(vers, "E:/����/loop2D.obj");
    //Eigen::MatrixXi edges = getLoopEdges(vers.rows());
    //debugWriteEdges("loopEdges", edges, vers);

    debugDisp("main() finished.");
 
    return 0;
}

 