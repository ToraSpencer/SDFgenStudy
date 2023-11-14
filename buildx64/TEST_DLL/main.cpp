#include "SDFGEN_DLL.h"
#include "triMesh.h"
#pragma comment(lib, "SDFGEN_DLL.lib")

#include <iostream>
#include <fstream>
#include <list>
#include <string>
  
// 已禁用vcpkg

  
/////////////////////////////////////////////////////////////////////////////////////////////////////////// DEBUG接口；
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

	static void debugWriteVers(const char* fileName, const std::vector<verF>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, vers);
	}

	static void debugWriteVers(const char* fileName, const std::vector<verD>& vers)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, vers);
	}

	static void debugWriteMesh(const char* fileName, const triMeshD& triMesh)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, triMesh);
	}

	static void debugWriteMesh(const char* fileName, const triMeshF& triMesh)
	{
		char path[512] = { 0 };
		sprintf_s(path, "%s%s.obj", g_debugPath.c_str(), fileName);
		writeOBJ(path, triMesh);
	}
}
using namespace MY_DEBUG;


////////////////////////////////////////////////////////////////////////////////////////////// 测试函数
namespace TEST_SDFGEN_3D 
{
	void test1()
    { 
        triMeshF mesh1;
		triMeshD mesh2;
        readOBJ(mesh1, "E:/材料/tooth.obj");
		readOBJ(mesh2, "E:/材料/twoTeeth.obj");
		debugWriteMesh("meshInput1", mesh1);
		debugWriteMesh("meshInput2", mesh2);

        // 1. 
        SDF_RESULT sdfData;
        genSDF3D(sdfData, mesh1, 0.5, 3);

        // 2. 生成距离场小于0的点云——距离场数据按照x优先y其次z最后的顺序存储，即按数组索引增大方向存储的是f(0, 0, 0), f(1, 0, 0), f(2, 0, 0)... f(0, 1, 0), f(1, 1, 0) ..f(0, 0, 1), f(1, 0, 1), ..
		{
			std::vector<verF> versOut;
			std::list<verF> versList;
			unsigned ni = sdfData.stepsCount[0];
			unsigned nj = sdfData.stepsCount[1];
			unsigned nk = sdfData.stepsCount[2];
			for (unsigned i = 0; i < ni; ++i)
				for (unsigned j = 0; j < nj; ++j)
					for (unsigned k = 0; k < nk; ++k)
						if (sdfData.SDFvalues[i + ni * (j + nj * k)] < 0)
							versList.push_back(verF{ sdfData.origin.x + i * sdfData.step, sdfData.origin.y + \
								 j * sdfData.step, sdfData.origin.z + k * sdfData.step });
			versOut.insert(versOut.end(), versList.begin(), versList.end());
			debugWriteVers("versOut1", versOut);
		}


		// 3. 
		genSDF3D(sdfData, mesh2, 0.5);
		{
			std::vector<verF> versOut;
			std::list<verF> versList;
			unsigned ni = sdfData.stepsCount[0];
			unsigned nj = sdfData.stepsCount[1];
			unsigned nk = sdfData.stepsCount[2];
			for (unsigned i = 0; i < ni; ++i)
				for (unsigned j = 0; j < nj; ++j)
					for (unsigned k = 0; k < nk; ++k)
						if (sdfData.SDFvalues[i + ni * (j + nj * k)] < 0)
							versList.push_back(verF{ sdfData.origin.x + i * sdfData.step, sdfData.origin.y + \
								 j * sdfData.step, sdfData.origin.z + k * sdfData.step });
			versOut.insert(versOut.end(), versList.begin(), versList.end());
			debugWriteVers("versOut2", versOut);
		}
		  
		debugDisp("finished.");
    }

}


int main(int argc, char** argv) 
{
	TEST_SDFGEN_3D::test1();

	return 0;
}