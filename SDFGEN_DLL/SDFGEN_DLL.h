#pragma once

#define WIN32_LEAN_AND_MEAN             // 从 Windows 头文件中排除极少使用的内容
#include <SDKDDKVer.h>	
#include <windows.h>	
#include <vector>
 
#include "triMesh.h"



#ifdef SDFGEN_DLL_EXPORTS
#define SDFGEN_DLL_API __declspec(dllexport) 
#else
#define SDFGEN_DLL_API __declspec(dllimport)
#endif


// 符号距离场结构体：
struct SDF_RESULT
{
	float step;
	int interCounts;
	verF origin;
	std::vector<unsigned> stepsCount;
	std::vector<float> SDFvalues;
};


// 输入三角网格，生成3D情形下的符号距离场； 
SDFGEN_DLL_API bool genSDF3D(SDF_RESULT& result, \
	const triMeshF& mesh,	const float step, const int interCounts = 5);

SDFGEN_DLL_API bool genSDF3D(SDF_RESULT& result, \
	const triMeshD& mesh, const float step, const int interCounts = 5);

