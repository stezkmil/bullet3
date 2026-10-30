#include <nvrtc.h>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace
{
struct Program
{
    nvrtcProgram handle = nullptr;
    ~Program()
    {
        if (handle)
            nvrtcDestroyProgram(&handle);
    }
};

bool check(nvrtcResult result)
{
    if (result == NVRTC_SUCCESS)
        return true;
    std::fprintf(stderr, "NVRTC: %s\n", nvrtcGetErrorString(result));
    return false;
}
}

int main(int argc, char** argv)
{
    if (argc != 3 && (argc != 4 || std::string(argv[3])!="--fast-math"))
    {
        std::fprintf(stderr, "Usage: XpbdCompilePtx source.cu output.h [--fast-math]\n");
        return 1;
    }
    std::ifstream input(argv[1], std::ios::binary);
    if (!input)
    {
        std::fprintf(stderr, "Cannot open kernel source: %s\n", argv[1]);
        return 1;
    }
    const std::string kernelSource((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const std::string sourcePath=argv[1];
    std::ifstream narrowphase(sourcePath.substr(0,sourcePath.find_last_of("/\\")+1)+"physx_gjk.cuh",std::ios::binary);
    if(!narrowphase) { std::fprintf(stderr,"Cannot read CUDA narrowphase source\n");return 1; }
    const std::string source=std::string((std::istreambuf_iterator<char>(narrowphase)),std::istreambuf_iterator<char>())+"\n"+kernelSource;
    if (input.bad())
        return 1;
    Program program;
    if (!check(nvrtcCreateProgram(&program.handle, source.c_str(), argv[1], 0, nullptr, nullptr)))
        return 1;
    const bool fast=argc==4;
    const char* options[] = {"--gpu-architecture=compute_75", fast?"--use_fast_math":"--fmad=false", "--std=c++14"};
    const auto compiled = nvrtcCompileProgram(program.handle, 3, options);
    size_t size = 0;
    if (!check(nvrtcGetProgramLogSize(program.handle, &size)))
        return 1;
    std::vector<char> log(size);
    if (size)
    {
        if (!check(nvrtcGetProgramLog(program.handle, log.data())))
            return 1;
        std::fputs(log.data(), stderr);
    }
    if (!check(compiled) || !check(nvrtcGetPTXSize(program.handle, &size)))
        return 1;
    std::vector<char> ptx(size);
    if (!check(nvrtcGetPTX(program.handle, ptx.data())))
        return 1;
    if (std::string(ptx.data()).find(".f64") != std::string::npos)
    {
        std::fprintf(stderr, "FP32 kernels unexpectedly contain FP64 PTX instructions\n");
        return 1;
    }
    // Do not create an output on compilation failure: a subsequent build must retry.
    std::ofstream output(argv[2], std::ios::binary | std::ios::trunc);
    // Numeric bytes avoid MSVC's string-literal size limit; NVRTC includes the final NUL.
    output << "static const char " << (fast?"xpbdFastPtx":"xpbdPtx") << "[] = {\n";
    for (size_t i = 0; i < ptx.size(); ++i)
    {
        output << static_cast<int>(ptx[i]) << ',';
        if ((i + 1) % 32 == 0)
            output << '\n';
    }
    output << "\n};\n";
    output.close();
    if (!output)
    {
        std::fprintf(stderr, "Cannot write embedded PTX: %s\n", argv[2]);
        std::remove(argv[2]);
        return 1;
    }
    return 0;
}
