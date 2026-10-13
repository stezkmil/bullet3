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
    if (argc != 3)
    {
        std::fprintf(stderr, "Usage: BulletVbdCompilePtx source.cu output.h\n");
        return 1;
    }
    std::ifstream input(argv[1], std::ios::binary);
    if (!input)
    {
        std::fprintf(stderr, "Cannot open kernel source: %s\n", argv[1]);
        return 1;
    }
    const std::string source((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    if (input.bad())
        return 1;
    Program program;
    if (!check(nvrtcCreateProgram(&program.handle, source.c_str(), argv[1], 0, nullptr, nullptr)))
        return 1;
    const char* options[] = {"--gpu-architecture=compute_75", "--fmad=false", "--std=c++14"};
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
    // Do not create an output on compilation failure: a subsequent build must retry.
    std::ofstream output(argv[2], std::ios::binary | std::ios::trunc);
    output << "static const char vbdGuardPtx[] = R\"vbdptx(" << ptx.data() << ")vbdptx\";\n";
    output.close();
    if (!output)
    {
        std::fprintf(stderr, "Cannot write embedded PTX: %s\n", argv[2]);
        std::remove(argv[2]);
        return 1;
    }
    return 0;
}
