#ifndef BT_GPU_XPBD_TGS_BACKEND_H
#define BT_GPU_XPBD_TGS_BACKEND_H

#include "btGpuXpbdTgsApi.h"
#include <string>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

// Owns the optional CUDA runtime and its context without linking Bullet to CUDA.
// The DLL is loaded from the executable directory.
class btGpuXpbdTgsBackend
{
    const btGpuXpbdTgsApi* m_api = nullptr;
    void* m_state = nullptr;
    bool m_attempted = false;
    std::string m_error;
#ifdef _WIN32
    HMODULE m_library = nullptr;
#endif
public:
    btGpuXpbdTgsBackend() = default;
    btGpuXpbdTgsBackend(const btGpuXpbdTgsBackend&) = delete;
    btGpuXpbdTgsBackend& operator=(const btGpuXpbdTgsBackend&) = delete;
    ~btGpuXpbdTgsBackend()
    {
        if (m_state) m_api->destroy(m_state);
#ifdef _WIN32
        if (m_library) FreeLibrary(m_library);
#endif
    }
    bool ready()
    {
        if (m_attempted) return m_state != nullptr;
        m_attempted = true;
#ifdef _WIN32
        wchar_t path[32768];
        const DWORD length = GetModuleFileNameW(nullptr, path, 32768);
        if (!length || length >= 32768)
        {
            m_error = "Cannot resolve GPU XPBD TGS executable directory";
            return false;
        }
        std::wstring filename(path, length);
        filename.resize(filename.find_last_of(L"\\/") + 1);
        filename += L"BulletGpuXpbdTgs.dll";
        m_library = LoadLibraryExW(filename.c_str(), nullptr, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_SYSTEM32);
        if (!m_library)
        {
            m_error = "BulletGpuXpbdTgs.dll or its CUDA driver dependency is unavailable";
            return false;
        }
        auto get = reinterpret_cast<btGpuXpbdTgsGetApi>(GetProcAddress(m_library, "btGetGpuXpbdTgsApi"));
        m_api = get ? get(BT_GPU_XPBD_ABI_VERSION) : nullptr;
        if (!m_api || !m_api->create || !m_api->destroy || !m_api->upload || !m_api->forces || !m_api->step ||
            !m_api->download || !m_api->error || !m_api->surface || !m_api->surfacePositions || !m_api->embed || !m_api->targets || !m_api->disabledPairs || !m_api->updateNodes || !m_api->friction || !m_api->springs)
        {
            m_api = nullptr;
            m_error = "Incompatible GPU XPBD TGS backend";
            return false;
        }
        m_state = m_api->create();
        if (!m_state)
        {
            const char* message = m_api->error(nullptr);
            m_error = message ? message : "GPU XPBD TGS context creation failed";
        }
        return m_state != nullptr;
#else
        m_error = "GPU XPBD TGS currently requires Windows";
        return false;
#endif
    }
    const btGpuXpbdTgsApi* api() const { return m_api; }
    void* state() const { return m_state; }
    std::string error() const
    {
        if (!m_state) return m_error;
        const char* message = m_api->error(m_state);
        return message ? message : "";
    }
};
#endif
