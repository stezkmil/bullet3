
import argparse,ctypes,os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--cuda',required=True);p.add_argument('--source',required=True);p.add_argument('--output',required=True);a=p.parse_args()
dll=next(Path(a.cuda).rglob('nvrtc64_*.dll'));directory=os.add_dll_directory(str(dll.parent));lib=ctypes.CDLL(str(dll))
program=ctypes.c_void_p()
def check(code):
 if code:raise RuntimeError('NVRTC error '+str(code))
check(lib.nvrtcCreateProgram(ctypes.byref(program),Path(a.source).read_bytes(),b'guard.cu',0,None,None))
options=(ctypes.c_char_p*3)(b'--gpu-architecture=compute_75',b'--fmad=false',b'--std=c++14')
code=lib.nvrtcCompileProgram(program,3,options)
size=ctypes.c_size_t();lib.nvrtcGetProgramLogSize(program,ctypes.byref(size));log=ctypes.create_string_buffer(size.value);lib.nvrtcGetProgramLog(program,log)
if log.value:print(log.value.decode())
check(code);check(lib.nvrtcGetPTXSize(program,ctypes.byref(size)));ptx=ctypes.create_string_buffer(size.value);check(lib.nvrtcGetPTX(program,ptx));check(lib.nvrtcDestroyProgram(ctypes.byref(program)))
Path(a.output).write_text('static const char vbdGuardPtx[] = R"vbdptx('+ptx.value.decode()+')vbdptx";\n')
