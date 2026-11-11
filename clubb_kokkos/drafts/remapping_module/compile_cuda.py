"""Compile isolated remapping TU/interfaces with root's actual CUDA toolchain; never execute."""
from pathlib import Path
import re,shlex,subprocess,json,hashlib
D=Path(__file__).resolve().parent;R=D.parents[2];lane=R/'build/parallel_cpu_lanes/mixing/source';B=R/'build/clubb_kokkos_nvhpc26_debug';out=R/'build/parallel_cpu_lanes/mixing/remapping_cuda_compile';out.mkdir(exist_ok=True)
F=(B/'CMakeFiles/clubb_kokkos_core_lib.dir/flags.make').read_text();cache=(B/'CMakeCache.txt').read_text()
compiler=re.search(r'^CMAKE_CXX_COMPILER:STRING=(.*)$',cache,re.M)[1]
def read_flags(key):return shlex.split(re.search(r'^'+key+r' = (.*)$',F,re.M)[1])
flags=read_flags('CXX_DEFINES')+read_flags('CXX_INCLUDES')+read_flags('CXX_FLAGS')
flags=[f.replace(str(R/'clubb_kokkos/src'),str(lane/'clubb_kokkos/src')) for f in flags]
flags+=['-I'+str(lane/'clubb_kokkos/src/CLUBB_core')]
records=[]
for name,source in [('remapping_module',lane/'clubb_kokkos/src/CLUBB_core/remapping_module.cpp'),('compile_interfaces',D/'compile_interfaces.cpp')]:
 obj=out/(name+'.o');cmd=[compiler,*flags,'-c',str(source),'-o',str(obj)];log=D/('cuda_'+name+'.log')
 with log.open('w') as stream:result=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT)
 row={'name':name,'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'command':cmd,'exit_code':result.returncode,'log':log.name,'object':str(obj),'object_sha256':hashlib.sha256(obj.read_bytes()).hexdigest() if result.returncode==0 else None};records.append(row)
 (D/'cuda_compile.json').write_text(json.dumps({'compile_only':True,'GPU_runtime_executed':False,'records':records},indent=2)+'\n')
 print(json.dumps({'name':name,'exit_code':result.returncode,'log':str(log)}),flush=True)
 if result.returncode:raise SystemExit(result.returncode)
