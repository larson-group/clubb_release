"""Generate only this isolated Native-order source draft, following SOURCE_OUTLINE."""
from pathlib import Path
from dataclasses import dataclass,field
import ast,contextlib,io,json,re,runpy,sys,textwrap
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,declaration
from generate_derived_types import split_top
with contextlib.redirect_stdout(io.StringIO()):utils=runpy.run_path('/tmp/clubb_clip_review.py')
split=utils['split'];balanced=utils['balanced']
source=(root/'src/CLUBB_core/mono_flux_limiter.F90').read_text()
names=re.findall(r'^  subroutine (\w+)',source,re.M)
blocks={n:re.search(r'^  subroutine '+n+r'\b.*?^  end subroutine '+n+r'\b',source,re.M|re.S)[0] for n in names}
public={'monotonic_turbulent_flux_limit','calc_turb_adv_range'}
metadata={};arguments={};initializers={}
for n,b in blocks.items():
    meta={};init={};ss=statements(b)
    for s in ss:
        if '::' not in s:continue
        meta.update(declaration(s))
        if 'parameter' in s.lower():
            for v in split_top(s.split('::',1)[1]):
                k,value=v.split('=',1);init[k.strip().lower()]=value.strip()
    metadata[n]=meta;initializers[n]=init
    arguments[n]=split_top(ss[0][ss[0].index('(')+1:ss[0].rindex(')')])

# Reuse the already reviewed Native interface/comment grouping logic, without
# importing its draft writer or mutating its owned folder.
env=globals()
for node in ast.parse((root/'clubb_kokkos/drafts/fill_holes/prepare.py').read_text()).body:
    if isinstance(node,ast.FunctionDef) and node.name in ('description','docs','signature'):
        code=ast.get_source_segment((root/'clubb_kokkos/drafts/fill_holes/prepare.py').read_text(),node)
        if node.name=='description':code=code.replace("'! Description:'","'! Description'")
        if node.name=='docs':
            code=code.replace('for i in range(start+1,len(lines)):',
              "limit=next((i for i in range(start+1,len(lines)) if 'Begin Code' in lines[i]),len(lines))\n        for i in range(start+1,limit):")
        if node.name=='signature':
            code=code.replace("('grid','hm_metadata_type','stats_type')","('grid','hm_metadata_type','stats_type','err_info_type')")
            code=code.replace("else 'Bool1D' if a.kind=='logical' else", "else ('Kokkos::View<int'+('*'*len(a.dims))+', Kokkos::LayoutLeft>') if a.kind=='integer' else 'Bool1D' if a.kind=='logical' else")
        exec(code,env)

@dataclass
class Node:
    kind:str
    text:str=''
    body:list=field(default_factory=list)
    other:list=field(default_factory=list)

def tokens(text):
    out=[];pending=''
    for raw in text.splitlines():
        line=raw.strip()
        if line.startswith(('!$acc','#')):continue
        if line.startswith('!'):out.append(Node('comment',line[1:]));continue
        if not line:continue
        quote=None;cut=len(raw)
        for i,ch in enumerate(raw):
            if quote:
                if ch==quote:quote=None
            elif ch in '\"\'':quote=ch
            elif ch=='!':cut=i;break
        code=raw[:cut].strip()
        if not code:continue
        pending+=' '+code.strip('& ')
        if code.endswith('&'):continue
        out.append(Node('stmt',pending.strip()));pending=''
    assert not pending
    return out

def parse(ts,pos=0,stop=(),single=False):
    out=[]
    while pos<len(ts):
        t=ts[pos];s=t.text
        if t.kind=='comment':out.append(t);pos+=1;continue
        if any(re.match(x,s,re.I) for x in stop):return out,pos
        if re.match(r'if\s*\(',s,re.I):
            condition,end=balanced(s,s.index('('));tail=s[end:].strip()
            if tail.lower()=='then':
                body,pos=parse(ts,pos+1,(r'else\b',r'elseif',r'end\s*if'))
                other=[]
                if re.match(r'else\s*if|elseif',ts[pos].text,re.I):
                    ts[pos].text=re.sub(r'^else\s*if|^elseif','if',ts[pos].text,flags=re.I)
                    other,pos=parse(ts,pos,single=True);pos-=1
                elif ts[pos].text.lower()=='else':other,pos=parse(ts,pos+1,(r'end\s*if',))
                assert re.match(r'end\s*if',ts[pos].text,re.I),(s,ts[pos].text)
                out.append(Node('if',condition,body,other));pos+=1
                if single:return out,pos
            else:out.append(Node('if',condition,[Node('stmt',tail)]));pos+=1
            continue
        if re.match(r'do\b',s,re.I):
            body,pos=parse(ts,pos+1,(r'end\s*do',));kind='for' if '=' in s else 'while'
            out.append(Node(kind,s,body));pos+=1;continue
        if re.match(r'select\s+case',s,re.I):
            cases=[];pos+=1
            while not re.match(r'end\s*select',ts[pos].text,re.I):
                if ts[pos].kind=='comment':cases.append(ts[pos]);pos+=1;continue
                label=ts[pos].text;assert re.match(r'case',label,re.I),label
                body,pos=parse(ts,pos+1,(r'case',r'end\s*select'));cases.append(Node('case',label,body))
            out.append(Node('select',s,cases));pos+=1;continue
        out.append(t);pos+=1
    return out,pos

def walk(nodes):
    for n in nodes:
        yield n
        yield from walk(n.body);yield from walk(n.other)

def type_of(a):
    if a.dims:return ('Real'+str(len(a.dims))+'D') if a.kind=='real' else 'Bool1D' if a.kind=='logical' else 'Kokkos::View<int'+('*'*len(a.dims))+',Kokkos::LayoutLeft>'
    return {'integer':'int','real':'real_t','logical':'bool','character':'std::string'}[a.kind]

grid_ranks={'zm':2,'zt':2,'dzm':2,'dzt':2,'invrs_dzm':2,'invrs_dzt':2,'weights_zt2zm':3,'weights_zm2zt':3}
all_bodies={};headers={}
for n,b in blocks.items():
    ts=tokens(b);last=max(i for i,x in enumerate(ts) if x.kind=='stmt' and '::' in x.text)
    all_bodies[n],_=parse(ts[last+1:-1]);headers[n]=[x for x in ts[1:last+1] if x.kind=='comment']
    # Full argument descriptions already follow their named signature. Remove
    # duplicate continuation fragments from relocated declaration-section prose.
    lines=b.splitlines();orphan=set();decl_start=next(i for i,s in enumerate(lines) if 'implicit none' in s.lower())
    decl_end=next(i for i in range(decl_start+1,len(lines)) if 'Begin Code' in lines[i])
    for i in range(decl_start+1,decl_end):
        if '!' not in lines[i] or lines[i].lstrip().startswith('!'):continue
        code,_=lines[i].split('!',1)
        if not any(re.search(r'\b'+re.escape(a)+r'\s*(?:,?\s*&)?\s*$',code) for a in arguments[n]):continue
        for line in lines[i+1:decl_end]:
            if not line.strip().startswith('!') or not line.strip(' !-'):break
            orphan.add(line.split('!',1)[1])
    headers[n]=[x for x in headers[n] if x.text not in orphan]

def power2(s):
    while '**' in s:
        pos=s.rindex('**');end=pos
        while s[end-1].isspace():end-=1
        start=end-1
        if s[start]==')':
            depth=1;start-=1
            while depth:depth+=(s[start]==')')-(s[start]=='(');start-=1
            start+=1;j=start-1
            while j>=0 and (s[j].isalnum() or s[j] in '_.%'):j-=1
            if j<start-1:start=j+1
        else:
            while start>=0 and (s[start].isalnum() or s[start] in '_.%'):start-=1
            start+=1
        exponent=pos+2
        while s[exponent].isspace():exponent+=1
        assert s[exponent]=='2',s
        value=s[start:end];s=s[:start]+'('+value+' * '+value+')'+s[exponent+1:]
    return s

class Emitter:
    def __init__(self,n):self.n=n;self.meta=metadata[n];self.counter=0
    def expr(self,s):
        strings=[]
        def store(m):strings.append(m[0]);return '@S'+str(len(strings)-1)+'@'
        s=re.sub(r'"[^"\n]*"|\'[^\'\n]*\'',store,s.strip())
        s=power2(s).replace('%','.')
        s=re.sub(r'any\(\s*err_info.err_code\s*==\s*clubb_fatal_error\s*\)',r'any_equal(err_info.err_code,clubb_fatal_error)',s,flags=re.I)
        for a,b in [('.and.','&&'),('.or.','||'),('.not.','!'),('.true.','true'),('.false.','false'),('/=','!=')]:s=re.sub(re.escape(a),b,s,flags=re.I)
        s=re.sub(r'((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd',lambda m:'real_t('+m[1]+')',s)
        for fun in ('sqrt','abs','exp','erf'):s=re.sub(r'(?<![\w:])'+fun+r'\s*\(',r'Kokkos::'+fun+'(',s,flags=re.I)
        for fun,op in [('min','<'),('max','>')]:
            while m:=re.search(r'(?<![\w:])'+fun+r'\s*\(',s,re.I):
                inside,end=balanced(s,m.end()-1);items=split(inside);assert len(items)==2,(self.n,s)
                a,b=items;s=s[:m.start()]+'(('+a+') '+op+' ('+b+') ? ('+a+') : ('+b+'))'+s[end:]
        for k in set(re.findall(r'gr\.(\w+)',s)):s=s.replace('gr.'+k,'gr_'+k)
        arrays={a.name:len(a.dims) for a in self.meta.values() if a.dims}|{'gr_'+k:v for k,v in grid_ranks.items()}
        for name,rank in sorted(arrays.items(),key=lambda v:-len(v[0])):
            while m:=re.search(r'(?<![\w.])'+name+r'\s*\(',s):
                inside,end=balanced(s,m.end()-1);indices=split(inside);assert len(indices)==rank,(self.n,name,inside)
                if all(a.strip()==':' for a in indices):value=name
                else:
                    assert not any(':' in a for a in indices),(self.n,s)
                    value=name+'('+','.join('('+self.expr(a)+')-1' for a in indices)+')'
                strings.append(value);s=s[:m.start()]+'@S'+str(len(strings)-1)+'@'+s[end:]
        s=s.replace('stats.l_sample','stats.flags().l_sample')
        s=re.sub(r'\bclubb_at_least_debug_level_api\s*\(','f2k::clubb_at_least_debug_level_api(',s)
        for i in reversed(range(len(strings))):s=s.replace('@S'+str(i)+'@',strings[i])
        return s
    def scalar_locals(self,nodes,exclude=()):
        text='\n'.join(x.text for x in walk(nodes));decl=[]
        for key,a in self.meta.items():
            if a.dims or a.name in arguments[self.n] or key in initializers[self.n] or a.name in exclude:continue
            if re.search(r'\b'+a.name+r'\s*=',text):decl.append(type_of(a)+' '+a.name+';')
        return '\n'.join(decl)
    def call(self,name,converted):
        # Preserve source physical continuation groups at ordinary science calls.
        match=re.search(r'(?im)^\s*call\s+'+re.escape(name)+r'\s*\(',blocks[self.n])
        if match:
            start=match.end()-1;_,end=balanced(blocks[self.n],start)
            # The Native source has ampersands/comments, so obtain physical group
            # lengths directly instead of parsing its raw multiline signature.
            raw=blocks[self.n][match.start():end].lstrip();counts=[]
            for i,line in enumerate(raw.splitlines()):
                code=line.split('!',1)[0].replace('&','').strip()
                if i==0:code=code[code.index('(')+1:]
                if i==len(raw.splitlines())-1:code=code.removesuffix(')')
                values=[x for x in split(code) if x.strip()]
                if values:counts.append(len(values))
        else:counts=[len(converted)]
        if sum(counts)!=len(converted):counts=[len(converted)]
        if name in metadata:intents=[metadata[name][k.lower()].intent for k in arguments[name]]
        elif name=='term_ma_zt_lhs':intents=['in']*9+['out']
        else:intents=['in']*(len(converted)-1)+['inout']
        groups=[];offset=0
        for count in counts:
            current=[];tag=None
            for j in range(offset,offset+count):
                intent={'in':'In','inout':'InOut','out':'Out'}[intents[j]]
                if current and (tag!=intent or len(','.join(current+[converted[j]]))>78):
                    groups.append((', '.join(current),tag));current=[]
                current.append(converted[j]);tag=intent
            if current:groups.append((', '.join(current),tag))
            offset+=count
        return name+'(\n'+'\n'.join('  '+line+(',' if i<len(groups)-1 else ');')+' // '+tag for i,(line,tag) in enumerate(groups))
    def emit(self,nodes,kernel=False):
        out=[]
        for x in nodes:
            kind,s=x.kind,x.text
            if kind=='comment':out.append('//'+s);continue
            if kind=='raw':out.append(s);continue
            if kind=='if':
                out+=['\nif ('+self.expr(s)+') {',self.emit(x.body,kernel),'}']
                if x.other:out[-1]+=' else {\n'+self.emit(x.other,kernel)+'\n}'
                continue
            if kind=='select':
                inside,_=balanced(s,s.index('('));out.append('\nswitch ('+self.expr(inside)+') {')
                for case in x.body:
                    if case.kind=='comment':out.append('//'+case.text);continue
                    label='default' if 'default' in case.text.lower() else 'case '+self.expr(case.text[case.text.index('(')+1:case.text.rindex(')')])
                    out += [label+': {',self.emit(case.body,kernel),'break;','}']
                out.append('}');continue
            if kind=='while':
                condition='true' if s.lower()=='do' else self.expr(s[s.index('(')+1:s.rindex(')')])
                out+=['\nwhile ('+condition+') {',self.emit(x.body,kernel),'}'];continue
            if kind=='for':
                m=re.match(r'do\s+(\w+)\s*=\s*(.*)',s,re.I);var=m[1];bounds=split(m[2]);lo,hi=bounds[:2];step=bounds[2] if len(bounds)>2 else '1'
                if not kernel:
                    # All outer scientific loops are kernels, except the Native
                    # top-spike diagnostic stage, replaced explicitly below.
                    loops=[(var,lo,hi,step)];body=x.body
                    while len([a for a in body if a.kind!='comment'])==1 and next(a for a in body if a.kind!='comment').kind=='for':
                        comments=[a for a in body if a.kind=='comment'];out.extend('//'+a.text for a in comments)
                        child=next(a for a in body if a.kind!='comment');m=re.match(r'do\s+(\w+)\s*=\s*(.*)',child.text,re.I);b=split(m[2]);loops.append((m[1],b[0],b[1],b[2] if len(b)>2 else '1'));body=child.body
                    count_exprs=['std::max(0,('+self.expr(b)+'-('+self.expr(a)+'))/('+self.expr(step)+')+1)' for _,a,b,step in loops]
                    self.counter+=1;reduction=any(a.kind=='stmt' and re.match(r'l_any_adjustment_needed\s*=',a.text) for a in walk(body))
                    counts=['loop_'+str(self.counter)+'_extent_'+str(j) for j in range(len(loops))]
                    out += ['const int '+name+'='+value+';' for name,value in zip(counts,count_exprs)]
                    prefix=('Kokkos::parallel_reduce' if reduction else 'Kokkos::parallel_for')+'("'+self.n+'::loop'+str(self.counter)+'",'+'*'.join(counts)+',KOKKOS_LAMBDA(const int cell'+(',int& adjustment_hits' if reduction else '')+') {'
                    out += ['\n'+prefix,'int remaining=cell;']
                    # Map contiguous columns to neighboring threads regardless
                    # of whether Native wrote i or k as the outer loop.
                    mapping=sorted(zip(loops,counts),key=lambda pair:0 if pair[0][0]=='i' else 1)
                    for (v,a,b,step),count in mapping:
                        out.append('const int '+v+'='+self.expr(a)+'+(remaining%('+count+'))*('+self.expr(step)+'); remaining/=('+count+');')
                    out.append(self.scalar_locals(body,{v for v,*_ in loops}|({'l_any_adjustment_needed'} if reduction else set())))
                    text=self.emit(body,True)
                    if reduction:text=text.replace('l_any_adjustment_needed = true;','adjustment_hits = 1;')
                    out.append(text)
                    if reduction:out += ['},adjustment_count);','l_any_adjustment_needed = adjustment_count != 0;']
                    else:out.append('});')
                    continue
                # Native DO bounds and direction remain explicit inside each kernel.
                condition=('('+self.expr(step)+')>0 ? '+var+'<='+self.expr(hi)+' : '+var+'>='+self.expr(hi))
                out += ['\nfor ('+var+'='+self.expr(lo)+'; '+condition+'; '+var+'+='+self.expr(step)+') {',self.emit(x.body,True),'}'];continue
            assert kind=='stmt',(kind,s)
            if s.lower()=='return':out.append('return;');continue
            if s.lower()=='exit':out.append('break;');continue
            if s.lower().startswith('call '):
                name=s.split()[1].split('(')[0];inside,_=balanced(s,s.index('('));actual=split(inside)
                converted=[self.expr(v) for v in actual]
                if name.startswith('stats_'):
                    if '/' in actual[1]:
                        arr=actual[1].split('/')[0].strip();self.counter+=1;values='stats_values_'+str(self.counter)
                        out += ['{','Real2D '+values+'(Kokkos::view_alloc(Kokkos::WithoutInitializing,"'+self.n+'::'+values+'"),'+arr+'.extent(0),'+arr+'.extent(1));',
                          'Kokkos::parallel_for("'+self.n+'::stats_values",'+arr+'.size(),KOKKOS_LAMBDA(int cell) {const int i=cell%'+arr+'.extent(0),k=cell/'+arr+'.extent(0);'+values+'(i,k)='+arr+'(i,k)/dt;});',
                          'f2k::'+name+'('+converted[0]+','+values+','+converted[2]+');','}']
                    else:out.append('f2k::'+name+'('+','.join(converted)+');')
                else:out.append(self.call(name,converted))
                continue
            if re.match(r'write\s*\(',s,re.I):
                inside,end=balanced(s,s.index('('));items=split(s[end:]);values=[self.expr(a) for a in items]
                out.append('std::cerr << '+' << " " << '.join(values)+' << "\\n";');continue
            m=re.match(r'(\w+(?:%\w+)?)(\(.*?\))?\s*=\s*(.*)',s,re.S)
            assert m,(self.n,s)
            lhs,index,rhs=m.groups()
            rank=len(self.meta[lhs.lower()].dims) if lhs.lower() in self.meta else 0
            if rank and (index is None or ':' in index):
                assert rhs.startswith(('zm2zt_api','zt2zm_api')),(self.n,s)
                out.append('Kokkos::deep_copy('+lhs+','+self.expr(rhs)+');')
            else:out.append(self.expr(lhs+(index or ''))+' = '+self.expr(rhs)+';')
        return '\n'.join(out)

def replace_spike(nodes):
    for index,x in enumerate(nodes):
        if x.kind=='for' and 'xm_vert_integral' in '\n'.join(a.text for a in walk(x.body)):
            nodes[index]=Node('raw','@SPIKE@');continue
        replace_spike(x.body);replace_spike(x.other)
replace_spike(all_bodies[names[0]])

def declarations(n,e):
    out=[]
    for key,a in metadata[n].items():
        if a.name in arguments[n]:continue
        if a.dims:
            dimensions=[]
            for d in a.dims:
                if ':' in d:
                    lo,hi=d.split(':');dimensions.append('('+e.expr(hi)+'-('+e.expr(lo)+')+1)')
                else:dimensions.append(e.expr(d))
            out.append(type_of(a)+' '+a.name+'(Kokkos::view_alloc(Kokkos::WithoutInitializing,"'+n+'::'+a.name+'"),'+','.join(dimensions)+');')
        else:out.append(('constexpr ' if key in initializers[n] else '')+type_of(a)+' '+a.name+('='+e.expr(initializers[n][key]) if key in initializers[n] else '')+';')
        # Attach actual Native declaration descriptions to their named owner.
        raw=blocks[n].split('implicit none',1)[1]
        for line in raw.splitlines():
            if '!' not in line or line.lstrip().startswith('!'):continue
            code,note=line.split('!',1)
            if re.search(r'\b'+re.escape(a.name)+r'\s*(?:,\s*&|&)?\s*$',code):
                out[-1]+=' // '+note.strip();break
    if n==names[0]:out.append('int adjustment_count=0;')
    executable='\n'.join(statements(blocks[n]))
    for key in sorted(set(re.findall(r'gr%(\w+)',executable))):out.append('const auto gr_'+key+'=gr.'+key+';')
    if n=='calc_mean_w_up_down_component':
        out=[s.replace('Variance of w (1st PDF component)            [m^2/s^2]',
                       'Standard deviation of w [m/s]') if s.startswith('real_t sigma_w_i_zm;') else s for s in out]
    return '\n'.join(out)

support='''// Bound turbulent fluxes and compensate the affected mean fields.
// Native derivations/routine order are retained; see mono_flux_limiter.hpp.
#include "mono_flux_limiter.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp" // l_force_descending_solves must come from its existing owner.
#include "mean_adv.hpp" // Reuse stable owner; promotion pending.
#include "matrix_solver_wrapper.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
namespace clubb_kokkos {
namespace {
// These selector values must match the Native advance_xm_wpxp_module owner.
constexpr int mono_flux_thlm=1,mono_flux_rtm=2,mono_flux_um=4,mono_flux_vm=5;
constexpr int ndiags3=3;
'''
for n in names:
    if n not in public:support+=signature(n)+'\n;\n'
support+='} // namespace\n'

def write():
    global support
    exec((folder/'special_bodies.py').read_text(),globals())
    out=support
    private_open=False
    for n in names:
        if n not in public and not private_open:out+='\nnamespace {\n';private_open=True
        if n in public and private_open:out+='} // namespace\n';private_open=False
        e=Emitter(n)
        out+='\n//=============================================================================\n'+signature(n)+'\n{\n'
        out+='\n'.join('//'+x.text for x in headers[n])+'\n'
        out+=declarations(n,e)+'\n'
        body=solver_body if n=='mfl_xm_solve' else e.emit(all_bodies[n])
        out+=body.replace('@SPIKE@',spike_body)+'\n}\n'
    if private_open:out+='} // namespace\n'
    out+='} // namespace clubb_kokkos\n'
    (folder/'mono_flux_limiter.cpp').write_text(out)
    h='''// Bound turbulent transport using the fields that can reach each vertical level.
//
// monotonic_turbulent_flux_limit limits wpxp and re-solves/adjusts xm, preserving
// Native budgets, spike removal and fatal ordering. calc_turb_adv_range computes
// Native1based lower/upper influencing-level indices from the PDF velocities.
// Fields and integer ranges are column-first LayoutLeft execution Views. Matrix
// scratch stays (band,column,level) until the actual solver boundary. No mutable
// scientific state survives a call. See INTEGRATION.md: draft only, uncompiled.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"
namespace clubb_kokkos {
'''
    h+='\n\n'.join(signature(n,True)+'\n;' for n in names if n in public)
    h+='\n} // namespace clubb_kokkos\n';(folder/'mono_flux_limiter.hpp').write_text(h)

if __name__=='__main__':write()
