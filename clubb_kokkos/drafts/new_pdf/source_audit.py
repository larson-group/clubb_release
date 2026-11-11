"""Compare Native and draft structure/math without compiling or running physics.

The audit preserves binary association when expanding integer squares/cubes,
normalizes explicit Native MIN/MAX adapters, and checks every assignment and
subroutine-call argument independently of the preparation script.
"""
from pathlib import Path
import ast,copy,hashlib,json,re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUTLINE=json.loads((HERE/'source_outline.json').read_text())
NATIVE={m:(ROOT/'src/CLUBB_core'/f'{m}.F90').read_text() for m in ('new_pdf_main','new_pdf')}
CPP={'new_pdf_main':(HERE/'new_pdf_main.cpp').read_text(),'new_pdf':(HERE/'new_pdf.hpp').read_text()}


def split_top(text):
    """Read comma-separated Native/target arguments at the current nesting depth."""
    result=[];level=0;begin=0
    for i,ch in enumerate(text):
        if ch=='(':level+=1
        elif ch==')':level-=1
        elif ch==',' and level==0:result.append(text[begin:i].strip());begin=i+1
    result.append(text[begin:].strip());return result


def f_statements(text):
    """Read Native logical statements, excluding comments/directives."""
    statements=[];pending=''
    for line in text.splitlines():
        code=line.split('!',1)[0].strip()
        if not code:continue
        pending+=' '+code.strip('& ')
        if not code.endswith('&'):statements.append(pending.strip());pending=''
    if pending:raise ValueError('unfinished Native continuation')
    return statements


def cpp_code(text):return re.sub(r'//[^\n]*','',text)


def definitions(text):
    text=cpp_code(text);result=[]
    for m in re.finditer(r'(?:KOKKOS_INLINE_FUNCTION\s+)?(?:real_t|void|auto|Kokkos::Array<real_t,\s*3>)\s+(\w+)\s*\(',text,re.S):
        depth=1;end=m.end()
        while depth:
            if text[end]=='(':depth+=1
            elif text[end]==')':depth-=1
            end+=1
        signature=text[m.end():end-1];tail=end
        while tail<len(text) and text[tail] not in '{;':tail+=1
        if tail==len(text) or text[tail]!='{':continue
        depth=1;finish=tail+1
        while depth:
            if text[finish]=='{':depth+=1
            elif text[finish]=='}':depth-=1
            finish+=1
        result.append({'name':m.group(1),'arguments':[a.strip().split()[-1].lstrip('&') for a in split_top(re.sub(r'Kokkos::Array<real_t,\s*3>', 'RealArray3', signature))],'signature':signature,'body':text[tail+1:finish-1]})
    return result


def ternaries(text):
    """Translate fully parenthesized C++ conditional expressions to Python AST syntax."""
    result='';i=0
    while i<len(text):
        if text[i]=='(':
            depth=1;end=i+1
            while depth:
                if text[end]=='(':depth+=1
                elif text[end]==')':depth-=1
                end+=1
            result+='('+ternaries(text[i+1:end-1])+')';i=end
        else:result+=text[i];i+=1
    depth=0;question=None;colon=None
    for i,ch in enumerate(result):
        if ch=='(':depth+=1
        elif ch==')':depth-=1
        elif depth==0 and ch=='?':question=i
        elif depth==0 and ch==':' and question is not None:colon=i;break
    if question is not None:
        if colon is None:raise ValueError(('ternary colon',result))
        result='('+result[question+1:colon]+' if '+result[:question]+' else '+result[colon+1:]+')'
    return result


class Normalize(ast.NodeTransformer):
    """Normalize only syntax/intrinsics and source integer powers, without reassociation."""
    def __init__(self,aliases=None):self.aliases=aliases or {}
    def visit_Name(self,node):
        if node.id in self.aliases:return copy.deepcopy(self.aliases[node.id])
        return node
    def visit_BinOp(self,node):
        self.generic_visit(node)
        if isinstance(node.op,ast.Pow):
            if isinstance(node.right,ast.Constant) and node.right.value in (2,3):
                base=node.left;value=base
                for _ in range(node.right.value-1):value=ast.BinOp(left=value,op=ast.Mult(),right=copy.deepcopy(base))
                return value
            return ast.Call(func=ast.Name(id='pow',ctx=ast.Load()),args=[node.left,node.right],keywords=[])
        return node
    def visit_Call(self,node):
        self.generic_visit(node)
        if isinstance(node.func,ast.Name) and node.func.id=='real_t':return node.args[0]
        if isinstance(node.func,ast.Name) and node.func.id in ('max','min') and len(node.args)>2:
            result=ast.Call(func=node.func,args=node.args[:2],keywords=[])
            for arg in node.args[2:]:result=ast.Call(func=node.func,args=[result,arg],keywords=[])
            return result
        return node
    def visit_IfExp(self,node):
        self.generic_visit(node);test=node.test
        if isinstance(test,ast.Compare) and len(test.ops)==1 and len(test.comparators)==1 and isinstance(test.ops[0],(ast.Gt,ast.Lt)) and ast.dump(node.body)==ast.dump(test.left) and ast.dump(node.orelse)==ast.dump(test.comparators[0]):
            return ast.Call(func=ast.Name(id='max' if isinstance(test.ops[0],ast.Gt) else 'min',ctx=ast.Load()),args=[test.left,test.comparators[0]],keywords=[])
        return node


def parse(expr,aliases=None):
    expr=expr.replace('Kokkos::','')
    expr=re.sub(r'\blambda_?\b','native_lambda',expr)
    expr=expr.replace('.true.','True').replace('.false.','False')
    expr=re.sub(r'\btrue\b','True',expr);expr=re.sub(r'\bfalse\b','False',expr)
    expr=re.sub(r'((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd\b',r'\1',expr)
    expr=expr.replace('.and.',' and ').replace('.or.',' or ').replace('.not.',' not ').replace('&&',' and ').replace('||',' or ')
    expr=re.sub(r'!(?!=)',' not ',expr)
    expr=ternaries(' '.join(expr.split()))
    tree=ast.parse(expr.strip(),mode='eval')
    return Normalize(aliases).visit(tree).body


def dump(expr,aliases=None):return ast.dump(parse(expr,aliases),include_attributes=False)



def block(info):
    source=NATIVE[info['module']];start=re.search(r'^  '+info['kind']+' '+info['name']+r'\b',source,re.M).start();end=source.index('  end '+info['kind']+' '+info['name'],start)
    return source[start:end]


def body_statements(info):
    records=f_statements(block(info));last=max(i for i,s in enumerate(records) if '::' in s);return records[last+1:]


def scalar_native(text,name):
    if name=='sort_roots':text=re.sub(r'(roots(?:_sorted)?)\(:,\s*([123])\)',lambda m:m.group(1)+'['+str(int(m.group(2))-1)+']',text)
    if name=='calc_setter_var_params':text=re.sub(r'(calc_mixture_fraction\s*\()\s*nz\s*,\s*',r'\1',text)
    return text


def assignments(native):
    result=[]
    for stmt in native:
        m=re.match(r'([A-Za-z_]\w*(?:%\w+)?(?:\([^)]*\)|\[[^]]*\])?)\s*=(?!=)\s*(.*)',stmt,re.S)
        if m:result.append(m.groups())
    return result


def c_assignments(body):
    return [(m.group(1),m.group(2)) for m in re.finditer(r'\b([A-Za-z_]\w*(?:\([^;={}]*?\)|\[[^]]*\])?)\s*=(?!=)\s*([^;{}]+);',body,re.S)]


def conditions(native):
    result=[]
    for stmt in native:
        m=re.match(r'(?:where|elsewhere|if|elseif)\s*\((.*)\)\s*(?:then)?$',stmt,re.I)
        if m and not m.group(1).strip().lower().startswith('any('):result.append(m.group(1))
    return result


def c_conditions(body):
    result=[]
    for m in re.finditer(r'\bif\s*\(',body):
        depth=1;end=m.end()
        while depth:
            if body[end]=='(':depth+=1
            elif body[end]==')':depth-=1
            end+=1
        result.append(body[m.end():end-1])
    return result


def norm_suffix(text):
    return re.sub(r'\b(\w+)_device\b',r'\1',text)


def calls(text,names):
    result=[]
    for m in re.finditer(r'\b('+ '|'.join(sorted(names,key=len,reverse=True))+r')\s*\(',text):
        depth=1;end=m.end()
        while depth:
            if text[end]=='(':depth+=1
            elif text[end]==')':depth-=1
            end+=1
        result.append((m.group(1),split_top(text[m.end():end-1])))
    return result


def view_normalize(text):
    text=norm_suffix(text).replace('clubb_params_host','clubb_params').replace('lambda_','lambda')
    text=re.sub(r'Kokkos::subview\(\s*(\w+)\s*,\s*i\s*,\s*Kokkos::ALL\(\)\s*\)',lambda m:'column('+m.group(1)+',i)',text)
    text=re.sub(r'\b(\w+)\(i,\s*:\)',lambda m:'column('+m.group(1)+',i)',text)
    return re.sub(r'\s+','',text)


def point_native(text,info):
    dims={}
    for stmt in f_statements(block(info)):
        if '::' not in stmt:continue
        attrs,names=stmt.split('::',1);dm=re.search(r'dimension\((.*?)\)',attrs,re.I)
        if dm:
            for n in split_top(names):dims[n.split('=',1)[0].strip()]=split_top(dm.group(1))
    text=re.sub(r'pdf_implicit_coefs_terms%(\w+)',lambda m:'pdf_'+m.group(1),text)
    for m in re.finditer(r'pdf_(\w+)\(i,:\)',text):dims['pdf_'+m.group(1)]=['ngrdcol','nz']
    out='';i=0
    while i<len(text):
        m=re.match(r'[A-Za-z_]\w*',text[i:])
        if not m:out+=text[i];i+=1;continue
        name=m.group();end=i+len(name);cursor=end
        while cursor<len(text) and text[cursor].isspace():cursor+=1
        if cursor<len(text) and text[cursor]=='(':
            level=1;finish=cursor+1
            while level:
                if text[finish]=='(':level+=1
                elif text[finish]==')':level-=1
                finish+=1
            args=split_top(text[cursor+1:finish-1]);shape=dims.get(name)
            out+=name+'('+', '.join('k' if a==':' else point_native(a,info) for a in args)+')'
            i=finish
        else:
            out+=name+'('+('k' if len(dims[name])==1 else 'i, k')+')' if name in dims else name;i=end
    return out


def cpp_limits_normalize(text):
    text=norm_suffix(text)
    return text


report={'status':'Source-only audit; uncompiled/untested','routines':[],'findings':[]}
DEFS={m:definitions(t) for m,t in CPP.items()}
for info in OUTLINE['routines']:
    name=info['name'];defs=[d for d in DEFS[info['module']] if d['name']==name]
    if name in ('new_pdf_driver','calc_responder_var','calc_limits_F_x_responder'):primary=defs[0]
    else:primary=next(d for d in defs if 'nz' not in d['arguments'])
    expected=[a['name'] for a in info['arguments']];profile=next((d for d in defs if d['arguments'] and d['arguments'][0]=='nz'),None)
    if profile is None or [a.removesuffix('_') if a=='lambda_' else a for a in profile['arguments']]!=expected:raise AssertionError(('profile signature',name,profile and profile['arguments'],expected))
    if primary is not profile:
        scalar=[a for a in expected if a!='nz'];scalar=['lambda_' if a=='lambda' else a for a in scalar]
        if primary['arguments']!=scalar:raise AssertionError(('scalar signature',name,primary['arguments'],scalar))
    native=body_statements(info)
    checked=0
    if name not in ('new_pdf_driver','calc_responder_var','calc_limits_F_x_responder'):
        na=assignments([scalar_native(s,name) for s in native]);ca=c_assignments(primary['body'])
        if len(na)!=len(ca):raise AssertionError(('assignment count',name,len(na),len(ca),na,ca))
        for (nl,nr),(cl,cr) in zip(na,ca):
            if dump(nl)!=dump(cl) or dump(nr)!=dump(cr):raise AssertionError(('assignment',name,nl,nr,cl,cr))
            checked+=1
        nc=[scalar_native(c,name) for c in conditions(native)];cc=c_conditions(primary['body'])
        if name=='calc_mixture_fraction':cc=cc[:-1]
        if [dump(c) for c in nc]!=[dump(c) for c in cc]:raise AssertionError(('conditions',name,nc,cc))
    elif name in ('new_pdf_driver','calc_responder_var'):
        na=assignments(native);ca=c_assignments(primary['body']);ca=[(norm_suffix(l),norm_suffix(r)) for l,r in ca if l!='i']
        ca=[(l,r) for l,r in ca if not any(x in r for x in ('to_execution_space','to_host_space','pdf_implicit_coefs_terms.')) and not l.startswith(('pdf_implicit_coefs_terms','clubb_params_host'))]
        na_math=[(l,r) for l,r in na if not re.search(r'\bcalc_\w+\(',r)]
        ca_math=[(l,r) for l,r in ca if not re.search(r'\bcalc_\w+\(',r)]
        if len(na_math)!=len(ca_math):raise AssertionError(('parent assignment count',name,len(na_math),len(ca_math),ca_math))
        for (nl,nr),(cl,cr) in zip(na_math,ca_math):
            if dump(point_native(nl,info))!=dump(cl) or dump(point_native(nr,info))!=dump(cr):raise AssertionError(('parent assignment',name,nl,nr,cl,cr))
            checked+=1
        nc=[point_native(c,info) for c in conditions(native)];cc=[norm_suffix(c) for c in c_conditions(primary['body'])]
        if [dump(c) for c in nc]!=[dump(c) for c in cc]:raise AssertionError(('parent conditions',name,nc,cc))
        names={x['name'] for x in OUTLINE['routines']};ncalls=calls('\n'.join(native),names);ccalls=calls(primary['body'],names)
        ncalls=[(n,[view_normalize(a) for a in aa]) for n,aa in ncalls];ccalls=[(n,[view_normalize(a) for a in aa]) for n,aa in ccalls]
        if ncalls!=ccalls:raise AssertionError(('call order/arguments',name,ncalls,ccalls))
    elif name=='calc_limits_F_x_responder':
        def limit_expr(text):
            text=re.sub(r'real\(\s*(sqrt_F_x_roots_[12])\(:,([123])\),\s*kind\s*=\s*core_rknd\s*\)',lambda m:m.group(1)+'(k,'+str(int(m.group(2))-1)+')',text)
            text=re.sub(r'(sqrt_F_x_roots_[12]_sorted)\(:,([123])\)',lambda m:m.group(1)+'(k,'+str(int(m.group(2))-1)+')',text)
            return point_native(text,info)
        na=[(l,r) for l,r in assignments(native) if not l.startswith('sqrt_F_x_roots')]
        ca=[(norm_suffix(l),norm_suffix(r)) for l,r in c_assignments(primary['body']) if re.match(r'(?:coef_[ABCD]_[12]|min_sqrt_F_x|max_sqrt_F_x|min_F_x_device|max_F_x_device)\(',l)]
        if len(na)!=len(ca):raise AssertionError(('limit assignment count',len(na),len(ca),na,ca))
        for (nl,nr),(cl,cr) in zip(na,ca):
            if dump(limit_expr(nl))!=dump(cl) or dump(limit_expr(nr))!=dump(cr):raise AssertionError(('limit assignment',nl,nr,cl,cr))
            checked+=1
        nc=[limit_expr(c) for c in conditions(native)];cc=[norm_suffix(c) for c in c_conditions(primary['body'])]
        if [dump(c) for c in nc]!=[dump(c) for c in cc]:raise AssertionError(('limit conditions',nc,cc))
        rootcalls=calls(primary['body'],{'cubic_solve'})
        expectedroots=[['nz']+['coef_'+a+'_'+str(c) for a in 'ABCD'] for c in (1,2)]
        for (_,actual),expectedroot in zip(rootcalls,expectedroots):
            normalized=[re.sub(r'to_host_space\(\s*(\w+)\s*\)',r'\1',a) for a in actual]
            if normalized!=expectedroot:raise AssertionError(('root call mapping',normalized,expectedroot))
        if len(rootcalls)!=2:raise AssertionError(('root crossing count',rootcalls))
        sortcalls=calls(primary['body'],{'sort_roots'})
        if [(n,[view_normalize(a) for a in aa]) for n,aa in sortcalls]!=[('sort_roots',['nz','sqrt_F_x_roots_1']),('sort_roots',['nz','sqrt_F_x_roots_2'])]:raise AssertionError(('sort call mapping',sortcalls))
    report['routines'].append({'name':name,'native_arguments':len(expected),'profile_contract':'exact Native argument order; no optionals','scalar_overload':primary is not profile,'assignment_trees_checked':checked,'source_branches':'independent AST comparison'})

# Keep relative procedure definition order, allowing adjacent overloads.
for module,defs in DEFS.items():
    names=[]
    for d in defs:
        if not names or names[-1]!=d['name']:names.append(d['name'])
    expected=[r['name'] for r in OUTLINE['routines'] if r['module']==module]
    if names!=expected:raise AssertionError(('routine order',module,names,expected))

for p in HERE.glob('*'):
    if p.suffix not in ('.hpp','.cpp'):continue
    for n,line in enumerate(p.read_text().splitlines(),1):
        if len(line)>100 or line.rstrip()!=line:raise AssertionError(('whitespace',p.name,n,len(line)))
report['hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [*(ROOT/'src/CLUBB_core'/f'{m}.F90' for m in NATIVE),*(HERE/n for n in ('new_pdf.hpp','new_pdf_main.cpp','new_pdf_main.hpp','f2k_dependencies.hpp'))]}
report['validation']='No compilation/tests/numerical execution. CPU Fortran versus CPU Kokkos BFB is required at root integration; GPU differences are milestone evidence.'
pairs={')':'(',']':'[','}':'{'}
for path in [p for p in HERE.iterdir() if p.suffix in ('.hpp','.cpp')]:
    code=cpp_code(path.read_text());code=re.sub(r'"(?:\\.|[^"\\])*"','""',code);stack=[]
    for index,ch in enumerate(code):
        if ch in '([{':stack.append((ch,index))
        elif ch in ')]}':
            if not stack or stack[-1][0]!=pairs[ch]:raise AssertionError(('delimiter',path.name,index,ch,stack[-1:]))
            stack.pop()
    if stack:raise AssertionError(('unclosed delimiter',path.name,stack))
for path in HERE.glob('*.md'):
    for link in re.findall(r'\]\(([^)]+)\)',path.read_text()):
        if not link.startswith(('http:','https:','#')) and not (path.parent/link.split('#',1)[0]).exists():raise AssertionError(('local link',path.name,link))
report['cleanup']={'cpp_delimiters':'balanced after comments/string removal','whitespace':'100 columns; no trailing whitespace','documentation_links':'all local destinations exist','manual_review':['all three advection flags and Native column order','same Fortran invalid-mixture condition including nonpositive/NaN spread branch','nine implicit coefficient/explicit-term packed fields','all scalar/profile array index adapters and caller-memory-space return paths','complex root real extraction, exact solve/sort/selection order','no statistics calls in Native; deprecated diagnostic F outputs retained'],'runtime':'not run; root owns promotion and CPU BFB gate'}
(HERE/'source_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print('Source-only audit passed for all13 routines; runtime and BFB remain unverified.')
