"""Prepare isolated Native-named hybrid PDF drafts; no builds or physics runs.

All writes stay in this draft directory. Source contracts and grouping remain
canonical; the generated draft needs separate source review and runtime checks.
"""
from pathlib import Path
import re,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUTLINE=json.loads((HERE/'source_outline.json').read_text())
CONTRACTS={r['name']:r for r in OUTLINE['routines']}
TEXTS={name:(ROOT/'src/CLUBB_core'/f'{name}.F90').read_text() for name in ('new_hybrid_pdf_main','new_hybrid_pdf')}


def split_top(text):
    """Split Native comma lists without splitting nested expressions."""
    result=[];depth=0;begin=0
    for i,c in enumerate(text):
        if c=='(':depth+=1
        elif c==')':depth-=1
        elif c==',' and depth==0:result.append(text[begin:i].strip());begin=i+1
    result.append(text[begin:].strip())
    return result


def records(text):
    """Join Native logical statements, retaining explanatory comments."""
    result=[];pending=''
    for raw in text.splitlines():
        stripped=raw.lstrip()
        if stripped.startswith('!'):
            if not stripped.startswith('!$'):result.append(('comment',stripped[1:]))
            continue
        code=raw.split('!',1)[0].strip()
        if not code:
            if not pending:result.append(('blank',''))
            continue
        pending+=' '+code.strip('& ')
        if not code.endswith('&'):result.append(('code',pending.strip()));pending=''
    if pending:raise ValueError('unfinished Native statement')
    return result


def routine(name):
    """Read one complete Native procedure and its local declaration inventory."""
    contract=CONTRACTS[name];text=TEXTS[contract['module']]
    start=re.search(r'^  (?:elemental )?(?:subroutine|function) '+name+r'\b',text,re.M).start()
    end=text.index('  end '+contract['kind']+' '+name,start)
    end=text.index('\n',end)
    block=text[start:end];rec=records(block)
    last=max(i for i,(kind,value) in enumerate(rec) if kind=='code' and '::' in value)
    decls={}
    for kind,stmt in rec[:last+1]:
        if kind!='code' or '::' not in stmt:continue
        attrs,vars=stmt.split('::',1)
        dm=re.search(r'dimension\s*\((.*?)\)',attrs,re.I)
        for var in split_top(vars):
            pair=var.split('=',1);var=pair[0].strip()
            decls[var]={'dims':split_top(dm.group(1)) if dm else [],'kind':'integer' if attrs.strip().startswith('integer') else 'implicit_coefs_terms' if 'type(' in attrs else 'real','parameter':'parameter' in attrs.lower(),'default':pair[1].strip() if len(pair)>1 else None}
    return dict(contract,block=block,records=rec[last+1:],decls=decls,start=start,end=end)

ROUTINES={name:routine(name) for name in CONTRACTS}


def comments(text,indent=0):
    """Translate a Native comment region without losing equations or citations."""
    result=[]
    for raw in text.splitlines():
        if raw.lstrip().startswith('!') and not raw.lstrip().startswith('!$'):
            result.append(' '*indent+'//'+raw.split('!',1)[1])
        elif not raw.strip():result.append('')
    return result


def description(info):
    block=info['block'];start=block.index('! Description:');end=block.index('\n    use ',start)
    return comments(block[start:end],2)


def signature(info,definition=True):
    typ='real_t' if info['kind']=='function' else 'void'
    prefix='KOKKOS_INLINE_FUNCTION ' if info['elemental'] else ''
    names={a['name']:a for a in info['arguments']}
    native_head=info['block'].split('(',1)[1].split(')',1)[0]
    groups=[[a.strip() for a in line.replace('&','').split('!',1)[0].split(',') if a.strip()] for line in native_head.splitlines()]
    groups=[g for g in groups if g]
    def value(name):
        a=names[name]
        typ='int' if a['kind']=='integer' else 'implicit_coefs_terms&' if a['kind']=='implicit_coefs_terms' else 'Real'+str(len(a['dims']))+'D' if a['dims'] else 'real_t'
        if a['kind']=='real' and not a['dims'] and a['intent']!='in':typ+='&'
        return typ+' '+name
    chunks=[]
    for group in groups:
        intents={names[a]['intent'] for a in group}
        annotation='InOut' if intents=={'inout'} else 'Out' if intents=={'out'} else 'In' if intents=={'in'} else 'In / InOut / Out'
        chunks.append(', '.join(value(a) for a in group))
        if group is not groups[-1]:chunks[-1]+=', // '+annotation
    result=prefix+typ+' '+info['name']+'('+ '\n    '.join(chunks)+')'
    return result+(' {' if definition else ';')


def powers(expr):
    """Replace complete integer/real power operands, preserving surrounding groups."""
    while '**' in expr:
        index=expr.rfind('**');end_base=index
        while end_base and expr[end_base-1].isspace():end_base-=1
        if expr[end_base-1]==')':
            depth=1;begin=end_base-2
            while depth:
                if expr[begin]==')':depth+=1
                elif expr[begin]=='(':depth-=1
                begin-=1
            begin+=1
            prefix=re.search(r'[A-Za-z_]\w*$',expr[:begin])
            if prefix:begin=prefix.start()
        else:
            match=re.search(r'(?:[A-Za-z_]\w*|(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)$',expr[:end_base])
            if not match:raise ValueError(('power base',expr,index))
            begin=match.start()
        match=re.match(r'\s*([A-Za-z_]\w*|\d+)',expr[index+2:])
        if not match:raise ValueError(('power exponent',expr,index))
        exponent=match.group(1);end=index+2+match.end();base=expr[begin:end_base].strip()
        if exponent.isdigit() and int(exponent) in (2,3):
            replacement='('+' * '.join('('+base+')' for _ in range(int(exponent)))+')'
        else:replacement='Kokkos::pow('+base+', '+exponent+')'
        expr=expr[:begin]+replacement+expr[end:]
    return expr


def minmax(expr):
    """Use ordered second-operand ties at the established Native intrinsic boundary."""
    pattern=re.compile(r'\b(min|max)\s*\(')
    while (match:=pattern.search(expr)):
        depth=1;end=match.end()
        while depth:
            if expr[end]=='(':depth+=1
            elif expr[end]==')':depth-=1
            end+=1
        args=[minmax(arg) for arg in split_top(expr[match.end():end-1])]
        value=args[0];op='>' if match.group(1)=='max' else '<'
        for other in args[1:]:value='(('+value+') '+op+' ('+other+') ? ('+value+') : ('+other+'))'
        expr=expr[:match.start()]+value+expr[end:]
    return expr


def expr(text):
    """Translate source syntax without changing arithmetic association."""
    text=powers(text)
    text=re.sub(r'(?<![\w.])((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd\b',lambda m:'real_t('+m.group(1).replace('D','e').replace('d','e')+')',text)
    text=re.sub(r'(?<!:)\b(sqrt|exp|abs)\s*\(',lambda m:'Kokkos::'+m.group(1)+'(',text)
    text=text.replace('.and.','&&').replace('.or.','||').replace('.not.','!')
    text=text.replace('.true.','true').replace('.false.','false')
    return minmax(text)


def scalar_body(info):
    """Translate each elemental routine's complete statement/comment order."""
    result=[];indent=2
    argnames={a['name'] for a in info['arguments']}
    result_name=None
    if info['kind']=='function':result_name=re.search(r'result\s*\(\s*(\w+)\s*\)',info['block'],re.I).group(1)
    for name,decl in info['decls'].items():
        if name in argnames:continue
        if decl['dims']:raise ValueError(('elemental local array',name))
        line=('constexpr ' if decl['parameter'] else '')+('int' if decl['kind']=='integer' else 'real_t')+' '+name
        if decl['default'] is not None:line+=' = '+expr(decl['default'])
        result.append('  '+line+';')
    result.append('')
    for kind,value in info['records']:
        if kind=='blank':result.append('');continue
        if kind=='comment':result.append(' '*indent+'//'+value);continue
        if value.lower().startswith('end '):continue
        if value.lower()=='return':result.append(' '*indent+'return'+(' '+result_name if result_name else '')+';');continue
        if re.match(r'end\s*if',value,re.I):indent-=2;result.append(' '*indent+'}');continue
        if re.match(r'elseif',value,re.I):
            indent-=2;condition=re.sub(r'^elseif\s*\((.*)\)\s*then$',r'\1',value,flags=re.I);result.append(' '*indent+'} else if ('+expr(condition)+') {');indent+=2;continue
        if value.lower()=='else':indent-=2;result.append(' '*indent+'} else {');indent+=2;continue
        if re.match(r'if\s*\(',value,re.I):
            condition=re.sub(r'^if\s*\((.*)\)\s*then$',r'\1',value,flags=re.I);result.append(' '*indent+'if ('+expr(condition)+') {');indent+=2;continue
        if re.match(r'call\s',value,re.I):result.append(' '*indent+expr(re.sub(r'^call\s+','',value,flags=re.I))+';');continue
        if info['name']=='calc_F_w_zeta_w' and re.match(r'max_corr_w_x_sqd\s*=\s*max\(',value):
            operands=split_top(value[value.index('max(')+4:-1])
            result.append(' '*indent+'// Fold the Native variadic MAX in order, selecting its second operand on ties.')
            result.append(' '*indent+'max_corr_w_x_sqd = '+expr('max('+', '.join(operands[:2])+')')+';')
            for operand in operands[2:]:
                result.append(' '*indent+'max_corr_w_x_sqd = '+expr('max(max_corr_w_x_sqd, '+operand+')')+';')
            continue
        result.append(' '*indent+expr(value)+';')
    if indent!=2:raise ValueError(('scalar nesting',info['name'],indent))
    return result


def scalar_definition(info):
    result=['// =============================================================================',signature(info)]+description(info)+['','  // Arguments retain Native order; scalar reference outputs implement elemental Out/InOut.']
    for a in info['arguments']:
        note = None
        for line in info['block'].splitlines():
            if '!' not in line:continue
            code,comment=line.split('!',1)
            if re.match(r'^\s*'+a['name']+r'\s*[, &]*$',code) and comment.strip() not in ('In','Out','In/Out','I/O'):
                note=comment.strip()
        if note:result.append('  // '+a['name']+': '+note)
    result+=['']+scalar_body(info)+['}']
    return '\n'.join(result)+'\n'

lower_header='''// Hybrid two-component normal PDF formulas from new_hybrid_pdf.F90.
//
// calculate_mixture_fraction and calculate_w_params close the setter variable w.
// calculate_responder_params closes rt, thl, wind or a scalar using that mixture.
// calculate_coef_wp4_implicit and calc_coef_wp2xp_implicit give implicit moment
// coefficients; calc_coefs_wpxp2_semiimpl and calc_coefs_wpxpyp_semiimpl return
// semi-implicit coefficients and explicit terms. Inputs retain their physical
// units, while mixture/spread/width parameters are dimensionless.
//
// All seven Native elemental interfaces are device-callable scalar forms with
// required reference outputs. No state is retained. The detailed Native method,
// equations and limiting cases remain beside their original procedure groups.
// This isolated draft is uncompiled and untested. Reference: Griffin and Larson (2020).
#pragma once
#include "array_types.hpp"
#include "constants_clubb.hpp"
namespace clubb_kokkos {
'''
lower=[lower_header];previous=TEXTS['new_hybrid_pdf'].index('  contains')+len('  contains')
for name,info in ROUTINES.items():
    if info['module']!='new_hybrid_pdf':continue
    gap=TEXTS['new_hybrid_pdf'][previous:info['start']]
    lower+=['\n'.join(comments(gap)),scalar_definition(info)];previous=info['end']
lower.append('} // namespace clubb_kokkos\n')
(HERE/'new_hybrid_pdf.hpp').write_text('\n'.join(lower))

main=ROUTINES['new_hybrid_pdf_driver']
main_header='''// Close hybrid PDF parameters and implicit moment terms for column/level fields.
//
// new_hybrid_pdf_driver computes the vertical-velocity mixture, then rt, thl,
// horizontal-wind and passive-scalar component means/variances. Responder
// skewnesses may be clipped. The caller-owned implicit_coefs_terms receives
// coefficients/explicit terms selected by the three Native advection flags.
// Model Views use column-first LayoutLeft; all scratch belongs to the call.
// This isolated draft has not been compiled, integrated or numerically tested.
#pragma once
#include "array_types.hpp"
#include "pdf_parameter_module.hpp"
namespace clubb_kokkos {
'''
(HERE/'new_hybrid_pdf_main.hpp').write_text(main_header+'\n'+signature(main,False)+'\n\n} // namespace clubb_kokkos\n')
# Main body translation is prepared separately so its array/domain adaptations
# remain explicit and reviewable rather than hidden in elemental wrappers.
(HERE/'main_scalar_helpers.tmp').write_text('\n'.join(scalar_definition(ROUTINES[n]) for n in ('calc_responder_driver','calc_F_w_zeta_w')))
print('Prepared seven scalar lower routines and two private scalar helper bodies; parent pending.')


def parse_nodes(rec,start=0,terminators=()):
    """Keep Native blocks while distinguishing host loops and profile kernels."""
    nodes=[];index=start
    while index<len(rec):
        kind,value=rec[index]
        if kind!='code':nodes.append((kind,value));index+=1;continue
        low=value.lower()
        if any(re.match(pattern,low) for pattern in terminators):return nodes,index
        if re.match(r'do\s+\w+\s*=',low):
            match=re.match(r'do\s+(\w+)\s*=\s*(.*)',value,re.I);bounds=split_top(match.group(2))
            body,end=parse_nodes(rec,index+1,(r'end\s*do',r'enddo'))
            nodes.append(('do',match.group(1),bounds,body));index=end+1;continue
        if re.match(r'if\s*\(',low):
            branches=[];condition=re.sub(r'^if\s*\((.*)\)\s*then$',r'\1',value,flags=re.I)
            body,end=parse_nodes(rec,index+1,(r'elseif',r'else$',r'end\s*if',r'endif'))
            branches.append((condition,body))
            while end<len(rec) and rec[end][0]=='code' and rec[end][1].lower().startswith('elseif'):
                condition=re.sub(r'^elseif\s*\((.*)\)\s*then$',r'\1',rec[end][1],flags=re.I)
                body,end=parse_nodes(rec,end+1,(r'elseif',r'else$',r'end\s*if',r'endif'))
                branches.append((condition,body))
            otherwise=None
            if end<len(rec) and rec[end][1].lower()=='else':otherwise,end=parse_nodes(rec,end+1,(r'end\s*if',r'endif'))
            nodes.append(('if',branches,otherwise));index=end+1;continue
        nodes.append(('code',value));index+=1
    return nodes,index

META=dict(main['decls'])
FIELDS=[]
for match in re.finditer(r'pdf_implicit_coefs_terms%(\w+)\((i,:[^)]*)\)',main['block']):
    name=match.group(1)
    if name not in FIELDS:FIELDS.append(name)
    META['pdf_'+name]={'dims':['ngrdcol','nz']+(['sclr_dim'] if match.group(2).count(':')==2 else [])}


def point(text):
    """Map elemental whole-profile operands to the current device cell."""
    text=re.sub(r'pdf_implicit_coefs_terms%(\w+)',lambda m:'pdf_'+m.group(1),text)
    result='';index=0
    while index<len(text):
        match=re.match(r'[A-Za-z_]\w*',text[index:])
        if not match:result+=text[index];index+=1;continue
        name=match.group();end=index+len(name);cursor=end
        while cursor<len(text) and text[cursor].isspace():cursor+=1
        if cursor<len(text) and text[cursor]=='(':
            depth=1;finish=cursor+1
            while depth:
                if text[finish]=='(':depth+=1
                elif text[finish]==')':depth-=1
                finish+=1
            args=split_top(text[cursor+1:finish-1])
            if name in META and META[name]['dims']:
                replacements=[]
                for axis,arg in enumerate(args):
                    if arg==':':
                        dim=META[name]['dims'][axis].strip()
                        replacements.append('k' if dim=='nz' else 'sclr' if dim=='sclr_dim' else 'i')
                    else:replacements.append(point(arg))
                result+=name+'('+', '.join(replacements)+')'
            else:result+=name+'('+', '.join(point(arg) for arg in args)+')'
            index=finish
        else:
            if name in META and META[name]['dims']:
                dims=META[name]['dims'];indices=['k' if d.strip()=='nz' else 'sclr' if d.strip()=='sclr_dim' else 'i' for d in dims]
                result+=name+'('+', '.join(indices)+')'
            else:result+=name
            index=end
    return result


def domain(code):
    """Find the parallel profile/species extent of a Native array assignment."""
    if '=' not in code:return None
    lhs=code.split('=',1)[0].strip()
    lhs=re.sub(r'pdf_implicit_coefs_terms%(\w+)',lambda m:'pdf_'+m.group(1),lhs)
    match=re.match(r'(\w+)(?:\((.*)\))?$',lhs)
    if not match or match.group(1) not in META or not META[match.group(1)]['dims']:return None
    dims=META[match.group(1)]['dims']
    if match.group(2) is None:return tuple(d.strip() for d in dims)
    args=split_top(match.group(2));return tuple(dims[j].strip() for j,arg in enumerate(args) if arg==':') or None

SEQUENCE=0


def first_destination(nodes):
    """Name a profile kernel after its first Native scientific destination."""
    for node in nodes:
        if node[0]=='code':
            match=re.match(r'(\w+(?:%\w+)?)\s*(?:\([^)]*\))?\s*=(?!=)',node[1])
            if match:return match.group(1).replace('%','_')
        if node[0]=='do':
            name=first_destination(node[3])
            if name:return name
        if node[0]=='if':
            for _,branch in node[1]:
                name=first_destination(branch)
                if name:return name
    return 'profile_values'


def kernel(lines,indent,body,axes,label):
    """Emit one explicit Kokkos launch for an independent Native array phase."""
    extent=' * '.join(axes)
    lines+=[' '*indent+'Kokkos::parallel_for("new_hybrid_pdf_driver::'+label+'", '+extent+',', ' '*indent+'  KOKKOS_LAMBDA(const int cell) {']
    inner=indent+4
    if axes==('nz',):lines.append(' '*inner+'const int k = cell;')
    elif axes==('nz','sclr_dim'):
        lines+=[' '*inner+'const int k = cell % nz;',' '*inner+'const int sclr = cell / nz;']
    else:raise ValueError(('kernel axes',axes))
    lines+=[' '*inner+line for line in body]
    lines+=[' '*indent+'  });']


def emit_nodes(nodes,indent=2,device=False):
    """Retain control flow/calls; map independent whole-array phases to kernels."""
    global SEQUENCE
    lines=[];index=0
    while index<len(nodes):
        node=nodes[index];kind=node[0]
        if kind=='comment':lines.append(' '*indent+'//'+node[1]);index+=1;continue
        if kind=='blank':lines.append('');index+=1;continue
        if kind=='do':
            _,var,bounds,body=node
            if bounds[0]!='1' or (len(bounds)==3 and bounds[2]!='1'):raise ValueError(('loop bounds',node))
            if var=='k' and not device:
                SEQUENCE+=1;bodylines=emit_nodes(body,0,True)
                # The Native profile loop writes each level independently.
                kernel(lines,indent,bodylines,('nz',),'profile_'+first_destination(body))
            else:
                lines.append(' '*indent+'for (int '+var+' = 0; '+var+' < '+bounds[1]+'; ++'+var+') {')
                lines+=emit_nodes(body,indent+2,device)
                lines.append(' '*indent+'}')
            index+=1;continue
        if kind=='if':
            _,branches,otherwise=node
            for branch,(condition,body) in enumerate(branches):
                if 'any(' in condition.lower():
                    if branch or device:raise ValueError(('ANY location',condition))
                    inside=re.search(r'any\s*\((.*)\)',condition,re.I).group(1)
                    lines.append(' '*indent+'int invalid_mixture = 0;')
                    lines+=[' '*indent+'Kokkos::parallel_reduce("new_hybrid_pdf_driver::mixture_check", nz,',' '*indent+'  KOKKOS_LAMBDA(const int k, int& count) {',' '*(indent+4)+'if ('+expr(point(inside))+') ++count;',' '*indent+'  }, invalid_mixture);']
                    lines.append(' '*indent+'if (invalid_mixture) {')
                else:lines.append(' '*indent+('if (' if branch==0 else '} else if (')+expr(point(condition) if device else condition)+') {')
                lines+=emit_nodes(body,indent+2,device)
            if otherwise is not None:
                lines.append(' '*indent+'} else {');lines+=emit_nodes(otherwise,indent+2,device)
            lines.append(' '*indent+'}');index+=1;continue
        code=node[1]
        if re.match(r'end\s+(subroutine|function)\b',code,re.I):index+=1;continue
        if code.lower()=='return':lines.append(' '*indent+'return;');index+=1;continue
        if code.lower()=='error stop':
            lines.append(' '*indent+'// Native ERROR STOP terminates this driver; the host exception reaches the runner.')
            lines.append(' '*indent+'throw std::runtime_error("new_hybrid_pdf_driver: invalid mixture fraction");');index+=1;continue
        if code.lower().startswith('write('):
            strings=re.findall(r'"([^\"]*)"',code)
            lines.append(' '*indent+'std::cerr << '+json.dumps(''.join(strings)+'\n')+';');index+=1;continue
        if code.lower().startswith('call '):
            call=re.sub(r'^call\s+','',code,flags=re.I)
            if call.startswith('compute_gamma_Skw'):
                lines.append(' '*indent+call+';')
            elif device:lines.append(' '*indent+expr(point(call))+';')
            else:
                SEQUENCE+=1
                first_arg=split_top(call[call.index('(')+1:-1])[0].split('(',1)[0]
                kernel(lines,indent,[expr(point(call))+';'],('nz',),call.split('(',1)[0]+'_'+first_arg)
            index+=1;continue
        if device and re.match(r'max_corr_w_sclr_sqd\(k\)\s*=\s*max\(',code):
            operand=split_top(code[code.index('max(')+4:-1])[0]
            lines.append(' '*indent+'// Evaluate this Native MAX operand once before applying its ordered comparison.')
            lines.append(' '*indent+'const real_t corr_w_sclr_sqd = '+expr(point(operand))+';')
            lines.append(' '*indent+'max_corr_w_sclr_sqd(k) = corr_w_sclr_sqd > max_corr_w_sclr_sqd(k) ? corr_w_sclr_sqd : max_corr_w_sclr_sqd(k);')
            index+=1;continue
        axes=domain(code) if not device else None
        if axes:
            group=[]
            while index<len(nodes) and nodes[index][0]=='code' and domain(nodes[index][1])==axes:
                group.append(expr(point(nodes[index][1]))+';');index+=1
            SEQUENCE+=1;kernel(lines,indent,group,axes,'set_'+code.split('=',1)[0].strip().split('(',1)[0].replace('%','_'));continue
        lines.append(' '*indent+expr(point(code) if device else code)+';');index+=1
    return lines

main_cpp=['''// Hybrid PDF orchestration; see the public header and draft integration notes.
#include "new_hybrid_pdf_main.hpp"
#include "new_hybrid_pdf.hpp"
#include "Skx_module.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp"
#include <iostream>
#include <stdexcept>
namespace clubb_kokkos {
namespace {
''']
for name in ('calc_responder_driver','calc_F_w_zeta_w'):main_cpp.append(signature(ROUTINES[name],False))
main_cpp+=['} // namespace','',signature(main)]+description(main)+['','  // Arguments']
for arg in main['arguments']:
    note=None
    for line in main['block'].splitlines():
        if '!' not in line:continue
        code,comment=line.split('!',1)
        if re.match(r'^\s*'+arg['name']+r'\s*[, &]*$',code) and comment.strip() not in ('In','Out','In/Out','I/O'):
            note=comment.strip()
    if note:main_cpp.append('  // '+arg['name']+': '+note)
main_cpp+=['','  //------------------------ Local Variables ------------------------']
argnames={a['name'] for a in main['arguments']}
for name,decl in main['decls'].items():
    if name in argnames or decl['kind']=='integer':continue
    dims=decl['dims']
    if not dims:raise ValueError(('driver local scalar',name))
    note=None
    for line in main['block'].splitlines():
        if '!' not in line:continue
        code,comment=line.split('!',1)
        if re.match(r'^\s*'+name+r'\s*[, &]*$',code):note=comment.strip()
    main_cpp.append('  Real'+str(len(dims))+'D '+name+'("new_hybrid_pdf_driver::'+name+'", '+', '.join(dims)+');'+(' // '+note if note else ''))
main_cpp+=['','  // Capture individual component Views, retaining the caller-owned UDT allocations.']
for field in FIELDS:main_cpp.append('  const auto pdf_'+field+' = pdf_implicit_coefs_terms.'+field+';')
main_cpp+=['','  //------------------------ Begin Code ------------------------', '  // Each independent profile phase uses a Kokkos kernel; coupled work within', '  // a cell and the Native column/call/branch order remain unchanged.']
nodes,_=parse_nodes(main['records'])
main_cpp+=emit_nodes(nodes)+['}','', 'namespace {', (HERE/'main_scalar_helpers.tmp').read_text(),'} // namespace','} // namespace clubb_kokkos','']
(HERE/'new_hybrid_pdf_main.cpp').write_text('\n'.join(main_cpp))
(HERE/'main_scalar_helpers.tmp').unlink()
print('Prepared complete56-argument parent with source scratch, control flow and device scalar calls.')
