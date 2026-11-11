"""Prepare isolated Native-named original PDF drafts; no builds or physics runs.

All writes stay in this draft directory. Source contracts and grouping remain
canonical; the generated draft needs separate source review and runtime checks.
"""
from pathlib import Path
import re,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUTLINE=json.loads((HERE/'source_outline.json').read_text())
CONTRACTS={r['name']:r for r in OUTLINE['routines']}
TEXTS={name:(ROOT/'src/CLUBB_core'/f'{name}.F90').read_text() for name in ('new_pdf_main','new_pdf')}


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
            decls[var]={'dims':split_top(dm.group(1)) if dm else [],'kind':'integer' if attrs.strip().startswith('integer') else 'logical' if attrs.strip().startswith('logical') else 'complex' if attrs.strip().startswith('complex') else 'implicit_coefs_terms' if 'type(' in attrs else 'real','parameter':'parameter' in attrs.lower(),'default':pair[1].strip() if len(pair)>1 else None}
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



def result_name(info):
    return re.search(r'result\s*\(\s*(\w+)\s*\)',info['block'],re.I).group(1) if info['kind']=='function' else None


def arg_notes(info):
    notes={}
    for name in info['decls']:
        for line in info['block'].splitlines():
            if '!' not in line:continue
            code,note=line.split('!',1)
            if re.match(r'^\s*'+name+r'\s*[, &]*$',code) and note.strip() not in ('In','Out','In/Out','In/Out.','I/O') and note.strip():notes[name]=note.strip()
    return notes


def scalar_signature(info,definition=True):
    args=[]
    for a in info['arguments']:
        if a['name']=='nz':continue
        typ='Kokkos::Array<real_t, 3>' if a['name']=='roots' else 'real_t'
        if a['intent']!='in':typ+='&'
        name='lambda_' if a['name']=='lambda' else a['name']
        args.append(typ+' '+name)
    typ='Kokkos::Array<real_t, 3>' if info['name']=='sort_roots' else 'real_t' if info['kind']=='function' else 'void'
    return 'KOKKOS_INLINE_FUNCTION '+typ+' '+info['name']+'('+', '.join(args)+')'+(' {' if definition else ';')


def profile_signature(info,definition=True):
    args=[];templates=[]
    for a in info['arguments']:
        if a['dims']:
            typ=a['name']+'View';templates.append('class '+typ)
        else:typ='int' if a['kind']=='integer' else 'real_t'
        name='lambda_' if a['name']=='lambda' else a['name']
        args.append(typ+' '+name)
    result = ''
    if info['kind']=='function':
        first = next(a['name'] for a in info['arguments'] if a['dims'])
        result = ' -> Kokkos::View<real_t'+('*'*len(info['result']['dims']))+', Kokkos::LayoutLeft, typename '+first+'View::memory_space>'
    return 'template <'+', '.join(templates)+'>\n'+('auto' if info['kind']=='function' else 'void')+' '+info['name']+'('+', '.join(args)+')'+result+(' {' if definition else ';')


def rewrite_scalar(code,info):
    code=re.sub(r'\blambda\b','lambda_',code)
    if info['name']=='sort_roots':
        code=re.sub(r'(roots(?:_sorted)?)\(:,\s*([123])\)',lambda m:m.group(1)+'['+str(int(m.group(2))-1)+']',code)
    code=expr(code)
    if re.search(r'\bcalc_mixture_fraction\s*\(',code):code=re.sub(r'(calc_mixture_fraction\s*\()\s*nz\s*,\s*',r'\1',code)
    return code


def cell_body(info,skip_mixture=False):
    lines=[];indent=2;args={a['name'] for a in info['arguments']};res=result_name(info);notes=arg_notes(info)
    for name,decl in info['decls'].items():
        if name in args:continue
        typ='bool' if decl['kind']=='logical' else 'Kokkos::Array<real_t, 3>' if name=='roots_sorted' else 'real_t'
        lines.append('  '+typ+' '+name+';'+(' // '+notes[name] if name in notes else ''))
    lines.append('')
    for kind,value in info['records']:
        low=value.lower()
        if kind=='comment':lines.append(' '*indent+'//'+value);continue
        if kind=='blank':lines.append('');continue
        if low.startswith('end subroutine') or low.startswith('end function'):continue
        if skip_mixture and value.startswith('mixt_frac = calc_mixture_fraction'):continue
        if re.match(r'if\s*\(\s*any',low):
            lines+=['  if (!l_calc_mixt_frac) {','    // Native ERROR STOP is fatal; device scalar calls use the Kokkos abort boundary.','    Kokkos::abort("Mixture fraction cannot be calculated.\\nThe value of F_x must be greater than 0 when | Skx | > 0.");','  }'];break
        if low=='return':lines.append(' '*indent+'return'+(' '+res if res else '')+';');continue
        if re.match(r'end\s*where',low):indent-=2;lines.append(' '*indent+'}');continue
        if low.startswith('elsewhere'):
            indent-=2;m=re.match(r'elsewhere\s*\((.*)\)',value,re.I)
            lines.append(' '*indent+'} else'+(' if ('+rewrite_scalar(m.group(1),info)+')' if m else '')+' {');indent+=2;continue
        if low.startswith('where'):
            m=re.match(r'where\s*\((.*)\)\s*$',value,re.I);lines.append(' '*indent+'if ('+rewrite_scalar(m.group(1),info)+') {');indent+=2;continue
        lines.append(' '*indent+rewrite_scalar(value,info)+';')
    if info['name']=='calc_mixture_fraction':lines.append('  return mixt_frac;')
    if indent!=2:raise ValueError((info['name'],indent))
    return lines


def scalar_definition(info):
    return '\n'.join([scalar_signature(info)]+description(info)+['','  // Independent Native profile cell; nz is represented by the enclosing View extent.']+['  // '+a['name']+': '+arg_notes(info)[a['name']] for a in info['arguments'] if a['name'] in arg_notes(info)]+['']+cell_body(info)+['}',''])


def stage_profile(info):
    lines=[]
    for a in info['arguments']:
        if a['dims']:lines.append('  auto '+a['name']+'_device = to_execution_space('+a['name']+(', false' if a['intent']=='out' else '')+');')
    return lines


def finish_profile(info):
    return ['  from_execution_space('+a['name']+', '+a['name']+'_device);' for a in info['arguments'] if a['dims'] and a['intent']!='in']


def profile_return(info,result='result'):
    first=next(a['name'] for a in info['arguments'] if a['dims']);rank=len(info['result']['dims']);alias='Real'+str(rank)+'D';dims='nz, 3' if rank==2 else 'nz'
    return ['  if constexpr (std::is_same_v<typename '+first+'View::memory_space, typename '+alias+'::memory_space>)','    return '+result+';','  else {','    Kokkos::View<real_t'+('*'*rank)+', Kokkos::LayoutLeft, typename '+first+'View::memory_space> output(','        Kokkos::view_alloc(Kokkos::WithoutInitializing, "'+info['name']+'"), '+dims+');','    from_execution_space(output, '+result+');','    return output;','  }']


def profile_description(info):
    description_lines=description(info)
    result=[]
    for line in description_lines:
        if result and line.strip() in ('', '//'):break
        result.append(line)
    return result+['  // Profile adaptation of the device scalar formula and Native method above.']


def profile_definition(info):
    name=info['name'];res=result_name(info);rank=len(info['result']['dims']) if res else 0
    lines=[profile_signature(info)]+profile_description(info)+['','  // Preserve Native profile argument order; stage Views at the public memory boundary.']
    lines+=stage_profile(info)
    if res:
        lines.append('  Real'+str(rank)+'D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "'+name+'"), nz'+(', 3' if rank==2 else '')+');')
    if name=='calc_setter_var_params':
        lines+=['  // Preserve the Native mixture validity check before computing any component outputs.','  const auto mixture = calc_mixture_fraction(nz, Skx, F_x, zeta_x, sgn_wpxp);','  const auto mixture_device = to_execution_space(mixture);']
    if name=='calc_mixture_fraction':
        lines+=['  int invalid_mixture = 0;','  Kokkos::parallel_reduce("calc_mixture_fraction::validity", nz,','      KOKKOS_LAMBDA(const int k, int& count) {','        if (!(F_x_device(k) > zero) && Kokkos::abs(Skx_device(k)) > zero) ++count;','      }, invalid_mixture);','  if (invalid_mixture) {','    std::cerr << "Mixture fraction cannot be calculated.\\n";','    std::cerr << "The value of F_x must be greater than 0 when | Skx | > 0.\\n";','    throw std::runtime_error("calc_mixture_fraction: invalid mixture fraction");','  }']
    lines+=['  Kokkos::parallel_for("'+name+'", nz, KOKKOS_LAMBDA(const int k) {']
    actual=[]
    for a in info['arguments']:
        if a['name']=='nz':continue
        actual.append(a['name']+'_device(k)' if a['dims'] else 'lambda_' if a['name']=='lambda' else a['name'])
    if name=='sort_roots':
        lines+=['    Kokkos::Array<real_t, 3> cell_roots = {roots_device(k, 0), roots_device(k, 1), roots_device(k, 2)};','    const auto sorted = new_pdf_detail::sort_roots(cell_roots);','    result(k, 0) = sorted[0];','    result(k, 1) = sorted[1];','    result(k, 2) = sorted[2];']
    elif name=='calc_setter_var_params':
        # The scalar body is emitted directly so vector validation remains before all outputs.
        argmap={a['name']:a['name']+'_device(k)' for a in info['arguments'] if a['dims']}
        body=cell_body(info,skip_mixture=True)
        # Out arguments are Views here; all local declarations from scalar body are absent.
        for raw in body:
            if raw.strip()=='return;':continue
            if not raw.lstrip().startswith('//'):
                for old,new in sorted(argmap.items(),key=lambda item:-len(item[0])):raw=re.sub(r'\b'+old+r'\b',new,raw)
            lines.append('  '+raw)
        # Assign mixture output before the other source statements, as the Native call does.
        insert=lines.index('  Kokkos::parallel_for("'+name+'", nz, KOKKOS_LAMBDA(const int k) {')+1
        lines.insert(insert,'    mixt_frac_device(k) = mixture_device(k);')
    else:
        lines.append('    '+('result(k) = ' if res else '')+name+'('+', '.join(actual)+');')
    lines+=['  });']+finish_profile(info)
    if res:lines+=profile_return(info)
    return '\n'.join(lines+['}',''])


def limits_definition(info):
    lines=[profile_signature(info)]+description(info)+['','  // Input dimensionless mixture, skewness/sign and squared-skewness bounds;', '  // write dimensionless min_F_x/max_F_x profiles of nz levels.', '  // Only the unported complex root backend crosses to owned HostSpace profiles.','  // All coefficient, selection and sorting arithmetic otherwise uses device Views.']+stage_profile(info)
    local=['coef_A_1','coef_B_1','coef_C_1','coef_D_1','coef_A_2','coef_B_2','coef_C_2','coef_D_2','min_sqrt_F_x','max_sqrt_F_x']
    for name in local:lines.append('  Real1D '+name+'("calc_limits_F_x_responder::'+name+'", nz);')
    lines+=['  Real2D sqrt_F_x_roots_1("calc_limits_F_x_responder::sqrt_F_x_roots_1", nz, 3);','  Real2D sqrt_F_x_roots_2("calc_limits_F_x_responder::sqrt_F_x_roots_2", nz, 3);']
    # Preserve first/second equation coefficient/solve/sort order. Each complex
    # result is copied to its real-valued View without changing root order.
    for comp in (1,2):
        one='one' if comp==1 else 'two';prefix='coef_';c=str(comp)
        lines+=['  // Set up the coefficients in the equation for the limit of sqrt(F_x) based', '  // on the '+('1st' if comp==1 else '2nd')+' PDF component standard deviation (sigma_x_'+c+') being greater than', '  // or equal to 0. This equation has the form:', '  // A * sqrt(F_x)^3 + B * sqrt(F_x)^2 + C * sqrt(F_x) + D = 0.','  Kokkos::parallel_for("calc_limits_F_x_responder::coefficients_'+c+'", nz, KOKKOS_LAMBDA(const int k) {','    coef_A_'+c+'(k) = -('+('one + mixt_frac_device(k)' if comp==1 else 'two - mixt_frac_device(k)')+');','    coef_B_'+c+'(k) = zero;','    coef_C_'+c+'(k) = three * '+('mixt_frac_device(k)' if comp==1 else '(one - mixt_frac_device(k))')+';','    coef_D_'+c+'(k) = '+('' if comp==1 else '-')+'Kokkos::sqrt(mixt_frac_device(k) * (one - mixt_frac_device(k))) * Skx_device(k) * sgn_wpxp_device(k);','  });','  // TODO: replace this transient public f2k call with calc_roots when promoted.','  const auto native_roots_'+c+' = f2k::cubic_solve(nz, to_host_space(coef_A_'+c+'), to_host_space(coef_B_'+c+'),','      to_host_space(coef_C_'+c+'), to_host_space(coef_D_'+c+'));','  auto real_roots_'+c+' = Kokkos::create_mirror_view(sqrt_F_x_roots_'+c+');','  for (int k = 0; k < nz; ++k)','    for (int root = 0; root < 3; ++root)','      real_roots_'+c+'(k, root) = native_roots_'+c+'(k, root).real();','  Kokkos::deep_copy(sqrt_F_x_roots_'+c+', real_roots_'+c+');','  // Sort the values of the roots (values of sqrt(F_x)) from smallest to', '  // largest. Ignore any complex component of the roots. The code below that', '  // uses sqrt_F_x_roots_'+c+'_sorted already factors the appropriate roots to use', '  // into account; the original first root remains separate for those branches.','  const auto sqrt_F_x_roots_'+c+'_sorted = new_pdf_detail::sort_roots(nz, sqrt_F_x_roots_'+c+');']
    # Translate only the source selection block; original unsorted roots become real Views.
    records=info['records'];begin=next(i for i,(k,v) in enumerate(records) if k=='code' and v.startswith('where ( Skx *'))
    lines+=['  // Find the minimum and maximum allowable values of sqrt(F_x) based on Skx', "  // and sgn( <w'x'> ).",'  Kokkos::parallel_for("calc_limits_F_x_responder::limits", nz, KOKKOS_LAMBDA(const int k) {']
    indent=4
    for kind,value in records[begin:]:
        low=value.lower()
        if kind=='comment':lines.append(' '*indent+'//'+value);continue
        if kind=='blank':lines.append('');continue
        if low in ('return',) or low.startswith('end subroutine'):continue
        if low.startswith('where'):
            v=re.match(r'where\s*\((.*)\)',value,re.I).group(1);lines.append(' '*indent+'if ('+limits_point(v,info)+') {');indent+=2
        elif low.startswith('elsewhere'):indent-=2;lines.append(' '*indent+'} else {');indent+=2
        elif low.startswith('endwhere'):indent-=2;lines.append(' '*indent+'}')
        else:lines.append(' '*indent+limits_point(value,info)+';')
    if indent!=4:raise ValueError(('limits',indent))
    lines+=['  });']+finish_profile(info)+['}',''];return '\n'.join(lines)


def limits_point(code,info):
    code=re.sub(r'real\(\s*(sqrt_F_x_roots_[12])\(:,([123])\),\s*kind\s*=\s*core_rknd\s*\)',lambda m:m.group(1)+'(k, '+str(int(m.group(2))-1)+')',code)
    code=re.sub(r'(sqrt_F_x_roots_[12]_sorted)\(:,([123])\)',lambda m:m.group(1)+'(k, '+str(int(m.group(2))-1)+')',code)
    for name in sorted([a['name'] for a in info['arguments'] if a['dims']],key=len,reverse=True):code=re.sub(r'\b'+name+r'\b',name+'_device(k)',code)
    for name in ('min_sqrt_F_x','max_sqrt_F_x'):code=re.sub(r'\b'+name+r'\b',name+'(k)',code)
    return expr(code)

HEADER='''// Original two-component normal PDF formulas for w, rt and thl.
//
// calc_mixture_fraction and calc_setter_var_params close the setter; responding
// variables use calc_responder_params and calc_limits_F_x_responder. The four
// calc_coef_wp4_implicit and calc_coef_wpxp2_implicit return implicit terms;
// calc_coefs_wp2xp_semiimpl and calc_coefs_wpxpyp_semiimpl return coefficient/
// explicit-term pairs.
// Public profile APIs retain Native nz and argument order. Independent-cell
// scalar overloads are device callable; the limits profile currently stages only
// its two complex root calls through the transient backend. Caller-owned outputs
// and returned fields retain the caller memory space. All scratch is call local.
// Reference: Griffin and Larson (2018). Isolated uncompiled/untested draft.
#pragma once
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "f2k_dependencies.hpp"
#include <Kokkos_Array.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <iostream>
#include <stdexcept>
#include <type_traits>
namespace clubb_kokkos {
namespace new_pdf_detail {
'''
lower=[HEADER,scalar_signature(ROUTINES['sort_roots'],False),profile_signature(ROUTINES['sort_roots'],False),'} // namespace new_pdf_detail\n']
previous=TEXTS['new_pdf'].index('  contains')+len('  contains')
for name,info in ROUTINES.items():
    if info['module']!='new_pdf':continue
    lower.append('\n'.join(comments(TEXTS['new_pdf'][previous:info['start']])));previous=info['end']
    if name=='sort_roots':lower.append('namespace new_pdf_detail {\n')
    if name=='calc_limits_F_x_responder':lower.append(limits_definition(info))
    else:lower.extend([scalar_definition(info),profile_definition(info)])
    if name=='sort_roots':lower.append('} // namespace new_pdf_detail\n')
lower.append('} // namespace clubb_kokkos\n')
(HERE/'new_pdf.hpp').write_text('\n'.join(lower))
(HERE/'f2k_dependencies.hpp').write_text('''// Pending public complex-root adapter; declaration only, with no fallback.
// Root owns promotion into src/f2k or replacement by the calc_roots C++ port.
#pragma once
#include "array_types.hpp"
#include <complex>
namespace f2k {
using HostComplex2D = Kokkos::View<std::complex<clubb_kokkos::real_t>**, Kokkos::LayoutLeft,
                                   Kokkos::HostSpace>;
using HostReal1D = Kokkos::View<clubb_kokkos::real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
HostComplex2D cubic_solve(int nz, HostReal1D a_coef, HostReal1D b_coef,
                          HostReal1D c_coef, HostReal1D d_coef);
} // namespace f2k
''')
print('Prepared all nine Native lower procedures and justified scalar/profile interfaces.')


def parse_nodes(rec,start=0,stops=()):
    nodes=[];index=start
    while index<len(rec):
        kind,value=rec[index];low=value.lower()
        if kind!='code':nodes.append((kind,value));index+=1;continue
        if any(re.match(p,low) for p in stops):return nodes,index
        if low.startswith('where') or re.match(r'if\s*\(',low):
            wh=low.startswith('where');term=(r'elsewhere',r'endwhere') if wh else (r'elseif',r'else$',r'end\s*if',r'endif')
            condition=re.match(r'(?:where|if)\s*\((.*)\)\s*(?:then)?$',value,re.I).group(1)
            body,end=parse_nodes(rec,index+1,term);branches=[(condition,body)];otherwise=None
            while end<len(rec) and rec[end][0]=='code' and re.match(r'(?:elsewhere|elseif)\s*\(',rec[end][1],re.I):
                condition=re.match(r'(?:elsewhere|elseif)\s*\((.*)\)\s*(?:then)?$',rec[end][1],re.I).group(1)
                body,end=parse_nodes(rec,end+1,term);branches.append((condition,body))
            if end<len(rec) and rec[end][1].lower() in ('elsewhere','else'):
                otherwise,end=parse_nodes(rec,end+1,(r'endwhere',) if wh else (r'end\s*if',r'endif'))
            nodes.append(('where' if wh else 'if',branches,otherwise));index=end+1;continue
        if re.match(r'do\s+\w+\s*=',low):
            m=re.match(r'do\s+(\w+)\s*=\s*(.*)',value,re.I);body,end=parse_nodes(rec,index+1,(r'end\s*do',r'enddo'))
            nodes.append(('do',m.group(1),split_top(m.group(2)),body));index=end+1;continue
        nodes.append(('code',value));index+=1
    return nodes,index


def main_signature(info,definition=True):
    args=[]
    for a in info['arguments']:
        typ='int' if a['kind']=='integer' else 'implicit_coefs_terms&' if a['kind']=='implicit_coefs_terms' else 'Real'+str(len(a['dims']))+'D' if a['dims'] else 'real_t'
        args.append(typ+' '+a['name'])
    return 'void '+info['name']+'('+', '.join(args)+')'+(' {' if definition else ';')


main=ROUTINES['new_pdf_driver'];FIELDS=[]
for m in re.finditer(r'pdf_implicit_coefs_terms%(\w+)\(i,:\)',main['block']):
    if m.group(1) not in FIELDS:FIELDS.append(m.group(1))


def point(code,info,profile=False):
    dims={n:d['dims'] for n,d in info['decls'].items() if d['dims']}
    code=re.sub(r'pdf_implicit_coefs_terms%(\w+)',lambda m:'pdf_'+m.group(1),code)
    dims.update({'pdf_'+n:['ngrdcol','nz'] for n in FIELDS})
    # Handle whole-column references before bare names; protected placeholders
    # prevent double indexing when mapping scientific field names.
    protected={}
    def keep(value):
        name='NATIVEPLACEHOLDER'+str(len(protected));protected[name]=value;return name
    for name,shape in sorted(dims.items(),key=lambda pair:-len(pair[0])):
        suffix='_device' if profile else ''
        code=re.sub(r'\b'+name+r'\s*\(i,\s*:\)',lambda _:keep(name+suffix+'(i, k)'),code)
        code=re.sub(r'\b'+name+r'\b(?!\s*\()',lambda _:keep(name+suffix+'('+('k' if len(shape)==1 else 'i, k')+')'),code)
    code=re.sub(r'\blambda\b','lambda_',code)
    for token,value in protected.items():code=code.replace(token,value)
    return expr(code)


def view_call(code,info,profile=False):
    code=re.sub(r'pdf_implicit_coefs_terms%(\w+)',lambda m:'pdf_'+m.group(1),code)
    if profile:
        for name,d in sorted(info['decls'].items(),key=lambda pair:-len(pair[0])):
            if d['dims']:code=re.sub(r'\b'+name+r'\b',name+'_device',code)
    else:
        code=re.sub(r'\b(\w+)\(i,\s*:\)',lambda m:'Kokkos::subview('+m.group(1)+', i, Kokkos::ALL())',code)
        # Native passes four tunable values as scalars to profile helpers. Copy
        # the tunable table once at this host orchestration boundary.
        code=code.replace('clubb_params(i,','clubb_params_host(i,')
    return code


def device_nodes(nodes,info,profile=False,indent=0):
    lines=[]
    for node in nodes:
        kind=node[0]
        if kind=='comment':lines.append(' '*indent+'//'+node[1]);continue
        if kind=='blank':lines.append('');continue
        if kind in ('where','if'):
            for j,(condition,body) in enumerate(node[1]):
                lines.append(' '*indent+('if (' if j==0 else '} else if (')+point(condition,info,profile)+') {');lines+=device_nodes(body,info,profile,indent+2)
            if node[2] is not None:lines.append(' '*indent+'} else {');lines+=device_nodes(node[2],info,profile,indent+2)
            lines.append(' '*indent+'}');continue
        if kind!='code':raise ValueError(('device node',node))
        if node[1].lower()=='return':continue
        lines.append(' '*indent+point(node[1],info,profile)+';')
    return lines


def profile_kernel(lines,info,nodes,label,indent=2,profile=False):
    lines+=[' '*indent+'Kokkos::parallel_for("'+info['name']+'::'+label+'", nz, KOKKOS_LAMBDA(const int k) {']
    lines+=device_nodes(nodes,info,profile,indent+2)
    lines+=[' '*indent+'});']


def responder_definition(info):
    nodes,_=parse_nodes(info['records']);first=next(i for i,n in enumerate(nodes) if n[0]=='where')
    lines=[profile_signature(info)]+description(info)+['','  // Native profile order is retained across the two complex-root solves.']+stage_profile(info)
    lines+=['  //'+n[1] if n[0]=='comment' else '' for n in nodes[:first] if n[0] in ('comment','blank')]
    profile_kernel(lines,info,[nodes[first]],'clip_Skx',profile=True)
    for node in nodes[first+1:]:
        if node[0]=='comment':lines.append('  //'+node[1]);continue
        if node[0]=='blank':lines.append('');continue
        code=node[1]
        if code.lower()=='return' or code.lower().startswith('end subroutine'):continue
        if code.lower().startswith('call '):code=code[5:]
        lines.append('  '+view_call(code,info,True)+';')
    lines+=finish_profile(info)+['}',''];return '\n'.join(lines)


MAINHEADER='''// Original two-component PDF driver for w, rt and thl fields.
//
// new_pdf_driver computes the setter mixture and component moments, then closes
// responder rt/thl parameters and packs implicit-coefficient/explicit-term fields.
// Skrt and Skthl may be clipped in place. Fields are column-first LayoutLeft;
// scratch belongs to the call and the caller owns implicit_coefs_terms storage.
// All three Native advection flags are retained. The only transient lower science
// dependency is the two cubic-root solves inside each responder's limit call.
// Isolated source draft: not compiled, integrated or numerically tested.
#pragma once
#include "array_types.hpp"
#include "pdf_parameter_module.hpp"
namespace clubb_kokkos {
'''
(HERE/'new_pdf_main.hpp').write_text(MAINHEADER+'\n'+main_signature(main,False)+'\n\n} // namespace clubb_kokkos\n')
cpp=['''// Setter/responder PDF orchestration; see the public header and integration notes.
#include "new_pdf_main.hpp"
#include "new_pdf.hpp"
#include "model_flags.hpp"
#include "parameter_indices.hpp"
namespace clubb_kokkos {
namespace {
''']
for name in ('calc_responder_var','calc_F_x_zeta_x_setter','calc_F_x_responder'):
    info=ROUTINES[name];cpp.append(profile_signature(info,False))
    if name!='calc_responder_var':cpp.append(scalar_signature(info,False))
cpp+=['} // namespace','',main_signature(main)]+description(main)+['','  // Arguments']
notes=arg_notes(main)
for a in main['arguments']:
    if a['name'] in notes:cpp.append('  // '+a['name']+': '+notes[a['name']])
cpp+=['','  //------------------------ Local Variables ------------------------', '  // Output only for recording statistics.', '  // Output functionality deprecated, see https://github.com/larson-group/clubb/issues/1176']
argnames={a['name'] for a in main['arguments']}
for name,d in main['decls'].items():
    if name in argnames or name=='i':continue
    if d['dims']:cpp.append('  Real'+str(len(d['dims']))+'D '+name+'("new_pdf_driver::'+name+'", '+', '.join(d['dims'])+');'+(' // '+notes[name] if name in notes else ''))
    else:cpp.append('  real_t '+name+';'+(' // '+notes[name] if name in notes else ''))
cpp+=['','  // Native profile helpers receive scalar tunables. Mirror once at this host boundary.','  const auto clubb_params_host = to_host_space(clubb_params);','  // Capture component Views rather than the owning derived type.']
for field in FIELDS:cpp.append('  const auto pdf_'+field+' = pdf_implicit_coefs_terms.'+field+';')
cpp+=['']
nodes,_=parse_nodes(main['records']);seq=0


def host_nodes(nodes,indent=2):
    global seq
    lines=[];index=0
    while index<len(nodes):
        node=nodes[index];kind=node[0]
        if kind=='comment':lines.append(' '*indent+'//'+node[1]);index+=1;continue
        if kind=='blank':lines.append('');index+=1;continue
        if kind=='do':
            if node[1]!='i' or node[2]!=['1','ngrdcol']:raise ValueError(node)
            lines.append(' '*indent+'for (int i = 0; i < ngrdcol; ++i) {');lines+=host_nodes(node[3],indent+2);lines.append(' '*indent+'}');index+=1;continue
        if kind=='if':
            for j,(condition,body) in enumerate(node[1]):lines.append(' '*indent+('if (' if j==0 else '} else if (')+expr(condition)+') {');lines+=host_nodes(body,indent+2)
            if node[2] is not None:lines.append(' '*indent+'} else {');lines+=host_nodes(node[2],indent+2)
            lines.append(' '*indent+'}');index+=1;continue
        if kind=='where':seq+=1;profile_kernel(lines,main,[node],'profile_'+str(seq),indent);index+=1;continue
        code=node[1];low=code.lower()
        if low.startswith('end subroutine'):index+=1;continue
        if low=='return':lines.append(' '*indent+'return;');index+=1;continue
        if low.startswith('call '):lines.append(' '*indent+view_call(code[5:],main)+';');index+=1;continue
        if ' = calc_' in code:
            lines.append(' '*indent+view_call(code,main)+';');index+=1;continue
        lhs=code.split('=',1)[0].strip();name=re.match(r'(\w+)',lhs).group(1)
        if '%' not in code and name in main['decls'] and not main['decls'][name]['dims']:
            lines.append(' '*indent+expr(code)+';');index+=1;continue
        # Fuse adjacent independent assignments in the same Native phase.
        group=[]
        while index<len(nodes) and nodes[index][0]=='code':
            candidate=nodes[index][1]
            if candidate.lower()=='return' or candidate.lower().startswith(('call ','end ')) or ' = calc_' in candidate:break
            lhs=candidate.split('=',1)[0].strip();nm=re.match(r'(\w+)',lhs).group(1)
            if '%' not in candidate and nm in main['decls'] and not main['decls'][nm]['dims']:break
            group.append(nodes[index]);index+=1
        if not group:raise ValueError(('host group',node))
        seq+=1;profile_kernel(lines,main,group,'profile_'+str(seq),indent)
    return lines

cpp+=host_nodes(nodes)+['}','', 'namespace {',responder_definition(ROUTINES['calc_responder_var'])]
for name in ('calc_F_x_zeta_x_setter','calc_F_x_responder'):
    cpp+=['',scalar_definition(ROUTINES[name]),profile_definition(ROUTINES[name])]
cpp+=['} // namespace','} // namespace clubb_kokkos','']
(HERE/'new_pdf_main.cpp').write_text('\n'.join(cpp))
print('Prepared the29-argument Native parent and all three private helpers.')
