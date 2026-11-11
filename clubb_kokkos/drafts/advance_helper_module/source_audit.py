"""Compare draft contracts and arithmetic structure; do not execute physics.

This is a source review aid, not a compiled/numerical test. Handwritten boundary
and ordered-reduction adaptations are reviewed separately in SOURCE_AUDIT.json.
"""
from pathlib import Path
import ast,hashlib,json,re,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,read_procedure,declaration
from generate_derived_types import split_top
NATIVE=ROOT/'src/CLUBB_core/advance_helper_module.F90'
TARGET=HERE/'advance_helper_module.hpp'
source=NATIVE.read_text();target=TARGET.read_text()
code=re.sub(r'//[^\n]*','',target)
pattern=re.compile(r'(?:KOKKOS_INLINE_FUNCTION|inline)\s+(?:real_t|void|Real[12]D)\s+(\w+)\s*\(([^{};]*?)\)\s*\{',re.S)
definitions=[]
for m in pattern.finditer(code):
 depth=1;end=m.end()
 while depth:
  depth+=(code[end]=='{')-(code[end]=='}');end+=1
 definitions.append({'name':m.group(1),'args':[x.strip().split()[-1].lstrip('&*') for x in split_top(m.group(2))],'body':code[m.end():end-1]})
names=re.findall(r'(?m)^  (?:subroutine|function) (\w+)\(',source)
assert [x['name'] for x in definitions[:26]]==names

# A small independent C++ expression reader retains every arithmetic association.
class Cpp:
 ops={'||':1,'&&':2,'==':3,'!=':3,'>':4,'<':4,'>=':4,'<=':4,'+':5,'-':5,'*':6,'/':6}
 def __init__(self,text,arrays):
  text=re.sub(r'static_cast<int>','as_int',text)
  self.tokens=re.findall(r'"(?:\\.|[^"\\])*"|(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?f?|[A-Za-z_]\w*(?:::\w+)*|->|&&|\|\||==|!=|>=|<=|[^\s]',text)
  self.i=0;self.arrays=arrays
 def pop(self):v=self.tokens[self.i];self.i+=1;return v
 def expr(self,minimum=0):
  v=self.pop()
  if v=='(':a=self.expr();assert self.pop()==')'
  elif v in ('-','+','!','*'):
   a=('unary',v,self.expr(7))
  elif v.startswith('"'):a=('str',json.loads(v))
  elif re.match(r'^(?:\d|\.)',v):a=('n',float(v.rstrip('f')))
  else:
   a=('v',v.lower().replace('kokkos::','').replace('f2k::',''))
   if a==('v','true'):a=('bool',True)
   if a==('v','false'):a=('bool',False)
  while self.i<len(self.tokens):
   op=self.tokens[self.i]
   if op in ('.','->'):
    self.pop();a=('field',a,self.pop().lower());continue
   if op=='(':
    self.pop();args=[]
    if self.tokens[self.i]!=')':
     while True:
      args.append(self.expr())
      if self.tokens[self.i]!=',':break
      self.pop()
    assert self.pop()==')'
    name=a[1] if a[0]=='v' else None
    if name=='real_t':a=args[0]
    elif name=='as_int':a=('call','int',tuple(args))
    else:
     if name in self.arrays:
      adjusted=[]
      for index,arg in enumerate(args):
       if name=='clubb_params' and index==1:adjusted.append(arg)
       else:
        assert arg[0:2]==('bin','-') and arg[3]==('n',1.0),(name,arg)
        adjusted.append(arg[2])
      args=adjusted
     a=('call',name or str(a),tuple(args))
    continue
   if op=='?' and minimum<=0:
    self.pop();yes=self.expr();assert self.pop()==':';no=self.expr();a=('if',a,yes,no);continue
   prec=self.ops.get(op,-1)
   if prec<minimum:break
   self.pop();a=('bin',op,a,self.expr(prec+1))
  return a
 def parse(self):
  a=self.expr();assert self.i==len(self.tokens),(self.tokens[self.i:],self.tokens);return a

def fexpr(text):
 text=re.sub(r'(\w+)%([\w]+)',r'\1_\2',text)
 text=re.sub(r'(\w+)\(\s*:\s*(?:,\s*:\s*)?\)',r'\1',text)
 text=re.sub(r'(?i)(\d+(?:\.\d*)?(?:e[+-]?\d+)?|\.\d+(?:e[+-]?\d+)?)_core_rknd',r'\1',text)
 text=text.lower().replace('.true.','True').replace('.false.','False')
 for a,b in [('.and.',' and '),('.or.',' or '),('.not.',' not '),('/=','!=')]:text=text.replace(a,b)
 def node(n):
  if isinstance(n,ast.Name):return ('v',n.id.lower())
  if isinstance(n,ast.Constant):return ('bool',n.value) if isinstance(n.value,bool) else ('str',n.value) if isinstance(n.value,str) else ('n',float(n.value))
  if isinstance(n,ast.UnaryOp):return ('unary',{ast.USub:'-',ast.UAdd:'+',ast.Not:'!'}[type(n.op)],node(n.operand))
  if isinstance(n,ast.BinOp):
   a,b=node(n.left),node(n.right)
   if isinstance(n.op,ast.Pow):
    assert b in [('n',2.),('n',3.)]
    return ('bin','*',a,a) if b==('n',2.) else ('bin','*',('bin','*',a,a),a)
   return ('bin',{ast.Add:'+',ast.Sub:'-',ast.Mult:'*',ast.Div:'/'}[type(n.op)],a,b)
  if isinstance(n,ast.BoolOp):
   values=[node(v) for v in n.values];a=values[0]
   for b in values[1:]:a=('bin','&&' if isinstance(n.op,ast.And) else '||',a,b)
   return a
  if isinstance(n,ast.Compare):return ('bin',{ast.Gt:'>',ast.Lt:'<',ast.GtE:'>=',ast.LtE:'<=',ast.Eq:'==',ast.NotEq:'!='}[type(n.ops[0])],node(n.left),node(n.comparators[0]))
  if isinstance(n,ast.Call):
   name=n.func.id.lower();args=tuple(node(a) for a in n.args)
   if name in {'max','min'}:return ('if',('bin','>' if name=='max' else '<',args[0],args[1]),args[0],args[1])
   return ('call',name,args)
  raise ValueError(ast.dump(n))
 return node(ast.parse(text.strip(),mode='eval').body)

reviews=[]
for d in definitions[:26]:
 name=d['name'];p=read_procedure(NATIVE,name)
 assert d['args']==[a.name for a in p.arguments],(name,'signature',d['args'])
 raw=re.search(r'(?ms)^  (?:subroutine|function) '+name+r'\(.*?^  end (?:subroutine|function) '+name+r'\b',source).group()
 decls={}
 for line in statements(raw):decls.update(declaration(line))
 arrays={key for key,a in decls.items() if a.dims}
 arrays.update({'gr_zm','gr_zt','gr_dzm','gr_dzt','gr_invrs_dzm','gr_invrs_dzt'})
 if name in {'set_boundary_conditions_lhs','set_boundary_conditions_rhs','vertical_avg','vertical_integral'}:
  reviews.append({'name':name,'arguments':len(p.arguments),'review':'Manual ordered boundary/serial reduction adaptation review'});continue
 native=[]
 for line in statements(raw):
  m=re.match(r'^(\w+(?:\([^()]*\))?)\s*=(?!=)\s*(.*)$',line)
  if m and not (name=='pvertinterp' and m.group(1)=='l_error'):
   native.append((fexpr(m.group(1)),fexpr(m.group(2))))
 actual=[];combined=[]
 for m in re.finditer(r'(?m)^\s*(\w+(?:\([^;]*?\))?)\s*=(?!=)\s*([^;]*);',d['body']):
  if name=='pvertinterp' and m.group(1)=='l_error':continue
  combined.append((m.start(),m.group(1),m.group(2)))
 for m in re.finditer(r'Kokkos::deep_copy\s*\(',d['body']):
  pos=m.end();depth=1;end=pos
  while depth:
   depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  args=split_top(d['body'][pos:end-1]);assert len(args)==2
  combined.append((m.start(),*args))
 for pos,lhs,rhs in sorted(combined):actual.append((Cpp(lhs,arrays).parse(),Cpp(rhs,arrays).parse()))
 assert native==actual,(name,'assignments',native,actual)
 native_conditions=[]
 for line in statements(raw):
  m=re.match(r'(?:if|elseif|else if)\s*\((.*)\)\s*then$',line,re.I)
  if m:native_conditions.append(fexpr(m.group(1)))
 actual_conditions=[]
 for m in re.finditer(r'\bif\s*\(',d['body']):
  pos=m.end();depth=1;end=pos
  while depth:
   depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  text=d['body'][pos:end-1]
  text=text.replace('stats != nullptr','present(stats)').replace('stats->flags().l_sample','stats_l_sample').replace('stats.flags().l_sample','stats_l_sample')
  actual_conditions.append(Cpp(text,arrays).parse())
 assert native_conditions==actual_conditions,(name,'branch predicates',native_conditions,actual_conditions)
 native_calls=[fexpr(re.sub(r'^call\s+','',line,flags=re.I)) for line in statements(raw) if re.match(r'call\s+',line,re.I)]
 actual_calls=[]
 for m in re.finditer(r'f2k::stats_update\s*\(',d['body']):
  pos=m.end();depth=1;end=pos
  while depth:
   depth+=(d['body'][end]=='(')-(d['body'][end]==')');end+=1
  text=d['body'][m.start():end].replace('*stats','stats')
  actual_calls.append(Cpp(text,arrays).parse())
 # calc_Ri_zm is the only non-stats CALL; its exact contract/order is visible.
 if name=='compute_Cx_fnc_Richardson':
  m=re.search(r'calc_Ri_zm\([^;]*\)',d['body']);actual_calls.insert(0,Cpp(m.group(),arrays).parse())
 assert native_calls==actual_calls,(name,'ordered calls',native_calls,actual_calls)
 reviews.append({'name':name,'arguments':len(p.arguments),'ordered_branch_predicates':len(native_conditions),'ordered_subroutine_calls':len(native_calls),'ordered_assignment_expression_trees':len(native),'review':'All arithmetic/grouping/index/argument trees exact after documented syntax/rank mappings'})
record={'canonical_source':str(NATIVE.relative_to(ROOT)),'canonical_sha256':hashlib.sha256(NATIVE.read_bytes()).hexdigest(),'target':str(TARGET.relative_to(ROOT)),'target_sha256':hashlib.sha256(TARGET.read_bytes()).hexdigest(),'native_definition_order':names,'routines':reviews,'status':'Source-only review aid. No build, execution or BFB claim.'}
(HERE/'expression_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print('26 Native signatures/order checked; scientific assignment trees matched for22 translated routines.')
