"""Verify final formatting and exact pre-cleanup scientific-token preservation."""
from pathlib import Path
import re,json,hashlib
D=Path(__file__).resolve().parent;B=json.loads((D/'cleanup_baseline.json').read_text());record={}
pattern=r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\''
for name,baseline in B.items():
 p=D/name;s=p.read_text();raw=s
 casts=re.findall(r'real_t\((zero|one_half)\)',s)
 if name.endswith('.cpp'):
  assert len(casts)==5
  s=re.sub(r'real_t\((zero|one_half)\)',r'\1',s)
 else:assert not casts
 clean=re.sub(pattern,lambda m:re.sub(r'[^\n]',' ',m[0]) if m[0].startswith('/') else m[0],s,flags=re.S)
 tokens=re.findall(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|\w+|[^\s]',clean)
 digest=hashlib.sha256(json.dumps(tokens,separators=(',',':')).encode()).hexdigest()
 assert digest==baseline['before_token_sha256'] and len(tokens)==baseline['token_count']
 assert all(len(l)<=100 and l.rstrip()==l for l in raw.splitlines())
 record[name]={**baseline,'final_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'final_token_sha256':digest,'exact_formatting_tokens_including_case_and_strings_preserved_after_only_approved_same_kind_value_materializations':True,'same_kind_value_materializations':len(casts),'separate_cast_proof':'device_prvalue_audit.json','maximum_columns':max(map(len,raw.splitlines())),'trailing_whitespace':False}
(D/'cleanup_audit.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
