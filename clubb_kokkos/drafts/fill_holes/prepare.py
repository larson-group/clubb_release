"""Build this isolated draft from the previously read Native source outline."""
from pathlib import Path
import ast,dataclasses,json,re,sys,textwrap
folder=Path(__file__).resolve().parent;root=folder.parents[2]
sys.path.insert(0,str(root/'clubb_kokkos/src/f2k'))
from fortran_procedures import statements,declaration
from generate_derived_types import split_top
source=(root/'src/CLUBB_core/fill_holes.F90').read_text()
outline=json.loads((folder/'SOURCE_OUTLINE.json').read_text())
names=[r['name'] for r in outline['routines'] if r['supported']]
blocks={n:re.search(r'^  subroutine '+n+r'\b.*?^  end subroutine '+n+r'\b',source,re.M|re.S)[0] for n in names}
metadata={};arguments={}
for n,b in blocks.items():
    ss=statements(b);meta={}
    for s in ss:
        if '::' in s:meta.update(declaration(s))
    metadata[n]=meta;arguments[n]=split_top(ss[0][ss[0].index('(')+1:ss[0].rindex(')')])

def description(n):
    lines=blocks[n].splitlines();start=next(i for i,s in enumerate(lines) if '! Description:' in s)
    end=next(i for i in range(start+1,len(lines)) if re.match(r'\s*use\b',lines[i],re.I))
    return '\n'.join('  //'+s.split('!',1)[1] if '!' in s else '' for s in lines[start:end]).rstrip()

def docs(n):
    # Only declaration comments after IMPLICIT NONE describe arguments.
    lines=blocks[n].splitlines();start=next(i for i,s in enumerate(lines) if 'implicit none' in s.lower());result={}
    for key,a in metadata[n].items():
        if key not in [x.lower() for x in arguments[n]]:continue
        for i in range(start+1,len(lines)):
            line=lines[i]
            if '!' not in line or line.lstrip().startswith('!'):continue
            code,note=line.split('!',1)
            if not re.search(r'\b'+re.escape(a.name)+r'\s*(?:,?\s*&)?\s*$',code):continue
            notes=[note.strip()]
            for tail in lines[i+1:]:
                if not tail.strip().startswith('!') or not tail.strip(' !-'):break
                notes.append(tail.split('!',1)[1].strip())
            result[key]=' '.join(notes);break
    defaults={'nz':'Number of grid levels [-]','nzm':'Number of momentum levels [-]',
      'nzt':'Number of thermodynamic levels [-]','ngrdcol':'Number of columns [-]',
      'hydromet_dim':'Number of hydrometeor species [-]',
      'hm_one_lev':'Same-phase species mixing ratios at one level [kg/kg]',
      'hm_one_lev_filled':'Same-phase mixing ratios after redistribution [kg/kg]',
      'hydromet':'Hydrometeor fields (level,species) [units vary]',
      'hydromet_filled':'Hydrometeor fields after frozen-species redistribution [units vary]',
      'hydromet_name':'Fixed-width species name [-]',
      'hydromet_list':'Fixed-width species names (character,species) [-]',
      'rvm_mc':'Water-vapor tendency [kg/kg/s]','thlm_mc':'Liquid-water potential-temperature tendency [K/s]',
      'gr':'Native grid coordinates and spacings [-,m]',
      'hm_metadata':'Hydrometeor species indices, flags and tolerances [-,units vary]',
      'l_fill_holes_hm':'Enable same-level frozen-species redistribution [-]',
      'stats':'Native statistics context, mutated by sampled budget calls [-]',
      'ihm':'Native1based species index [-]',
      'l_use_wp3_lim_with_smth_Heaviside':'Select smooth skewness limiter [-]'}
    for key in arguments[n]:
        if key.lower() not in result and key.lower() in defaults:result[key.lower()]=defaults[key.lower()]
    if n=='fill_holes_wp2_from_horz_tke':result['up2']='Zonal TKE component supplying the hole [m^2/s^2]'
    return result

def signature(n,prototype=False):
    typed={};meta=metadata[n]
    for key in arguments[n]:
        a=meta[key.lower()]
        if a.dims:
            typ=('Real'+str(len(a.dims))+'D') if a.kind=='real' else 'Bool1D' if a.kind=='logical' else 'Kokkos::View<char**, Kokkos::LayoutLeft>'
        elif a.kind in ('grid','hm_metadata_type','stats_type'):
            typ=('f2k::' if a.kind=='stats_type' else '')+a.kind+'&'
            if a.intent=='in':typ='const '+typ
        elif a.kind=='character':typ='std::string*' if a.optional else 'const std::string&'
        else:typ={'real':'real_t','integer':'int','logical':'bool'}[a.kind]+('&' if a.intent in ('inout','out') else '')
        typed[key]=typ+' '+key+(' = nullptr' if prototype and a.optional else '')
    # Use Native continuation groups, splitting only for C++ width/intent.
    groups=[];pending=[]
    for raw in blocks[n].splitlines():
        code=raw.split('!',1)[0]
        if '(' in code and not pending and not groups:code=code.split('(',1)[1]
        code=code.replace('&','').rstrip().removesuffix(')')
        keys=[s.strip() for s in code.split(',') if s.strip()]
        if keys:groups.append(keys)
        if ')' in raw.split('!',1)[0]:break
    flattened=[x for g in groups for x in g];assert flattened==arguments[n],(n,flattened,arguments[n])
    lines=[]
    for group in groups:
        current=[];intent=None
        for key in group:
            a=meta[key.lower()];tag=('Optional ' if a.optional else '')+{'in':'In','inout':'InOut','out':'Out'}[a.intent]
            if current and (intent!=tag or len(', '.join(current+[typed[key]]))>76):
                lines.append((', '.join(current),intent));current=[]
            current.append(typed[key]);intent=tag
        if current:lines.append((', '.join(current),intent))
    d=docs(n);out=[]
    for key in arguments[n]:
        if key.lower() in d:out += ['// '+x for x in textwrap.wrap(key+': '+d[key.lower()],width=94)]
    out.append('void '+n+'(')
    for i,(text,tag) in enumerate(lines):out.append('    '+text+(',' if i<len(lines)-1 else ')')+' // '+tag)
    return '\n'.join(out)

# A small source translator is restricted to scalar cell calculations below.
# Arrays and loops retain the original Native symbols and association.
tree=ast.parse((root/'clubb_kokkos/drafts/advance_xm_wpxp_module/prepare.py').read_text())
namespace={'dataclass':dataclasses.dataclass,'field':dataclasses.field,'re':re,'split_top':split_top}
for obj in tree.body:
    if isinstance(obj,(ast.FunctionDef,ast.ClassDef)) and obj.name in ('Node','tokens','parse','powers','split_call'):
        exec(compile(ast.Module(body=[obj],type_ignores=[]),'<source-parser>','exec'),namespace)
tokens=namespace['tokens'];parse=namespace['parse'];powers=namespace['powers'];split_call=namespace['split_call']

def expr(s,n):
    strings=[]
    def store(m):strings.append(m[0]);return '@S'+str(len(strings)-1)+'@'
    s=re.sub(r'"[^"\n]*"|\'[^\'\n]*\'',store,s)
    s=powers(s)
    s=re.sub(r'((?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)_core_rknd',lambda m:'real_t('+m[1]+')',s)
    for a,b in [('.and.','&&'),('.or.','||'),('.not.','!'),('.true.','true'),('.false.','false'),('/=','!=')]:s=re.sub(re.escape(a),b,s,flags=re.I)
    for fun in ('abs','sqrt'):s=re.sub(r'\b'+fun+r'\s*\(',r'Kokkos::'+fun+'(',s,flags=re.I)
    s=re.sub(r'epsilon\s*\([^()]*\)',r'std::numeric_limits<real_t>::epsilon()',s,flags=re.I)
    for fun,op in [('min','<'),('max','>')]:
        while m:=re.search(r'(?<![\w:])'+fun+r'\s*\(',s,re.I):
            args,end=split_call(s,m.end());assert len(args)==2,(n,s)
            a,b=args;s=s[:m.start()]+'(('+a+') '+op+' ('+b+') ? ('+a+') : ('+b+'))'+s[end:]
    for a in metadata[n].values():
        if not a.dims:continue
        while m:=re.search(r'\b'+a.name+r'\s*\(',s):
            indexes,end=split_call(s,m.end());assert all(':' not in x for x in indexes),(n,s)
            value=a.name+'('+', '.join('('+expr(x,n)+') - 1' for x in indexes)+')'
            strings.append(value);s=s[:m.start()]+'@S'+str(len(strings)-1)+'@'+s[end:]
    for i in reversed(range(len(strings))):s=s.replace('@S'+str(i)+'@',strings[i])
    return s.strip()

def emit(nodes,n):
    out=[]
    for x in nodes:
        if x.kind=='comment':
            if not x.text.lstrip().startswith(('print','if ( .true.')):out.append('//'+x.text)
        elif x.kind=='if':
            out += ['if ('+expr(x.text,n)+') {']+['  '+s for s in emit(x.body,n)]+['}']
            if x.other:out += ['else {']+['  '+s for s in emit(x.other,n)]+['}']
        elif x.kind=='for':
            m=re.match(r'do\s+(\w+)\s*=\s*(.*)',x.text,re.I);idx=m[1];bounds=split_top(m[2]);lo,hi=bounds[:2];step=bounds[2] if len(bounds)>2 else '1'
            out += ['for (int '+idx+' = '+expr(lo,n)+'; '+idx+' <= '+expr(hi,n)+'; '+idx+' += '+expr(step,n)+') {']+['  '+s for s in emit(x.body,n)]+['}']
        else:
            assert not re.match(r'write|print|call',x.text,re.I),x.text
            out.append(expr(x.text,n)+';')
    return out

def native_wp2_science():
    b=blocks['fill_holes_wp2_from_horz_tke'];part=b[b.index('        tke_x2_old ='):b.index('        if ( l_clubb_at_least_debug_level_3 ) then')]
    nodes,_=parse(tokens(part));return '\n'.join(emit(nodes,'fill_holes_wp2_from_horz_tke'))

def wp2_warnings():
    b=blocks['fill_holes_wp2_from_horz_tke'];part=b[b.index('          ! Check for conservation'):b.index('        end if\n\n      end do')]
    # Host diagnostics are rendered directly from their Native logical statements.
    out=[];level=0
    for s in statements(part):
        if s.lower().startswith('if '):
            cond=s[s.index('(')+1:s.rindex(')')];out.append('  '*level+'if ('+expr(cond,'fill_holes_wp2_from_horz_tke')+') {');level+=1
        elif s.lower().startswith('end'):level-=1;out.append('  '*level+'}')
        elif s.lower().startswith('write'):
            # Quoted diagnostic labels contain commas; the bridge declaration
            # splitter intentionally handles declarations only, so do not use it here.
            payload=[];start=0;depth=0;quote=None;raw=s[s.index(')')+1:]
            for i,ch in enumerate(raw):
                if quote:
                    if ch==quote:quote=None
                elif ch in '\"\'':quote=ch
                elif ch=='(':depth+=1
                elif ch==')':depth-=1
                elif ch==',' and depth==0:payload.append(raw[start:i].strip());start=i+1
            payload.append(raw[start:].strip())
            items=[p if p.startswith(('"',"'")) else expr(p,'fill_holes_wp2_from_horz_tke') for p in payload]
            out.append('  '*level+'std::cerr << '+ ' << " " << '.join(items)+' << "\\n";')
        else:out.append('  '*level+expr(s,'fill_holes_wp2_from_horz_tke')+';')
    return '\n'.join(out)

bodies={}
bodies['fill_holes_vertical_api']=r'''
  Real2D field_initial;
  if (l_debug || l_print_before_after) {
    field_initial = Real2D("fill_holes_vertical_api::field_initial", ngrdcol, nz);
    Kokkos::deep_copy(field_initial, field);
    std::cerr << "field threshold: " << threshold << "\nfield: ";
    print_field(field); std::cerr << " -- over\n";
  }
  // Only bother with a fill call if there are values below threshold.
  bool l_field_below_threshold = any_below(field, threshold);
  if (l_field_below_threshold) {
    if (fill_holes_type == global_fill) {
      // Modify the entire range: maximally effective, computationally cheap, minimally local.
      fill_holes_global(nz, ngrdcol, threshold,
                        lower_hf_level, upper_hf_level,
                        dz, rho_ds, grid_dir_indx,
                        field);
    } else if (fill_holes_type == sliding_window) {
      // Consecutive overlapping windows are serial; global fill is the Native fallback.
      fill_holes_sliding_window(nz, ngrdcol, threshold,
                                lower_hf_level, upper_hf_level,
                                dz, rho_ds, grid_dir_indx,
                                field);
    } else {
      // Same supported scope as JAX: the four other Native methods are pending.
      // Preserve the Native no-hole return before method dispatch/termination.
      std::cerr << "ERROR: unknown fill_holes_type: " << fill_holes_type << "\n";
      std::exit(EXIT_SUCCESS);
    }
  }
  if (l_debug || l_print_before_after) {
    auto values = to_host_space(field);
    auto initial = to_host_space(field_initial);
    auto density = to_host_space(rho_ds);
    auto spacing = to_host_space(dz);
    if (any_below(field, threshold)) {
      std::cerr << "THERE ARE BELOW THRESHOLD VALUES IN THE FIELD AFTER THE FILL\n"
                   "EVEN THOUGH THE FIELD AVERAGE IS BELOW THRESHOLD. This is not\n"
                   "a gauranteed error, since smoothing contraints may prevent \n"
                   "some above threshold levels from having all their mass exhausted.\n";
      for (int i = 1; i <= ngrdcol; ++i) {
        real_t numer = zero, denom = zero;
        for (int k = 1; k <= nz; ++k) {
          numer += values(i-1,k-1) * density(i-1,k-1) * spacing(i-1,k-1);
          denom += density(i-1,k-1) * spacing(i-1,k-1);
        }
        if (numer / denom > threshold) {
          for (int k = 1; k <= nz; ++k) {
            if (values(i-1,k-1) < threshold)
              std::cout << "field(" << i << " " << k << ") = " << values(i-1,k-1)
                        << " < " << threshold << "\n";
          }
          std::cerr << "column " << i << " field average = " << numer / denom
                    << " -- threshold = " << threshold << "\n";
        } else {
          std::cerr << "FIELD AVERAGE BELOW THRESHOLD - FILLING INCOMPLETE AS EXPECTED\n";
        }
      }
    }
    real_t before = zero, after = zero;
    for (int k = 1; k <= nz; ++k) for (int i = 1; i <= ngrdcol; ++i) {
      before += initial(i-1,k-1) * density(i-1,k-1) * spacing(i-1,k-1);
      after += values(i-1,k-1) * density(i-1,k-1) * spacing(i-1,k-1);
    }
    std::cerr << "difference in mass (before-after): " << before-after << "\nfield after: ";
    print_field(field); std::cerr << " -- over\n";
  }
'''
bodies['fill_holes_global']=r'''
  if (l_debug) std::cout << "fill_holes_type: GLOBAL FILL\n";
  Real2D rho_ds_dz(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_global::rho_ds_dz"),ngrdcol,nz);
  Real2D field_clipped(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_global::field_clipped"),ngrdcol,nz);
  Real1D denom_integral_global(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_global::denom_integral_global"),ngrdcol);
  Real1D numer_integral_global(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_global::numer_integral_global"),ngrdcol);
  Real1D field_avg_global(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_global::field_avg_global"),ngrdcol);
  Real1D mass_fraction_global(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_global::mass_fraction_global"),ngrdcol);
  // Independent columns share a kernel; each vertical integral stays sequential in Native order.
  Kokkos::parallel_for("fill_holes_global::columns",ngrdcol,KOKKOS_LAMBDA(int column) {
    const int i = column + 1;
    for (int k = 1; k <= nz; ++k) rho_ds_dz(i-1,k-1) = rho_ds(i-1,k-1) * dz(i-1,k-1);
    // Compute the numerator and denominator integrals.
    numer_integral_global(i-1) = real_t(0.0);
    denom_integral_global(i-1) = real_t(0.0);
    for (int k = lower_hf_level; grid_dir_indx*(k-upper_hf_level) <= 0; k += grid_dir_indx) {
      numer_integral_global(i-1) = numer_integral_global(i-1) + rho_ds_dz(i-1,k-1)*field(i-1,k-1);
      denom_integral_global(i-1) = denom_integral_global(i-1) + rho_ds_dz(i-1,k-1);
    }
    // Vertical average: see advance_helper_module::vertical_avg.
    field_avg_global(i-1) = numer_integral_global(i-1) / denom_integral_global(i-1);
    // Clip upward if all holes can be filled; otherwise remove all above-threshold mass.
    for (int k = lower_hf_level; grid_dir_indx*(k-upper_hf_level) <= 0; k += grid_dir_indx) {
      if (field_avg_global(i-1) >= threshold)
        field_clipped(i-1,k-1) = threshold > field(i-1,k-1) ? threshold : field(i-1,k-1);
      else
        field_clipped(i-1,k-1) = threshold < field(i-1,k-1) ? threshold : field(i-1,k-1);
    }
    // Only the numerator changes for the clipped vertical integral.
    numer_integral_global(i-1) = real_t(0.0);
    for (int k = lower_hf_level; grid_dir_indx*(k-upper_hf_level) <= 0; k += grid_dir_indx)
      numer_integral_global(i-1) = numer_integral_global(i-1) + rho_ds_dz(i-1,k-1)*field_clipped(i-1,k-1);
    bool l_column_below_threshold = false;
    for (int k = lower_hf_level; grid_dir_indx*(k-upper_hf_level) <= 0; k += grid_dir_indx)
      if (field(i-1,k-1) < threshold) l_column_below_threshold = true;
    // Do not update a column without holes inside the filling range.
    if (l_column_below_threshold) {
      const real_t field_clipped_avg = numer_integral_global(i-1) / denom_integral_global(i-1);
      // Avoid division by zero when the clipped average is approximately threshold.
      if (Kokkos::abs(field_clipped_avg-threshold) > Kokkos::abs(field_clipped_avg+threshold)*eps/2) {
        // Multiplicative normalization conserves the original integrated mass.
        mass_fraction_global(i-1) = (field_avg_global(i-1)-threshold)/(field_clipped_avg-threshold);
        for (int k = lower_hf_level; grid_dir_indx*(k-upper_hf_level) <= 0; k += grid_dir_indx)
          field(i-1,k-1) = threshold + mass_fraction_global(i-1)*(field_clipped(i-1,k-1)-threshold);
      }
    }
  });
'''
bodies['fill_holes_sliding_window']=r'''
  if (l_debug) std::cout << "fill_holes_type: SLIDING WINDOW FILL WITH GLOBAL FALLBACK\n";
  bool l_field_below_threshold = any_below(field,threshold);
  if (l_field_below_threshold) {
    Real2D rho_ds_dz(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_sliding_window::rho_ds_dz"),ngrdcol,nz);
    Real2D invrs_denom_integral(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_sliding_window::invrs_denom_integral"),ngrdcol,nz);
    Real2D field_clipped(Kokkos::view_alloc(Kokkos::WithoutInitializing,"fill_holes_sliding_window::field_clipped"),ngrdcol,nz);
    Kokkos::parallel_for("fill_holes_sliding_window::rho_ds_dz",ngrdcol*nz,KOKKOS_LAMBDA(int cell) {
      const int i = cell%ngrdcol, k = cell/ngrdcol;
      rho_ds_dz(i,k) = rho_ds(i,k)*dz(i,k);
    });
    const int start_indx = lower_hf_level + grid_dir_indx*num_hf_draw_points;
    const int stop_indx = upper_hf_level - grid_dir_indx*num_hf_draw_points;
    // The denominator does not change. Explicit source-order sums retain the Native
    // OpenMP-offload improvement: https://github.com/larson-group/clubb/issues/1138#issuecomment-1974918151
    Kokkos::parallel_for("fill_holes_sliding_window::denominator",ngrdcol,KOKKOS_LAMBDA(int column) {
      const int i = column+1;
      for (int k = start_indx; grid_dir_indx*(k-stop_indx) <= 0; k += grid_dir_indx) {
        real_t rho_k_sum = real_t(0.0);
        const int start_indx_j = k-grid_dir_indx*num_hf_draw_points;
        const int stop_indx_j = k+grid_dir_indx*num_hf_draw_points;
        for (int j = start_indx_j; grid_dir_indx*(j-stop_indx_j) <= 0; j += grid_dir_indx)
          rho_k_sum = rho_k_sum + rho_ds_dz(i-1,j-1);
        invrs_denom_integral(i-1,k-1) = one/rho_k_sum;
      }
    });
    // Make one pass through each profile. Overlapping windows read prior updates;
    // the k loop must stay serial, including for the descending grid.
    Kokkos::parallel_for("fill_holes_sliding_window::columns",ngrdcol,KOKKOS_LAMBDA(int column) {
      const int i = column+1;
      for (int k = start_indx; grid_dir_indx*(k-stop_indx) <= 0; k += grid_dir_indx) {
        const int k_start = k-grid_dir_indx*num_hf_draw_points;
        const int k_end = k+grid_dir_indx*num_hf_draw_points;
        bool l_window_below_threshold = false;
        for (int j = k_start; grid_dir_indx*(j-k_end) <= 0; j += grid_dir_indx)
          if (field(i-1,j-1) < threshold) l_window_below_threshold = true;
        if (l_window_below_threshold) {
          // Compute the original field average, which must be conserved.
          real_t field_avg = real_t(0.0);
          for (int j = k_start; grid_dir_indx*(j-k_end) <= 0; j += grid_dir_indx)
            field_avg = field_avg + rho_ds_dz(i-1,j-1)*field(i-1,j-1);
          field_avg = field_avg*invrs_denom_integral(i-1,k-1);
          // Remove either the holes or the available above-threshold material.
          for (int j = k_start; grid_dir_indx*(j-k_end) <= 0; j += grid_dir_indx) {
            if (field_avg >= threshold)
              field_clipped(i-1,j-1) = threshold > field(i-1,j-1) ? threshold : field(i-1,j-1);
            else
              field_clipped(i-1,j-1) = threshold < field(i-1,j-1) ? threshold : field(i-1,j-1);
          }
          // The clipped total mass is no smaller when filling all holes is possible.
          real_t field_clipped_avg = real_t(0.0);
          for (int j = k_start; grid_dir_indx*(j-k_end) <= 0; j += grid_dir_indx)
            field_clipped_avg = field_clipped_avg + rho_ds_dz(i-1,j-1)*field_clipped(i-1,j-1);
          field_clipped_avg = field_clipped_avg*invrs_denom_integral(i-1,k-1);
          if (Kokkos::abs(field_clipped_avg-threshold) > Kokkos::abs(field_clipped_avg+threshold)*eps/2) {
            const real_t mass_fraction = (field_avg-threshold)/(field_clipped_avg-threshold);
            for (int j = k_start; grid_dir_indx*(j-k_end) <= 0; j += grid_dir_indx)
              field(i-1,j-1) = threshold + mass_fraction*(field_clipped(i-1,j-1)-threshold);
          }
        }
      }
    });
  }
  // Check the entire field, as Native does, before the global fallback.
  l_field_below_threshold = any_below(field,threshold);
  if (l_field_below_threshold) {
    if (l_debug) std::cout << "first pass insufficient -- doing global fill\n";
    fill_holes_global(nz,ngrdcol,threshold,
                      lower_hf_level,upper_hf_level,
                      dz,rho_ds,grid_dir_indx,
                      field);
  }
'''

# Native cell science and diagnostics are independently selected from the source.
bodies['fill_holes_wp2_from_horz_tke']=r'''
  bool l_field_below_threshold = any_below(wp2,threshold);
  if (!l_field_below_threshold) return;
  const bool l_clubb_at_least_debug_level_3 = f2k::clubb_at_least_debug_level_api(3);
  Real2D wp2_initial,up2_initial,vp2_initial;
  if (l_clubb_at_least_debug_level_3) {
    wp2_initial = Real2D("fill_holes_wp2_from_horz_tke::wp2_initial",ngrdcol,nz);
    up2_initial = Real2D("fill_holes_wp2_from_horz_tke::up2_initial",ngrdcol,nz);
    vp2_initial = Real2D("fill_holes_wp2_from_horz_tke::vp2_initial",ngrdcol,nz);
    Kokkos::deep_copy(wp2_initial,wp2);Kokkos::deep_copy(up2_initial,up2);Kokkos::deep_copy(vp2_initial,vp2);
  }
  // Each level transfers TKE independently; numeric lower..upper bounds are unchanged.
  Kokkos::parallel_for("fill_holes_wp2_from_horz_tke::cells",ngrdcol*std::max(0,upper_hf_level-lower_hf_level+1),KOKKOS_LAMBDA(int cell) {
    const int i = cell%ngrdcol+1,k = cell/ngrdcol+lower_hf_level;
    real_t up2_old,vp2_old,wp2_old,tke_x2_old,missing_wp2,up2_avail,vp2_avail,up2_vp2_avail,ratio;
@SCIENCE@
  });
  if (l_clubb_at_least_debug_level_3) {
    // Mirror old/new values for Native ordered diagnostics; no I/O inside kernels.
    auto wp2_new = to_host_space(wp2),up2_new = to_host_space(up2),vp2_new = to_host_space(vp2);
    auto wp2_before = to_host_space(wp2_initial),up2_before = to_host_space(up2_initial),vp2_before = to_host_space(vp2_initial);
    for (int i = 1;i <= ngrdcol;++i) for (int k = lower_hf_level;k <= upper_hf_level;++k) {
      const real_t wp2_old=wp2_before(i-1,k-1),up2_old=up2_before(i-1,k-1),vp2_old=vp2_before(i-1,k-1);
      const real_t tke_x2_old=up2_old+vp2_old+wp2_old;
      real_t tke_diff;
@WARNINGS@
    }
  }
'''.replace('@SCIENCE@','\n'.join('    '+s for s in native_wp2_science().splitlines()))
warnings=wp2_warnings()
for a in ('wp2','up2','vp2'):warnings=re.sub(r'\b'+a+r'\(',a+'_new(',warnings)
bodies['fill_holes_wp2_from_horz_tke']=bodies['fill_holes_wp2_from_horz_tke'].replace('@WARNINGS@','\n'.join('      '+s for s in warnings.splitlines()))

header_top='''// Fill negative or sub-threshold fields while conserving the available mass or TKE.
//
// fill_holes_vertical_api dispatches global/sliding-window vertical filling;
// fill_holes_wp2_from_horz_tke transfers horizontal TKE at the same level.
// hole_filling_hm_one_lev and fill_holes_hydromet_api redistribute frozen species;
// fill_holes_wv transfers residual holes to vapor/latent-heat tendencies.
// fill_holes_driver_api coordinates those stages and Native budget calls.
// clip_hydromet_conc_mvr enforces particle-radius limits; setup_stats_names
// supplies species budgets/fall speeds. Arrays use LayoutLeft; batched physical
// fields are column first, while single-column APIs retain Native level/species
// axes. Grid/metadata indices remain Native1based. Other vertical methods are
// outside actual JAX-supported scope. Draft only; see INTEGRATION.md.
#pragma once
#include "array_types.hpp"
#include "grid_class.hpp"
#include "corr_varnce_module.hpp"
#include "f2k/stats_bridge.hpp"
#include "constants_clubb.hpp"
#include <string>
#include <limits>
namespace clubb_kokkos {
'''

def write():
    assert set(bodies)==set(names),set(names)-set(bodies)
    public=[n for n in names if n not in ('fill_holes_global','fill_holes_sliding_window')]
    header=header_top+'\n\n'.join(signature(n,True)+'\n;' for n in public)+'\n} // namespace clubb_kokkos\n'
    cpp='''// Native-order hole filling and hydrometeor clipping.
// See fill_holes.hpp for public entries and INTEGRATION.md for scope/dependencies.
#include "fill_holes.hpp"
#include "model_flags.hpp"
#include "f2k/procedure_bridge.hpp"
#include "f2k_dependencies.hpp" // Pending public lower contracts only.
#include <Kokkos_MathematicalFunctions.hpp>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
namespace clubb_kokkos {
namespace {
constexpr bool l_debug = false;
constexpr bool l_print_before_after = false;
// Native ANY(field<threshold) is shared by three host dispatch boundaries.
template<class View> bool any_below(View field,real_t threshold) {
  int below = 0;
  Kokkos::parallel_reduce("fill_holes::any_below",field.size(),KOKKOS_LAMBDA(int cell,int& count) {
    if (field.data()[cell] < threshold) ++count;
  },below);
  return below != 0;
}
template<class View> void print_field(View field) {
  auto values=to_host_space(field);
  std::cerr << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (std::size_t k=0;k<values.extent(1);++k) for (std::size_t i=0;i<values.extent(0);++i)
    std::cerr << values(i,k) << " ";
  std::cerr << "\\n";
}
'''
    cpp += '\n\n'.join(signature(n)+'\n;' for n in ('fill_holes_global','fill_holes_sliding_window'))+'\n} // namespace\n'
    for n in names:
        private=n in ('fill_holes_global','fill_holes_sliding_window')
        if private:cpp+='\nnamespace {\n'
        cpp+='\n//=============================================================================\n'+signature(n)+'\n{\n'+description(n)+'\n'+bodies[n].strip('\n')+'\n}\n'
        if private:cpp+='} // namespace\n'
    cpp+='} // namespace clubb_kokkos\n'
    (folder/'fill_holes.hpp').write_text(header);(folder/'fill_holes.cpp').write_text(cpp)

if __name__=='__main__':
    # Remaining manually reviewed Native bodies are held separately for readability.
    exec((folder/'remaining_bodies.py').read_text(),globals())
    write()
