"""Morrison two-moment microphysics: module_mp_graupel.F90 (CLUBB build).

Public routine names and input/inout argument order follow Fortran. Output
buffers become a returned mapping with the Fortran names; columns have an extra
leading batch axis and run independently through lax.map. Bounds are static,
one-based full-column bounds, as supplied by the CLUBB caller.

GRAUPEL_INIT returns the source module constants as a pure parameter mapping,
constructed at tracing time from parameters_microphys and module-owned host
configuration scalars rather than stored as mutable JAX globals. Internal
scalar-column control flow is intentionally kept separate from this public
boundary. Default REAL state is float32; imported CLUBB constants retain float64.
Native arithmetic allows compiler contraction.

The unused SAM microphysics.F90 wrapper is excluded from the Fortran build.
WRF-only MP_GRAUPEL and radar diagnostics are outside this CLUBB port. Atmospheric
validation currently covers warm-cloud/rain LBA, not active ice or graupel.
"""
import jax
import jax.numpy as jnp
from jax import lax
from clubb_jax.src.CLUBB_core import constants_clubb as C
from clubb_jax.src.Microphys import parameters_microphys as P

F = jnp.float32
D = jnp.float64
I = jnp.int32
B = jnp.bool_

# Adding coefficient term for clex9_oct14 case. This will reduce NNUCCD and NNUCCC
# by some factor to allow cloud to persist at realistic time intervals.
NNUCCD_REDUCE_COEF = 1.0
NNUCCC_REDUCE_COEF = 1.0



_DERF1_A = (
    # block 0
    0.00000000005958930743, -0.00000000113739022964, 0.00000001466005199839,
    -0.00000016350354461960, 0.00000164610044809620, -0.00001492559551950604,
    0.00012055331122299265, -0.00085483269811296660, 0.00522397762482322257,
    -0.02686617064507733420, 0.11283791670954881569, -0.37612638903183748117,
    1.12837916709551257377,
    # block 1
    0.00000000002372510631, -0.00000000045493253732, 0.00000000590362766598,
    -0.00000006642090827576, 0.00000067595634268133, -0.00000621188515924000,
    0.00005103883009709690, -0.00037015410692956173, 0.00233307631218880978,
    -0.01254988477182192210, 0.05657061146827041994, -0.21379664776456006580,
    0.84270079294971486929,
    # block 2
    0.00000000000949905026, -0.00000000018310229805, 0.00000000239463074000,
    -0.00000002721444369609, 0.00000028045522331686, -0.00000261830022482897,
    0.00002195455056768781, -0.00016358986921372656, 0.00107052153564110318,
    -0.00608284718113590151, 0.02986978465246258244, -0.13055593046562267625,
    0.67493323603965504676,
    # block 3
    0.00000000000382722073, -0.00000000007421598602, 0.00000000097930574080,
    -0.00000001126008898854, 0.00000011775134830784, -0.00000111992758382650,
    0.00000962023443095201, -0.00007404402135070773, 0.00050689993654144881,
    -0.00307553051439272889, 0.01668977892553165586, -0.08548534594781312114,
    0.56909076642393639985,
    # block 4
    0.00000000000155296588, -0.00000000003032205868, 0.00000000040424830707,
    -0.00000000471135111493, 0.00000005011915876293, -0.00000048722516178974,
    0.00000430683284629395, -0.00003445026145385764, 0.00024879276133931664,
    -0.00162940941748079288, 0.00988786373932350462, -0.05962426839442303805,
    0.49766113250947636708,
)
_DERF1_B = (
    # block 0 (B 0-12)
    -0.00000000029734388465, 0.00000000269776334046, -0.00000000640788827665,
    -0.00000001667820132100, -0.00000021854388148686, 0.00000266246030457984,
    0.00001612722157047886, -0.00025616361025506629, 0.00015380842432375365,
    0.00815533022524927908, -0.01402283663896319337, -0.19746892495383021487,
    0.71511720328842845913,
    # block 1 (B 13-25)
    -0.00000000001951073787, -0.00000000032302692214, 0.00000000522461866919,
    0.00000000342940918551, -0.00000035772874310272, 0.00000019999935792654,
    0.00002687044575042908, -0.00011843240273775776, -0.00080991728956032271,
    0.00661062970502241174, 0.00909530922354827295, -0.20160072778491013140,
    0.51169696718727644908,
    # block 2 (B 26-38)
    0.00000000003147682272, -0.00000000048465972408, 0.00000000063675740242,
    0.00000003377623323271, -0.00000015451139637086, -0.00000203340624738438,
    0.00001947204525295057, 0.00002854147231653228, -0.00101565063152200272,
    0.00271187003520095655, 0.02328095035422810727, -0.16725021123116877197,
    0.32490054966649436974,
    # block 3 (B 39-51)
    0.00000000002319363370, -0.00000000006303206648, -0.00000000264888267434,
    0.00000002050708040581, 0.00000011371857327578, -0.00000211211337219663,
    0.00000368797328322935, 0.00009823686253424796, -0.00065860243990455368,
    -0.00075285814895230877, 0.02585434424202960464, -0.11637092784486193258,
    0.18267336775296612024,
    # block 4 (B 52-64)
    -0.00000000000367789363, 0.00000000020876046746, -0.00000000193319027226,
    -0.00000000435953392472, 0.00000018006992266137, -0.00000078441223763969,
    -0.00000675407647949153, 0.00008428418334440096, -0.00017604388937031815,
    -0.00239729611435071610, 0.02064129023876022970, -0.06905562880005864105,
    0.09084526782065478489,
)
_DERF1_A_BLOCKS = jnp.asarray(_DERF1_A).reshape(5, 13)
_DERF1_B_BLOCKS = jnp.asarray(_DERF1_B).reshape(5, 13)


def GRAUPEL_INIT():
    """Initialize physical constants and switches from CLUBB module parameters."""
    s = {}
    s['nnuccd_reduce_coef'] = jnp.zeros((), dtype=jnp.float32)
    s['nnuccc_reduce_coef'] = jnp.zeros((), dtype=jnp.float32)
    s['pi'] = jnp.zeros((), dtype=jnp.float32)
    s['sqrtpi'] = jnp.zeros((), dtype=jnp.float32)
    s['doicemicro'] = jnp.zeros((), dtype=jnp.bool_)
    s['dograupel'] = jnp.zeros((), dtype=jnp.bool_)
    s['dohail'] = jnp.zeros((), dtype=jnp.bool_)
    s['dosb_warm_rain'] = jnp.zeros((), dtype=jnp.bool_)
    s['dopredictnc'] = jnp.zeros((), dtype=jnp.bool_)
    s['dosubgridw'] = jnp.zeros((), dtype=jnp.bool_)
    s['doarcticicenucl'] = jnp.zeros((), dtype=jnp.bool_)
    s['docloudedgeactivation'] = jnp.zeros((), dtype=jnp.bool_)
    s['dofix_pgam'] = jnp.zeros((), dtype=jnp.bool_)
    s['aerosol_mode'] = jnp.zeros((), dtype=jnp.int32)
    s['nc0'] = jnp.zeros((), dtype=jnp.float32)
    s['ccnconst'] = jnp.zeros((), dtype=jnp.float32)
    s['ccnexpnt'] = jnp.zeros((), dtype=jnp.float32)
    s['aer_rm1'] = jnp.zeros((), dtype=jnp.float32)
    s['aer_rm2'] = jnp.zeros((), dtype=jnp.float32)
    s['aer_n1'] = jnp.zeros((), dtype=jnp.float32)
    s['aer_n2'] = jnp.zeros((), dtype=jnp.float32)
    s['aer_sig1'] = jnp.zeros((), dtype=jnp.float32)
    s['aer_sig2'] = jnp.zeros((), dtype=jnp.float32)
    s['pgam_fixed'] = jnp.zeros((), dtype=jnp.float32)
    s['iact'] = jnp.zeros((), dtype=jnp.int32)
    s['inum'] = jnp.zeros((), dtype=jnp.int32)
    s['ndcnst'] = jnp.zeros((), dtype=jnp.float32)
    s['iliq'] = jnp.zeros((), dtype=jnp.int32)
    s['inuc'] = jnp.zeros((), dtype=jnp.int32)
    s['ibase'] = jnp.zeros((), dtype=jnp.int32)
    s['isub'] = jnp.zeros((), dtype=jnp.int32)
    s['igraup'] = jnp.zeros((), dtype=jnp.int32)
    s['ihail'] = jnp.zeros((), dtype=jnp.int32)
    s['irain'] = jnp.zeros((), dtype=jnp.int32)
    s['isatadj'] = jnp.zeros((), dtype=jnp.int32)
    s['ai'] = jnp.zeros((), dtype=jnp.float32)
    s['ac'] = jnp.zeros((), dtype=jnp.float32)
    s['as'] = jnp.zeros((), dtype=jnp.float32)
    s['ar'] = jnp.zeros((), dtype=jnp.float32)
    s['ag'] = jnp.zeros((), dtype=jnp.float32)
    s['bi'] = jnp.zeros((), dtype=jnp.float32)
    s['bc'] = jnp.zeros((), dtype=jnp.float32)
    s['bs'] = jnp.zeros((), dtype=jnp.float32)
    s['br'] = jnp.zeros((), dtype=jnp.float32)
    s['bg'] = jnp.zeros((), dtype=jnp.float32)
    s['r'] = jnp.zeros((), dtype=jnp.float32)
    s['rhosu'] = jnp.zeros((), dtype=jnp.float32)
    s['rhow'] = jnp.zeros((), dtype=jnp.float32)
    s['rhoi'] = jnp.zeros((), dtype=jnp.float32)
    s['rhosn'] = jnp.zeros((), dtype=jnp.float32)
    s['rhog'] = jnp.zeros((), dtype=jnp.float32)
    s['aimm'] = jnp.zeros((), dtype=jnp.float32)
    s['bimm'] = jnp.zeros((), dtype=jnp.float32)
    s['ecr'] = jnp.zeros((), dtype=jnp.float32)
    s['dcs'] = jnp.zeros((), dtype=jnp.float32)
    s['mi0'] = jnp.zeros((), dtype=jnp.float32)
    s['mg0'] = jnp.zeros((), dtype=jnp.float32)
    s['f1s'] = jnp.zeros((), dtype=jnp.float32)
    s['f2s'] = jnp.zeros((), dtype=jnp.float32)
    s['f1r'] = jnp.zeros((), dtype=jnp.float32)
    s['f2r'] = jnp.zeros((), dtype=jnp.float32)
    s['g'] = jnp.zeros((), dtype=jnp.float32)
    s['qsmall'] = jnp.zeros((), dtype=jnp.float32)
    s['ci'] = jnp.zeros((), dtype=jnp.float32)
    s['di'] = jnp.zeros((), dtype=jnp.float32)
    s['cs'] = jnp.zeros((), dtype=jnp.float32)
    s['ds'] = jnp.zeros((), dtype=jnp.float32)
    s['cg'] = jnp.zeros((), dtype=jnp.float32)
    s['dg'] = jnp.zeros((), dtype=jnp.float32)
    s['eii'] = jnp.zeros((), dtype=jnp.float32)
    s['eci'] = jnp.zeros((), dtype=jnp.float32)
    s['rin'] = jnp.zeros((), dtype=jnp.float32)
    s['tmelt'] = jnp.zeros((), dtype=jnp.float32)
    s['cpw'] = jnp.zeros((), dtype=jnp.float32)
    s['c1'] = jnp.zeros((), dtype=jnp.float32)
    s['k1'] = jnp.zeros((), dtype=jnp.float32)
    s['mw'] = jnp.zeros((), dtype=jnp.float32)
    s['osm'] = jnp.zeros((), dtype=jnp.float32)
    s['vi'] = jnp.zeros((), dtype=jnp.float32)
    s['epsm'] = jnp.zeros((), dtype=jnp.float32)
    s['rhoa'] = jnp.zeros((), dtype=jnp.float32)
    s['map'] = jnp.zeros((), dtype=jnp.float32)
    s['ma'] = jnp.zeros((), dtype=jnp.float32)
    s['rr'] = jnp.zeros((), dtype=jnp.float32)
    s['bact'] = jnp.zeros((), dtype=jnp.float32)
    s['rm1'] = jnp.zeros((), dtype=jnp.float32)
    s['rm2'] = jnp.zeros((), dtype=jnp.float32)
    s['nanew1'] = jnp.zeros((), dtype=jnp.float32)
    s['nanew2'] = jnp.zeros((), dtype=jnp.float32)
    s['sig1'] = jnp.zeros((), dtype=jnp.float32)
    s['sig2'] = jnp.zeros((), dtype=jnp.float32)
    s['f11'] = jnp.zeros((), dtype=jnp.float32)
    s['f12'] = jnp.zeros((), dtype=jnp.float32)
    s['f21'] = jnp.zeros((), dtype=jnp.float32)
    s['f22'] = jnp.zeros((), dtype=jnp.float32)
    s['mmult'] = jnp.zeros((), dtype=jnp.float32)
    s['lammaxi'] = jnp.zeros((), dtype=jnp.float32)
    s['lammini'] = jnp.zeros((), dtype=jnp.float32)
    s['lammaxr'] = jnp.zeros((), dtype=jnp.float32)
    s['lamminr'] = jnp.zeros((), dtype=jnp.float32)
    s['lammaxs'] = jnp.zeros((), dtype=jnp.float32)
    s['lammins'] = jnp.zeros((), dtype=jnp.float32)
    s['lammaxg'] = jnp.zeros((), dtype=jnp.float32)
    s['lamming'] = jnp.zeros((), dtype=jnp.float32)
    s['cons1'] = jnp.zeros((), dtype=jnp.float32)
    s['cons2'] = jnp.zeros((), dtype=jnp.float32)
    s['cons3'] = jnp.zeros((), dtype=jnp.float32)
    s['cons4'] = jnp.zeros((), dtype=jnp.float32)
    s['cons5'] = jnp.zeros((), dtype=jnp.float32)
    s['cons6'] = jnp.zeros((), dtype=jnp.float32)
    s['cons7'] = jnp.zeros((), dtype=jnp.float32)
    s['cons8'] = jnp.zeros((), dtype=jnp.float32)
    s['cons9'] = jnp.zeros((), dtype=jnp.float32)
    s['cons10'] = jnp.zeros((), dtype=jnp.float32)
    s['cons11'] = jnp.zeros((), dtype=jnp.float32)
    s['cons12'] = jnp.zeros((), dtype=jnp.float32)
    s['cons13'] = jnp.zeros((), dtype=jnp.float32)
    s['cons14'] = jnp.zeros((), dtype=jnp.float32)
    s['cons15'] = jnp.zeros((), dtype=jnp.float32)
    s['cons16'] = jnp.zeros((), dtype=jnp.float32)
    s['cons17'] = jnp.zeros((), dtype=jnp.float32)
    s['cons18'] = jnp.zeros((), dtype=jnp.float32)
    s['cons19'] = jnp.zeros((), dtype=jnp.float32)
    s['cons20'] = jnp.zeros((), dtype=jnp.float32)
    s['cons21'] = jnp.zeros((), dtype=jnp.float32)
    s['cons22'] = jnp.zeros((), dtype=jnp.float32)
    s['cons23'] = jnp.zeros((), dtype=jnp.float32)
    s['cons24'] = jnp.zeros((), dtype=jnp.float32)
    s['cons25'] = jnp.zeros((), dtype=jnp.float32)
    s['cons26'] = jnp.zeros((), dtype=jnp.float32)
    s['cons27'] = jnp.zeros((), dtype=jnp.float32)
    s['cons28'] = jnp.zeros((), dtype=jnp.float32)
    s['cons29'] = jnp.zeros((), dtype=jnp.float32)
    s['cons30'] = jnp.zeros((), dtype=jnp.float32)
    s['cons31'] = jnp.zeros((), dtype=jnp.float32)
    s['cons32'] = jnp.zeros((), dtype=jnp.float32)
    s['cons33'] = jnp.zeros((), dtype=jnp.float32)
    s['cons34'] = jnp.zeros((), dtype=jnp.float32)
    s['cons35'] = jnp.zeros((), dtype=jnp.float32)
    s['cons36'] = jnp.zeros((), dtype=jnp.float32)
    s['cons37'] = jnp.zeros((), dtype=jnp.float32)
    s['cons38'] = jnp.zeros((), dtype=jnp.float32)
    s['cons39'] = jnp.zeros((), dtype=jnp.float32)
    s['cons40'] = jnp.zeros((), dtype=jnp.float32)
    s['cons41'] = jnp.zeros((), dtype=jnp.float32)
    s['dnu'] = jnp.zeros((16,), dtype=jnp.float32)
    s['cloud_frac_thresh'] = jnp.zeros((), dtype=jnp.float32)
    s['_pc'] = jnp.zeros((), dtype=jnp.int32)
    s['n'] = jnp.zeros((), dtype=jnp.int32)
    s['i'] = jnp.zeros((), dtype=jnp.int32)
    s['lv'] = jnp.zeros((), dtype=jnp.float64)
    s['ls'] = jnp.zeros((), dtype=jnp.float64)
    s['cp'] = jnp.zeros((), dtype=jnp.float64)
    s['rv'] = jnp.zeros((), dtype=jnp.float64)
    s['rd'] = jnp.zeros((), dtype=jnp.float64)
    s['t_freeze_k'] = jnp.zeros((), dtype=jnp.float64)
    s['rho_lw'] = jnp.zeros((), dtype=jnp.float64)
    s['grav'] = jnp.zeros((), dtype=jnp.float64)
    s['ep_2'] = jnp.zeros((), dtype=jnp.float64)
    s['nnuccd_reduce_coef'] = F(NNUCCD_REDUCE_COEF)
    s['nnuccc_reduce_coef'] = F(NNUCCC_REDUCE_COEF)
    s['pi'] = F(F(3.141592653589793))
    s['sqrtpi'] = F(F(0.9189385332046728))
    s['cloud_frac_thresh'] = F(F(0.005))
    s['lv'] = D(C.Lv)
    s['ls'] = D(C.Ls)
    s['cp'] = D(C.Cp)
    s['rv'] = D(C.Rv)
    s['rd'] = D(C.Rd)
    s['t_freeze_k'] = D(C.T_freeze_K)
    s['rho_lw'] = D(C.rho_lw)
    s['grav'] = D(C.grav)
    s['ep_2'] = D(C.ep)
    s['doicemicro'] = B(P.l_ice_microphys)
    s['dograupel'] = B(P.l_ice_microphys and P.l_graupel)
    s['dohail'] = B(P.l_hail)
    s['dosb_warm_rain'] = B(P.l_seifert_beheng)
    s['dopredictnc'] = B(P.l_predict_Nc)
    s['dosubgridw'] = B(P.l_subgrid_w)
    s['doarcticicenucl'] = B(P.l_arctic_nucl)
    s['docloudedgeactivation'] = B(P.l_cloud_edge_activation)
    s['dofix_pgam'] = B(P.l_fix_pgam)
    s["nc0"] = F(P.Nc0_in_cloud / 1.e6)
    s["aerosol_mode"] = I({"morrison_no_aerosol":0,"morrison_power_law":1,"morrison_lognormal":2}[P.specify_aerosol])
    s['ccnconst'] = F(120.0)
    s['ccnexpnt'] = F(0.4)
    s['aer_rm1'] = F(1.1e-08)
    s['aer_rm2'] = F(6e-08)
    s['aer_n1'] = F(125000000.0)
    s['aer_n2'] = F(65000000.0)
    s['aer_sig1'] = F(1.2)
    s['aer_sig2'] = F(1.7)
    s['pgam_fixed'] = F(5.0)
    old_1 = s['inum']
    s['inum'] = I(1)
    s['inum'] = jnp.where(s["_pc"] == 0, s['inum'], old_1)
    def yes_2(s):
        s = dict(s)
        old_3 = s['inum']
        s['inum'] = I(0)
        s['inum'] = jnp.where(s["_pc"] == 0, s['inum'], old_3)
        return s
    def no_2(s):
        s = dict(s)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['dopredictnc']), yes_2, no_2, s)
    old_4 = s['ndcnst']
    s['ndcnst'] = F(s['nc0'])
    s['ndcnst'] = jnp.where(s["_pc"] == 0, s['ndcnst'], old_4)
    def yes_5(s):
        s = dict(s)
        old_6 = s['iact']
        s['iact'] = I(2)
        s['iact'] = jnp.where(s["_pc"] == 0, s['iact'], old_6)
        return s
    def no_5(s):
        s = dict(s)
        def yes_7(s):
            s = dict(s)
            old_8 = s['iact']
            s['iact'] = I(1)
            s['iact'] = jnp.where(s["_pc"] == 0, s['iact'], old_8)
            return s
        def no_7(s):
            s = dict(s)
            old_9 = s['iact']
            s['iact'] = I(0)
            s['iact'] = jnp.where(s["_pc"] == 0, s['iact'], old_9)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['aerosol_mode'] == 1)), yes_7, no_7, s)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['aerosol_mode'] == 2)), yes_5, no_5, s)
    def yes_10(s):
        s = dict(s)
        old_11 = s['ibase']
        s['ibase'] = I(2)
        s['ibase'] = jnp.where(s["_pc"] == 0, s['ibase'], old_11)
        return s
    def no_10(s):
        s = dict(s)
        old_12 = s['ibase']
        s['ibase'] = I(1)
        s['ibase'] = jnp.where(s["_pc"] == 0, s['ibase'], old_12)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['docloudedgeactivation']), yes_10, no_10, s)
    def yes_13(s):
        s = dict(s)
        old_14 = s['isub']
        s['isub'] = I(0)
        s['isub'] = jnp.where(s["_pc"] == 0, s['isub'], old_14)
        return s
    def no_13(s):
        s = dict(s)
        old_15 = s['isub']
        s['isub'] = I(1)
        s['isub'] = jnp.where(s["_pc"] == 0, s['isub'], old_15)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['dosubgridw']), yes_13, no_13, s)
    def yes_16(s):
        s = dict(s)
        old_17 = s['iliq']
        s['iliq'] = I(0)
        s['iliq'] = jnp.where(s["_pc"] == 0, s['iliq'], old_17)
        return s
    def no_16(s):
        s = dict(s)
        old_18 = s['iliq']
        s['iliq'] = I(1)
        s['iliq'] = jnp.where(s["_pc"] == 0, s['iliq'], old_18)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['doicemicro']), yes_16, no_16, s)
    def yes_19(s):
        s = dict(s)
        old_20 = s['inuc']
        s['inuc'] = I(1)
        s['inuc'] = jnp.where(s["_pc"] == 0, s['inuc'], old_20)
        return s
    def no_19(s):
        s = dict(s)
        old_21 = s['inuc']
        s['inuc'] = I(0)
        s['inuc'] = jnp.where(s["_pc"] == 0, s['inuc'], old_21)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['doarcticicenucl']), yes_19, no_19, s)
    def yes_22(s):
        s = dict(s)
        old_23 = s['igraup']
        s['igraup'] = I(0)
        s['igraup'] = jnp.where(s["_pc"] == 0, s['igraup'], old_23)
        return s
    def no_22(s):
        s = dict(s)
        old_24 = s['igraup']
        s['igraup'] = I(1)
        s['igraup'] = jnp.where(s["_pc"] == 0, s['igraup'], old_24)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['dograupel']), yes_22, no_22, s)
    def yes_25(s):
        s = dict(s)
        old_26 = s['ihail']
        s['ihail'] = I(1)
        s['ihail'] = jnp.where(s["_pc"] == 0, s['ihail'], old_26)
        return s
    def no_25(s):
        s = dict(s)
        old_27 = s['ihail']
        s['ihail'] = I(0)
        s['ihail'] = jnp.where(s["_pc"] == 0, s['ihail'], old_27)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['dohail']), yes_25, no_25, s)
    def yes_28(s):
        s = dict(s)
        old_29 = s['irain']
        s['irain'] = I(1)
        s['irain'] = jnp.where(s["_pc"] == 0, s['irain'], old_29)
        return s
    def no_28(s):
        s = dict(s)
        old_30 = s['irain']
        s['irain'] = I(0)
        s['irain'] = jnp.where(s["_pc"] == 0, s['irain'], old_30)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['dosb_warm_rain']), yes_28, no_28, s)
    old_31 = s['isatadj']
    s['isatadj'] = I(0)
    s['isatadj'] = jnp.where(s["_pc"] == 0, s['isatadj'], old_31)
    old_32 = s['ai']
    s['ai'] = F(F(700.0))
    s['ai'] = jnp.where(s["_pc"] == 0, s['ai'], old_32)
    old_33 = s['ac']
    s['ac'] = F(F(30000000.0))
    s['ac'] = jnp.where(s["_pc"] == 0, s['ac'], old_33)
    old_34 = s['as']
    s['as'] = F(F(11.72))
    s['as'] = jnp.where(s["_pc"] == 0, s['as'], old_34)
    old_35 = s['ar']
    s['ar'] = F(F(841.99667))
    s['ar'] = jnp.where(s["_pc"] == 0, s['ar'], old_35)
    old_36 = s['bi']
    s['bi'] = F(F(1.0))
    s['bi'] = jnp.where(s["_pc"] == 0, s['bi'], old_36)
    old_37 = s['bc']
    s['bc'] = F(F(2.0))
    s['bc'] = jnp.where(s["_pc"] == 0, s['bc'], old_37)
    old_38 = s['bs']
    s['bs'] = F(F(0.41))
    s['bs'] = jnp.where(s["_pc"] == 0, s['bs'], old_38)
    old_39 = s['br']
    s['br'] = F(F(0.8))
    s['br'] = jnp.where(s["_pc"] == 0, s['br'], old_39)
    def yes_40(s):
        s = dict(s)
        old_41 = s['ag']
        s['ag'] = F(F(19.3))
        s['ag'] = jnp.where(s["_pc"] == 0, s['ag'], old_41)
        old_42 = s['bg']
        s['bg'] = F(F(0.37))
        s['bg'] = jnp.where(s["_pc"] == 0, s['bg'], old_42)
        return s
    def no_40(s):
        s = dict(s)
        old_43 = s['ag']
        s['ag'] = F(F(114.5))
        s['ag'] = jnp.where(s["_pc"] == 0, s['ag'], old_43)
        old_44 = s['bg']
        s['bg'] = F(F(0.5))
        s['bg'] = jnp.where(s["_pc"] == 0, s['bg'], old_44)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['ihail'] == 0)), yes_40, no_40, s)
    old_45 = s['r']
    s['r'] = F(s['rd'])
    s['r'] = jnp.where(s["_pc"] == 0, s['r'], old_45)
    old_46 = s['rhow']
    s['rhow'] = F(s['rho_lw'])
    s['rhow'] = jnp.where(s["_pc"] == 0, s['rhow'], old_46)
    old_47 = s['tmelt']
    s['tmelt'] = F(s['t_freeze_k'])
    s['tmelt'] = jnp.where(s["_pc"] == 0, s['tmelt'], old_47)
    old_48 = s['rhosu']
    s['rhosu'] = F(_div(F(85000.0), _arith(s['r'], s['tmelt'], 'mul')))
    s['rhosu'] = jnp.where(s["_pc"] == 0, s['rhosu'], old_48)
    old_49 = s['rhosu']
    s['rhosu'] = F(_div(F(85000.0), _arith(s['r'], s['tmelt'], 'mul')))
    s['rhosu'] = jnp.where(s["_pc"] == 0, s['rhosu'], old_49)
    old_50 = s['rhow']
    s['rhow'] = F(F(997.0))
    s['rhow'] = jnp.where(s["_pc"] == 0, s['rhow'], old_50)
    old_51 = s['rhoi']
    s['rhoi'] = F(F(500.0))
    s['rhoi'] = jnp.where(s["_pc"] == 0, s['rhoi'], old_51)
    old_52 = s['rhosn']
    s['rhosn'] = F(F(100.0))
    s['rhosn'] = jnp.where(s["_pc"] == 0, s['rhosn'], old_52)
    def yes_53(s):
        s = dict(s)
        old_54 = s['rhog']
        s['rhog'] = F(F(400.0))
        s['rhog'] = jnp.where(s["_pc"] == 0, s['rhog'], old_54)
        return s
    def no_53(s):
        s = dict(s)
        old_55 = s['rhog']
        s['rhog'] = F(F(900.0))
        s['rhog'] = jnp.where(s["_pc"] == 0, s['rhog'], old_55)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['ihail'] == 0)), yes_53, no_53, s)
    old_56 = s['aimm']
    s['aimm'] = F(F(0.66))
    s['aimm'] = jnp.where(s["_pc"] == 0, s['aimm'], old_56)
    old_57 = s['bimm']
    s['bimm'] = F(F(100.0))
    s['bimm'] = jnp.where(s["_pc"] == 0, s['bimm'], old_57)
    old_58 = s['ecr']
    s['ecr'] = F(F(1.0))
    s['ecr'] = jnp.where(s["_pc"] == 0, s['ecr'], old_58)
    old_59 = s['dcs']
    s['dcs'] = F(F(0.000125))
    s['dcs'] = jnp.where(s["_pc"] == 0, s['dcs'], old_59)
    old_60 = s['mi0']
    s['mi0'] = F(_arith(_arith(_arith(_div(F(4.0), F(3.0)), s['pi'], 'mul'), s['rhoi'], 'mul'), _arith(F(1e-05), 3, 'pow'), 'mul'))
    s['mi0'] = jnp.where(s["_pc"] == 0, s['mi0'], old_60)
    old_61 = s['mg0']
    s['mg0'] = F(F(1.6e-10))
    s['mg0'] = jnp.where(s["_pc"] == 0, s['mg0'], old_61)
    old_62 = s['f1s']
    s['f1s'] = F(F(0.86))
    s['f1s'] = jnp.where(s["_pc"] == 0, s['f1s'], old_62)
    old_63 = s['f2s']
    s['f2s'] = F(F(0.28))
    s['f2s'] = jnp.where(s["_pc"] == 0, s['f2s'], old_63)
    old_64 = s['f1r']
    s['f1r'] = F(F(0.78))
    s['f1r'] = jnp.where(s["_pc"] == 0, s['f1r'], old_64)
    old_65 = s['f2r']
    s['f2r'] = F(F(0.308))
    s['f2r'] = jnp.where(s["_pc"] == 0, s['f2r'], old_65)
    old_66 = s['g']
    s['g'] = F(s['grav'])
    s['g'] = jnp.where(s["_pc"] == 0, s['g'], old_66)
    old_67 = s['qsmall']
    s['qsmall'] = F(F(1e-14))
    s['qsmall'] = jnp.where(s["_pc"] == 0, s['qsmall'], old_67)
    old_68 = s['eii']
    s['eii'] = F(F(0.1))
    s['eii'] = jnp.where(s["_pc"] == 0, s['eii'], old_68)
    old_69 = s['eci']
    s['eci'] = F(F(0.7))
    s['eci'] = jnp.where(s["_pc"] == 0, s['eci'], old_69)
    old_70 = s['cpw']
    s['cpw'] = F(F(4218.0))
    s['cpw'] = jnp.where(s["_pc"] == 0, s['cpw'], old_70)
    old_71 = s['ci']
    s['ci'] = F(_div(_arith(s['rhoi'], s['pi'], 'mul'), F(6.0)))
    s['ci'] = jnp.where(s["_pc"] == 0, s['ci'], old_71)
    old_72 = s['di']
    s['di'] = F(F(3.0))
    s['di'] = jnp.where(s["_pc"] == 0, s['di'], old_72)
    old_73 = s['cs']
    s['cs'] = F(_div(_arith(s['rhosn'], s['pi'], 'mul'), F(6.0)))
    s['cs'] = jnp.where(s["_pc"] == 0, s['cs'], old_73)
    old_74 = s['ds']
    s['ds'] = F(F(3.0))
    s['ds'] = jnp.where(s["_pc"] == 0, s['ds'], old_74)
    old_75 = s['cg']
    s['cg'] = F(_div(_arith(s['rhog'], s['pi'], 'mul'), F(6.0)))
    s['cg'] = jnp.where(s["_pc"] == 0, s['cg'], old_75)
    old_76 = s['dg']
    s['dg'] = F(F(3.0))
    s['dg'] = jnp.where(s["_pc"] == 0, s['dg'], old_76)
    old_77 = s['rin']
    s['rin'] = F(F(1e-07))
    s['rin'] = jnp.where(s["_pc"] == 0, s['rin'], old_77)
    old_78 = s['mmult']
    s['mmult'] = F(_arith(_arith(_arith(_div(F(4.0), F(3.0)), s['pi'], 'mul'), s['rhoi'], 'mul'), _arith(F(5e-06), 3, 'pow'), 'mul'))
    s['mmult'] = jnp.where(s["_pc"] == 0, s['mmult'], old_78)
    old_79 = s['lammaxi']
    s['lammaxi'] = F(_div(F(1.0), F(1e-06)))
    s['lammaxi'] = jnp.where(s["_pc"] == 0, s['lammaxi'], old_79)
    old_80 = s['lammini']
    s['lammini'] = F(_div(F(1.0), _arith(_arith(F(2.0), s['dcs'], 'mul'), F(0.0001), 'add')))
    s['lammini'] = jnp.where(s["_pc"] == 0, s['lammini'], old_80)
    old_81 = s['lammaxr']
    s['lammaxr'] = F(_div(F(1.0), F(2e-05)))
    s['lammaxr'] = jnp.where(s["_pc"] == 0, s['lammaxr'], old_81)
    old_82 = s['lamminr']
    s['lamminr'] = F(_div(F(1.0), F(0.0028)))
    s['lamminr'] = jnp.where(s["_pc"] == 0, s['lamminr'], old_82)
    old_83 = s['lammaxs']
    s['lammaxs'] = F(_div(F(1.0), F(1e-05)))
    s['lammaxs'] = jnp.where(s["_pc"] == 0, s['lammaxs'], old_83)
    old_84 = s['lammins']
    s['lammins'] = F(_div(F(1.0), F(0.002)))
    s['lammins'] = jnp.where(s["_pc"] == 0, s['lammins'], old_84)
    old_85 = s['lammaxg']
    s['lammaxg'] = F(_div(F(1.0), F(2e-05)))
    s['lammaxg'] = jnp.where(s["_pc"] == 0, s['lammaxg'], old_85)
    old_86 = s['lamming']
    s['lamming'] = F(_div(F(1.0), F(0.002)))
    s['lamming'] = jnp.where(s["_pc"] == 0, s['lamming'], old_86)
    old_87 = s['k1']
    s['k1'] = F(s['ccnexpnt'])
    s['k1'] = jnp.where(s["_pc"] == 0, s['k1'], old_87)
    old_88 = s['c1']
    s['c1'] = F(s['ccnconst'])
    s['c1'] = jnp.where(s["_pc"] == 0, s['c1'], old_88)
    old_89 = s['mw']
    s['mw'] = F(F(0.018))
    s['mw'] = jnp.where(s["_pc"] == 0, s['mw'], old_89)
    old_90 = s['osm']
    s['osm'] = F(F(1.0))
    s['osm'] = jnp.where(s["_pc"] == 0, s['osm'], old_90)
    old_91 = s['vi']
    s['vi'] = F(F(3.0))
    s['vi'] = jnp.where(s["_pc"] == 0, s['vi'], old_91)
    old_92 = s['epsm']
    s['epsm'] = F(F(0.7))
    s['epsm'] = jnp.where(s["_pc"] == 0, s['epsm'], old_92)
    old_93 = s['rhoa']
    s['rhoa'] = F(F(1777.0))
    s['rhoa'] = jnp.where(s["_pc"] == 0, s['rhoa'], old_93)
    old_94 = s['map']
    s['map'] = F(F(0.132))
    s['map'] = jnp.where(s["_pc"] == 0, s['map'], old_94)
    old_95 = s['ma']
    s['ma'] = F(F(0.0284))
    s['ma'] = jnp.where(s["_pc"] == 0, s['ma'], old_95)
    old_96 = s['rr']
    s['rr'] = F(F(8.3187))
    s['rr'] = jnp.where(s["_pc"] == 0, s['rr'], old_96)
    old_97 = s['bact']
    s['bact'] = F(_div(_arith(_arith(_arith(_arith(s['vi'], s['osm'], 'mul'), s['epsm'], 'mul'), s['mw'], 'mul'), s['rhoa'], 'mul'), _arith(s['map'], s['rhow'], 'mul')))
    s['bact'] = jnp.where(s["_pc"] == 0, s['bact'], old_97)
    old_98 = s['rm1']
    s['rm1'] = F(s['aer_rm1'])
    s['rm1'] = jnp.where(s["_pc"] == 0, s['rm1'], old_98)
    old_99 = s['sig1']
    s['sig1'] = F(s['aer_sig1'])
    s['sig1'] = jnp.where(s["_pc"] == 0, s['sig1'], old_99)
    old_100 = s['nanew1']
    s['nanew1'] = F(s['aer_n1'])
    s['nanew1'] = jnp.where(s["_pc"] == 0, s['nanew1'], old_100)
    old_101 = s['f11']
    s['f11'] = F(_arith(F(0.5), _intrinsic('exp', _arith(F(2.5), _arith(_intrinsic('log', s['sig1']), 2, 'pow'), 'mul')), 'mul'))
    s['f11'] = jnp.where(s["_pc"] == 0, s['f11'], old_101)
    old_102 = s['f21']
    s['f21'] = F(_arith(F(1.0), _arith(F(0.25), _intrinsic('log', s['sig1']), 'mul'), 'add'))
    s['f21'] = jnp.where(s["_pc"] == 0, s['f21'], old_102)
    old_103 = s['rm2']
    s['rm2'] = F(s['aer_rm2'])
    s['rm2'] = jnp.where(s["_pc"] == 0, s['rm2'], old_103)
    old_104 = s['sig2']
    s['sig2'] = F(s['aer_sig2'])
    s['sig2'] = jnp.where(s["_pc"] == 0, s['sig2'], old_104)
    old_105 = s['nanew2']
    s['nanew2'] = F(s['aer_n2'])
    s['nanew2'] = jnp.where(s["_pc"] == 0, s['nanew2'], old_105)
    old_106 = s['f12']
    s['f12'] = F(_arith(F(0.5), _intrinsic('exp', _arith(F(2.5), _arith(_intrinsic('log', s['sig2']), 2, 'pow'), 'mul')), 'mul'))
    s['f12'] = jnp.where(s["_pc"] == 0, s['f12'], old_106)
    old_107 = s['f22']
    s['f22'] = F(_arith(F(1.0), _arith(F(0.25), _intrinsic('log', s['sig2']), 'mul'), 'add'))
    s['f22'] = jnp.where(s["_pc"] == 0, s['f22'], old_107)
    old_108 = s['cons1']
    s['cons1'] = F(_arith(GAMMA(_arith(F(1.0), s['ds'], 'add')), s['cs'], 'mul'))
    s['cons1'] = jnp.where(s["_pc"] == 0, s['cons1'], old_108)
    old_109 = s['cons2']
    s['cons2'] = F(_arith(GAMMA(_arith(F(1.0), s['dg'], 'add')), s['cg'], 'mul'))
    s['cons2'] = jnp.where(s["_pc"] == 0, s['cons2'], old_109)
    old_110 = s['cons3']
    s['cons3'] = F(_div(GAMMA(_arith(F(4.0), s['bs'], 'add')), F(6.0)))
    s['cons3'] = jnp.where(s["_pc"] == 0, s['cons3'], old_110)
    old_111 = s['cons4']
    s['cons4'] = F(_div(GAMMA(_arith(F(4.0), s['br'], 'add')), F(6.0)))
    s['cons4'] = jnp.where(s["_pc"] == 0, s['cons4'], old_111)
    old_112 = s['cons5']
    s['cons5'] = F(GAMMA(_arith(F(1.0), s['bs'], 'add')))
    s['cons5'] = jnp.where(s["_pc"] == 0, s['cons5'], old_112)
    old_113 = s['cons6']
    s['cons6'] = F(GAMMA(_arith(F(1.0), s['br'], 'add')))
    s['cons6'] = jnp.where(s["_pc"] == 0, s['cons6'], old_113)
    old_114 = s['cons7']
    s['cons7'] = F(_div(GAMMA(_arith(F(4.0), s['bg'], 'add')), F(6.0)))
    s['cons7'] = jnp.where(s["_pc"] == 0, s['cons7'], old_114)
    old_115 = s['cons8']
    s['cons8'] = F(GAMMA(_arith(F(1.0), s['bg'], 'add')))
    s['cons8'] = jnp.where(s["_pc"] == 0, s['cons8'], old_115)
    old_116 = s['cons9']
    s['cons9'] = F(GAMMA(_arith(_div(F(5.0), F(2.0)), _div(s['br'], F(2.0)), 'add')))
    s['cons9'] = jnp.where(s["_pc"] == 0, s['cons9'], old_116)
    old_117 = s['cons10']
    s['cons10'] = F(GAMMA(_arith(_div(F(5.0), F(2.0)), _div(s['bs'], F(2.0)), 'add')))
    s['cons10'] = jnp.where(s["_pc"] == 0, s['cons10'], old_117)
    old_118 = s['cons11']
    s['cons11'] = F(GAMMA(_arith(_div(F(5.0), F(2.0)), _div(s['bg'], F(2.0)), 'add')))
    s['cons11'] = jnp.where(s["_pc"] == 0, s['cons11'], old_118)
    old_119 = s['cons12']
    s['cons12'] = F(_arith(GAMMA(_arith(F(1.0), s['di'], 'add')), s['ci'], 'mul'))
    s['cons12'] = jnp.where(s["_pc"] == 0, s['cons12'], old_119)
    old_120 = s['cons13']
    s['cons13'] = F(_arith(_div(_arith(GAMMA(_arith(s['bs'], F(3.0), 'add')), s['pi'], 'mul'), F(4.0)), s['eci'], 'mul'))
    s['cons13'] = jnp.where(s["_pc"] == 0, s['cons13'], old_120)
    old_121 = s['cons14']
    s['cons14'] = F(_arith(_div(_arith(GAMMA(_arith(s['bg'], F(3.0), 'add')), s['pi'], 'mul'), F(4.0)), s['eci'], 'mul'))
    s['cons14'] = jnp.where(s["_pc"] == 0, s['cons14'], old_121)
    old_122 = s['cons15']
    s['cons15'] = F(_div(_arith(_arith(_arith((-F(1108.0)), s['eii'], 'mul'), _arith(s['pi'], _div(_arith(F(1.0), s['bs'], 'sub'), F(3.0)), 'pow'), 'mul'), _arith(s['rhosn'], _div(_arith((-F(2.0)), s['bs'], 'sub'), F(3.0)), 'pow'), 'mul'), _arith(F(4.0), F(720.0), 'mul')))
    s['cons15'] = jnp.where(s["_pc"] == 0, s['cons15'], old_122)
    old_123 = s['cons16']
    s['cons16'] = F(_arith(_div(_arith(GAMMA(_arith(s['bi'], F(3.0), 'add')), s['pi'], 'mul'), F(4.0)), s['eci'], 'mul'))
    s['cons16'] = jnp.where(s["_pc"] == 0, s['cons16'], old_123)
    old_124 = s['cons17']
    s['cons17'] = F(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(F(4.0), F(2.0), 'mul'), F(3.0), 'mul'), s['rhosu'], 'mul'), s['pi'], 'mul'), s['eci'], 'mul'), s['eci'], 'mul'), GAMMA(_arith(_arith(F(2.0), s['bs'], 'mul'), F(2.0), 'add')), 'mul'), _arith(F(8.0), _arith(s['rhog'], s['rhosn'], 'sub'), 'mul')))
    s['cons17'] = jnp.where(s["_pc"] == 0, s['cons17'], old_124)
    old_125 = s['cons18']
    s['cons18'] = F(_arith(s['rhosn'], s['rhosn'], 'mul'))
    s['cons18'] = jnp.where(s["_pc"] == 0, s['cons18'], old_125)
    old_126 = s['cons19']
    s['cons19'] = F(_arith(s['rhow'], s['rhow'], 'mul'))
    s['cons19'] = jnp.where(s["_pc"] == 0, s['cons19'], old_126)
    old_127 = s['cons20']
    s['cons20'] = F(_arith(_arith(_arith(_arith(F(20.0), s['pi'], 'mul'), s['pi'], 'mul'), s['rhow'], 'mul'), s['bimm'], 'mul'))
    s['cons20'] = jnp.where(s["_pc"] == 0, s['cons20'], old_127)
    old_128 = s['cons21']
    s['cons21'] = F(_div(F(4.0), _arith(s['dcs'], s['rhoi'], 'mul')))
    s['cons21'] = jnp.where(s["_pc"] == 0, s['cons21'], old_128)
    old_129 = s['cons22']
    s['cons22'] = F(_div(_arith(_arith(s['pi'], s['rhoi'], 'mul'), _arith(s['dcs'], 3, 'pow'), 'mul'), F(6.0)))
    s['cons22'] = jnp.where(s["_pc"] == 0, s['cons22'], old_129)
    old_130 = s['cons23']
    s['cons23'] = F(_arith(_arith(_div(s['pi'], F(4.0)), s['eii'], 'mul'), GAMMA(_arith(s['bs'], F(3.0), 'add')), 'mul'))
    s['cons23'] = jnp.where(s["_pc"] == 0, s['cons23'], old_130)
    old_131 = s['cons24']
    s['cons24'] = F(_arith(_arith(_div(s['pi'], F(4.0)), s['ecr'], 'mul'), GAMMA(_arith(s['br'], F(3.0), 'add')), 'mul'))
    s['cons24'] = jnp.where(s["_pc"] == 0, s['cons24'], old_131)
    old_132 = s['cons25']
    s['cons25'] = F(_arith(_arith(_arith(_div(_arith(s['pi'], s['pi'], 'mul'), F(24.0)), s['rhow'], 'mul'), s['ecr'], 'mul'), GAMMA(_arith(s['br'], F(6.0), 'add')), 'mul'))
    s['cons25'] = jnp.where(s["_pc"] == 0, s['cons25'], old_132)
    old_133 = s['cons26']
    s['cons26'] = F(_arith(_div(s['pi'], F(6.0)), s['rhow'], 'mul'))
    s['cons26'] = jnp.where(s["_pc"] == 0, s['cons26'], old_133)
    old_134 = s['cons27']
    s['cons27'] = F(GAMMA(_arith(F(1.0), s['bi'], 'add')))
    s['cons27'] = jnp.where(s["_pc"] == 0, s['cons27'], old_134)
    old_135 = s['cons28']
    s['cons28'] = F(_div(GAMMA(_arith(F(4.0), s['bi'], 'add')), F(6.0)))
    s['cons28'] = jnp.where(s["_pc"] == 0, s['cons28'], old_135)
    old_136 = s['cons29']
    s['cons29'] = F(_arith(_arith(_arith(_div(F(4.0), F(3.0)), s['pi'], 'mul'), s['rhow'], 'mul'), _arith(F(2.5e-05), 3, 'pow'), 'mul'))
    s['cons29'] = jnp.where(s["_pc"] == 0, s['cons29'], old_136)
    old_137 = s['cons30']
    s['cons30'] = F(_arith(_arith(_div(F(4.0), F(3.0)), s['pi'], 'mul'), s['rhow'], 'mul'))
    s['cons30'] = jnp.where(s["_pc"] == 0, s['cons30'], old_137)
    old_138 = s['cons31']
    s['cons31'] = F(_arith(_arith(_arith(s['pi'], s['pi'], 'mul'), s['ecr'], 'mul'), s['rhosn'], 'mul'))
    s['cons31'] = jnp.where(s["_pc"] == 0, s['cons31'], old_138)
    old_139 = s['cons32']
    s['cons32'] = F(_arith(_div(s['pi'], F(2.0)), s['ecr'], 'mul'))
    s['cons32'] = jnp.where(s["_pc"] == 0, s['cons32'], old_139)
    old_140 = s['cons33']
    s['cons33'] = F(_arith(_arith(_arith(s['pi'], s['pi'], 'mul'), s['ecr'], 'mul'), s['rhog'], 'mul'))
    s['cons33'] = jnp.where(s["_pc"] == 0, s['cons33'], old_140)
    old_141 = s['cons34']
    s['cons34'] = F(_arith(_div(F(5.0), F(2.0)), _div(s['br'], F(2.0)), 'add'))
    s['cons34'] = jnp.where(s["_pc"] == 0, s['cons34'], old_141)
    old_142 = s['cons35']
    s['cons35'] = F(_arith(_div(F(5.0), F(2.0)), _div(s['bs'], F(2.0)), 'add'))
    s['cons35'] = jnp.where(s["_pc"] == 0, s['cons35'], old_142)
    old_143 = s['cons36']
    s['cons36'] = F(_arith(_div(F(5.0), F(2.0)), _div(s['bg'], F(2.0)), 'add'))
    s['cons36'] = jnp.where(s["_pc"] == 0, s['cons36'], old_143)
    old_144 = s['cons37']
    s['cons37'] = F(_div(_arith(_arith(F(4.0), s['pi'], 'mul'), F(1.38e-23), 'mul'), _arith(_arith(F(6.0), s['pi'], 'mul'), s['rin'], 'mul')))
    s['cons37'] = jnp.where(s["_pc"] == 0, s['cons37'], old_144)
    old_145 = s['cons38']
    s['cons38'] = F(_arith(_div(_arith(s['pi'], s['pi'], 'mul'), F(3.0)), s['rhow'], 'mul'))
    s['cons38'] = jnp.where(s["_pc"] == 0, s['cons38'], old_145)
    old_146 = s['cons39']
    s['cons39'] = F(_arith(_arith(_div(_arith(s['pi'], s['pi'], 'mul'), F(36.0)), s['rhow'], 'mul'), s['bimm'], 'mul'))
    s['cons39'] = jnp.where(s["_pc"] == 0, s['cons39'], old_146)
    old_147 = s['cons40']
    s['cons40'] = F(_arith(_div(s['pi'], F(6.0)), s['bimm'], 'mul'))
    s['cons40'] = jnp.where(s["_pc"] == 0, s['cons40'], old_147)
    old_148 = s['cons41']
    s['cons41'] = F(_arith(_arith(_arith(s['pi'], s['pi'], 'mul'), s['ecr'], 'mul'), s['rhow'], 'mul'))
    s['cons41'] = jnp.where(s["_pc"] == 0, s['cons41'], old_148)
    old_149 = s['dnu']
    s['dnu'] = s['dnu'].at[1 - 1].set(F((-F(0.557))))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_149)
    old_150 = s['dnu']
    s['dnu'] = s['dnu'].at[2 - 1].set(F((-F(0.557))))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_150)
    old_151 = s['dnu']
    s['dnu'] = s['dnu'].at[3 - 1].set(F((-F(0.43))))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_151)
    old_152 = s['dnu']
    s['dnu'] = s['dnu'].at[4 - 1].set(F((-F(0.307))))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_152)
    old_153 = s['dnu']
    s['dnu'] = s['dnu'].at[5 - 1].set(F((-F(0.186))))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_153)
    old_154 = s['dnu']
    s['dnu'] = s['dnu'].at[6 - 1].set(F((-F(0.067))))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_154)
    old_155 = s['dnu']
    s['dnu'] = s['dnu'].at[7 - 1].set(F(F(0.05)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_155)
    old_156 = s['dnu']
    s['dnu'] = s['dnu'].at[8 - 1].set(F(F(0.167)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_156)
    old_157 = s['dnu']
    s['dnu'] = s['dnu'].at[9 - 1].set(F(F(0.282)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_157)
    old_158 = s['dnu']
    s['dnu'] = s['dnu'].at[10 - 1].set(F(F(0.397)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_158)
    old_159 = s['dnu']
    s['dnu'] = s['dnu'].at[11 - 1].set(F(F(0.512)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_159)
    old_160 = s['dnu']
    s['dnu'] = s['dnu'].at[12 - 1].set(F(F(0.626)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_160)
    old_161 = s['dnu']
    s['dnu'] = s['dnu'].at[13 - 1].set(F(F(0.739)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_161)
    old_162 = s['dnu']
    s['dnu'] = s['dnu'].at[14 - 1].set(F(F(0.853)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_162)
    old_163 = s['dnu']
    s['dnu'] = s['dnu'].at[15 - 1].set(F(F(0.966)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_163)
    old_164 = s['dnu']
    s['dnu'] = s['dnu'].at[16 - 1].set(F(F(0.966)))
    s['dnu'] = jnp.where(s["_pc"] == 0, s['dnu'], old_164)
    return s


def M2005MICRO_GRAUPEL(
        QC3D, QI3D, QNI3D, QR3D, NC3D, NI3D, NS3D, NR3D,
        T3D, QV3D, PRES, RHO, DZQ, W3D, WVAR,
        DT,
        IMS, IME, JMS, JME, KMS, KME,
        ITS, ITE, JTS, JTE, KTS, KTE,
        QG3D, NG3D,
        CF3D):
    """Advance independent columns with the source's process-split scheme.

    Input/inout arguments retain their Fortran names and relative order.
    Output-only arguments (tendencies, fallout, radii and process diagnostics)
    are returned with updated inout fields, keyed by uppercase Fortran names.
    The source initializes its tendency buffers internally, so none are inputs.
    The leading column axis replaces the caller's Fortran loop over columns.
    """
    if (IMS, IME, JMS, JME, KMS, ITS, ITE, JTS, JTE, KTS) != (1,)*10:
        raise ValueError('Morrison requires CLUBB single-column horizontal bounds and KMS=KTS=1')
    if QC3D.ndim != 2 or KME != QC3D.shape[-1] or KTE != KME:
        raise ValueError('Morrison requires batched full columns with KME=KTE=nzt')
    inputs = dict(qc3d=QC3D, qi3d=QI3D, qni3d=QNI3D, qr3d=QR3D,
                  nc3d=NC3D, ni3d=NI3D, ns3d=NS3D, nr3d=NR3D,
                  t3d=T3D, qv3d=QV3D, pres=PRES, rho=RHO,
                  dzq=DZQ, w3d=W3D, wvar=WVAR,
                  qg3d=QG3D, ng3d=NG3D, cf3d=CF3D)
    inputs = {key: jnp.asarray(value, dtype=jnp.float32) for key, value in inputs.items()}
    params = GRAUPEL_INIT()
    output = lax.map(lambda col: _column({**col, 'dt': jnp.float32(DT)}, params), inputs)
    return {name.upper(): value for name, value in output.items()}



def _column(inputs, params):
    nz = inputs["qc3d"].shape[0]
    s = dict(params)
    s['ims'] = jnp.zeros((), dtype=jnp.int32)
    s['ime'] = jnp.zeros((), dtype=jnp.int32)
    s['jms'] = jnp.zeros((), dtype=jnp.int32)
    s['jme'] = jnp.zeros((), dtype=jnp.int32)
    s['kms'] = jnp.zeros((), dtype=jnp.int32)
    s['kme'] = jnp.zeros((), dtype=jnp.int32)
    s['its'] = jnp.zeros((), dtype=jnp.int32)
    s['ite'] = jnp.zeros((), dtype=jnp.int32)
    s['jts'] = jnp.zeros((), dtype=jnp.int32)
    s['jte'] = jnp.zeros((), dtype=jnp.int32)
    s['kts'] = jnp.zeros((), dtype=jnp.int32)
    s['kte'] = jnp.zeros((), dtype=jnp.int32)
    s['qc3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qi3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qni3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qr3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nc3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ni3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ns3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nr3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qc3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qi3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qni3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qr3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nc3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ni3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ns3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nr3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['t3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qv3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['t3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qv3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pres'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['rho'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dzq'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['w3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['wvar'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qg3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ng3dten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qg3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ng3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qgsten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qrsten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qisten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qnisten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qcsten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ngsten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nrsten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nisten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nssten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ncsten'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['cf3d'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['precrt'] = jnp.zeros((), dtype=jnp.float32)
    s['snowrt'] = jnp.zeros((), dtype=jnp.float32)
    s['effc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['effi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['effs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['effr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['effg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dt'] = jnp.zeros((), dtype=jnp.float32)
    s['lamc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['lami'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['lams'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['lamr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['lamg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['cdist1'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['n0i'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['n0s'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['n0rr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['n0g'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pgam'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsubc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsubi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsubs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsubr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['prd'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pre'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['prds'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nnuccc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['mnuccc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pra'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['prc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pcc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nnuccd'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['mnuccd'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['mnuccr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nnuccr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['npra'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nragg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsagg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nprc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nprc1'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['prai'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['prci'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['psacws'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['npsacws'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['psacwi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['npsacwi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nprci'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nprai'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nmults'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nmultr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qmults'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qmultr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pracs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['npracs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pccn'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['psmlt'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['evpms'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsmlts'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsmltr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['piacr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['niacr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['praci'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['piacrs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['niacrs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pracis'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['eprd'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['eprds'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pracg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['psacwg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pgsacw'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pgracs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['prdg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['eprdg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['evpmg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['pgmlt'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['npracg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['npsacwg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nscng'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ngracs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ngmltg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ngmltr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nsubg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['psacr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nmultg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nmultrg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qmultg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qmultrg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nact'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['sizefix_nr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['sizefix_nc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['sizefix_ni'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['sizefix_ns'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['sizefix_ng'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['negfix_ni'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['negfix_ns'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['negfix_nc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['negfix_nr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['negfix_ng'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nim_morr_cl'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qc_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qr_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qi_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qs_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qg_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nc_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nr_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ni_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ns_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ng_inst'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['kap'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['evs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['eis'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qvs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qvi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qvqvs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qvqvsi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dv'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['xxls'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['xxlv'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['cpm'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['mu'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['sc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['xlf'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ab'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['abi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dap'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nacnt'] = jnp.zeros((), dtype=jnp.float32)
    s['fmult'] = jnp.zeros((), dtype=jnp.float32)
    s['coffi'] = jnp.zeros((), dtype=jnp.float32)
    s['dumi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumfni'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumfng'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['uni'] = jnp.zeros((), dtype=jnp.float32)
    s['umi'] = jnp.zeros((), dtype=jnp.float32)
    s['umr'] = jnp.zeros((), dtype=jnp.float32)
    s['fr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['fi'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['fni'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['fg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['fng'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['rgvm'] = jnp.zeros((), dtype=jnp.float32)
    s['faloutr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['falouti'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faloutni'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faltndr'] = jnp.zeros((), dtype=jnp.float32)
    s['faltndi'] = jnp.zeros((), dtype=jnp.float32)
    s['faltndni'] = jnp.zeros((), dtype=jnp.float32)
    s['rho2'] = jnp.zeros((), dtype=jnp.float32)
    s['dumqs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumfns'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ums'] = jnp.zeros((), dtype=jnp.float32)
    s['uns'] = jnp.zeros((), dtype=jnp.float32)
    s['fs'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['fns'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['falouts'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faloutns'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faloutg'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faloutng'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faltnds'] = jnp.zeros((), dtype=jnp.float32)
    s['faltndns'] = jnp.zeros((), dtype=jnp.float32)
    s['unr'] = jnp.zeros((), dtype=jnp.float32)
    s['faltndg'] = jnp.zeros((), dtype=jnp.float32)
    s['faltndng'] = jnp.zeros((), dtype=jnp.float32)
    s['dumc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumfnc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['unc'] = jnp.zeros((), dtype=jnp.float32)
    s['umc'] = jnp.zeros((), dtype=jnp.float32)
    s['ung'] = jnp.zeros((), dtype=jnp.float32)
    s['umg'] = jnp.zeros((), dtype=jnp.float32)
    s['fc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faloutc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faloutnc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faltndc'] = jnp.zeros((), dtype=jnp.float32)
    s['faltndnc'] = jnp.zeros((), dtype=jnp.float32)
    s['fnc'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumfnr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faloutnr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['faltndnr'] = jnp.zeros((), dtype=jnp.float32)
    s['fnr'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ain'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['arn'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['asn'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['acn'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['agn'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dum'] = jnp.zeros((), dtype=jnp.float32)
    s['dum1'] = jnp.zeros((), dtype=jnp.float32)
    s['dum2'] = jnp.zeros((), dtype=jnp.float32)
    s['dumt'] = jnp.zeros((), dtype=jnp.float32)
    s['dumqv'] = jnp.zeros((), dtype=jnp.float32)
    s['dumqss'] = jnp.zeros((), dtype=jnp.float32)
    s['dumqsi'] = jnp.zeros((), dtype=jnp.float32)
    s['dums'] = jnp.zeros((), dtype=jnp.float32)
    s['tmpnum'] = jnp.zeros((), dtype=jnp.float32)
    s['dqsdt'] = jnp.zeros((), dtype=jnp.float32)
    s['dqsidt'] = jnp.zeros((), dtype=jnp.float32)
    s['epsi'] = jnp.zeros((), dtype=jnp.float32)
    s['epss'] = jnp.zeros((), dtype=jnp.float32)
    s['epsr'] = jnp.zeros((), dtype=jnp.float32)
    s['epsg'] = jnp.zeros((), dtype=jnp.float32)
    s['tauc'] = jnp.zeros((), dtype=jnp.float32)
    s['taur'] = jnp.zeros((), dtype=jnp.float32)
    s['taui'] = jnp.zeros((), dtype=jnp.float32)
    s['taus'] = jnp.zeros((), dtype=jnp.float32)
    s['taug'] = jnp.zeros((), dtype=jnp.float32)
    s['dumact'] = jnp.zeros((), dtype=jnp.float32)
    s['dum3'] = jnp.zeros((), dtype=jnp.float32)
    s['k'] = jnp.zeros((), dtype=jnp.int32)
    s['nstep'] = jnp.zeros((), dtype=jnp.int32)
    s['ltrue'] = jnp.zeros((), dtype=jnp.int32)
    s['ct'] = jnp.zeros((), dtype=jnp.float32)
    s['temp1'] = jnp.zeros((), dtype=jnp.float32)
    s['sat1'] = jnp.zeros((), dtype=jnp.float32)
    s['sigvl'] = jnp.zeros((), dtype=jnp.float32)
    s['kel'] = jnp.zeros((), dtype=jnp.float32)
    s['kc2'] = jnp.zeros((), dtype=jnp.float32)
    s['cry'] = jnp.zeros((), dtype=jnp.float32)
    s['kry'] = jnp.zeros((), dtype=jnp.float32)
    s['dumqi'] = jnp.zeros((), dtype=jnp.float32)
    s['dumni'] = jnp.zeros((), dtype=jnp.float32)
    s['dc0'] = jnp.zeros((), dtype=jnp.float32)
    s['ds0'] = jnp.zeros((), dtype=jnp.float32)
    s['dg0'] = jnp.zeros((), dtype=jnp.float32)
    s['dumqc'] = jnp.zeros((), dtype=jnp.float32)
    s['dumqr'] = jnp.zeros((), dtype=jnp.float32)
    s['ratio'] = jnp.zeros((), dtype=jnp.float32)
    s['sum_dep'] = jnp.zeros((), dtype=jnp.float32)
    s['fudgef'] = jnp.zeros((), dtype=jnp.float32)
    s['wef'] = jnp.zeros((), dtype=jnp.float32)
    s['anuc'] = jnp.zeros((), dtype=jnp.float32)
    s['bnuc'] = jnp.zeros((), dtype=jnp.float32)
    s['aact'] = jnp.zeros((), dtype=jnp.float32)
    s['gamm'] = jnp.zeros((), dtype=jnp.float32)
    s['gg'] = jnp.zeros((), dtype=jnp.float32)
    s['psi'] = jnp.zeros((), dtype=jnp.float32)
    s['eta1'] = jnp.zeros((), dtype=jnp.float32)
    s['eta2'] = jnp.zeros((), dtype=jnp.float32)
    s['sm1'] = jnp.zeros((), dtype=jnp.float32)
    s['sm2'] = jnp.zeros((), dtype=jnp.float32)
    s['smax'] = jnp.zeros((), dtype=jnp.float32)
    s['uu1'] = jnp.zeros((), dtype=jnp.float32)
    s['uu2'] = jnp.zeros((), dtype=jnp.float32)
    s['alpha'] = jnp.zeros((), dtype=jnp.float32)
    s['dlams'] = jnp.zeros((), dtype=jnp.float32)
    s['dlamr'] = jnp.zeros((), dtype=jnp.float32)
    s['dlami'] = jnp.zeros((), dtype=jnp.float32)
    s['dlamc'] = jnp.zeros((), dtype=jnp.float32)
    s['dlamg'] = jnp.zeros((), dtype=jnp.float32)
    s['lammax'] = jnp.zeros((), dtype=jnp.float32)
    s['lammin'] = jnp.zeros((), dtype=jnp.float32)
    s['idrop'] = jnp.zeros((), dtype=jnp.int32)
    s['nu'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['dumii'] = jnp.zeros((), dtype=jnp.int32)
    s['qc3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qi3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qni3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qr3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nc3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ni3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ns3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['nr3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['t3d_temp'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qv3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qg3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['ng3d_init'] = jnp.zeros((nz,), dtype=jnp.float32)
    s['qv_init'] = jnp.zeros((), dtype=jnp.float32)
    s['qsat_init'] = jnp.zeros((), dtype=jnp.float32)
    s['tmpqsmall'] = jnp.zeros((), dtype=jnp.float32)
    s['t3d_init'] = jnp.zeros((), dtype=jnp.float32)
    s['ims'] = I(1)
    s['ime'] = I(1)
    s['jms'] = I(1)
    s['jme'] = I(1)
    s['kms'] = I(1)
    s['its'] = I(1)
    s['ite'] = I(1)
    s['jts'] = I(1)
    s['jte'] = I(1)
    s['kts'] = I(1)
    s['kme'] = I(nz)
    s['kte'] = I(nz)
    s.update(inputs)
    old_1 = s['ltrue']
    s['ltrue'] = I(0)
    s['ltrue'] = jnp.where(s["_pc"] == 0, s['ltrue'], old_1)
    old_2 = s['effc']
    s['effc'] = jnp.broadcast_to(F(F(25.0)), s['effc'].shape)
    s['effc'] = jnp.where(s["_pc"] == 0, s['effc'], old_2)
    old_3 = s['effi']
    s['effi'] = jnp.broadcast_to(F(F(25.0)), s['effi'].shape)
    s['effi'] = jnp.where(s["_pc"] == 0, s['effi'], old_3)
    old_4 = s['effs']
    s['effs'] = jnp.broadcast_to(F(F(25.0)), s['effs'].shape)
    s['effs'] = jnp.where(s["_pc"] == 0, s['effs'], old_4)
    old_5 = s['effr']
    s['effr'] = jnp.broadcast_to(F(F(25.0)), s['effr'].shape)
    s['effr'] = jnp.where(s["_pc"] == 0, s['effr'], old_5)
    old_6 = s['effg']
    s['effg'] = jnp.broadcast_to(F(F(25.0)), s['effg'].shape)
    s['effg'] = jnp.where(s["_pc"] == 0, s['effg'], old_6)
    old_7 = s['prc']
    s['prc'] = jnp.broadcast_to(F(F(0.0)), s['prc'].shape)
    s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_7)
    old_8 = s['nprc']
    s['nprc'] = jnp.broadcast_to(F(F(0.0)), s['nprc'].shape)
    s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_8)
    old_9 = s['nprc1']
    s['nprc1'] = jnp.broadcast_to(F(F(0.0)), s['nprc1'].shape)
    s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_9)
    old_10 = s['pra']
    s['pra'] = jnp.broadcast_to(F(F(0.0)), s['pra'].shape)
    s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_10)
    old_11 = s['npra']
    s['npra'] = jnp.broadcast_to(F(F(0.0)), s['npra'].shape)
    s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_11)
    old_12 = s['nragg']
    s['nragg'] = jnp.broadcast_to(F(F(0.0)), s['nragg'].shape)
    s['nragg'] = jnp.where(s["_pc"] == 0, s['nragg'], old_12)
    old_13 = s['psmlt']
    s['psmlt'] = jnp.broadcast_to(F(F(0.0)), s['psmlt'].shape)
    s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_13)
    old_14 = s['nsmlts']
    s['nsmlts'] = jnp.broadcast_to(F(F(0.0)), s['nsmlts'].shape)
    s['nsmlts'] = jnp.where(s["_pc"] == 0, s['nsmlts'], old_14)
    old_15 = s['nsmltr']
    s['nsmltr'] = jnp.broadcast_to(F(F(0.0)), s['nsmltr'].shape)
    s['nsmltr'] = jnp.where(s["_pc"] == 0, s['nsmltr'], old_15)
    old_16 = s['pcc']
    s['pcc'] = jnp.broadcast_to(F(F(0.0)), s['pcc'].shape)
    s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_16)
    old_17 = s['pre']
    s['pre'] = jnp.broadcast_to(F(F(0.0)), s['pre'].shape)
    s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_17)
    old_18 = s['nsubc']
    s['nsubc'] = jnp.broadcast_to(F(F(0.0)), s['nsubc'].shape)
    s['nsubc'] = jnp.where(s["_pc"] == 0, s['nsubc'], old_18)
    old_19 = s['nsubr']
    s['nsubr'] = jnp.broadcast_to(F(F(0.0)), s['nsubr'].shape)
    s['nsubr'] = jnp.where(s["_pc"] == 0, s['nsubr'], old_19)
    old_20 = s['pracg']
    s['pracg'] = jnp.broadcast_to(F(F(0.0)), s['pracg'].shape)
    s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_20)
    old_21 = s['npracg']
    s['npracg'] = jnp.broadcast_to(F(F(0.0)), s['npracg'].shape)
    s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_21)
    old_22 = s['psmlt']
    s['psmlt'] = jnp.broadcast_to(F(F(0.0)), s['psmlt'].shape)
    s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_22)
    old_23 = s['evpms']
    s['evpms'] = jnp.broadcast_to(F(F(0.0)), s['evpms'].shape)
    s['evpms'] = jnp.where(s["_pc"] == 0, s['evpms'], old_23)
    old_24 = s['pgmlt']
    s['pgmlt'] = jnp.broadcast_to(F(F(0.0)), s['pgmlt'].shape)
    s['pgmlt'] = jnp.where(s["_pc"] == 0, s['pgmlt'], old_24)
    old_25 = s['evpmg']
    s['evpmg'] = jnp.broadcast_to(F(F(0.0)), s['evpmg'].shape)
    s['evpmg'] = jnp.where(s["_pc"] == 0, s['evpmg'], old_25)
    old_26 = s['pracs']
    s['pracs'] = jnp.broadcast_to(F(F(0.0)), s['pracs'].shape)
    s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_26)
    old_27 = s['npracs']
    s['npracs'] = jnp.broadcast_to(F(F(0.0)), s['npracs'].shape)
    s['npracs'] = jnp.where(s["_pc"] == 0, s['npracs'], old_27)
    old_28 = s['ngmltg']
    s['ngmltg'] = jnp.broadcast_to(F(F(0.0)), s['ngmltg'].shape)
    s['ngmltg'] = jnp.where(s["_pc"] == 0, s['ngmltg'], old_28)
    old_29 = s['ngmltr']
    s['ngmltr'] = jnp.broadcast_to(F(F(0.0)), s['ngmltr'].shape)
    s['ngmltr'] = jnp.where(s["_pc"] == 0, s['ngmltr'], old_29)
    old_30 = s['fr']
    s['fr'] = jnp.broadcast_to(F(F(0.0)), s['fr'].shape)
    s['fr'] = jnp.where(s["_pc"] == 0, s['fr'], old_30)
    old_31 = s['nprc1']
    s['nprc1'] = jnp.broadcast_to(F(F(0.0)), s['nprc1'].shape)
    s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_31)
    old_32 = s['nragg']
    s['nragg'] = jnp.broadcast_to(F(F(0.0)), s['nragg'].shape)
    s['nragg'] = jnp.where(s["_pc"] == 0, s['nragg'], old_32)
    old_33 = s['npracg']
    s['npracg'] = jnp.broadcast_to(F(F(0.0)), s['npracg'].shape)
    s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_33)
    old_34 = s['nsubr']
    s['nsubr'] = jnp.broadcast_to(F(F(0.0)), s['nsubr'].shape)
    s['nsubr'] = jnp.where(s["_pc"] == 0, s['nsubr'], old_34)
    old_35 = s['nsmltr']
    s['nsmltr'] = jnp.broadcast_to(F(F(0.0)), s['nsmltr'].shape)
    s['nsmltr'] = jnp.where(s["_pc"] == 0, s['nsmltr'], old_35)
    old_36 = s['ngmltr']
    s['ngmltr'] = jnp.broadcast_to(F(F(0.0)), s['ngmltr'].shape)
    s['ngmltr'] = jnp.where(s["_pc"] == 0, s['ngmltr'], old_36)
    old_37 = s['npracs']
    s['npracs'] = jnp.broadcast_to(F(F(0.0)), s['npracs'].shape)
    s['npracs'] = jnp.where(s["_pc"] == 0, s['npracs'], old_37)
    old_38 = s['nnuccr']
    s['nnuccr'] = jnp.broadcast_to(F(F(0.0)), s['nnuccr'].shape)
    s['nnuccr'] = jnp.where(s["_pc"] == 0, s['nnuccr'], old_38)
    old_39 = s['niacr']
    s['niacr'] = jnp.broadcast_to(F(F(0.0)), s['niacr'].shape)
    s['niacr'] = jnp.where(s["_pc"] == 0, s['niacr'], old_39)
    old_40 = s['niacrs']
    s['niacrs'] = jnp.broadcast_to(F(F(0.0)), s['niacrs'].shape)
    s['niacrs'] = jnp.where(s["_pc"] == 0, s['niacrs'], old_40)
    old_41 = s['ngracs']
    s['ngracs'] = jnp.broadcast_to(F(F(0.0)), s['ngracs'].shape)
    s['ngracs'] = jnp.where(s["_pc"] == 0, s['ngracs'], old_41)
    old_42 = s['nact']
    s['nact'] = jnp.broadcast_to(F(F(0.0)), s['nact'].shape)
    s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_42)
    old_43 = s['sizefix_nr']
    s['sizefix_nr'] = jnp.broadcast_to(F(F(0.0)), s['sizefix_nr'].shape)
    s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_43)
    old_44 = s['sizefix_nc']
    s['sizefix_nc'] = jnp.broadcast_to(F(F(0.0)), s['sizefix_nc'].shape)
    s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_44)
    old_45 = s['sizefix_ni']
    s['sizefix_ni'] = jnp.broadcast_to(F(F(0.0)), s['sizefix_ni'].shape)
    s['sizefix_ni'] = jnp.where(s["_pc"] == 0, s['sizefix_ni'], old_45)
    old_46 = s['sizefix_ns']
    s['sizefix_ns'] = jnp.broadcast_to(F(F(0.0)), s['sizefix_ns'].shape)
    s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_46)
    old_47 = s['sizefix_ng']
    s['sizefix_ng'] = jnp.broadcast_to(F(F(0.0)), s['sizefix_ng'].shape)
    s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_47)
    old_48 = s['negfix_ni']
    s['negfix_ni'] = jnp.broadcast_to(F(F(0.0)), s['negfix_ni'].shape)
    s['negfix_ni'] = jnp.where(s["_pc"] == 0, s['negfix_ni'], old_48)
    old_49 = s['negfix_ns']
    s['negfix_ns'] = jnp.broadcast_to(F(F(0.0)), s['negfix_ns'].shape)
    s['negfix_ns'] = jnp.where(s["_pc"] == 0, s['negfix_ns'], old_49)
    old_50 = s['negfix_nc']
    s['negfix_nc'] = jnp.broadcast_to(F(F(0.0)), s['negfix_nc'].shape)
    s['negfix_nc'] = jnp.where(s["_pc"] == 0, s['negfix_nc'], old_50)
    old_51 = s['negfix_nr']
    s['negfix_nr'] = jnp.broadcast_to(F(F(0.0)), s['negfix_nr'].shape)
    s['negfix_nr'] = jnp.where(s["_pc"] == 0, s['negfix_nr'], old_51)
    old_52 = s['negfix_ng']
    s['negfix_ng'] = jnp.broadcast_to(F(F(0.0)), s['negfix_ng'].shape)
    s['negfix_ng'] = jnp.where(s["_pc"] == 0, s['negfix_ng'], old_52)
    old_53 = s['nim_morr_cl']
    s['nim_morr_cl'] = jnp.broadcast_to(F(F(0.0)), s['nim_morr_cl'].shape)
    s['nim_morr_cl'] = jnp.where(s["_pc"] == 0, s['nim_morr_cl'], old_53)
    old_54 = s['qc_inst']
    s['qc_inst'] = jnp.broadcast_to(F(F(0.0)), s['qc_inst'].shape)
    s['qc_inst'] = jnp.where(s["_pc"] == 0, s['qc_inst'], old_54)
    old_55 = s['qr_inst']
    s['qr_inst'] = jnp.broadcast_to(F(F(0.0)), s['qr_inst'].shape)
    s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_55)
    old_56 = s['qi_inst']
    s['qi_inst'] = jnp.broadcast_to(F(F(0.0)), s['qi_inst'].shape)
    s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_56)
    old_57 = s['qs_inst']
    s['qs_inst'] = jnp.broadcast_to(F(F(0.0)), s['qs_inst'].shape)
    s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_57)
    old_58 = s['qg_inst']
    s['qg_inst'] = jnp.broadcast_to(F(F(0.0)), s['qg_inst'].shape)
    s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_58)
    old_59 = s['nc_inst']
    s['nc_inst'] = jnp.broadcast_to(F(F(0.0)), s['nc_inst'].shape)
    s['nc_inst'] = jnp.where(s["_pc"] == 0, s['nc_inst'], old_59)
    old_60 = s['nr_inst']
    s['nr_inst'] = jnp.broadcast_to(F(F(0.0)), s['nr_inst'].shape)
    s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_60)
    old_61 = s['ni_inst']
    s['ni_inst'] = jnp.broadcast_to(F(F(0.0)), s['ni_inst'].shape)
    s['ni_inst'] = jnp.where(s["_pc"] == 0, s['ni_inst'], old_61)
    old_62 = s['ns_inst']
    s['ns_inst'] = jnp.broadcast_to(F(F(0.0)), s['ns_inst'].shape)
    s['ns_inst'] = jnp.where(s["_pc"] == 0, s['ns_inst'], old_62)
    old_63 = s['ng_inst']
    s['ng_inst'] = jnp.broadcast_to(F(F(0.0)), s['ng_inst'].shape)
    s['ng_inst'] = jnp.where(s["_pc"] == 0, s['ng_inst'], old_63)
    old_64 = s['qc3d_init']
    s['qc3d_init'] = jnp.broadcast_to(F(s['qc3d']), s['qc3d_init'].shape)
    s['qc3d_init'] = jnp.where(s["_pc"] == 0, s['qc3d_init'], old_64)
    old_65 = s['qi3d_init']
    s['qi3d_init'] = jnp.broadcast_to(F(s['qi3d']), s['qi3d_init'].shape)
    s['qi3d_init'] = jnp.where(s["_pc"] == 0, s['qi3d_init'], old_65)
    old_66 = s['qni3d_init']
    s['qni3d_init'] = jnp.broadcast_to(F(s['qni3d']), s['qni3d_init'].shape)
    s['qni3d_init'] = jnp.where(s["_pc"] == 0, s['qni3d_init'], old_66)
    old_67 = s['qr3d_init']
    s['qr3d_init'] = jnp.broadcast_to(F(s['qr3d']), s['qr3d_init'].shape)
    s['qr3d_init'] = jnp.where(s["_pc"] == 0, s['qr3d_init'], old_67)
    old_68 = s['nc3d_init']
    s['nc3d_init'] = jnp.broadcast_to(F(s['nc3d']), s['nc3d_init'].shape)
    s['nc3d_init'] = jnp.where(s["_pc"] == 0, s['nc3d_init'], old_68)
    old_69 = s['ni3d_init']
    s['ni3d_init'] = jnp.broadcast_to(F(s['ni3d']), s['ni3d_init'].shape)
    s['ni3d_init'] = jnp.where(s["_pc"] == 0, s['ni3d_init'], old_69)
    old_70 = s['ns3d_init']
    s['ns3d_init'] = jnp.broadcast_to(F(s['ns3d']), s['ns3d_init'].shape)
    s['ns3d_init'] = jnp.where(s["_pc"] == 0, s['ns3d_init'], old_70)
    old_71 = s['nr3d_init']
    s['nr3d_init'] = jnp.broadcast_to(F(s['nr3d']), s['nr3d_init'].shape)
    s['nr3d_init'] = jnp.where(s["_pc"] == 0, s['nr3d_init'], old_71)
    old_72 = s['qg3d_init']
    s['qg3d_init'] = jnp.broadcast_to(F(s['qg3d']), s['qg3d_init'].shape)
    s['qg3d_init'] = jnp.where(s["_pc"] == 0, s['qg3d_init'], old_72)
    old_73 = s['ng3d_init']
    s['ng3d_init'] = jnp.broadcast_to(F(s['ng3d']), s['ng3d_init'].shape)
    s['ng3d_init'] = jnp.where(s["_pc"] == 0, s['ng3d_init'], old_73)
    old_74 = s['t3d_temp']
    s['t3d_temp'] = jnp.broadcast_to(F(s['t3d']), s['t3d_temp'].shape)
    s['t3d_temp'] = jnp.where(s["_pc"] == 0, s['t3d_temp'], old_74)
    old_75 = s['qv3d_init']
    s['qv3d_init'] = jnp.broadcast_to(F(s['qv3d']), s['qv3d_init'].shape)
    s['qv3d_init'] = jnp.where(s["_pc"] == 0, s['qv3d_init'], old_75)
    def loop_76(iteration, s):
        s = dict(s)
        s['k'] = I(s['kts'] + iteration * (1))
        old_77 = s['xxlv']
        s['xxlv'] = jnp.broadcast_to(F(s['lv']), s['xxlv'].shape)
        s['xxlv'] = jnp.where(s["_pc"] == 0, s['xxlv'], old_77)
        old_78 = s['xxls']
        s['xxls'] = s['xxls'].at[s['k'] - 1].set(F(s['ls']))
        s['xxls'] = jnp.where(s["_pc"] == 0, s['xxls'], old_78)
        old_79 = s['cpm']
        s['cpm'] = s['cpm'].at[s['k'] - 1].set(F(s['cp']))
        s['cpm'] = jnp.where(s["_pc"] == 0, s['cpm'], old_79)
        old_80 = s['evs']
        s['evs'] = s['evs'].at[s['k'] - 1].set(F(jnp.minimum(_arith(F(0.99), s['pres'][s['k'] - 1], 'mul'), POLYSVP(s['t3d'][s['k'] - 1], 0))))
        s['evs'] = jnp.where(s["_pc"] == 0, s['evs'], old_80)
        old_81 = s['eis']
        s['eis'] = s['eis'].at[s['k'] - 1].set(F(jnp.minimum(_arith(F(0.99), s['pres'][s['k'] - 1], 'mul'), POLYSVP(s['t3d'][s['k'] - 1], 1))))
        s['eis'] = jnp.where(s["_pc"] == 0, s['eis'], old_81)
        def yes_82(s):
            s = dict(s)
            old_83 = s['eis']
            s['eis'] = s['eis'].at[s['k'] - 1].set(F(s['evs'][s['k'] - 1]))
            s['eis'] = jnp.where(s["_pc"] == 0, s['eis'], old_83)
            return s
        def no_82(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['eis'][s['k'] - 1] > s['evs'][s['k'] - 1])), yes_82, no_82, s)
        old_84 = s['qvs']
        s['qvs'] = s['qvs'].at[s['k'] - 1].set(F(_div(_arith(s['ep_2'], s['evs'][s['k'] - 1], 'mul'), _arith(s['pres'][s['k'] - 1], s['evs'][s['k'] - 1], 'sub'))))
        s['qvs'] = jnp.where(s["_pc"] == 0, s['qvs'], old_84)
        old_85 = s['qvi']
        s['qvi'] = s['qvi'].at[s['k'] - 1].set(F(_div(_arith(s['ep_2'], s['eis'][s['k'] - 1], 'mul'), _arith(s['pres'][s['k'] - 1], s['eis'][s['k'] - 1], 'sub'))))
        s['qvi'] = jnp.where(s["_pc"] == 0, s['qvi'], old_85)
        def yes_86(s):
            s = dict(s)
            old_87 = s['t3d_init']
            s['t3d_init'] = F(s['t3d'][s['k'] - 1])
            s['t3d_init'] = jnp.where(s["_pc"] == 0, s['t3d_init'], old_87)
            old_88 = s['qv_init']
            s['qv_init'] = F(s['qv3d'][s['k'] - 1])
            s['qv_init'] = jnp.where(s["_pc"] == 0, s['qv_init'], old_88)
            old_89 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(s['qvs'][s['k'] - 1]))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_89)
            old_90 = s['qsat_init']
            s['qsat_init'] = F(s['qvs'][s['k'] - 1])
            s['qsat_init'] = jnp.where(s["_pc"] == 0, s['qsat_init'], old_90)
            old_91 = s['qc3d']
            s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(_div(s['qc3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
            s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_91)
            def yes_92(s):
                s = dict(s)
                old_93 = s['nc3d']
                s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(s['nc3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
                s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_93)
                return s
            def no_92(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['inum'] == 0)), yes_92, no_92, s)
            old_94 = s['qr3d']
            s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(_div(s['qr3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
            s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_94)
            old_95 = s['nr3d']
            s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_div(s['nr3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
            s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_95)
            def yes_96(s):
                s = dict(s)
                old_97 = s['qi3d']
                s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(_div(s['qi3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
                s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_97)
                old_98 = s['ni3d']
                s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_div(s['ni3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
                s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_98)
                old_99 = s['qni3d']
                s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(_div(s['qni3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
                s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_99)
                old_100 = s['ns3d']
                s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_div(s['ns3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
                s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_100)
                return s
            def no_96(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['iliq'] == 0)), yes_96, no_96, s)
            def yes_101(s):
                s = dict(s)
                old_102 = s['qg3d']
                s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(_div(s['qg3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
                s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_102)
                old_103 = s['ng3d']
                s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_div(s['ng3d'][s['k'] - 1], s['cf3d'][s['k'] - 1])))
                s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_103)
                return s
            def no_101(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 0)), yes_101, no_101, s)
            return s
        def no_86(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['cf3d'][s['k'] - 1] > s['cloud_frac_thresh'])), yes_86, no_86, s)
        old_104 = s['qvqvs']
        s['qvqvs'] = s['qvqvs'].at[s['k'] - 1].set(F(_div(s['qv3d'][s['k'] - 1], s['qvs'][s['k'] - 1])))
        s['qvqvs'] = jnp.where(s["_pc"] == 0, s['qvqvs'], old_104)
        old_105 = s['qvqvsi']
        s['qvqvsi'] = s['qvqvsi'].at[s['k'] - 1].set(F(_div(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1])))
        s['qvqvsi'] = jnp.where(s["_pc"] == 0, s['qvqvsi'], old_105)
        def yes_106(s):
            s = dict(s)
            def yes_107(s):
                s = dict(s)
                old_108 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_108)
                old_109 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qr3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_109)
                old_110 = s['qr_inst']
                s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_110)
                old_111 = s['qr3d']
                s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_111)
                return s
            def no_107(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] < F(1e-08))), yes_107, no_107, s)
            def yes_112(s):
                s = dict(s)
                old_113 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qc3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_113)
                old_114 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qc3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_114)
                old_115 = s['qc_inst']
                s['qc_inst'] = s['qc_inst'].at[s['k'] - 1].set(F(_arith(s['qc_inst'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qc_inst'] = jnp.where(s["_pc"] == 0, s['qc_inst'], old_115)
                old_116 = s['qc3d']
                s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_116)
                return s
            def no_112(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] < F(1e-08))), yes_112, no_112, s)
            return s
        def no_106(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qvqvs'][s['k'] - 1] < F(0.9))), yes_106, no_106, s)
        def yes_117(s):
            s = dict(s)
            def yes_118(s):
                s = dict(s)
                old_119 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qi3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_119)
                old_120 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qi3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_120)
                old_121 = s['qi_inst']
                s['qi_inst'] = s['qi_inst'].at[s['k'] - 1].set(F(_arith(s['qi_inst'][s['k'] - 1], _div(s['qi3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_121)
                old_122 = s['qi3d']
                s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_122)
                return s
            def no_118(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] < F(1e-08))), yes_118, no_118, s)
            def yes_123(s):
                s = dict(s)
                old_124 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qni3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_124)
                old_125 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qni3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_125)
                old_126 = s['qs_inst']
                s['qs_inst'] = s['qs_inst'].at[s['k'] - 1].set(F(_arith(s['qs_inst'][s['k'] - 1], _div(s['qni3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_126)
                old_127 = s['qni3d']
                s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_127)
                return s
            def no_123(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] < F(1e-08))), yes_123, no_123, s)
            def yes_128(s):
                s = dict(s)
                old_129 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qg3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_129)
                old_130 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qg3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_130)
                old_131 = s['qg_inst']
                s['qg_inst'] = s['qg_inst'].at[s['k'] - 1].set(F(_arith(s['qg_inst'][s['k'] - 1], _div(s['qg3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_131)
                old_132 = s['qg3d']
                s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_132)
                return s
            def no_128(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] < F(1e-08))), yes_128, no_128, s)
            return s
        def no_117(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qvqvsi'][s['k'] - 1] < F(0.9))), yes_117, no_117, s)
        old_133 = s['xlf']
        s['xlf'] = s['xlf'].at[s['k'] - 1].set(F(_arith(s['xxls'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'sub')))
        s['xlf'] = jnp.where(s["_pc"] == 0, s['xlf'], old_133)
        def yes_134(s):
            s = dict(s)
            old_135 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qc3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_135)
            old_136 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qc3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_136)
            old_137 = s['nc_inst']
            s['nc_inst'] = s['nc_inst'].at[s['k'] - 1].set(F(_arith(s['nc_inst'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']), 'sub')))
            s['nc_inst'] = jnp.where(s["_pc"] == 0, s['nc_inst'], old_137)
            old_138 = s['qc_inst']
            s['qc_inst'] = s['qc_inst'].at[s['k'] - 1].set(F(_arith(s['qc_inst'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qc_inst'] = jnp.where(s["_pc"] == 0, s['qc_inst'], old_138)
            old_139 = s['qc3d']
            s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_139)
            old_140 = s['nc3d']
            s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_140)
            old_141 = s['effc']
            s['effc'] = s['effc'].at[s['k'] - 1].set(F(F(0.0)))
            s['effc'] = jnp.where(s["_pc"] == 0, s['effc'], old_141)
            return s
        def no_134(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] < s['qsmall'])), yes_134, no_134, s)
        def yes_142(s):
            s = dict(s)
            old_143 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_143)
            old_144 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qr3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_144)
            old_145 = s['nr_inst']
            s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'sub')))
            s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_145)
            old_146 = s['qr_inst']
            s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_146)
            old_147 = s['qr3d']
            s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_147)
            old_148 = s['nr3d']
            s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_148)
            old_149 = s['effr']
            s['effr'] = s['effr'].at[s['k'] - 1].set(F(F(0.0)))
            s['effr'] = jnp.where(s["_pc"] == 0, s['effr'], old_149)
            return s
        def no_142(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] < s['qsmall'])), yes_142, no_142, s)
        def yes_150(s):
            s = dict(s)
            old_151 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qi3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_151)
            old_152 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qi3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_152)
            old_153 = s['ni_inst']
            s['ni_inst'] = s['ni_inst'].at[s['k'] - 1].set(F(_arith(s['ni_inst'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']), 'sub')))
            s['ni_inst'] = jnp.where(s["_pc"] == 0, s['ni_inst'], old_153)
            old_154 = s['qi_inst']
            s['qi_inst'] = s['qi_inst'].at[s['k'] - 1].set(F(_arith(s['qi_inst'][s['k'] - 1], _div(s['qi3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_154)
            old_155 = s['qi3d']
            s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_155)
            old_156 = s['ni3d']
            s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_156)
            old_157 = s['effi']
            s['effi'] = s['effi'].at[s['k'] - 1].set(F(F(0.0)))
            s['effi'] = jnp.where(s["_pc"] == 0, s['effi'], old_157)
            return s
        def no_150(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] < s['qsmall'])), yes_150, no_150, s)
        def yes_158(s):
            s = dict(s)
            old_159 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qni3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_159)
            old_160 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qni3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_160)
            old_161 = s['ns_inst']
            s['ns_inst'] = s['ns_inst'].at[s['k'] - 1].set(F(_arith(s['ns_inst'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']), 'sub')))
            s['ns_inst'] = jnp.where(s["_pc"] == 0, s['ns_inst'], old_161)
            old_162 = s['qs_inst']
            s['qs_inst'] = s['qs_inst'].at[s['k'] - 1].set(F(_arith(s['qs_inst'][s['k'] - 1], _div(s['qni3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_162)
            old_163 = s['qni3d']
            s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_163)
            old_164 = s['ns3d']
            s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_164)
            old_165 = s['effs']
            s['effs'] = s['effs'].at[s['k'] - 1].set(F(F(0.0)))
            s['effs'] = jnp.where(s["_pc"] == 0, s['effs'], old_165)
            return s
        def no_158(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] < s['qsmall'])), yes_158, no_158, s)
        def yes_166(s):
            s = dict(s)
            old_167 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qg3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_167)
            old_168 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qg3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_168)
            old_169 = s['ng_inst']
            s['ng_inst'] = s['ng_inst'].at[s['k'] - 1].set(F(_arith(s['ng_inst'][s['k'] - 1], _div(s['ng3d'][s['k'] - 1], s['dt']), 'sub')))
            s['ng_inst'] = jnp.where(s["_pc"] == 0, s['ng_inst'], old_169)
            old_170 = s['qg_inst']
            s['qg_inst'] = s['qg_inst'].at[s['k'] - 1].set(F(_arith(s['qg_inst'][s['k'] - 1], _div(s['qg3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_170)
            old_171 = s['qg3d']
            s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_171)
            old_172 = s['ng3d']
            s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_172)
            old_173 = s['effg']
            s['effg'] = s['effg'].at[s['k'] - 1].set(F(F(0.0)))
            s['effg'] = jnp.where(s["_pc"] == 0, s['effg'], old_173)
            return s
        def no_166(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] < s['qsmall'])), yes_166, no_166, s)
        old_174 = s['qrsten']
        s['qrsten'] = s['qrsten'].at[s['k'] - 1].set(F(F(0.0)))
        s['qrsten'] = jnp.where(s["_pc"] == 0, s['qrsten'], old_174)
        old_175 = s['qisten']
        s['qisten'] = s['qisten'].at[s['k'] - 1].set(F(F(0.0)))
        s['qisten'] = jnp.where(s["_pc"] == 0, s['qisten'], old_175)
        old_176 = s['qnisten']
        s['qnisten'] = s['qnisten'].at[s['k'] - 1].set(F(F(0.0)))
        s['qnisten'] = jnp.where(s["_pc"] == 0, s['qnisten'], old_176)
        old_177 = s['qcsten']
        s['qcsten'] = s['qcsten'].at[s['k'] - 1].set(F(F(0.0)))
        s['qcsten'] = jnp.where(s["_pc"] == 0, s['qcsten'], old_177)
        old_178 = s['qgsten']
        s['qgsten'] = s['qgsten'].at[s['k'] - 1].set(F(F(0.0)))
        s['qgsten'] = jnp.where(s["_pc"] == 0, s['qgsten'], old_178)
        old_179 = s['nrsten']
        s['nrsten'] = s['nrsten'].at[s['k'] - 1].set(F(F(0.0)))
        s['nrsten'] = jnp.where(s["_pc"] == 0, s['nrsten'], old_179)
        old_180 = s['nisten']
        s['nisten'] = s['nisten'].at[s['k'] - 1].set(F(F(0.0)))
        s['nisten'] = jnp.where(s["_pc"] == 0, s['nisten'], old_180)
        old_181 = s['nssten']
        s['nssten'] = s['nssten'].at[s['k'] - 1].set(F(F(0.0)))
        s['nssten'] = jnp.where(s["_pc"] == 0, s['nssten'], old_181)
        old_182 = s['ncsten']
        s['ncsten'] = s['ncsten'].at[s['k'] - 1].set(F(F(0.0)))
        s['ncsten'] = jnp.where(s["_pc"] == 0, s['ncsten'], old_182)
        old_183 = s['ngsten']
        s['ngsten'] = s['ngsten'].at[s['k'] - 1].set(F(F(0.0)))
        s['ngsten'] = jnp.where(s["_pc"] == 0, s['ngsten'], old_183)
        old_184 = s['mu']
        s['mu'] = s['mu'].at[s['k'] - 1].set(F(_div(_arith(F(1.496e-06), _arith(s['t3d'][s['k'] - 1], F(1.5), 'pow'), 'mul'), _arith(s['t3d'][s['k'] - 1], F(120.0), 'add'))))
        s['mu'] = jnp.where(s["_pc"] == 0, s['mu'], old_184)
        old_185 = s['dum']
        s['dum'] = F(_arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.54), 'pow'))
        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_185)
        old_186 = s['ain']
        s['ain'] = s['ain'].at[s['k'] - 1].set(F(_arith(_arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.35), 'pow'), s['ai'], 'mul')))
        s['ain'] = jnp.where(s["_pc"] == 0, s['ain'], old_186)
        old_187 = s['arn']
        s['arn'] = s['arn'].at[s['k'] - 1].set(F(_arith(s['dum'], s['ar'], 'mul')))
        s['arn'] = jnp.where(s["_pc"] == 0, s['arn'], old_187)
        old_188 = s['asn']
        s['asn'] = s['asn'].at[s['k'] - 1].set(F(_arith(s['dum'], s['as'], 'mul')))
        s['asn'] = jnp.where(s["_pc"] == 0, s['asn'], old_188)
        old_189 = s['acn']
        s['acn'] = s['acn'].at[s['k'] - 1].set(F(_div(_arith(s['g'], s['rhow'], 'mul'), _arith(F(18.0), s['mu'][s['k'] - 1], 'mul'))))
        s['acn'] = jnp.where(s["_pc"] == 0, s['acn'], old_189)
        old_190 = s['agn']
        s['agn'] = s['agn'].at[s['k'] - 1].set(F(_arith(s['dum'], s['ag'], 'mul')))
        s['agn'] = jnp.where(s["_pc"] == 0, s['agn'], old_190)
        old_191 = s['lami']
        s['lami'] = s['lami'].at[s['k'] - 1].set(F(F(0.0)))
        s['lami'] = jnp.where(s["_pc"] == 0, s['lami'], old_191)
        def yes_192(s):
            s = dict(s)
            def yes_193(s):
                s = dict(s)
                s["_pc"] = jnp.where(s["_pc"] == 0, I(200), s["_pc"])
                return s
            def no_193(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] < s['tmelt'])) & ((s['qvqvsi'][s['k'] - 1] < F(0.999))))), yes_193, no_193, s)
            def yes_194(s):
                s = dict(s)
                s["_pc"] = jnp.where(s["_pc"] == 0, I(200), s["_pc"])
                return s
            def no_194(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] >= s['tmelt'])) & ((s['qvqvs'][s['k'] - 1] < F(0.999))))), yes_194, no_194, s)
            return s
        def no_192(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] < s['qsmall'])) & ((s['qi3d'][s['k'] - 1] < s['qsmall'])) & ((s['qni3d'][s['k'] - 1] < s['qsmall'])) & ((s['qr3d'][s['k'] - 1] < s['qsmall'])) & ((s['qg3d'][s['k'] - 1] < s['qsmall'])))), yes_192, no_192, s)
        old_195 = s['kap']
        s['kap'] = s['kap'].at[s['k'] - 1].set(F(_arith(F(1414.0), s['mu'][s['k'] - 1], 'mul')))
        s['kap'] = jnp.where(s["_pc"] == 0, s['kap'], old_195)
        old_196 = s['dv']
        s['dv'] = s['dv'].at[s['k'] - 1].set(F(_div(_arith(F(8.794e-05), _arith(s['t3d'][s['k'] - 1], F(1.81), 'pow'), 'mul'), s['pres'][s['k'] - 1])))
        s['dv'] = jnp.where(s["_pc"] == 0, s['dv'], old_196)
        old_197 = s['sc']
        s['sc'] = s['sc'].at[s['k'] - 1].set(F(_div(s['mu'][s['k'] - 1], _arith(s['rho'][s['k'] - 1], s['dv'][s['k'] - 1], 'mul'))))
        s['sc'] = jnp.where(s["_pc"] == 0, s['sc'], old_197)
        old_198 = s['dum']
        s['dum'] = F(_arith(s['rv'], _arith(s['t3d'][s['k'] - 1], 2, 'pow'), 'mul'))
        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_198)
        old_199 = s['dqsdt']
        s['dqsdt'] = F(_div(_arith(s['xxlv'][s['k'] - 1], s['qvs'][s['k'] - 1], 'mul'), s['dum']))
        s['dqsdt'] = jnp.where(s["_pc"] == 0, s['dqsdt'], old_199)
        old_200 = s['dqsidt']
        s['dqsidt'] = F(_div(_arith(s['xxls'][s['k'] - 1], s['qvi'][s['k'] - 1], 'mul'), s['dum']))
        s['dqsidt'] = jnp.where(s["_pc"] == 0, s['dqsidt'], old_200)
        old_201 = s['abi']
        s['abi'] = s['abi'].at[s['k'] - 1].set(F(_arith(F(1.0), _div(_arith(s['dqsidt'], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'add')))
        s['abi'] = jnp.where(s["_pc"] == 0, s['abi'], old_201)
        old_202 = s['ab']
        s['ab'] = s['ab'].at[s['k'] - 1].set(F(_arith(F(1.0), _div(_arith(s['dqsdt'], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'add')))
        s['ab'] = jnp.where(s["_pc"] == 0, s['ab'], old_202)
        def yes_203(s):
            s = dict(s)
            def yes_204(s):
                s = dict(s)
                old_205 = s['nc3d']
                s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_arith(s['ndcnst'], F(1000000.0), 'mul'), s['rho'][s['k'] - 1])))
                s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_205)
                return s
            def no_204(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['inum'] == 1)), yes_204, no_204, s)
            def yes_206(s):
                s = dict(s)
                old_207 = s['qr_inst']
                s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qni3d'][s['k'] - 1], s['dt']), 'add')))
                s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_207)
                old_208 = s['qs_inst']
                s['qs_inst'] = s['qs_inst'].at[s['k'] - 1].set(F(_arith(s['qs_inst'][s['k'] - 1], _div(s['qni3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_208)
                old_209 = s['nr_inst']
                s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']), 'add')))
                s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_209)
                old_210 = s['ns_inst']
                s['ns_inst'] = s['ns_inst'].at[s['k'] - 1].set(F(_arith(s['ns_inst'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']), 'sub')))
                s['ns_inst'] = jnp.where(s["_pc"] == 0, s['ns_inst'], old_210)
                old_211 = s['qr3d']
                s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(_arith(s['qr3d'][s['k'] - 1], s['qni3d'][s['k'] - 1], 'add')))
                s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_211)
                old_212 = s['nr3d']
                s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], s['ns3d'][s['k'] - 1], 'add')))
                s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_212)
                old_213 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qni3d'][s['k'] - 1], s['xlf'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_213)
                old_214 = s['qni3d']
                s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_214)
                old_215 = s['ns3d']
                s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_215)
                return s
            def no_206(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] < F(1e-06))), yes_206, no_206, s)
            def yes_216(s):
                s = dict(s)
                old_217 = s['qr_inst']
                s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qg3d'][s['k'] - 1], s['dt']), 'add')))
                s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_217)
                old_218 = s['qg_inst']
                s['qg_inst'] = s['qg_inst'].at[s['k'] - 1].set(F(_arith(s['qg_inst'][s['k'] - 1], _div(s['qg3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_218)
                old_219 = s['nr_inst']
                s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], _div(s['ng3d'][s['k'] - 1], s['dt']), 'add')))
                s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_219)
                old_220 = s['ng_inst']
                s['ng_inst'] = s['ng_inst'].at[s['k'] - 1].set(F(_arith(s['ng_inst'][s['k'] - 1], _div(s['ng3d'][s['k'] - 1], s['dt']), 'sub')))
                s['ng_inst'] = jnp.where(s["_pc"] == 0, s['ng_inst'], old_220)
                old_221 = s['qr3d']
                s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(_arith(s['qr3d'][s['k'] - 1], s['qg3d'][s['k'] - 1], 'add')))
                s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_221)
                old_222 = s['nr3d']
                s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], s['ng3d'][s['k'] - 1], 'add')))
                s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_222)
                old_223 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qg3d'][s['k'] - 1], s['xlf'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_223)
                old_224 = s['qg3d']
                s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_224)
                old_225 = s['ng3d']
                s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_225)
                return s
            def no_216(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] < F(1e-06))), yes_216, no_216, s)
            def yes_226(s):
                s = dict(s)
                s["_pc"] = jnp.where(s["_pc"] == 0, I(300), s["_pc"])
                return s
            def no_226(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] < s['qsmall'])) & ((s['qni3d'][s['k'] - 1] < F(1e-08))) & ((s['qr3d'][s['k'] - 1] < s['qsmall'])) & ((s['qg3d'][s['k'] - 1] < F(1e-08))))), yes_226, no_226, s)
            def yes_227(s):
                s = dict(s)
                old_228 = s['negfix_ni']
                s['negfix_ni'] = s['negfix_ni'].at[s['k'] - 1].set(F(_arith(s['negfix_ni'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_ni'] = jnp.where(s["_pc"] == 0, s['negfix_ni'], old_228)
                return s
            def no_227(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['ni3d'][s['k'] - 1] < F(0.0))), yes_227, no_227, s)
            def yes_229(s):
                s = dict(s)
                old_230 = s['negfix_ns']
                s['negfix_ns'] = s['negfix_ns'].at[s['k'] - 1].set(F(_arith(s['negfix_ns'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_ns'] = jnp.where(s["_pc"] == 0, s['negfix_ns'], old_230)
                return s
            def no_229(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['ns3d'][s['k'] - 1] < F(0.0))), yes_229, no_229, s)
            def yes_231(s):
                s = dict(s)
                old_232 = s['negfix_nc']
                s['negfix_nc'] = s['negfix_nc'].at[s['k'] - 1].set(F(_arith(s['negfix_nc'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_nc'] = jnp.where(s["_pc"] == 0, s['negfix_nc'], old_232)
                return s
            def no_231(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['nc3d'][s['k'] - 1] < F(0.0))), yes_231, no_231, s)
            def yes_233(s):
                s = dict(s)
                old_234 = s['negfix_nr']
                s['negfix_nr'] = s['negfix_nr'].at[s['k'] - 1].set(F(_arith(s['negfix_nr'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_nr'] = jnp.where(s["_pc"] == 0, s['negfix_nr'], old_234)
                return s
            def no_233(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['nr3d'][s['k'] - 1] < F(0.0))), yes_233, no_233, s)
            def yes_235(s):
                s = dict(s)
                old_236 = s['negfix_ng']
                s['negfix_ng'] = s['negfix_ng'].at[s['k'] - 1].set(F(_arith(s['negfix_ng'][s['k'] - 1], _div(s['ng3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_ng'] = jnp.where(s["_pc"] == 0, s['negfix_ng'], old_236)
                return s
            def no_235(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['ng3d'][s['k'] - 1] < F(0.0))), yes_235, no_235, s)
            old_237 = s['ns3d']
            s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ns3d'][s['k'] - 1])))
            s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_237)
            old_238 = s['nc3d']
            s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['nc3d'][s['k'] - 1])))
            s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_238)
            old_239 = s['nr3d']
            s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['nr3d'][s['k'] - 1])))
            s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_239)
            old_240 = s['ng3d']
            s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ng3d'][s['k'] - 1])))
            s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_240)
            def yes_241(s):
                s = dict(s)
                old_242 = s['lamr']
                s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(s['pi'], s['rhow'], 'mul'), s['nr3d'][s['k'] - 1], 'mul'), s['qr3d'][s['k'] - 1]), _div(F(1.0), F(3.0)), 'pow')))
                s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_242)
                old_243 = s['n0rr']
                s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], s['lamr'][s['k'] - 1], 'mul')))
                s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_243)
                old_244 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_244)
                def yes_245(s):
                    s = dict(s)
                    old_246 = s['tmpnum']
                    s['tmpnum'] = F(s['nr3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_246)
                    old_247 = s['lamr']
                    s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(s['lamminr']))
                    s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_247)
                    old_248 = s['n0rr']
                    s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamr'][s['k'] - 1], 4, 'pow'), s['qr3d'][s['k'] - 1], 'mul'), _arith(s['pi'], s['rhow'], 'mul'))))
                    s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_248)
                    old_249 = s['nr3d']
                    s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_div(s['n0rr'][s['k'] - 1], s['lamr'][s['k'] - 1])))
                    s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_249)
                    old_250 = s['sizefix_nr']
                    s['sizefix_nr'] = s['sizefix_nr'].at[s['k'] - 1].set(F(_arith(s['sizefix_nr'][s['k'] - 1], _div(_arith(s['nr3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_250)
                    return s
                def no_245(s):
                    s = dict(s)
                    def yes_251(s):
                        s = dict(s)
                        old_252 = s['tmpnum']
                        s['tmpnum'] = F(s['nr3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_252)
                        old_253 = s['lamr']
                        s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(s['lammaxr']))
                        s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_253)
                        old_254 = s['n0rr']
                        s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamr'][s['k'] - 1], 4, 'pow'), s['qr3d'][s['k'] - 1], 'mul'), _arith(s['pi'], s['rhow'], 'mul'))))
                        s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_254)
                        old_255 = s['nr3d']
                        s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_div(s['n0rr'][s['k'] - 1], s['lamr'][s['k'] - 1])))
                        s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_255)
                        old_256 = s['sizefix_nr']
                        s['sizefix_nr'] = s['sizefix_nr'].at[s['k'] - 1].set(F(_arith(s['sizefix_nr'][s['k'] - 1], _div(_arith(s['nr3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_256)
                        return s
                    def no_251(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lamr'][s['k'] - 1] > s['lammaxr'])), yes_251, no_251, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamr'][s['k'] - 1] < s['lamminr'])), yes_245, no_245, s)
                return s
            def no_241(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])), yes_241, no_241, s)
            def yes_257(s):
                s = dict(s)
                def yes_258(s):
                    s = dict(s)
                    old_259 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(s['pgam_fixed']))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_259)
                    return s
                def no_258(s):
                    s = dict(s)
                    old_260 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_arith(F(0.0005714), _arith(_div(s['nc3d'][s['k'] - 1], F(1000000.0)), s['rho'][s['k'] - 1], 'mul'), 'mul'), F(0.2714), 'add')))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_260)
                    old_261 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_div(F(1.0), _arith(s['pgam'][s['k'] - 1], 2, 'pow')), F(1.0), 'sub')))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_261)
                    old_262 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.maximum(s['pgam'][s['k'] - 1], F(2.0))))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_262)
                    old_263 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.minimum(s['pgam'][s['k'] - 1], F(10.0))))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_263)
                    return s
                s = lax.cond((s["_pc"] == 0) & (s['dofix_pgam']), yes_258, no_258, s)
                old_264 = s['dumii']
                s['dumii'] = I(I(s['pgam'][s['k'] - 1]))
                s['dumii'] = jnp.where(s["_pc"] == 0, s['dumii'], old_264)
                old_265 = s['nu']
                s['nu'] = s['nu'].at[s['k'] - 1].set(F(_arith(s['dnu'][s['dumii'] - 1], _arith(_arith(s['dnu'][_arith(s['dumii'], 1, 'add') - 1], s['dnu'][s['dumii'] - 1], 'sub'), _arith(s['pgam'][s['k'] - 1], F(s['dumii']), 'sub'), 'mul'), 'add')))
                s['nu'] = jnp.where(s["_pc"] == 0, s['nu'], old_265)
                old_266 = s['lamc']
                s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(s['cons26'], s['nc3d'][s['k'] - 1], 'mul'), GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add')), 'mul'), _arith(s['qc3d'][s['k'] - 1], GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add')), 'mul')), _div(F(1.0), F(3.0)), 'pow')))
                s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_266)
                old_267 = s['lammin']
                s['lammin'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(6e-05)))
                s['lammin'] = jnp.where(s["_pc"] == 0, s['lammin'], old_267)
                old_268 = s['lammax']
                s['lammax'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(1e-06)))
                s['lammax'] = jnp.where(s["_pc"] == 0, s['lammax'], old_268)
                old_269 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_269)
                def yes_270(s):
                    s = dict(s)
                    old_271 = s['lamc']
                    s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(s['lammin']))
                    s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_271)
                    old_272 = s['tmpnum']
                    s['tmpnum'] = F(s['nc3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_272)
                    old_273 = s['nc3d']
                    s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_intrinsic('exp', _arith(_arith(_arith(_arith(F(3.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), _intrinsic('log', s['qc3d'][s['k'] - 1]), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'))), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add'))), 'sub')), s['cons26'])))
                    s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_273)
                    old_274 = s['sizefix_nc']
                    s['sizefix_nc'] = s['sizefix_nc'].at[s['k'] - 1].set(F(_arith(s['sizefix_nc'][s['k'] - 1], _div(_arith(s['nc3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_274)
                    return s
                def no_270(s):
                    s = dict(s)
                    def yes_275(s):
                        s = dict(s)
                        old_276 = s['lamc']
                        s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(s['lammax']))
                        s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_276)
                        old_277 = s['tmpnum']
                        s['tmpnum'] = F(s['nc3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_277)
                        old_278 = s['nc3d']
                        s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_intrinsic('exp', _arith(_arith(_arith(_arith(F(3.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), _intrinsic('log', s['qc3d'][s['k'] - 1]), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'))), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add'))), 'sub')), s['cons26'])))
                        s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_278)
                        old_279 = s['sizefix_nc']
                        s['sizefix_nc'] = s['sizefix_nc'].at[s['k'] - 1].set(F(_arith(s['sizefix_nc'][s['k'] - 1], _div(_arith(s['nc3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_279)
                        return s
                    def no_275(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lamc'][s['k'] - 1] > s['lammax'])), yes_275, no_275, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamc'][s['k'] - 1] < s['lammin'])), yes_270, no_270, s)
                return s
            def no_257(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])), yes_257, no_257, s)
            def yes_280(s):
                s = dict(s)
                old_281 = s['lams']
                s['lams'] = s['lams'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons1'], s['ns3d'][s['k'] - 1], 'mul'), s['qni3d'][s['k'] - 1]), _div(F(1.0), s['ds']), 'pow')))
                s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_281)
                old_282 = s['n0s']
                s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_arith(s['ns3d'][s['k'] - 1], s['lams'][s['k'] - 1], 'mul')))
                s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_282)
                old_283 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_283)
                def yes_284(s):
                    s = dict(s)
                    old_285 = s['lams']
                    s['lams'] = s['lams'].at[s['k'] - 1].set(F(s['lammins']))
                    s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_285)
                    old_286 = s['n0s']
                    s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lams'][s['k'] - 1], _arith(s['ds'], F(1.0), 'add'), 'pow'), s['qni3d'][s['k'] - 1], 'mul'), s['cons1'])))
                    s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_286)
                    old_287 = s['tmpnum']
                    s['tmpnum'] = F(s['ns3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_287)
                    old_288 = s['ns3d']
                    s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_div(s['n0s'][s['k'] - 1], s['lams'][s['k'] - 1])))
                    s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_288)
                    old_289 = s['sizefix_ns']
                    s['sizefix_ns'] = s['sizefix_ns'].at[s['k'] - 1].set(F(_arith(s['sizefix_ns'][s['k'] - 1], _div(_arith(s['ns3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_289)
                    return s
                def no_284(s):
                    s = dict(s)
                    def yes_290(s):
                        s = dict(s)
                        old_291 = s['lams']
                        s['lams'] = s['lams'].at[s['k'] - 1].set(F(s['lammaxs']))
                        s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_291)
                        old_292 = s['n0s']
                        s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lams'][s['k'] - 1], _arith(s['ds'], F(1.0), 'add'), 'pow'), s['qni3d'][s['k'] - 1], 'mul'), s['cons1'])))
                        s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_292)
                        old_293 = s['tmpnum']
                        s['tmpnum'] = F(s['ns3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_293)
                        old_294 = s['ns3d']
                        s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_div(s['n0s'][s['k'] - 1], s['lams'][s['k'] - 1])))
                        s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_294)
                        old_295 = s['sizefix_ns']
                        s['sizefix_ns'] = s['sizefix_ns'].at[s['k'] - 1].set(F(_arith(s['sizefix_ns'][s['k'] - 1], _div(_arith(s['ns3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_295)
                        return s
                    def no_290(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lams'][s['k'] - 1] > s['lammaxs'])), yes_290, no_290, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lams'][s['k'] - 1] < s['lammins'])), yes_284, no_284, s)
                return s
            def no_280(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])), yes_280, no_280, s)
            def yes_296(s):
                s = dict(s)
                old_297 = s['lamg']
                s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons2'], s['ng3d'][s['k'] - 1], 'mul'), s['qg3d'][s['k'] - 1]), _div(F(1.0), s['dg']), 'pow')))
                s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_297)
                old_298 = s['n0g']
                s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_arith(s['ng3d'][s['k'] - 1], s['lamg'][s['k'] - 1], 'mul')))
                s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_298)
                old_299 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_299)
                def yes_300(s):
                    s = dict(s)
                    old_301 = s['lamg']
                    s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(s['lamming']))
                    s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_301)
                    old_302 = s['n0g']
                    s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamg'][s['k'] - 1], _arith(s['dg'], F(1.0), 'add'), 'pow'), s['qg3d'][s['k'] - 1], 'mul'), s['cons2'])))
                    s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_302)
                    old_303 = s['tmpnum']
                    s['tmpnum'] = F(s['ng3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_303)
                    old_304 = s['ng3d']
                    s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_div(s['n0g'][s['k'] - 1], s['lamg'][s['k'] - 1])))
                    s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_304)
                    old_305 = s['sizefix_ng']
                    s['sizefix_ng'] = s['sizefix_ng'].at[s['k'] - 1].set(F(_arith(s['sizefix_ng'][s['k'] - 1], _div(_arith(s['ng3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_305)
                    return s
                def no_300(s):
                    s = dict(s)
                    def yes_306(s):
                        s = dict(s)
                        old_307 = s['lamg']
                        s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(s['lammaxg']))
                        s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_307)
                        old_308 = s['n0g']
                        s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamg'][s['k'] - 1], _arith(s['dg'], F(1.0), 'add'), 'pow'), s['qg3d'][s['k'] - 1], 'mul'), s['cons2'])))
                        s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_308)
                        old_309 = s['tmpnum']
                        s['tmpnum'] = F(s['ng3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_309)
                        old_310 = s['ng3d']
                        s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_div(s['n0g'][s['k'] - 1], s['lamg'][s['k'] - 1])))
                        s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_310)
                        old_311 = s['sizefix_ng']
                        s['sizefix_ng'] = s['sizefix_ng'].at[s['k'] - 1].set(F(_arith(s['sizefix_ng'][s['k'] - 1], _div(_arith(s['ng3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_311)
                        return s
                    def no_306(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lamg'][s['k'] - 1] > s['lammaxg'])), yes_306, no_306, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamg'][s['k'] - 1] < s['lamming'])), yes_300, no_300, s)
                return s
            def no_296(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] >= s['qsmall'])), yes_296, no_296, s)
            old_312 = s['prc']
            s['prc'] = s['prc'].at[s['k'] - 1].set(F(F(0.0)))
            s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_312)
            old_313 = s['nprc']
            s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(F(0.0)))
            s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_313)
            old_314 = s['nprc1']
            s['nprc1'] = s['nprc1'].at[s['k'] - 1].set(F(F(0.0)))
            s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_314)
            old_315 = s['pra']
            s['pra'] = s['pra'].at[s['k'] - 1].set(F(F(0.0)))
            s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_315)
            old_316 = s['npra']
            s['npra'] = s['npra'].at[s['k'] - 1].set(F(F(0.0)))
            s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_316)
            old_317 = s['nragg']
            s['nragg'] = s['nragg'].at[s['k'] - 1].set(F(F(0.0)))
            s['nragg'] = jnp.where(s["_pc"] == 0, s['nragg'], old_317)
            old_318 = s['psmlt']
            s['psmlt'] = s['psmlt'].at[s['k'] - 1].set(F(F(0.0)))
            s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_318)
            old_319 = s['nsmlts']
            s['nsmlts'] = s['nsmlts'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsmlts'] = jnp.where(s["_pc"] == 0, s['nsmlts'], old_319)
            old_320 = s['nsmltr']
            s['nsmltr'] = s['nsmltr'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsmltr'] = jnp.where(s["_pc"] == 0, s['nsmltr'], old_320)
            old_321 = s['pcc']
            s['pcc'] = s['pcc'].at[s['k'] - 1].set(F(F(0.0)))
            s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_321)
            old_322 = s['pre']
            s['pre'] = s['pre'].at[s['k'] - 1].set(F(F(0.0)))
            s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_322)
            old_323 = s['nsubc']
            s['nsubc'] = s['nsubc'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsubc'] = jnp.where(s["_pc"] == 0, s['nsubc'], old_323)
            old_324 = s['nsubr']
            s['nsubr'] = s['nsubr'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsubr'] = jnp.where(s["_pc"] == 0, s['nsubr'], old_324)
            old_325 = s['pracg']
            s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(F(0.0)))
            s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_325)
            old_326 = s['npracg']
            s['npracg'] = s['npracg'].at[s['k'] - 1].set(F(F(0.0)))
            s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_326)
            old_327 = s['psmlt']
            s['psmlt'] = s['psmlt'].at[s['k'] - 1].set(F(F(0.0)))
            s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_327)
            old_328 = s['evpms']
            s['evpms'] = s['evpms'].at[s['k'] - 1].set(F(F(0.0)))
            s['evpms'] = jnp.where(s["_pc"] == 0, s['evpms'], old_328)
            old_329 = s['pgmlt']
            s['pgmlt'] = s['pgmlt'].at[s['k'] - 1].set(F(F(0.0)))
            s['pgmlt'] = jnp.where(s["_pc"] == 0, s['pgmlt'], old_329)
            old_330 = s['evpmg']
            s['evpmg'] = s['evpmg'].at[s['k'] - 1].set(F(F(0.0)))
            s['evpmg'] = jnp.where(s["_pc"] == 0, s['evpmg'], old_330)
            old_331 = s['pracs']
            s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(F(0.0)))
            s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_331)
            old_332 = s['npracs']
            s['npracs'] = s['npracs'].at[s['k'] - 1].set(F(F(0.0)))
            s['npracs'] = jnp.where(s["_pc"] == 0, s['npracs'], old_332)
            old_333 = s['ngmltg']
            s['ngmltg'] = s['ngmltg'].at[s['k'] - 1].set(F(F(0.0)))
            s['ngmltg'] = jnp.where(s["_pc"] == 0, s['ngmltg'], old_333)
            old_334 = s['ngmltr']
            s['ngmltr'] = s['ngmltr'].at[s['k'] - 1].set(F(F(0.0)))
            s['ngmltr'] = jnp.where(s["_pc"] == 0, s['ngmltr'], old_334)
            def yes_335(s):
                s = dict(s)
                def yes_336(s):
                    s = dict(s)
                    old_337 = s['prc']
                    s['prc'] = s['prc'].at[s['k'] - 1].set(F(_arith(_arith(F(1350.0), _arith(s['qc3d'][s['k'] - 1], F(2.47), 'pow'), 'mul'), _arith(_arith(_div(s['nc3d'][s['k'] - 1], F(1000000.0)), s['rho'][s['k'] - 1], 'mul'), (-F(1.79)), 'pow'), 'mul')))
                    s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_337)
                    old_338 = s['nprc1']
                    s['nprc1'] = s['nprc1'].at[s['k'] - 1].set(F(_div(s['prc'][s['k'] - 1], s['cons29'])))
                    s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_338)
                    old_339 = s['nprc']
                    s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(_div(s['prc'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['nc3d'][s['k'] - 1]))))
                    s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_339)
                    old_340 = s['nprc']
                    s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(jnp.minimum(s['nprc'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']))))
                    s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_340)
                    return s
                def no_336(s):
                    s = dict(s)
                    def yes_341(s):
                        s = dict(s)
                        old_342 = s['dum']
                        s['dum'] = F(_arith(F(1.0), _div(s['qc3d'][s['k'] - 1], _arith(s['qc3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')), 'sub'))
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_342)
                        old_343 = s['dum1']
                        s['dum1'] = F(_arith(_arith(F(600.0), _arith(s['dum'], F(0.68), 'pow'), 'mul'), _arith(_arith(F(1.0), _arith(s['dum'], F(0.68), 'pow'), 'sub'), 3, 'pow'), 'mul'))
                        s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_343)
                        old_344 = s['prc']
                        s['prc'] = s['prc'].at[s['k'] - 1].set(F(_div(_arith(_arith(_div(_arith(_div(_arith(_arith(_div(F(9440000000.0), _arith(F(20.0), F(2.6e-07), 'mul')), _arith(s['nu'][s['k'] - 1], F(2.0), 'add'), 'mul'), _arith(s['nu'][s['k'] - 1], F(4.0), 'add'), 'mul'), _arith(_arith(s['nu'][s['k'] - 1], F(1.0), 'add'), 2, 'pow')), _arith(_div(_arith(s['rho'][s['k'] - 1], s['qc3d'][s['k'] - 1], 'mul'), F(1000.0)), 4, 'pow'), 'mul'), _arith(_div(_arith(s['rho'][s['k'] - 1], s['nc3d'][s['k'] - 1], 'mul'), F(1000000.0)), 2, 'pow')), _arith(F(1.0), _div(s['dum1'], _arith(_arith(F(1.0), s['dum'], 'sub'), 2, 'pow')), 'add'), 'mul'), F(1000.0), 'mul'), s['rho'][s['k'] - 1])))
                        s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_344)
                        old_345 = s['nprc']
                        s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['prc'][s['k'] - 1], F(2.0), 'mul'), F(2.6e-07)), F(1000.0), 'mul')))
                        s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_345)
                        old_346 = s['nprc1']
                        s['nprc1'] = s['nprc1'].at[s['k'] - 1].set(F(_arith(F(0.5), s['nprc'][s['k'] - 1], 'mul')))
                        s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_346)
                        return s
                    def no_341(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 1)), yes_341, no_341, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 0)), yes_336, no_336, s)
                return s
            def no_335(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] >= F(1e-06))), yes_335, no_335, s)
            def yes_347(s):
                s = dict(s)
                old_348 = s['ums']
                s['ums'] = F(_div(_arith(s['asn'][s['k'] - 1], s['cons3'], 'mul'), _arith(s['lams'][s['k'] - 1], s['bs'], 'pow')))
                s['ums'] = jnp.where(s["_pc"] == 0, s['ums'], old_348)
                old_349 = s['umr']
                s['umr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons4'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_349)
                old_350 = s['uns']
                s['uns'] = F(_div(_arith(s['asn'][s['k'] - 1], s['cons5'], 'mul'), _arith(s['lams'][s['k'] - 1], s['bs'], 'pow')))
                s['uns'] = jnp.where(s["_pc"] == 0, s['uns'], old_350)
                old_351 = s['unr']
                s['unr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons6'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_351)
                old_352 = s['dum']
                s['dum'] = F(_arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.54), 'pow'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_352)
                old_353 = s['ums']
                s['ums'] = F(jnp.minimum(s['ums'], _arith(F(1.2), s['dum'], 'mul')))
                s['ums'] = jnp.where(s["_pc"] == 0, s['ums'], old_353)
                old_354 = s['uns']
                s['uns'] = F(jnp.minimum(s['uns'], _arith(F(1.2), s['dum'], 'mul')))
                s['uns'] = jnp.where(s["_pc"] == 0, s['uns'], old_354)
                old_355 = s['umr']
                s['umr'] = F(jnp.minimum(s['umr'], _arith(F(9.1), s['dum'], 'mul')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_355)
                old_356 = s['unr']
                s['unr'] = F(jnp.minimum(s['unr'], _arith(F(9.1), s['dum'], 'mul')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_356)
                old_357 = s['pracs']
                s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(_arith(s['cons31'], _arith(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(F(1.2), s['umr'], 'mul'), _arith(F(0.95), s['ums'], 'mul'), 'sub'), 2, 'pow'), _arith(_arith(F(0.08), s['ums'], 'mul'), s['umr'], 'mul'), 'add'), F(0.5), 'pow'), s['rho'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(s['lams'][s['k'] - 1], 3, 'pow')), _arith(_arith(_div(F(5.0), _arith(_arith(s['lams'][s['k'] - 1], 3, 'pow'), s['lamr'][s['k'] - 1], 'mul')), _div(F(2.0), _arith(_arith(s['lams'][s['k'] - 1], 2, 'pow'), _arith(s['lamr'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(0.5), _arith(s['lams'][s['k'] - 1], _arith(s['lamr'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul'), 'mul')))
                s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_357)
                return s
            def no_347(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qr3d'][s['k'] - 1] >= F(1e-08))) & ((s['qni3d'][s['k'] - 1] >= F(1e-08))))), yes_347, no_347, s)
            def yes_358(s):
                s = dict(s)
                old_359 = s['umg']
                s['umg'] = F(_div(_arith(s['agn'][s['k'] - 1], s['cons7'], 'mul'), _arith(s['lamg'][s['k'] - 1], s['bg'], 'pow')))
                s['umg'] = jnp.where(s["_pc"] == 0, s['umg'], old_359)
                old_360 = s['umr']
                s['umr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons4'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_360)
                old_361 = s['ung']
                s['ung'] = F(_div(_arith(s['agn'][s['k'] - 1], s['cons8'], 'mul'), _arith(s['lamg'][s['k'] - 1], s['bg'], 'pow')))
                s['ung'] = jnp.where(s["_pc"] == 0, s['ung'], old_361)
                old_362 = s['unr']
                s['unr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons6'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_362)
                old_363 = s['dum']
                s['dum'] = F(_arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.54), 'pow'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_363)
                old_364 = s['umg']
                s['umg'] = F(jnp.minimum(s['umg'], _arith(F(20.0), s['dum'], 'mul')))
                s['umg'] = jnp.where(s["_pc"] == 0, s['umg'], old_364)
                old_365 = s['ung']
                s['ung'] = F(jnp.minimum(s['ung'], _arith(F(20.0), s['dum'], 'mul')))
                s['ung'] = jnp.where(s["_pc"] == 0, s['ung'], old_365)
                old_366 = s['umr']
                s['umr'] = F(jnp.minimum(s['umr'], _arith(F(9.1), s['dum'], 'mul')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_366)
                old_367 = s['unr']
                s['unr'] = F(jnp.minimum(s['unr'], _arith(F(9.1), s['dum'], 'mul')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_367)
                old_368 = s['pracg']
                s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(_arith(s['cons41'], _arith(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(F(1.2), s['umr'], 'mul'), _arith(F(0.95), s['umg'], 'mul'), 'sub'), 2, 'pow'), _arith(_arith(F(0.08), s['umg'], 'mul'), s['umr'], 'mul'), 'add'), F(0.5), 'pow'), s['rho'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0g'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], 3, 'pow')), _arith(_arith(_div(F(5.0), _arith(_arith(s['lamr'][s['k'] - 1], 3, 'pow'), s['lamg'][s['k'] - 1], 'mul')), _div(F(2.0), _arith(_arith(s['lamr'][s['k'] - 1], 2, 'pow'), _arith(s['lamg'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(0.5), _arith(s['lamr'][s['k'] - 1], _arith(s['lamg'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul'), 'mul')))
                s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_368)
                old_369 = s['dum']
                s['dum'] = F(_div(s['pracg'][s['k'] - 1], F(5.2e-07)))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_369)
                old_370 = s['npracg']
                s['npracg'] = s['npracg'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_arith(_arith(s['cons32'], s['rho'][s['k'] - 1], 'mul'), _arith(_arith(_arith(F(1.7), _arith(_arith(s['unr'], s['ung'], 'sub'), 2, 'pow'), 'mul'), _arith(_arith(F(0.3), s['unr'], 'mul'), s['ung'], 'mul'), 'add'), F(0.5), 'pow'), 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0g'][s['k'] - 1], 'mul'), _arith(_arith(_div(F(1.0), _arith(_arith(s['lamr'][s['k'] - 1], 3, 'pow'), s['lamg'][s['k'] - 1], 'mul')), _div(F(1.0), _arith(_arith(s['lamr'][s['k'] - 1], 2, 'pow'), _arith(s['lamg'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(1.0), _arith(s['lamr'][s['k'] - 1], _arith(s['lamg'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul')))
                s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_370)
                old_371 = s['npracg']
                s['npracg'] = s['npracg'].at[s['k'] - 1].set(F(jnp.maximum(_arith(s['npracg'][s['k'] - 1], s['dum'], 'sub'), F(0.0))))
                s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_371)
                return s
            def no_358(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qr3d'][s['k'] - 1] >= F(1e-08))) & ((s['qg3d'][s['k'] - 1] >= F(1e-08))))), yes_358, no_358, s)
            def yes_372(s):
                s = dict(s)
                def yes_373(s):
                    s = dict(s)
                    old_374 = s['dum']
                    s['dum'] = F(_arith(s['qc3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'mul'))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_374)
                    old_375 = s['pra']
                    s['pra'] = s['pra'].at[s['k'] - 1].set(F(_arith(F(67.0), _arith(s['dum'], F(1.15), 'pow'), 'mul')))
                    s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_375)
                    old_376 = s['npra']
                    s['npra'] = s['npra'].at[s['k'] - 1].set(F(_div(s['pra'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['nc3d'][s['k'] - 1]))))
                    s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_376)
                    return s
                def no_373(s):
                    s = dict(s)
                    def yes_377(s):
                        s = dict(s)
                        old_378 = s['dum']
                        s['dum'] = F(_arith(F(1.0), _div(s['qc3d'][s['k'] - 1], _arith(s['qc3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')), 'sub'))
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_378)
                        old_379 = s['dum1']
                        s['dum1'] = F(_arith(_div(s['dum'], _arith(s['dum'], F(0.0005), 'add')), 4, 'pow'))
                        s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_379)
                        old_380 = s['pra']
                        s['pra'] = s['pra'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_div(_arith(F(5780.0), s['rho'][s['k'] - 1], 'mul'), F(1000.0)), s['qc3d'][s['k'] - 1], 'mul'), s['qr3d'][s['k'] - 1], 'mul'), s['dum1'], 'mul')))
                        s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_380)
                        old_381 = s['npra']
                        s['npra'] = s['npra'].at[s['k'] - 1].set(F(_div(_arith(_div(_arith(_div(_arith(s['pra'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), F(1000.0)), _div(_arith(s['nc3d'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), F(1000000.0)), 'mul'), _div(_arith(s['qc3d'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), F(1000.0))), F(1000000.0), 'mul'), s['rho'][s['k'] - 1])))
                        s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_381)
                        return s
                    def no_377(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 1)), yes_377, no_377, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 0)), yes_373, no_373, s)
                return s
            def no_372(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qr3d'][s['k'] - 1] >= F(1e-08))) & ((s['qc3d'][s['k'] - 1] >= F(1e-08))))), yes_372, no_372, s)
            def yes_382(s):
                s = dict(s)
                old_383 = s['dum1']
                s['dum1'] = F(F(0.0003))
                s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_383)
                def yes_384(s):
                    s = dict(s)
                    old_385 = s['dum']
                    s['dum'] = F(F(1.0))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_385)
                    return s
                def no_384(s):
                    s = dict(s)
                    def yes_386(s):
                        s = dict(s)
                        old_387 = s['dum']
                        s['dum'] = F(_arith(F(2.0), _intrinsic('exp', _arith(F(2300.0), _arith(_div(F(1.0), s['lamr'][s['k'] - 1]), s['dum1'], 'sub'), 'mul')), 'sub'))
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_387)
                        return s
                    def no_386(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((_div(F(1.0), s['lamr'][s['k'] - 1]) >= s['dum1'])), yes_386, no_386, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((_div(F(1.0), s['lamr'][s['k'] - 1]) < s['dum1'])), yes_384, no_384, s)
                old_388 = s['nragg']
                s['nragg'] = s['nragg'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_arith((-F(5.78)), s['dum'], 'mul'), s['nr3d'][s['k'] - 1], 'mul'), s['qr3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul')))
                s['nragg'] = jnp.where(s["_pc"] == 0, s['nragg'], old_388)
                return s
            def no_382(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= F(1e-08))), yes_382, no_382, s)
            def yes_389(s):
                s = dict(s)
                old_390 = s['epsr']
                s['epsr'] = F(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), _arith(_div(s['f1r'], _arith(s['lamr'][s['k'] - 1], s['lamr'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2r'], _arith(_div(_arith(s['arn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons9'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['cons34'], 'pow')), 'add'), 'mul'))
                s['epsr'] = jnp.where(s["_pc"] == 0, s['epsr'], old_390)
                return s
            def no_389(s):
                s = dict(s)
                old_391 = s['epsr']
                s['epsr'] = F(F(0.0))
                s['epsr'] = jnp.where(s["_pc"] == 0, s['epsr'], old_391)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])), yes_389, no_389, s)
            def yes_392(s):
                s = dict(s)
                old_393 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(_div(_arith(s['epsr'], _arith(s['qv3d'][s['k'] - 1], s['qvs'][s['k'] - 1], 'sub'), 'mul'), s['ab'][s['k'] - 1])))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_393)
                old_394 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(jnp.minimum(s['pre'][s['k'] - 1], F(0.0))))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_394)
                return s
            def no_392(s):
                s = dict(s)
                old_395 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(F(0.0)))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_395)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qv3d'][s['k'] - 1] < s['qvs'][s['k'] - 1])), yes_392, no_392, s)
            def yes_396(s):
                s = dict(s)
                old_397 = s['dum']
                s['dum'] = F(_arith(_arith(_div((-s['cpw']), s['xlf'][s['k'] - 1]), jnp.maximum(_arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), F(0.0)), 'mul'), s['pracs'][s['k'] - 1], 'mul'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_397)
                old_398 = s['psmlt']
                s['psmlt'] = s['psmlt'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_div(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0s'][s['k'] - 1], 'mul'), s['kap'][s['k'] - 1], 'mul'), jnp.minimum(_arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), F(0.0)), 'mul'), s['xlf'][s['k'] - 1]), s['rho'][s['k'] - 1], 'mul'), _arith(_div(s['f1s'], _arith(s['lams'][s['k'] - 1], s['lams'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2s'], _arith(_div(_arith(s['asn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons10'], 'mul'), _arith(s['lams'][s['k'] - 1], s['cons35'], 'pow')), 'add'), 'mul'), s['dum'], 'add')))
                s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_398)
                def yes_399(s):
                    s = dict(s)
                    old_400 = s['epss']
                    s['epss'] = F(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0s'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), _arith(_div(s['f1s'], _arith(s['lams'][s['k'] - 1], s['lams'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2s'], _arith(_div(_arith(s['asn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons10'], 'mul'), _arith(s['lams'][s['k'] - 1], s['cons35'], 'pow')), 'add'), 'mul'))
                    s['epss'] = jnp.where(s["_pc"] == 0, s['epss'], old_400)
                    old_401 = s['evpms']
                    s['evpms'] = s['evpms'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['qv3d'][s['k'] - 1], s['qvs'][s['k'] - 1], 'sub'), s['epss'], 'mul'), s['ab'][s['k'] - 1])))
                    s['evpms'] = jnp.where(s["_pc"] == 0, s['evpms'], old_401)
                    old_402 = s['evpms']
                    s['evpms'] = s['evpms'].at[s['k'] - 1].set(F(jnp.maximum(s['evpms'][s['k'] - 1], s['psmlt'][s['k'] - 1])))
                    s['evpms'] = jnp.where(s["_pc"] == 0, s['evpms'], old_402)
                    old_403 = s['psmlt']
                    s['psmlt'] = s['psmlt'].at[s['k'] - 1].set(F(_arith(s['psmlt'][s['k'] - 1], s['evpms'][s['k'] - 1], 'sub')))
                    s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_403)
                    return s
                def no_399(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['qvqvs'][s['k'] - 1] < F(1.0))), yes_399, no_399, s)
                return s
            def no_396(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= F(1e-08))), yes_396, no_396, s)
            def yes_404(s):
                s = dict(s)
                old_405 = s['dum']
                s['dum'] = F(_arith(_arith(_div((-s['cpw']), s['xlf'][s['k'] - 1]), jnp.maximum(_arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), F(0.0)), 'mul'), s['pracg'][s['k'] - 1], 'mul'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_405)
                old_406 = s['pgmlt']
                s['pgmlt'] = s['pgmlt'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_div(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0g'][s['k'] - 1], 'mul'), s['kap'][s['k'] - 1], 'mul'), _arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), 'mul'), s['xlf'][s['k'] - 1]), s['rho'][s['k'] - 1], 'mul'), _arith(_div(s['f1s'], _arith(s['lamg'][s['k'] - 1], s['lamg'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2s'], _arith(_div(_arith(s['agn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons11'], 'mul'), _arith(s['lamg'][s['k'] - 1], s['cons36'], 'pow')), 'add'), 'mul'), s['dum'], 'add')))
                s['pgmlt'] = jnp.where(s["_pc"] == 0, s['pgmlt'], old_406)
                def yes_407(s):
                    s = dict(s)
                    old_408 = s['epsg']
                    s['epsg'] = F(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0g'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), _arith(_div(s['f1s'], _arith(s['lamg'][s['k'] - 1], s['lamg'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2s'], _arith(_div(_arith(s['agn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons11'], 'mul'), _arith(s['lamg'][s['k'] - 1], s['cons36'], 'pow')), 'add'), 'mul'))
                    s['epsg'] = jnp.where(s["_pc"] == 0, s['epsg'], old_408)
                    old_409 = s['evpmg']
                    s['evpmg'] = s['evpmg'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['qv3d'][s['k'] - 1], s['qvs'][s['k'] - 1], 'sub'), s['epsg'], 'mul'), s['ab'][s['k'] - 1])))
                    s['evpmg'] = jnp.where(s["_pc"] == 0, s['evpmg'], old_409)
                    old_410 = s['evpmg']
                    s['evpmg'] = s['evpmg'].at[s['k'] - 1].set(F(jnp.maximum(s['evpmg'][s['k'] - 1], s['pgmlt'][s['k'] - 1])))
                    s['evpmg'] = jnp.where(s["_pc"] == 0, s['evpmg'], old_410)
                    old_411 = s['pgmlt']
                    s['pgmlt'] = s['pgmlt'].at[s['k'] - 1].set(F(_arith(s['pgmlt'][s['k'] - 1], s['evpmg'][s['k'] - 1], 'sub')))
                    s['pgmlt'] = jnp.where(s["_pc"] == 0, s['pgmlt'], old_411)
                    return s
                def no_407(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['qvqvs'][s['k'] - 1] < F(1.0))), yes_407, no_407, s)
                return s
            def no_404(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] >= F(1e-08))), yes_404, no_404, s)
            old_412 = s['pracg']
            s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(F(0.0)))
            s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_412)
            old_413 = s['pracs']
            s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(F(0.0)))
            s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_413)
            old_414 = s['dum']
            s['dum'] = F(_arith(_arith(s['prc'][s['k'] - 1], s['pra'][s['k'] - 1], 'add'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_414)
            def yes_415(s):
                s = dict(s)
                old_416 = s['ratio']
                s['ratio'] = F(_div(s['qc3d'][s['k'] - 1], s['dum']))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_416)
                old_417 = s['prc']
                s['prc'] = s['prc'].at[s['k'] - 1].set(F(_arith(s['prc'][s['k'] - 1], s['ratio'], 'mul')))
                s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_417)
                old_418 = s['pra']
                s['pra'] = s['pra'].at[s['k'] - 1].set(F(_arith(s['pra'][s['k'] - 1], s['ratio'], 'mul')))
                s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_418)
                return s
            def no_415(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qc3d'][s['k'] - 1])) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])))), yes_415, no_415, s)
            old_419 = s['dum']
            s['dum'] = F(_arith(_arith(_arith((-s['psmlt'][s['k'] - 1]), s['evpms'][s['k'] - 1], 'sub'), s['pracs'][s['k'] - 1], 'add'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_419)
            def yes_420(s):
                s = dict(s)
                old_421 = s['ratio']
                s['ratio'] = F(_div(s['qni3d'][s['k'] - 1], s['dum']))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_421)
                old_422 = s['psmlt']
                s['psmlt'] = s['psmlt'].at[s['k'] - 1].set(F(_arith(s['psmlt'][s['k'] - 1], s['ratio'], 'mul')))
                s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_422)
                old_423 = s['evpms']
                s['evpms'] = s['evpms'].at[s['k'] - 1].set(F(_arith(s['evpms'][s['k'] - 1], s['ratio'], 'mul')))
                s['evpms'] = jnp.where(s["_pc"] == 0, s['evpms'], old_423)
                old_424 = s['pracs']
                s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(_arith(s['pracs'][s['k'] - 1], s['ratio'], 'mul')))
                s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_424)
                return s
            def no_420(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qni3d'][s['k'] - 1])) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])))), yes_420, no_420, s)
            old_425 = s['dum']
            s['dum'] = F(_arith(_arith(_arith((-s['pgmlt'][s['k'] - 1]), s['evpmg'][s['k'] - 1], 'sub'), s['pracg'][s['k'] - 1], 'add'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_425)
            def yes_426(s):
                s = dict(s)
                old_427 = s['ratio']
                s['ratio'] = F(_div(s['qg3d'][s['k'] - 1], s['dum']))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_427)
                old_428 = s['pgmlt']
                s['pgmlt'] = s['pgmlt'].at[s['k'] - 1].set(F(_arith(s['pgmlt'][s['k'] - 1], s['ratio'], 'mul')))
                s['pgmlt'] = jnp.where(s["_pc"] == 0, s['pgmlt'], old_428)
                old_429 = s['evpmg']
                s['evpmg'] = s['evpmg'].at[s['k'] - 1].set(F(_arith(s['evpmg'][s['k'] - 1], s['ratio'], 'mul')))
                s['evpmg'] = jnp.where(s["_pc"] == 0, s['evpmg'], old_429)
                old_430 = s['pracg']
                s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(_arith(s['pracg'][s['k'] - 1], s['ratio'], 'mul')))
                s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_430)
                return s
            def no_426(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qg3d'][s['k'] - 1])) & ((s['qg3d'][s['k'] - 1] >= s['qsmall'])))), yes_426, no_426, s)
            old_431 = s['dum']
            s['dum'] = F(_arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['pracs'][s['k'] - 1]), s['pracg'][s['k'] - 1], 'sub'), s['pre'][s['k'] - 1], 'sub'), s['pra'][s['k'] - 1], 'sub'), s['prc'][s['k'] - 1], 'sub'), s['psmlt'][s['k'] - 1], 'add'), s['pgmlt'][s['k'] - 1], 'add'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_431)
            def yes_432(s):
                s = dict(s)
                old_433 = s['ratio']
                s['ratio'] = F(_div(_arith(_arith(_arith(_arith(_arith(_arith(_div(s['qr3d'][s['k'] - 1], s['dt']), s['pracs'][s['k'] - 1], 'add'), s['pracg'][s['k'] - 1], 'add'), s['pra'][s['k'] - 1], 'add'), s['prc'][s['k'] - 1], 'add'), s['psmlt'][s['k'] - 1], 'sub'), s['pgmlt'][s['k'] - 1], 'sub'), (-s['pre'][s['k'] - 1])))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_433)
                old_434 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(_arith(s['pre'][s['k'] - 1], s['ratio'], 'mul')))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_434)
                return s
            def no_432(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qr3d'][s['k'] - 1])) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])))), yes_432, no_432, s)
            def yes_435(s):
                s = dict(s)
                return s
            def no_435(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['pre'][s['k'] - 1] < F(0.0))) & ((s['cf3d'][s['k'] - 1] > s['cloud_frac_thresh'])))), yes_435, no_435, s)
            def yes_436(s):
                s = dict(s)
                lax.cond(s["_pc"] == 0, lambda _: jax.debug.callback(_fatal, s["_pc"], ordered=True), lambda _: None, None)
                return s
            def no_436(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['pracs'][s['k'] - 1] != F(0.0))) | ((s['pracg'][s['k'] - 1] != F(0.0))))), yes_436, no_436, s)
            old_437 = s['qv3dten']
            s['qv3dten'] = s['qv3dten'].at[s['k'] - 1].set(F(_arith(s['qv3dten'][s['k'] - 1], _arith(_arith((-s['pre'][s['k'] - 1]), s['evpms'][s['k'] - 1], 'sub'), s['evpmg'][s['k'] - 1], 'sub'), 'add')))
            s['qv3dten'] = jnp.where(s["_pc"] == 0, s['qv3dten'], old_437)
            old_438 = s['t3dten']
            s['t3dten'] = s['t3dten'].at[s['k'] - 1].set(F(_arith(s['t3dten'][s['k'] - 1], _div(_arith(_arith(_arith(s['pre'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(s['evpms'][s['k'] - 1], s['evpmg'][s['k'] - 1], 'add'), s['xxls'][s['k'] - 1], 'mul'), 'add'), _arith(_arith(_arith(_arith(s['psmlt'][s['k'] - 1], s['pgmlt'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'sub'), s['pracg'][s['k'] - 1], 'sub'), s['xlf'][s['k'] - 1], 'mul'), 'add'), s['cpm'][s['k'] - 1]), 'add')))
            s['t3dten'] = jnp.where(s["_pc"] == 0, s['t3dten'], old_438)
            old_439 = s['qc3dten']
            s['qc3dten'] = s['qc3dten'].at[s['k'] - 1].set(F(_arith(s['qc3dten'][s['k'] - 1], _arith((-s['pra'][s['k'] - 1]), s['prc'][s['k'] - 1], 'sub'), 'add')))
            s['qc3dten'] = jnp.where(s["_pc"] == 0, s['qc3dten'], old_439)
            old_440 = s['qr3dten']
            s['qr3dten'] = s['qr3dten'].at[s['k'] - 1].set(F(_arith(s['qr3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(s['pre'][s['k'] - 1], s['pra'][s['k'] - 1], 'add'), s['prc'][s['k'] - 1], 'add'), s['psmlt'][s['k'] - 1], 'sub'), s['pgmlt'][s['k'] - 1], 'sub'), s['pracs'][s['k'] - 1], 'add'), s['pracg'][s['k'] - 1], 'add'), 'add')))
            s['qr3dten'] = jnp.where(s["_pc"] == 0, s['qr3dten'], old_440)
            old_441 = s['qni3dten']
            s['qni3dten'] = s['qni3dten'].at[s['k'] - 1].set(F(_arith(s['qni3dten'][s['k'] - 1], _arith(_arith(s['psmlt'][s['k'] - 1], s['evpms'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'sub'), 'add')))
            s['qni3dten'] = jnp.where(s["_pc"] == 0, s['qni3dten'], old_441)
            old_442 = s['qg3dten']
            s['qg3dten'] = s['qg3dten'].at[s['k'] - 1].set(F(_arith(s['qg3dten'][s['k'] - 1], _arith(_arith(s['pgmlt'][s['k'] - 1], s['evpmg'][s['k'] - 1], 'add'), s['pracg'][s['k'] - 1], 'sub'), 'add')))
            s['qg3dten'] = jnp.where(s["_pc"] == 0, s['qg3dten'], old_442)
            old_443 = s['nc3dten']
            s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], _arith((-s['npra'][s['k'] - 1]), s['nprc'][s['k'] - 1], 'sub'), 'add')))
            s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_443)
            old_444 = s['nr3dten']
            s['nr3dten'] = s['nr3dten'].at[s['k'] - 1].set(F(_arith(s['nr3dten'][s['k'] - 1], _arith(_arith(s['nprc1'][s['k'] - 1], s['nragg'][s['k'] - 1], 'add'), s['npracg'][s['k'] - 1], 'sub'), 'add')))
            s['nr3dten'] = jnp.where(s["_pc"] == 0, s['nr3dten'], old_444)
            def yes_445(s):
                s = dict(s)
                old_446 = s['dum']
                s['dum'] = F(_div(_arith(s['pre'][s['k'] - 1], s['dt'], 'mul'), s['qr3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_446)
                old_447 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_447)
                old_448 = s['nsubr']
                s['nsubr'] = s['nsubr'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['nr3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['nsubr'] = jnp.where(s["_pc"] == 0, s['nsubr'], old_448)
                return s
            def no_445(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['pre'][s['k'] - 1] < F(0.0))), yes_445, no_445, s)
            def yes_449(s):
                s = dict(s)
                old_450 = s['dum']
                s['dum'] = F(_div(_arith(_arith(s['evpms'][s['k'] - 1], s['psmlt'][s['k'] - 1], 'add'), s['dt'], 'mul'), s['qni3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_450)
                old_451 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_451)
                old_452 = s['nsmlts']
                s['nsmlts'] = s['nsmlts'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['ns3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['nsmlts'] = jnp.where(s["_pc"] == 0, s['nsmlts'], old_452)
                return s
            def no_449(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((_arith(s['evpms'][s['k'] - 1], s['psmlt'][s['k'] - 1], 'add') < F(0.0))), yes_449, no_449, s)
            def yes_453(s):
                s = dict(s)
                old_454 = s['dum']
                s['dum'] = F(_div(_arith(s['psmlt'][s['k'] - 1], s['dt'], 'mul'), s['qni3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_454)
                old_455 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_455)
                old_456 = s['nsmltr']
                s['nsmltr'] = s['nsmltr'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['ns3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['nsmltr'] = jnp.where(s["_pc"] == 0, s['nsmltr'], old_456)
                return s
            def no_453(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['psmlt'][s['k'] - 1] < F(0.0))), yes_453, no_453, s)
            def yes_457(s):
                s = dict(s)
                old_458 = s['dum']
                s['dum'] = F(_div(_arith(_arith(s['evpmg'][s['k'] - 1], s['pgmlt'][s['k'] - 1], 'add'), s['dt'], 'mul'), s['qg3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_458)
                old_459 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_459)
                old_460 = s['ngmltg']
                s['ngmltg'] = s['ngmltg'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['ng3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['ngmltg'] = jnp.where(s["_pc"] == 0, s['ngmltg'], old_460)
                return s
            def no_457(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((_arith(s['evpmg'][s['k'] - 1], s['pgmlt'][s['k'] - 1], 'add') < F(0.0))), yes_457, no_457, s)
            def yes_461(s):
                s = dict(s)
                old_462 = s['dum']
                s['dum'] = F(_div(_arith(s['pgmlt'][s['k'] - 1], s['dt'], 'mul'), s['qg3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_462)
                old_463 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_463)
                old_464 = s['ngmltr']
                s['ngmltr'] = s['ngmltr'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['ng3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['ngmltr'] = jnp.where(s["_pc"] == 0, s['ngmltr'], old_464)
                return s
            def no_461(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['pgmlt'][s['k'] - 1] < F(0.0))), yes_461, no_461, s)
            old_465 = s['ns3dten']
            s['ns3dten'] = s['ns3dten'].at[s['k'] - 1].set(F(_arith(s['ns3dten'][s['k'] - 1], s['nsmlts'][s['k'] - 1], 'add')))
            s['ns3dten'] = jnp.where(s["_pc"] == 0, s['ns3dten'], old_465)
            old_466 = s['ng3dten']
            s['ng3dten'] = s['ng3dten'].at[s['k'] - 1].set(F(_arith(s['ng3dten'][s['k'] - 1], s['ngmltg'][s['k'] - 1], 'add')))
            s['ng3dten'] = jnp.where(s["_pc"] == 0, s['ng3dten'], old_466)
            old_467 = s['nr3dten']
            s['nr3dten'] = s['nr3dten'].at[s['k'] - 1].set(F(_arith(s['nr3dten'][s['k'] - 1], _arith(_arith(s['nsubr'][s['k'] - 1], s['nsmltr'][s['k'] - 1], 'sub'), s['ngmltr'][s['k'] - 1], 'sub'), 'add')))
            s['nr3dten'] = jnp.where(s["_pc"] == 0, s['nr3dten'], old_467)
            s["_pc"] = jnp.where(s["_pc"] == 300, I(0), s["_pc"])
            def yes_468(s):
                s = dict(s)
                old_469 = s['dumt']
                s['dumt'] = F(_arith(s['t3d'][s['k'] - 1], _arith(s['dt'], s['t3dten'][s['k'] - 1], 'mul'), 'add'))
                s['dumt'] = jnp.where(s["_pc"] == 0, s['dumt'], old_469)
                old_470 = s['dumqv']
                s['dumqv'] = F(_arith(s['qv3d'][s['k'] - 1], _arith(s['dt'], s['qv3dten'][s['k'] - 1], 'mul'), 'add'))
                s['dumqv'] = jnp.where(s["_pc"] == 0, s['dumqv'], old_470)
                old_471 = s['dum']
                s['dum'] = F(jnp.minimum(_arith(F(0.99), s['pres'][s['k'] - 1], 'mul'), POLYSVP(s['dumt'], 0)))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_471)
                old_472 = s['dumqss']
                s['dumqss'] = F(_div(_arith(s['ep_2'], s['dum'], 'mul'), _arith(s['pres'][s['k'] - 1], s['dum'], 'sub')))
                s['dumqss'] = jnp.where(s["_pc"] == 0, s['dumqss'], old_472)
                old_473 = s['dumqc']
                s['dumqc'] = F(_arith(s['qc3d'][s['k'] - 1], _arith(s['dt'], s['qc3dten'][s['k'] - 1], 'mul'), 'add'))
                s['dumqc'] = jnp.where(s["_pc"] == 0, s['dumqc'], old_473)
                old_474 = s['dumqc']
                s['dumqc'] = F(jnp.maximum(s['dumqc'], F(0.0)))
                s['dumqc'] = jnp.where(s["_pc"] == 0, s['dumqc'], old_474)
                old_475 = s['dums']
                s['dums'] = F(_arith(s['dumqv'], s['dumqss'], 'sub'))
                s['dums'] = jnp.where(s["_pc"] == 0, s['dums'], old_475)
                old_476 = s['pcc']
                s['pcc'] = s['pcc'].at[s['k'] - 1].set(F(_div(_div(s['dums'], _arith(F(1.0), _div(_arith(_arith(s['xxlv'][s['k'] - 1], 2, 'pow'), s['dumqss'], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rv'], 'mul'), _arith(s['dumt'], 2, 'pow'), 'mul')), 'add')), s['dt'])))
                s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_476)
                def yes_477(s):
                    s = dict(s)
                    old_478 = s['pcc']
                    s['pcc'] = s['pcc'].at[s['k'] - 1].set(F(_div((-_arith(s['qc3d'][s['k'] - 1], _arith(s['dt'], s['qc3dten'][s['k'] - 1], 'mul'), 'add')), s['dt'])))
                    s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_478)
                    return s
                def no_477(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((_arith(_arith(_arith(s['pcc'][s['k'] - 1], s['dt'], 'mul'), s['qc3d'][s['k'] - 1], 'add'), _arith(s['dt'], s['qc3dten'][s['k'] - 1], 'mul'), 'add') < F(0.0))), yes_477, no_477, s)
                old_479 = s['qv3dten']
                s['qv3dten'] = s['qv3dten'].at[s['k'] - 1].set(F(_arith(s['qv3dten'][s['k'] - 1], s['pcc'][s['k'] - 1], 'sub')))
                s['qv3dten'] = jnp.where(s["_pc"] == 0, s['qv3dten'], old_479)
                old_480 = s['t3dten']
                s['t3dten'] = s['t3dten'].at[s['k'] - 1].set(F(_arith(s['t3dten'][s['k'] - 1], _div(_arith(s['pcc'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'add')))
                s['t3dten'] = jnp.where(s["_pc"] == 0, s['t3dten'], old_480)
                old_481 = s['qc3dten']
                s['qc3dten'] = s['qc3dten'].at[s['k'] - 1].set(F(_arith(s['qc3dten'][s['k'] - 1], s['pcc'][s['k'] - 1], 'add')))
                s['qc3dten'] = jnp.where(s["_pc"] == 0, s['qc3dten'], old_481)
                return s
            def no_468(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['isatadj'] == 0)), yes_468, no_468, s)
            def yes_482(s):
                s = dict(s)
                def yes_483(s):
                    s = dict(s)
                    old_484 = s['dum']
                    s['dum'] = F(_arith(s['w3d'][s['k'] - 1], s['wvar'][s['k'] - 1], 'add'))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_484)
                    old_485 = s['dum']
                    s['dum'] = F(jnp.maximum(s['dum'], F(0.1)))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_485)
                    return s
                def no_483(s):
                    s = dict(s)
                    def yes_486(s):
                        s = dict(s)
                        old_487 = s['dum']
                        s['dum'] = F(s['w3d'][s['k'] - 1])
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_487)
                        return s
                    def no_486(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['isub'] == 1)), yes_486, no_486, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['isub'] == 0)), yes_483, no_483, s)
                def yes_488(s):
                    s = dict(s)
                    def yes_489(s):
                        s = dict(s)
                        old_490 = s['idrop']
                        s['idrop'] = I(0)
                        s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_490)
                        def yes_491(s):
                            s = dict(s)
                            old_492 = s['idrop']
                            s['idrop'] = I(1)
                            s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_492)
                            return s
                        def no_491(s):
                            s = dict(s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] <= _div(F(5e-05), s['rho'][s['k'] - 1]))), yes_491, no_491, s)
                        def yes_493(s):
                            s = dict(s)
                            old_494 = s['idrop']
                            s['idrop'] = I(1)
                            s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_494)
                            return s
                        def no_493(s):
                            s = dict(s)
                            def yes_495(s):
                                s = dict(s)
                                def yes_496(s):
                                    s = dict(s)
                                    old_497 = s['idrop']
                                    s['idrop'] = I(1)
                                    s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_497)
                                    return s
                                def no_496(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] > _div(F(5e-05), s['rho'][s['k'] - 1]))) & ((s['qc3d'][_arith(s['k'], 1, 'sub') - 1] <= _div(F(5e-05), s['rho'][_arith(s['k'], 1, 'sub') - 1]))))), yes_496, no_496, s)
                                return s
                            def no_495(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['k'] >= 2)), yes_495, no_495, s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['k'] == 1)), yes_493, no_493, s)
                        def yes_498(s):
                            s = dict(s)
                            def yes_499(s):
                                s = dict(s)
                                old_500 = s['dum']
                                s['dum'] = F(_arith(s['dum'], F(100.0), 'mul'))
                                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_500)
                                old_501 = s['dum2']
                                s['dum2'] = F(_arith(_arith(F(0.88), _arith(s['c1'], _div(F(2.0), _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'), _arith(_arith(F(0.07), _arith(s['dum'], F(1.5), 'pow'), 'mul'), _div(s['k1'], _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_501)
                                old_502 = s['dum2']
                                s['dum2'] = F(_arith(s['dum2'], F(1000000.0), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_502)
                                old_503 = s['dum2']
                                s['dum2'] = F(_div(s['dum2'], s['rho'][s['k'] - 1]))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_503)
                                old_504 = s['dum2']
                                s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_504)
                                old_505 = s['dum2']
                                s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_505)
                                old_506 = s['nc3dten']
                                s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_506)
                                old_507 = s['nact']
                                s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_507)
                                return s
                            def no_499(s):
                                s = dict(s)
                                def yes_508(s):
                                    s = dict(s)
                                    old_509 = s['sigvl']
                                    s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                    s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_509)
                                    old_510 = s['aact']
                                    s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                    s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_510)
                                    old_511 = s['alpha']
                                    s['alpha'] = F(_arith(_div(_arith(_arith(s['g'], s['mw'], 'mul'), s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rr'], 'mul'), _arith(s['t3d'][s['k'] - 1], 2, 'pow'), 'mul')), _div(_arith(s['g'], s['ma'], 'mul'), _arith(s['rr'], s['t3d'][s['k'] - 1], 'mul')), 'sub'))
                                    s['alpha'] = jnp.where(s["_pc"] == 0, s['alpha'], old_511)
                                    old_512 = s['gamm']
                                    s['gamm'] = F(_arith(_div(_arith(s['rr'], s['t3d'][s['k'] - 1], 'mul'), _arith(s['evs'][s['k'] - 1], s['mw'], 'mul')), _div(_arith(s['mw'], _arith(s['xxlv'][s['k'] - 1], 2, 'pow'), 'mul'), _arith(_arith(_arith(s['cpm'][s['k'] - 1], s['pres'][s['k'] - 1], 'mul'), s['ma'], 'mul'), s['t3d'][s['k'] - 1], 'mul')), 'add'))
                                    s['gamm'] = jnp.where(s["_pc"] == 0, s['gamm'], old_512)
                                    old_513 = s['gg']
                                    s['gg'] = F(_div(F(1.0), _arith(_div(_arith(_arith(s['rhow'], s['rr'], 'mul'), s['t3d'][s['k'] - 1], 'mul'), _arith(_arith(s['evs'][s['k'] - 1], s['dv'][s['k'] - 1], 'mul'), s['mw'], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['rhow'], 'mul'), _arith(s['kap'][s['k'] - 1], s['t3d'][s['k'] - 1], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['mw'], 'mul'), _arith(s['t3d'][s['k'] - 1], s['rr'], 'mul')), F(1.0), 'sub'), 'mul'), 'add')))
                                    s['gg'] = jnp.where(s["_pc"] == 0, s['gg'], old_513)
                                    old_514 = s['psi']
                                    s['psi'] = F(_arith(_arith(_div(F(2.0), F(3.0)), _arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(0.5), 'pow'), 'mul'), s['aact'], 'mul'))
                                    s['psi'] = jnp.where(s["_pc"] == 0, s['psi'], old_514)
                                    old_515 = s['eta1']
                                    s['eta1'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew1'], 'mul')))
                                    s['eta1'] = jnp.where(s["_pc"] == 0, s['eta1'], old_515)
                                    old_516 = s['eta2']
                                    s['eta2'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew2'], 'mul')))
                                    s['eta2'] = jnp.where(s["_pc"] == 0, s['eta2'], old_516)
                                    old_517 = s['sm1']
                                    s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_517)
                                    old_518 = s['sm2']
                                    s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_518)
                                    old_519 = s['dum1']
                                    s['dum1'] = F(_arith(_div(F(1.0), _arith(s['sm1'], 2, 'pow')), _arith(_arith(s['f11'], _arith(_div(s['psi'], s['eta1']), F(1.5), 'pow'), 'mul'), _arith(s['f21'], _arith(_div(_arith(s['sm1'], 2, 'pow'), _arith(s['eta1'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_519)
                                    old_520 = s['dum2']
                                    s['dum2'] = F(_arith(_div(F(1.0), _arith(s['sm2'], 2, 'pow')), _arith(_arith(s['f12'], _arith(_div(s['psi'], s['eta2']), F(1.5), 'pow'), 'mul'), _arith(s['f22'], _arith(_div(_arith(s['sm2'], 2, 'pow'), _arith(s['eta2'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_520)
                                    old_521 = s['smax']
                                    s['smax'] = F(_div(F(1.0), _arith(_arith(s['dum1'], s['dum2'], 'add'), F(0.5), 'pow')))
                                    s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_521)
                                    old_522 = s['uu1']
                                    s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                    s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_522)
                                    old_523 = s['uu2']
                                    s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                    s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_523)
                                    old_524 = s['dum1']
                                    s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_524)
                                    old_525 = s['dum2']
                                    s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_525)
                                    old_526 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_526)
                                    old_527 = s['dum2']
                                    s['dum2'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_527)
                                    old_528 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_528)
                                    old_529 = s['dum2']
                                    s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_529)
                                    old_530 = s['nc3dten']
                                    s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_530)
                                    old_531 = s['nact']
                                    s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_531)
                                    return s
                                def no_508(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 2)), yes_508, no_508, s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 1)), yes_499, no_499, s)
                            return s
                        def no_498(s):
                            s = dict(s)
                            def yes_532(s):
                                s = dict(s)
                                old_533 = s['tauc']
                                s['tauc'] = F(_div(F(1.0), _div(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), s['nc3d'][s['k'] - 1], 'mul'), _arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), 'mul'), s['lamc'][s['k'] - 1])))
                                s['tauc'] = jnp.where(s["_pc"] == 0, s['tauc'], old_533)
                                def yes_534(s):
                                    s = dict(s)
                                    old_535 = s['taur']
                                    s['taur'] = F(_div(F(1.0), s['epsr']))
                                    s['taur'] = jnp.where(s["_pc"] == 0, s['taur'], old_535)
                                    return s
                                def no_534(s):
                                    s = dict(s)
                                    old_536 = s['taur']
                                    s['taur'] = F(F(100000000.0))
                                    s['taur'] = jnp.where(s["_pc"] == 0, s['taur'], old_536)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['epsr'] > F(1e-08))), yes_534, no_534, s)
                                old_537 = s['dum3']
                                s['dum3'] = F(_arith(_arith(_arith(_div(_arith(s['qvs'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), _arith(s['pres'][s['k'] - 1], s['evs'][s['k'] - 1], 'sub')), _div(s['dqsdt'], s['cp']), 'add'), s['g'], 'mul'), s['dum'], 'mul'))
                                s['dum3'] = jnp.where(s["_pc"] == 0, s['dum3'], old_537)
                                old_538 = s['dum3']
                                s['dum3'] = F(_div(_arith(_arith(s['dum3'], s['tauc'], 'mul'), s['taur'], 'mul'), _arith(s['tauc'], s['taur'], 'add')))
                                s['dum3'] = jnp.where(s["_pc"] == 0, s['dum3'], old_538)
                                def yes_539(s):
                                    s = dict(s)
                                    def yes_540(s):
                                        s = dict(s)
                                        old_541 = s['dum']
                                        s['dum'] = F(_arith(s['dum'], F(100.0), 'mul'))
                                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_541)
                                        old_542 = s['dumact']
                                        s['dumact'] = F(_arith(_arith(F(0.88), _arith(s['c1'], _div(F(2.0), _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'), _arith(_arith(F(0.07), _arith(s['dum'], F(1.5), 'pow'), 'mul'), _div(s['k1'], _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'))
                                        s['dumact'] = jnp.where(s["_pc"] == 0, s['dumact'], old_542)
                                        old_543 = s['dum3']
                                        s['dum3'] = F(_arith(_div(s['dum3'], s['qvs'][s['k'] - 1]), F(100.0), 'mul'))
                                        s['dum3'] = jnp.where(s["_pc"] == 0, s['dum3'], old_543)
                                        old_544 = s['dum2']
                                        s['dum2'] = F(_arith(s['c1'], _arith(s['dum3'], s['k1'], 'pow'), 'mul'))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_544)
                                        old_545 = s['dum2']
                                        s['dum2'] = F(jnp.minimum(s['dum2'], s['dumact']))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_545)
                                        old_546 = s['dum2']
                                        s['dum2'] = F(_arith(s['dum2'], F(1000000.0), 'mul'))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_546)
                                        old_547 = s['dum2']
                                        s['dum2'] = F(_div(s['dum2'], s['rho'][s['k'] - 1]))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_547)
                                        old_548 = s['dum2']
                                        s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_548)
                                        old_549 = s['dum2']
                                        s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_549)
                                        old_550 = s['nc3dten']
                                        s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                        s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_550)
                                        old_551 = s['nact']
                                        s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                        s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_551)
                                        return s
                                    def no_540(s):
                                        s = dict(s)
                                        def yes_552(s):
                                            s = dict(s)
                                            old_553 = s['sigvl']
                                            s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                            s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_553)
                                            old_554 = s['aact']
                                            s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                            s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_554)
                                            old_555 = s['alpha']
                                            s['alpha'] = F(_arith(_div(_arith(_arith(s['g'], s['mw'], 'mul'), s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rr'], 'mul'), _arith(s['t3d'][s['k'] - 1], 2, 'pow'), 'mul')), _div(_arith(s['g'], s['ma'], 'mul'), _arith(s['rr'], s['t3d'][s['k'] - 1], 'mul')), 'sub'))
                                            s['alpha'] = jnp.where(s["_pc"] == 0, s['alpha'], old_555)
                                            old_556 = s['gamm']
                                            s['gamm'] = F(_arith(_div(_arith(s['rr'], s['t3d'][s['k'] - 1], 'mul'), _arith(s['evs'][s['k'] - 1], s['mw'], 'mul')), _div(_arith(s['mw'], _arith(s['xxlv'][s['k'] - 1], 2, 'pow'), 'mul'), _arith(_arith(_arith(s['cpm'][s['k'] - 1], s['pres'][s['k'] - 1], 'mul'), s['ma'], 'mul'), s['t3d'][s['k'] - 1], 'mul')), 'add'))
                                            s['gamm'] = jnp.where(s["_pc"] == 0, s['gamm'], old_556)
                                            old_557 = s['gg']
                                            s['gg'] = F(_div(F(1.0), _arith(_div(_arith(_arith(s['rhow'], s['rr'], 'mul'), s['t3d'][s['k'] - 1], 'mul'), _arith(_arith(s['evs'][s['k'] - 1], s['dv'][s['k'] - 1], 'mul'), s['mw'], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['rhow'], 'mul'), _arith(s['kap'][s['k'] - 1], s['t3d'][s['k'] - 1], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['mw'], 'mul'), _arith(s['t3d'][s['k'] - 1], s['rr'], 'mul')), F(1.0), 'sub'), 'mul'), 'add')))
                                            s['gg'] = jnp.where(s["_pc"] == 0, s['gg'], old_557)
                                            old_558 = s['psi']
                                            s['psi'] = F(_arith(_arith(_div(F(2.0), F(3.0)), _arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(0.5), 'pow'), 'mul'), s['aact'], 'mul'))
                                            s['psi'] = jnp.where(s["_pc"] == 0, s['psi'], old_558)
                                            old_559 = s['eta1']
                                            s['eta1'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew1'], 'mul')))
                                            s['eta1'] = jnp.where(s["_pc"] == 0, s['eta1'], old_559)
                                            old_560 = s['eta2']
                                            s['eta2'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew2'], 'mul')))
                                            s['eta2'] = jnp.where(s["_pc"] == 0, s['eta2'], old_560)
                                            old_561 = s['sm1']
                                            s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_561)
                                            old_562 = s['sm2']
                                            s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_562)
                                            old_563 = s['dum1']
                                            s['dum1'] = F(_arith(_div(F(1.0), _arith(s['sm1'], 2, 'pow')), _arith(_arith(s['f11'], _arith(_div(s['psi'], s['eta1']), F(1.5), 'pow'), 'mul'), _arith(s['f21'], _arith(_div(_arith(s['sm1'], 2, 'pow'), _arith(s['eta1'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                            s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_563)
                                            old_564 = s['dum2']
                                            s['dum2'] = F(_arith(_div(F(1.0), _arith(s['sm2'], 2, 'pow')), _arith(_arith(s['f12'], _arith(_div(s['psi'], s['eta2']), F(1.5), 'pow'), 'mul'), _arith(s['f22'], _arith(_div(_arith(s['sm2'], 2, 'pow'), _arith(s['eta2'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_564)
                                            old_565 = s['smax']
                                            s['smax'] = F(_div(F(1.0), _arith(_arith(s['dum1'], s['dum2'], 'add'), F(0.5), 'pow')))
                                            s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_565)
                                            old_566 = s['uu1']
                                            s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                            s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_566)
                                            old_567 = s['uu2']
                                            s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                            s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_567)
                                            old_568 = s['dum1']
                                            s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                            s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_568)
                                            old_569 = s['dum2']
                                            s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_569)
                                            old_570 = s['dum2']
                                            s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_570)
                                            old_571 = s['dumact']
                                            s['dumact'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                            s['dumact'] = jnp.where(s["_pc"] == 0, s['dumact'], old_571)
                                            old_572 = s['sigvl']
                                            s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                            s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_572)
                                            old_573 = s['aact']
                                            s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                            s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_573)
                                            old_574 = s['sm1']
                                            s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_574)
                                            old_575 = s['sm2']
                                            s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_575)
                                            old_576 = s['smax']
                                            s['smax'] = F(_div(s['dum3'], s['qvs'][s['k'] - 1]))
                                            s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_576)
                                            old_577 = s['uu1']
                                            s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                            s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_577)
                                            old_578 = s['uu2']
                                            s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                            s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_578)
                                            old_579 = s['dum1']
                                            s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                            s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_579)
                                            old_580 = s['dum2']
                                            s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_580)
                                            old_581 = s['dum2']
                                            s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_581)
                                            old_582 = s['dum2']
                                            s['dum2'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_582)
                                            old_583 = s['dum2']
                                            s['dum2'] = F(jnp.minimum(s['dum2'], s['dumact']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_583)
                                            old_584 = s['dum2']
                                            s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_584)
                                            old_585 = s['dum2']
                                            s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_585)
                                            old_586 = s['nc3dten']
                                            s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                            s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_586)
                                            old_587 = s['nact']
                                            s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                            s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_587)
                                            return s
                                        def no_552(s):
                                            s = dict(s)
                                            return s
                                        s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 2)), yes_552, no_552, s)
                                        return s
                                    s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 1)), yes_540, no_540, s)
                                    return s
                                def no_539(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((_div(s['dum3'], s['qvs'][s['k'] - 1]) >= F(1e-06))), yes_539, no_539, s)
                                return s
                            def no_532(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['idrop'] == 0)), yes_532, no_532, s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['idrop'] == 1)), yes_498, no_498, s)
                        return s
                    def no_489(s):
                        s = dict(s)
                        def yes_588(s):
                            s = dict(s)
                            def yes_589(s):
                                s = dict(s)
                                old_590 = s['dum']
                                s['dum'] = F(_arith(s['dum'], F(100.0), 'mul'))
                                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_590)
                                old_591 = s['dum2']
                                s['dum2'] = F(_arith(_arith(F(0.88), _arith(s['c1'], _div(F(2.0), _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'), _arith(_arith(F(0.07), _arith(s['dum'], F(1.5), 'pow'), 'mul'), _div(s['k1'], _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_591)
                                old_592 = s['dum2']
                                s['dum2'] = F(_arith(s['dum2'], F(1000000.0), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_592)
                                old_593 = s['dum2']
                                s['dum2'] = F(_div(s['dum2'], s['rho'][s['k'] - 1]))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_593)
                                old_594 = s['dum2']
                                s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_594)
                                old_595 = s['dum2']
                                s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_595)
                                old_596 = s['nc3dten']
                                s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_596)
                                old_597 = s['nact']
                                s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_597)
                                return s
                            def no_589(s):
                                s = dict(s)
                                def yes_598(s):
                                    s = dict(s)
                                    old_599 = s['sigvl']
                                    s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                    s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_599)
                                    old_600 = s['aact']
                                    s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                    s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_600)
                                    old_601 = s['alpha']
                                    s['alpha'] = F(_arith(_div(_arith(_arith(s['g'], s['mw'], 'mul'), s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rr'], 'mul'), _arith(s['t3d'][s['k'] - 1], 2, 'pow'), 'mul')), _div(_arith(s['g'], s['ma'], 'mul'), _arith(s['rr'], s['t3d'][s['k'] - 1], 'mul')), 'sub'))
                                    s['alpha'] = jnp.where(s["_pc"] == 0, s['alpha'], old_601)
                                    old_602 = s['gamm']
                                    s['gamm'] = F(_arith(_div(_arith(s['rr'], s['t3d'][s['k'] - 1], 'mul'), _arith(s['evs'][s['k'] - 1], s['mw'], 'mul')), _div(_arith(s['mw'], _arith(s['xxlv'][s['k'] - 1], 2, 'pow'), 'mul'), _arith(_arith(_arith(s['cpm'][s['k'] - 1], s['pres'][s['k'] - 1], 'mul'), s['ma'], 'mul'), s['t3d'][s['k'] - 1], 'mul')), 'add'))
                                    s['gamm'] = jnp.where(s["_pc"] == 0, s['gamm'], old_602)
                                    old_603 = s['gg']
                                    s['gg'] = F(_div(F(1.0), _arith(_div(_arith(_arith(s['rhow'], s['rr'], 'mul'), s['t3d'][s['k'] - 1], 'mul'), _arith(_arith(s['evs'][s['k'] - 1], s['dv'][s['k'] - 1], 'mul'), s['mw'], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['rhow'], 'mul'), _arith(s['kap'][s['k'] - 1], s['t3d'][s['k'] - 1], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['mw'], 'mul'), _arith(s['t3d'][s['k'] - 1], s['rr'], 'mul')), F(1.0), 'sub'), 'mul'), 'add')))
                                    s['gg'] = jnp.where(s["_pc"] == 0, s['gg'], old_603)
                                    old_604 = s['psi']
                                    s['psi'] = F(_arith(_arith(_div(F(2.0), F(3.0)), _arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(0.5), 'pow'), 'mul'), s['aact'], 'mul'))
                                    s['psi'] = jnp.where(s["_pc"] == 0, s['psi'], old_604)
                                    old_605 = s['eta1']
                                    s['eta1'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew1'], 'mul')))
                                    s['eta1'] = jnp.where(s["_pc"] == 0, s['eta1'], old_605)
                                    old_606 = s['eta2']
                                    s['eta2'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew2'], 'mul')))
                                    s['eta2'] = jnp.where(s["_pc"] == 0, s['eta2'], old_606)
                                    old_607 = s['sm1']
                                    s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_607)
                                    old_608 = s['sm2']
                                    s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_608)
                                    old_609 = s['dum1']
                                    s['dum1'] = F(_arith(_div(F(1.0), _arith(s['sm1'], 2, 'pow')), _arith(_arith(s['f11'], _arith(_div(s['psi'], s['eta1']), F(1.5), 'pow'), 'mul'), _arith(s['f21'], _arith(_div(_arith(s['sm1'], 2, 'pow'), _arith(s['eta1'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_609)
                                    old_610 = s['dum2']
                                    s['dum2'] = F(_arith(_div(F(1.0), _arith(s['sm2'], 2, 'pow')), _arith(_arith(s['f12'], _arith(_div(s['psi'], s['eta2']), F(1.5), 'pow'), 'mul'), _arith(s['f22'], _arith(_div(_arith(s['sm2'], 2, 'pow'), _arith(s['eta2'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_610)
                                    old_611 = s['smax']
                                    s['smax'] = F(_div(F(1.0), _arith(_arith(s['dum1'], s['dum2'], 'add'), F(0.5), 'pow')))
                                    s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_611)
                                    old_612 = s['uu1']
                                    s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                    s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_612)
                                    old_613 = s['uu2']
                                    s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                    s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_613)
                                    old_614 = s['dum1']
                                    s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_614)
                                    old_615 = s['dum2']
                                    s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_615)
                                    old_616 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_616)
                                    old_617 = s['dum2']
                                    s['dum2'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_617)
                                    old_618 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_618)
                                    old_619 = s['dum2']
                                    s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_619)
                                    old_620 = s['nc3dten']
                                    s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_620)
                                    old_621 = s['nact']
                                    s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_621)
                                    return s
                                def no_598(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 2)), yes_598, no_598, s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 1)), yes_589, no_589, s)
                            return s
                        def no_588(s):
                            s = dict(s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['ibase'] == 2)), yes_588, no_588, s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['ibase'] == 1)), yes_489, no_489, s)
                    return s
                def no_488(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['dum'] >= F(0.001))), yes_488, no_488, s)
                return s
            def no_482(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((_arith(s['qc3d'][s['k'] - 1], _arith(s['qc3dten'][s['k'] - 1], s['dt'], 'mul'), 'add') >= s['qsmall'])) & ((s['inum'] == 0)))), yes_482, no_482, s)
            return s
        def no_203(s):
            s = dict(s)
            def yes_622(s):
                s = dict(s)
                old_623 = s['nc3d']
                s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_arith(s['ndcnst'], F(1000000.0), 'mul'), s['rho'][s['k'] - 1])))
                s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_623)
                return s
            def no_622(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['inum'] == 1)), yes_622, no_622, s)
            def yes_624(s):
                s = dict(s)
                old_625 = s['negfix_ni']
                s['negfix_ni'] = s['negfix_ni'].at[s['k'] - 1].set(F(_arith(s['negfix_ni'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_ni'] = jnp.where(s["_pc"] == 0, s['negfix_ni'], old_625)
                return s
            def no_624(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['ni3d'][s['k'] - 1] < F(0.0))), yes_624, no_624, s)
            def yes_626(s):
                s = dict(s)
                old_627 = s['negfix_ns']
                s['negfix_ns'] = s['negfix_ns'].at[s['k'] - 1].set(F(_arith(s['negfix_ns'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_ns'] = jnp.where(s["_pc"] == 0, s['negfix_ns'], old_627)
                return s
            def no_626(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['ns3d'][s['k'] - 1] < F(0.0))), yes_626, no_626, s)
            def yes_628(s):
                s = dict(s)
                old_629 = s['negfix_nc']
                s['negfix_nc'] = s['negfix_nc'].at[s['k'] - 1].set(F(_arith(s['negfix_nc'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_nc'] = jnp.where(s["_pc"] == 0, s['negfix_nc'], old_629)
                return s
            def no_628(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['nc3d'][s['k'] - 1] < F(0.0))), yes_628, no_628, s)
            def yes_630(s):
                s = dict(s)
                old_631 = s['negfix_nr']
                s['negfix_nr'] = s['negfix_nr'].at[s['k'] - 1].set(F(_arith(s['negfix_nr'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_nr'] = jnp.where(s["_pc"] == 0, s['negfix_nr'], old_631)
                return s
            def no_630(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['nr3d'][s['k'] - 1] < F(0.0))), yes_630, no_630, s)
            def yes_632(s):
                s = dict(s)
                old_633 = s['negfix_ng']
                s['negfix_ng'] = s['negfix_ng'].at[s['k'] - 1].set(F(_arith(s['negfix_ng'][s['k'] - 1], _div(s['ng3d'][s['k'] - 1], s['dt']), 'add')))
                s['negfix_ng'] = jnp.where(s["_pc"] == 0, s['negfix_ng'], old_633)
                return s
            def no_632(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['ng3d'][s['k'] - 1] < F(0.0))), yes_632, no_632, s)
            old_634 = s['ni3d']
            s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ni3d'][s['k'] - 1])))
            s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_634)
            old_635 = s['ns3d']
            s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ns3d'][s['k'] - 1])))
            s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_635)
            old_636 = s['nc3d']
            s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['nc3d'][s['k'] - 1])))
            s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_636)
            old_637 = s['nr3d']
            s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['nr3d'][s['k'] - 1])))
            s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_637)
            old_638 = s['ng3d']
            s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ng3d'][s['k'] - 1])))
            s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_638)
            def yes_639(s):
                s = dict(s)
                old_640 = s['lami']
                s['lami'] = s['lami'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons12'], s['ni3d'][s['k'] - 1], 'mul'), s['qi3d'][s['k'] - 1]), _div(F(1.0), s['di']), 'pow')))
                s['lami'] = jnp.where(s["_pc"] == 0, s['lami'], old_640)
                old_641 = s['n0i']
                s['n0i'] = s['n0i'].at[s['k'] - 1].set(F(_arith(s['ni3d'][s['k'] - 1], s['lami'][s['k'] - 1], 'mul')))
                s['n0i'] = jnp.where(s["_pc"] == 0, s['n0i'], old_641)
                old_642 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_642)
                def yes_643(s):
                    s = dict(s)
                    old_644 = s['lami']
                    s['lami'] = s['lami'].at[s['k'] - 1].set(F(s['lammini']))
                    s['lami'] = jnp.where(s["_pc"] == 0, s['lami'], old_644)
                    old_645 = s['tmpnum']
                    s['tmpnum'] = F(s['ni3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_645)
                    old_646 = s['n0i']
                    s['n0i'] = s['n0i'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lami'][s['k'] - 1], _arith(s['di'], F(1.0), 'add'), 'pow'), s['qi3d'][s['k'] - 1], 'mul'), s['cons12'])))
                    s['n0i'] = jnp.where(s["_pc"] == 0, s['n0i'], old_646)
                    old_647 = s['ni3d']
                    s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_div(s['n0i'][s['k'] - 1], s['lami'][s['k'] - 1])))
                    s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_647)
                    old_648 = s['sizefix_ni']
                    s['sizefix_ni'] = s['sizefix_ni'].at[s['k'] - 1].set(F(_arith(s['sizefix_ni'][s['k'] - 1], _div(_arith(s['ni3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ni'] = jnp.where(s["_pc"] == 0, s['sizefix_ni'], old_648)
                    return s
                def no_643(s):
                    s = dict(s)
                    def yes_649(s):
                        s = dict(s)
                        old_650 = s['lami']
                        s['lami'] = s['lami'].at[s['k'] - 1].set(F(s['lammaxi']))
                        s['lami'] = jnp.where(s["_pc"] == 0, s['lami'], old_650)
                        old_651 = s['n0i']
                        s['n0i'] = s['n0i'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lami'][s['k'] - 1], _arith(s['di'], F(1.0), 'add'), 'pow'), s['qi3d'][s['k'] - 1], 'mul'), s['cons12'])))
                        s['n0i'] = jnp.where(s["_pc"] == 0, s['n0i'], old_651)
                        old_652 = s['tmpnum']
                        s['tmpnum'] = F(s['ni3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_652)
                        old_653 = s['ni3d']
                        s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_div(s['n0i'][s['k'] - 1], s['lami'][s['k'] - 1])))
                        s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_653)
                        old_654 = s['sizefix_ni']
                        s['sizefix_ni'] = s['sizefix_ni'].at[s['k'] - 1].set(F(_arith(s['sizefix_ni'][s['k'] - 1], _div(_arith(s['ni3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_ni'] = jnp.where(s["_pc"] == 0, s['sizefix_ni'], old_654)
                        return s
                    def no_649(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lami'][s['k'] - 1] > s['lammaxi'])), yes_649, no_649, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lami'][s['k'] - 1] < s['lammini'])), yes_643, no_643, s)
                return s
            def no_639(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] >= s['qsmall'])), yes_639, no_639, s)
            def yes_655(s):
                s = dict(s)
                old_656 = s['lamr']
                s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(s['pi'], s['rhow'], 'mul'), s['nr3d'][s['k'] - 1], 'mul'), s['qr3d'][s['k'] - 1]), _div(F(1.0), F(3.0)), 'pow')))
                s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_656)
                old_657 = s['n0rr']
                s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], s['lamr'][s['k'] - 1], 'mul')))
                s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_657)
                old_658 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_658)
                def yes_659(s):
                    s = dict(s)
                    old_660 = s['tmpnum']
                    s['tmpnum'] = F(s['nr3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_660)
                    old_661 = s['lamr']
                    s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(s['lamminr']))
                    s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_661)
                    old_662 = s['n0rr']
                    s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamr'][s['k'] - 1], 4, 'pow'), s['qr3d'][s['k'] - 1], 'mul'), _arith(s['pi'], s['rhow'], 'mul'))))
                    s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_662)
                    old_663 = s['nr3d']
                    s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_div(s['n0rr'][s['k'] - 1], s['lamr'][s['k'] - 1])))
                    s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_663)
                    old_664 = s['sizefix_nr']
                    s['sizefix_nr'] = s['sizefix_nr'].at[s['k'] - 1].set(F(_arith(s['sizefix_nr'][s['k'] - 1], _div(_arith(s['nr3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_664)
                    return s
                def no_659(s):
                    s = dict(s)
                    def yes_665(s):
                        s = dict(s)
                        old_666 = s['tmpnum']
                        s['tmpnum'] = F(s['nr3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_666)
                        old_667 = s['lamr']
                        s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(s['lammaxr']))
                        s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_667)
                        old_668 = s['n0rr']
                        s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamr'][s['k'] - 1], 4, 'pow'), s['qr3d'][s['k'] - 1], 'mul'), _arith(s['pi'], s['rhow'], 'mul'))))
                        s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_668)
                        old_669 = s['nr3d']
                        s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_div(s['n0rr'][s['k'] - 1], s['lamr'][s['k'] - 1])))
                        s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_669)
                        old_670 = s['sizefix_nr']
                        s['sizefix_nr'] = s['sizefix_nr'].at[s['k'] - 1].set(F(_arith(s['sizefix_nr'][s['k'] - 1], _div(_arith(s['nr3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_670)
                        return s
                    def no_665(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lamr'][s['k'] - 1] > s['lammaxr'])), yes_665, no_665, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamr'][s['k'] - 1] < s['lamminr'])), yes_659, no_659, s)
                return s
            def no_655(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])), yes_655, no_655, s)
            def yes_671(s):
                s = dict(s)
                def yes_672(s):
                    s = dict(s)
                    old_673 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(s['pgam_fixed']))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_673)
                    return s
                def no_672(s):
                    s = dict(s)
                    old_674 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_arith(F(0.0005714), _arith(_div(s['nc3d'][s['k'] - 1], F(1000000.0)), s['rho'][s['k'] - 1], 'mul'), 'mul'), F(0.2714), 'add')))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_674)
                    old_675 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_div(F(1.0), _arith(s['pgam'][s['k'] - 1], 2, 'pow')), F(1.0), 'sub')))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_675)
                    old_676 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.maximum(s['pgam'][s['k'] - 1], F(2.0))))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_676)
                    old_677 = s['pgam']
                    s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.minimum(s['pgam'][s['k'] - 1], F(10.0))))
                    s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_677)
                    return s
                s = lax.cond((s["_pc"] == 0) & (s['dofix_pgam']), yes_672, no_672, s)
                old_678 = s['dumii']
                s['dumii'] = I(I(s['pgam'][s['k'] - 1]))
                s['dumii'] = jnp.where(s["_pc"] == 0, s['dumii'], old_678)
                old_679 = s['nu']
                s['nu'] = s['nu'].at[s['k'] - 1].set(F(_arith(s['dnu'][s['dumii'] - 1], _arith(_arith(s['dnu'][_arith(s['dumii'], 1, 'add') - 1], s['dnu'][s['dumii'] - 1], 'sub'), _arith(s['pgam'][s['k'] - 1], F(s['dumii']), 'sub'), 'mul'), 'add')))
                s['nu'] = jnp.where(s["_pc"] == 0, s['nu'], old_679)
                old_680 = s['lamc']
                s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(s['cons26'], s['nc3d'][s['k'] - 1], 'mul'), GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add')), 'mul'), _arith(s['qc3d'][s['k'] - 1], GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add')), 'mul')), _div(F(1.0), F(3.0)), 'pow')))
                s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_680)
                old_681 = s['lammin']
                s['lammin'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(6e-05)))
                s['lammin'] = jnp.where(s["_pc"] == 0, s['lammin'], old_681)
                old_682 = s['lammax']
                s['lammax'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(1e-06)))
                s['lammax'] = jnp.where(s["_pc"] == 0, s['lammax'], old_682)
                old_683 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_683)
                def yes_684(s):
                    s = dict(s)
                    old_685 = s['lamc']
                    s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(s['lammin']))
                    s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_685)
                    old_686 = s['tmpnum']
                    s['tmpnum'] = F(s['nc3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_686)
                    old_687 = s['nc3d']
                    s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_intrinsic('exp', _arith(_arith(_arith(_arith(F(3.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), _intrinsic('log', s['qc3d'][s['k'] - 1]), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'))), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add'))), 'sub')), s['cons26'])))
                    s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_687)
                    old_688 = s['sizefix_nc']
                    s['sizefix_nc'] = s['sizefix_nc'].at[s['k'] - 1].set(F(_arith(s['sizefix_nc'][s['k'] - 1], _div(_arith(s['nc3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_688)
                    return s
                def no_684(s):
                    s = dict(s)
                    def yes_689(s):
                        s = dict(s)
                        old_690 = s['lamc']
                        s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(s['lammax']))
                        s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_690)
                        old_691 = s['tmpnum']
                        s['tmpnum'] = F(s['nc3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_691)
                        old_692 = s['nc3d']
                        s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_intrinsic('exp', _arith(_arith(_arith(_arith(F(3.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), _intrinsic('log', s['qc3d'][s['k'] - 1]), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'))), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add'))), 'sub')), s['cons26'])))
                        s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_692)
                        old_693 = s['sizefix_nc']
                        s['sizefix_nc'] = s['sizefix_nc'].at[s['k'] - 1].set(F(_arith(s['sizefix_nc'][s['k'] - 1], _div(_arith(s['nc3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_693)
                        return s
                    def no_689(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lamc'][s['k'] - 1] > s['lammax'])), yes_689, no_689, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamc'][s['k'] - 1] < s['lammin'])), yes_684, no_684, s)
                old_694 = s['cdist1']
                s['cdist1'] = s['cdist1'].at[s['k'] - 1].set(F(_div(s['nc3d'][s['k'] - 1], GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add')))))
                s['cdist1'] = jnp.where(s["_pc"] == 0, s['cdist1'], old_694)
                return s
            def no_671(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])), yes_671, no_671, s)
            def yes_695(s):
                s = dict(s)
                old_696 = s['lams']
                s['lams'] = s['lams'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons1'], s['ns3d'][s['k'] - 1], 'mul'), s['qni3d'][s['k'] - 1]), _div(F(1.0), s['ds']), 'pow')))
                s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_696)
                old_697 = s['n0s']
                s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_arith(s['ns3d'][s['k'] - 1], s['lams'][s['k'] - 1], 'mul')))
                s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_697)
                old_698 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_698)
                def yes_699(s):
                    s = dict(s)
                    old_700 = s['lams']
                    s['lams'] = s['lams'].at[s['k'] - 1].set(F(s['lammins']))
                    s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_700)
                    old_701 = s['n0s']
                    s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lams'][s['k'] - 1], _arith(s['ds'], F(1.0), 'add'), 'pow'), s['qni3d'][s['k'] - 1], 'mul'), s['cons1'])))
                    s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_701)
                    old_702 = s['tmpnum']
                    s['tmpnum'] = F(s['ns3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_702)
                    old_703 = s['ns3d']
                    s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_div(s['n0s'][s['k'] - 1], s['lams'][s['k'] - 1])))
                    s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_703)
                    old_704 = s['sizefix_ns']
                    s['sizefix_ns'] = s['sizefix_ns'].at[s['k'] - 1].set(F(_arith(s['sizefix_ns'][s['k'] - 1], _div(_arith(s['ns3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_704)
                    return s
                def no_699(s):
                    s = dict(s)
                    def yes_705(s):
                        s = dict(s)
                        old_706 = s['lams']
                        s['lams'] = s['lams'].at[s['k'] - 1].set(F(s['lammaxs']))
                        s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_706)
                        old_707 = s['n0s']
                        s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lams'][s['k'] - 1], _arith(s['ds'], F(1.0), 'add'), 'pow'), s['qni3d'][s['k'] - 1], 'mul'), s['cons1'])))
                        s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_707)
                        old_708 = s['tmpnum']
                        s['tmpnum'] = F(s['ns3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_708)
                        old_709 = s['ns3d']
                        s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_div(s['n0s'][s['k'] - 1], s['lams'][s['k'] - 1])))
                        s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_709)
                        old_710 = s['sizefix_ns']
                        s['sizefix_ns'] = s['sizefix_ns'].at[s['k'] - 1].set(F(_arith(s['sizefix_ns'][s['k'] - 1], _div(_arith(s['ns3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_710)
                        return s
                    def no_705(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lams'][s['k'] - 1] > s['lammaxs'])), yes_705, no_705, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lams'][s['k'] - 1] < s['lammins'])), yes_699, no_699, s)
                return s
            def no_695(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])), yes_695, no_695, s)
            def yes_711(s):
                s = dict(s)
                old_712 = s['lamg']
                s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons2'], s['ng3d'][s['k'] - 1], 'mul'), s['qg3d'][s['k'] - 1]), _div(F(1.0), s['dg']), 'pow')))
                s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_712)
                old_713 = s['n0g']
                s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_arith(s['ng3d'][s['k'] - 1], s['lamg'][s['k'] - 1], 'mul')))
                s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_713)
                old_714 = s['tmpnum']
                s['tmpnum'] = F(F(0.0))
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_714)
                def yes_715(s):
                    s = dict(s)
                    old_716 = s['lamg']
                    s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(s['lamming']))
                    s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_716)
                    old_717 = s['n0g']
                    s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamg'][s['k'] - 1], _arith(s['dg'], F(1.0), 'add'), 'pow'), s['qg3d'][s['k'] - 1], 'mul'), s['cons2'])))
                    s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_717)
                    old_718 = s['tmpnum']
                    s['tmpnum'] = F(s['ng3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_718)
                    old_719 = s['ng3d']
                    s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_div(s['n0g'][s['k'] - 1], s['lamg'][s['k'] - 1])))
                    s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_719)
                    old_720 = s['sizefix_ng']
                    s['sizefix_ng'] = s['sizefix_ng'].at[s['k'] - 1].set(F(_arith(s['sizefix_ng'][s['k'] - 1], _div(_arith(s['ng3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_720)
                    return s
                def no_715(s):
                    s = dict(s)
                    def yes_721(s):
                        s = dict(s)
                        old_722 = s['lamg']
                        s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(s['lammaxg']))
                        s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_722)
                        old_723 = s['n0g']
                        s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamg'][s['k'] - 1], _arith(s['dg'], F(1.0), 'add'), 'pow'), s['qg3d'][s['k'] - 1], 'mul'), s['cons2'])))
                        s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_723)
                        old_724 = s['tmpnum']
                        s['tmpnum'] = F(s['ng3d'][s['k'] - 1])
                        s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_724)
                        old_725 = s['ng3d']
                        s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_div(s['n0g'][s['k'] - 1], s['lamg'][s['k'] - 1])))
                        s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_725)
                        old_726 = s['sizefix_ng']
                        s['sizefix_ng'] = s['sizefix_ng'].at[s['k'] - 1].set(F(_arith(s['sizefix_ng'][s['k'] - 1], _div(_arith(s['ng3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                        s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_726)
                        return s
                    def no_721(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['lamg'][s['k'] - 1] > s['lammaxg'])), yes_721, no_721, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamg'][s['k'] - 1] < s['lamming'])), yes_715, no_715, s)
                return s
            def no_711(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] >= s['qsmall'])), yes_711, no_711, s)
            old_727 = s['mnuccc']
            s['mnuccc'] = s['mnuccc'].at[s['k'] - 1].set(F(F(0.0)))
            s['mnuccc'] = jnp.where(s["_pc"] == 0, s['mnuccc'], old_727)
            old_728 = s['nnuccc']
            s['nnuccc'] = s['nnuccc'].at[s['k'] - 1].set(F(F(0.0)))
            s['nnuccc'] = jnp.where(s["_pc"] == 0, s['nnuccc'], old_728)
            old_729 = s['prc']
            s['prc'] = s['prc'].at[s['k'] - 1].set(F(F(0.0)))
            s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_729)
            old_730 = s['nprc']
            s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(F(0.0)))
            s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_730)
            old_731 = s['nprc1']
            s['nprc1'] = s['nprc1'].at[s['k'] - 1].set(F(F(0.0)))
            s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_731)
            old_732 = s['nsagg']
            s['nsagg'] = s['nsagg'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsagg'] = jnp.where(s["_pc"] == 0, s['nsagg'], old_732)
            old_733 = s['psacws']
            s['psacws'] = s['psacws'].at[s['k'] - 1].set(F(F(0.0)))
            s['psacws'] = jnp.where(s["_pc"] == 0, s['psacws'], old_733)
            old_734 = s['npsacws']
            s['npsacws'] = s['npsacws'].at[s['k'] - 1].set(F(F(0.0)))
            s['npsacws'] = jnp.where(s["_pc"] == 0, s['npsacws'], old_734)
            old_735 = s['psacwi']
            s['psacwi'] = s['psacwi'].at[s['k'] - 1].set(F(F(0.0)))
            s['psacwi'] = jnp.where(s["_pc"] == 0, s['psacwi'], old_735)
            old_736 = s['npsacwi']
            s['npsacwi'] = s['npsacwi'].at[s['k'] - 1].set(F(F(0.0)))
            s['npsacwi'] = jnp.where(s["_pc"] == 0, s['npsacwi'], old_736)
            old_737 = s['pracs']
            s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(F(0.0)))
            s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_737)
            old_738 = s['npracs']
            s['npracs'] = s['npracs'].at[s['k'] - 1].set(F(F(0.0)))
            s['npracs'] = jnp.where(s["_pc"] == 0, s['npracs'], old_738)
            old_739 = s['nmults']
            s['nmults'] = s['nmults'].at[s['k'] - 1].set(F(F(0.0)))
            s['nmults'] = jnp.where(s["_pc"] == 0, s['nmults'], old_739)
            old_740 = s['qmults']
            s['qmults'] = s['qmults'].at[s['k'] - 1].set(F(F(0.0)))
            s['qmults'] = jnp.where(s["_pc"] == 0, s['qmults'], old_740)
            old_741 = s['nmultr']
            s['nmultr'] = s['nmultr'].at[s['k'] - 1].set(F(F(0.0)))
            s['nmultr'] = jnp.where(s["_pc"] == 0, s['nmultr'], old_741)
            old_742 = s['qmultr']
            s['qmultr'] = s['qmultr'].at[s['k'] - 1].set(F(F(0.0)))
            s['qmultr'] = jnp.where(s["_pc"] == 0, s['qmultr'], old_742)
            old_743 = s['nmultg']
            s['nmultg'] = s['nmultg'].at[s['k'] - 1].set(F(F(0.0)))
            s['nmultg'] = jnp.where(s["_pc"] == 0, s['nmultg'], old_743)
            old_744 = s['qmultg']
            s['qmultg'] = s['qmultg'].at[s['k'] - 1].set(F(F(0.0)))
            s['qmultg'] = jnp.where(s["_pc"] == 0, s['qmultg'], old_744)
            old_745 = s['nmultrg']
            s['nmultrg'] = s['nmultrg'].at[s['k'] - 1].set(F(F(0.0)))
            s['nmultrg'] = jnp.where(s["_pc"] == 0, s['nmultrg'], old_745)
            old_746 = s['qmultrg']
            s['qmultrg'] = s['qmultrg'].at[s['k'] - 1].set(F(F(0.0)))
            s['qmultrg'] = jnp.where(s["_pc"] == 0, s['qmultrg'], old_746)
            old_747 = s['mnuccr']
            s['mnuccr'] = s['mnuccr'].at[s['k'] - 1].set(F(F(0.0)))
            s['mnuccr'] = jnp.where(s["_pc"] == 0, s['mnuccr'], old_747)
            old_748 = s['nnuccr']
            s['nnuccr'] = s['nnuccr'].at[s['k'] - 1].set(F(F(0.0)))
            s['nnuccr'] = jnp.where(s["_pc"] == 0, s['nnuccr'], old_748)
            old_749 = s['pra']
            s['pra'] = s['pra'].at[s['k'] - 1].set(F(F(0.0)))
            s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_749)
            old_750 = s['npra']
            s['npra'] = s['npra'].at[s['k'] - 1].set(F(F(0.0)))
            s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_750)
            old_751 = s['nragg']
            s['nragg'] = s['nragg'].at[s['k'] - 1].set(F(F(0.0)))
            s['nragg'] = jnp.where(s["_pc"] == 0, s['nragg'], old_751)
            old_752 = s['prci']
            s['prci'] = s['prci'].at[s['k'] - 1].set(F(F(0.0)))
            s['prci'] = jnp.where(s["_pc"] == 0, s['prci'], old_752)
            old_753 = s['nprci']
            s['nprci'] = s['nprci'].at[s['k'] - 1].set(F(F(0.0)))
            s['nprci'] = jnp.where(s["_pc"] == 0, s['nprci'], old_753)
            old_754 = s['prai']
            s['prai'] = s['prai'].at[s['k'] - 1].set(F(F(0.0)))
            s['prai'] = jnp.where(s["_pc"] == 0, s['prai'], old_754)
            old_755 = s['nprai']
            s['nprai'] = s['nprai'].at[s['k'] - 1].set(F(F(0.0)))
            s['nprai'] = jnp.where(s["_pc"] == 0, s['nprai'], old_755)
            old_756 = s['nnuccd']
            s['nnuccd'] = s['nnuccd'].at[s['k'] - 1].set(F(F(0.0)))
            s['nnuccd'] = jnp.where(s["_pc"] == 0, s['nnuccd'], old_756)
            old_757 = s['mnuccd']
            s['mnuccd'] = s['mnuccd'].at[s['k'] - 1].set(F(F(0.0)))
            s['mnuccd'] = jnp.where(s["_pc"] == 0, s['mnuccd'], old_757)
            old_758 = s['pcc']
            s['pcc'] = s['pcc'].at[s['k'] - 1].set(F(F(0.0)))
            s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_758)
            old_759 = s['pre']
            s['pre'] = s['pre'].at[s['k'] - 1].set(F(F(0.0)))
            s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_759)
            old_760 = s['prd']
            s['prd'] = s['prd'].at[s['k'] - 1].set(F(F(0.0)))
            s['prd'] = jnp.where(s["_pc"] == 0, s['prd'], old_760)
            old_761 = s['prds']
            s['prds'] = s['prds'].at[s['k'] - 1].set(F(F(0.0)))
            s['prds'] = jnp.where(s["_pc"] == 0, s['prds'], old_761)
            old_762 = s['eprd']
            s['eprd'] = s['eprd'].at[s['k'] - 1].set(F(F(0.0)))
            s['eprd'] = jnp.where(s["_pc"] == 0, s['eprd'], old_762)
            old_763 = s['eprds']
            s['eprds'] = s['eprds'].at[s['k'] - 1].set(F(F(0.0)))
            s['eprds'] = jnp.where(s["_pc"] == 0, s['eprds'], old_763)
            old_764 = s['nsubc']
            s['nsubc'] = s['nsubc'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsubc'] = jnp.where(s["_pc"] == 0, s['nsubc'], old_764)
            old_765 = s['nsubi']
            s['nsubi'] = s['nsubi'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsubi'] = jnp.where(s["_pc"] == 0, s['nsubi'], old_765)
            old_766 = s['nsubs']
            s['nsubs'] = s['nsubs'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsubs'] = jnp.where(s["_pc"] == 0, s['nsubs'], old_766)
            old_767 = s['nsubr']
            s['nsubr'] = s['nsubr'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsubr'] = jnp.where(s["_pc"] == 0, s['nsubr'], old_767)
            old_768 = s['piacr']
            s['piacr'] = s['piacr'].at[s['k'] - 1].set(F(F(0.0)))
            s['piacr'] = jnp.where(s["_pc"] == 0, s['piacr'], old_768)
            old_769 = s['niacr']
            s['niacr'] = s['niacr'].at[s['k'] - 1].set(F(F(0.0)))
            s['niacr'] = jnp.where(s["_pc"] == 0, s['niacr'], old_769)
            old_770 = s['praci']
            s['praci'] = s['praci'].at[s['k'] - 1].set(F(F(0.0)))
            s['praci'] = jnp.where(s["_pc"] == 0, s['praci'], old_770)
            old_771 = s['piacrs']
            s['piacrs'] = s['piacrs'].at[s['k'] - 1].set(F(F(0.0)))
            s['piacrs'] = jnp.where(s["_pc"] == 0, s['piacrs'], old_771)
            old_772 = s['niacrs']
            s['niacrs'] = s['niacrs'].at[s['k'] - 1].set(F(F(0.0)))
            s['niacrs'] = jnp.where(s["_pc"] == 0, s['niacrs'], old_772)
            old_773 = s['pracis']
            s['pracis'] = s['pracis'].at[s['k'] - 1].set(F(F(0.0)))
            s['pracis'] = jnp.where(s["_pc"] == 0, s['pracis'], old_773)
            old_774 = s['pracg']
            s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(F(0.0)))
            s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_774)
            old_775 = s['psacr']
            s['psacr'] = s['psacr'].at[s['k'] - 1].set(F(F(0.0)))
            s['psacr'] = jnp.where(s["_pc"] == 0, s['psacr'], old_775)
            old_776 = s['psacwg']
            s['psacwg'] = s['psacwg'].at[s['k'] - 1].set(F(F(0.0)))
            s['psacwg'] = jnp.where(s["_pc"] == 0, s['psacwg'], old_776)
            old_777 = s['pgsacw']
            s['pgsacw'] = s['pgsacw'].at[s['k'] - 1].set(F(F(0.0)))
            s['pgsacw'] = jnp.where(s["_pc"] == 0, s['pgsacw'], old_777)
            old_778 = s['pgracs']
            s['pgracs'] = s['pgracs'].at[s['k'] - 1].set(F(F(0.0)))
            s['pgracs'] = jnp.where(s["_pc"] == 0, s['pgracs'], old_778)
            old_779 = s['prdg']
            s['prdg'] = s['prdg'].at[s['k'] - 1].set(F(F(0.0)))
            s['prdg'] = jnp.where(s["_pc"] == 0, s['prdg'], old_779)
            old_780 = s['eprdg']
            s['eprdg'] = s['eprdg'].at[s['k'] - 1].set(F(F(0.0)))
            s['eprdg'] = jnp.where(s["_pc"] == 0, s['eprdg'], old_780)
            old_781 = s['npracg']
            s['npracg'] = s['npracg'].at[s['k'] - 1].set(F(F(0.0)))
            s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_781)
            old_782 = s['npsacwg']
            s['npsacwg'] = s['npsacwg'].at[s['k'] - 1].set(F(F(0.0)))
            s['npsacwg'] = jnp.where(s["_pc"] == 0, s['npsacwg'], old_782)
            old_783 = s['nscng']
            s['nscng'] = s['nscng'].at[s['k'] - 1].set(F(F(0.0)))
            s['nscng'] = jnp.where(s["_pc"] == 0, s['nscng'], old_783)
            old_784 = s['ngracs']
            s['ngracs'] = s['ngracs'].at[s['k'] - 1].set(F(F(0.0)))
            s['ngracs'] = jnp.where(s["_pc"] == 0, s['ngracs'], old_784)
            old_785 = s['nsubg']
            s['nsubg'] = s['nsubg'].at[s['k'] - 1].set(F(F(0.0)))
            s['nsubg'] = jnp.where(s["_pc"] == 0, s['nsubg'], old_785)
            def yes_786(s):
                s = dict(s)
                old_787 = s['nacnt']
                s['nacnt'] = F(_arith(_intrinsic('exp', _arith((-F(2.8)), _arith(F(0.262), _arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), 'mul'), 'add')), F(1000.0), 'mul'))
                s['nacnt'] = jnp.where(s["_pc"] == 0, s['nacnt'], old_787)
                old_788 = s['dum']
                s['dum'] = F(_div(_div(_arith(F(7.37), s['t3d'][s['k'] - 1], 'mul'), _arith(_arith(F(288.0), F(10.0), 'mul'), s['pres'][s['k'] - 1], 'mul')), F(100.0)))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_788)
                old_789 = s['dap']
                s['dap'] = s['dap'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['cons37'], s['t3d'][s['k'] - 1], 'mul'), _arith(F(1.0), _div(s['dum'], s['rin']), 'add'), 'mul'), s['mu'][s['k'] - 1])))
                s['dap'] = jnp.where(s["_pc"] == 0, s['dap'], old_789)
                old_790 = s['mnuccc']
                s['mnuccc'] = s['mnuccc'].at[s['k'] - 1].set(F(_arith(_arith(_arith(s['cons38'], s['dap'][s['k'] - 1], 'mul'), s['nacnt'], 'mul'), _intrinsic('exp', _arith(_arith(_intrinsic('log', s['cdist1'][s['k'] - 1]), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(5.0), 'add'))), 'add'), _arith(F(4.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), 'sub')), 'mul')))
                s['mnuccc'] = jnp.where(s["_pc"] == 0, s['mnuccc'], old_790)
                old_791 = s['nnuccc']
                s['nnuccc'] = s['nnuccc'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['dap'][s['k'] - 1], 'mul'), s['nacnt'], 'mul'), s['cdist1'][s['k'] - 1], 'mul'), GAMMA(_arith(s['pgam'][s['k'] - 1], F(2.0), 'add')), 'mul'), s['lamc'][s['k'] - 1])))
                s['nnuccc'] = jnp.where(s["_pc"] == 0, s['nnuccc'], old_791)
                old_792 = s['mnuccc']
                s['mnuccc'] = s['mnuccc'].at[s['k'] - 1].set(F(_arith(s['mnuccc'][s['k'] - 1], _arith(_arith(s['cons39'], _intrinsic('exp', _arith(_arith(_intrinsic('log', s['cdist1'][s['k'] - 1]), _intrinsic('log', GAMMA(_arith(F(7.0), s['pgam'][s['k'] - 1], 'add'))), 'add'), _arith(F(6.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), 'sub')), 'mul'), _intrinsic('exp', _arith(s['aimm'], _arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), 'mul')), 'mul'), 'add')))
                s['mnuccc'] = jnp.where(s["_pc"] == 0, s['mnuccc'], old_792)
                old_793 = s['nnuccc']
                s['nnuccc'] = s['nnuccc'].at[s['k'] - 1].set(F(_arith(s['nnuccc'][s['k'] - 1], _arith(_arith(s['cons40'], _intrinsic('exp', _arith(_arith(_intrinsic('log', s['cdist1'][s['k'] - 1]), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add'))), 'add'), _arith(F(3.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), 'sub')), 'mul'), _intrinsic('exp', _arith(s['aimm'], _arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), 'mul')), 'mul'), 'add')))
                s['nnuccc'] = jnp.where(s["_pc"] == 0, s['nnuccc'], old_793)
                old_794 = s['nnuccc']
                s['nnuccc'] = s['nnuccc'].at[s['k'] - 1].set(F(jnp.minimum(s['nnuccc'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']))))
                s['nnuccc'] = jnp.where(s["_pc"] == 0, s['nnuccc'], old_794)
                return s
            def no_786(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] >= s['qsmall'])) & ((s['t3d'][s['k'] - 1] < F(269.15))))), yes_786, no_786, s)
            old_795 = s['nnuccc']
            s['nnuccc'] = s['nnuccc'].at[s['k'] - 1].set(F(_arith(s['nnuccc'][s['k'] - 1], s['nnuccc_reduce_coef'], 'mul')))
            s['nnuccc'] = jnp.where(s["_pc"] == 0, s['nnuccc'], old_795)
            def yes_796(s):
                s = dict(s)
                def yes_797(s):
                    s = dict(s)
                    old_798 = s['prc']
                    s['prc'] = s['prc'].at[s['k'] - 1].set(F(_arith(_arith(F(1350.0), _arith(s['qc3d'][s['k'] - 1], F(2.47), 'pow'), 'mul'), _arith(_arith(_div(s['nc3d'][s['k'] - 1], F(1000000.0)), s['rho'][s['k'] - 1], 'mul'), (-F(1.79)), 'pow'), 'mul')))
                    s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_798)
                    old_799 = s['nprc1']
                    s['nprc1'] = s['nprc1'].at[s['k'] - 1].set(F(_div(s['prc'][s['k'] - 1], s['cons29'])))
                    s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_799)
                    old_800 = s['nprc']
                    s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(_div(s['prc'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['nc3d'][s['k'] - 1]))))
                    s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_800)
                    old_801 = s['nprc']
                    s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(jnp.minimum(s['nprc'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']))))
                    s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_801)
                    return s
                def no_797(s):
                    s = dict(s)
                    def yes_802(s):
                        s = dict(s)
                        old_803 = s['dum']
                        s['dum'] = F(_arith(F(1.0), _div(s['qc3d'][s['k'] - 1], _arith(s['qc3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')), 'sub'))
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_803)
                        old_804 = s['dum1']
                        s['dum1'] = F(_arith(_arith(F(600.0), _arith(s['dum'], F(0.68), 'pow'), 'mul'), _arith(_arith(F(1.0), _arith(s['dum'], F(0.68), 'pow'), 'sub'), 3, 'pow'), 'mul'))
                        s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_804)
                        old_805 = s['prc']
                        s['prc'] = s['prc'].at[s['k'] - 1].set(F(_div(_arith(_arith(_div(_arith(_div(_arith(_arith(_div(F(9440000000.0), _arith(F(20.0), F(2.6e-07), 'mul')), _arith(s['nu'][s['k'] - 1], F(2.0), 'add'), 'mul'), _arith(s['nu'][s['k'] - 1], F(4.0), 'add'), 'mul'), _arith(_arith(s['nu'][s['k'] - 1], F(1.0), 'add'), 2, 'pow')), _arith(_div(_arith(s['rho'][s['k'] - 1], s['qc3d'][s['k'] - 1], 'mul'), F(1000.0)), 4, 'pow'), 'mul'), _arith(_div(_arith(s['rho'][s['k'] - 1], s['nc3d'][s['k'] - 1], 'mul'), F(1000000.0)), 2, 'pow')), _arith(F(1.0), _div(s['dum1'], _arith(_arith(F(1.0), s['dum'], 'sub'), 2, 'pow')), 'add'), 'mul'), F(1000.0), 'mul'), s['rho'][s['k'] - 1])))
                        s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_805)
                        old_806 = s['nprc']
                        s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['prc'][s['k'] - 1], F(2.0), 'mul'), F(2.6e-07)), F(1000.0), 'mul')))
                        s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_806)
                        old_807 = s['nprc1']
                        s['nprc1'] = s['nprc1'].at[s['k'] - 1].set(F(_arith(F(0.5), s['nprc'][s['k'] - 1], 'mul')))
                        s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_807)
                        return s
                    def no_802(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 1)), yes_802, no_802, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 0)), yes_797, no_797, s)
                return s
            def no_796(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] >= F(1e-06))), yes_796, no_796, s)
            def yes_808(s):
                s = dict(s)
                old_809 = s['nsagg']
                s['nsagg'] = s['nsagg'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons15'], s['asn'][s['k'] - 1], 'mul'), _arith(s['rho'][s['k'] - 1], _div(_arith(F(2.0), s['bs'], 'add'), F(3.0)), 'pow'), 'mul'), _arith(s['qni3d'][s['k'] - 1], _div(_arith(F(2.0), s['bs'], 'add'), F(3.0)), 'pow'), 'mul'), _arith(_arith(s['ns3d'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), _div(_arith(F(4.0), s['bs'], 'sub'), F(3.0)), 'pow'), 'mul'), s['rho'][s['k'] - 1])))
                s['nsagg'] = jnp.where(s["_pc"] == 0, s['nsagg'], old_809)
                return s
            def no_808(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= F(1e-08))), yes_808, no_808, s)
            def yes_810(s):
                s = dict(s)
                old_811 = s['psacws']
                s['psacws'] = s['psacws'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons13'], s['asn'][s['k'] - 1], 'mul'), s['qc3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(s['lams'][s['k'] - 1], _arith(s['bs'], F(3.0), 'add'), 'pow'))))
                s['psacws'] = jnp.where(s["_pc"] == 0, s['psacws'], old_811)
                old_812 = s['npsacws']
                s['npsacws'] = s['npsacws'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons13'], s['asn'][s['k'] - 1], 'mul'), s['nc3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(s['lams'][s['k'] - 1], _arith(s['bs'], F(3.0), 'add'), 'pow'))))
                s['npsacws'] = jnp.where(s["_pc"] == 0, s['npsacws'], old_812)
                return s
            def no_810(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qni3d'][s['k'] - 1] >= F(1e-08))) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])))), yes_810, no_810, s)
            def yes_813(s):
                s = dict(s)
                old_814 = s['psacwg']
                s['psacwg'] = s['psacwg'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons14'], s['agn'][s['k'] - 1], 'mul'), s['qc3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0g'][s['k'] - 1], 'mul'), _arith(s['lamg'][s['k'] - 1], _arith(s['bg'], F(3.0), 'add'), 'pow'))))
                s['psacwg'] = jnp.where(s["_pc"] == 0, s['psacwg'], old_814)
                old_815 = s['npsacwg']
                s['npsacwg'] = s['npsacwg'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons14'], s['agn'][s['k'] - 1], 'mul'), s['nc3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0g'][s['k'] - 1], 'mul'), _arith(s['lamg'][s['k'] - 1], _arith(s['bg'], F(3.0), 'add'), 'pow'))))
                s['npsacwg'] = jnp.where(s["_pc"] == 0, s['npsacwg'], old_815)
                return s
            def no_813(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qg3d'][s['k'] - 1] >= F(1e-08))) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])))), yes_813, no_813, s)
            def yes_816(s):
                s = dict(s)
                def yes_817(s):
                    s = dict(s)
                    old_818 = s['psacwi']
                    s['psacwi'] = s['psacwi'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons16'], s['ain'][s['k'] - 1], 'mul'), s['qc3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0i'][s['k'] - 1], 'mul'), _arith(s['lami'][s['k'] - 1], _arith(s['bi'], F(3.0), 'add'), 'pow'))))
                    s['psacwi'] = jnp.where(s["_pc"] == 0, s['psacwi'], old_818)
                    old_819 = s['npsacwi']
                    s['npsacwi'] = s['npsacwi'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons16'], s['ain'][s['k'] - 1], 'mul'), s['nc3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0i'][s['k'] - 1], 'mul'), _arith(s['lami'][s['k'] - 1], _arith(s['bi'], F(3.0), 'add'), 'pow'))))
                    s['npsacwi'] = jnp.where(s["_pc"] == 0, s['npsacwi'], old_819)
                    return s
                def no_817(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((_div(F(1.0), s['lami'][s['k'] - 1]) >= F(0.0001))), yes_817, no_817, s)
                return s
            def no_816(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qi3d'][s['k'] - 1] >= F(1e-08))) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])))), yes_816, no_816, s)
            def yes_820(s):
                s = dict(s)
                old_821 = s['ums']
                s['ums'] = F(_div(_arith(s['asn'][s['k'] - 1], s['cons3'], 'mul'), _arith(s['lams'][s['k'] - 1], s['bs'], 'pow')))
                s['ums'] = jnp.where(s["_pc"] == 0, s['ums'], old_821)
                old_822 = s['umr']
                s['umr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons4'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_822)
                old_823 = s['uns']
                s['uns'] = F(_div(_arith(s['asn'][s['k'] - 1], s['cons5'], 'mul'), _arith(s['lams'][s['k'] - 1], s['bs'], 'pow')))
                s['uns'] = jnp.where(s["_pc"] == 0, s['uns'], old_823)
                old_824 = s['unr']
                s['unr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons6'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_824)
                old_825 = s['dum']
                s['dum'] = F(_arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.54), 'pow'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_825)
                old_826 = s['ums']
                s['ums'] = F(jnp.minimum(s['ums'], _arith(F(1.2), s['dum'], 'mul')))
                s['ums'] = jnp.where(s["_pc"] == 0, s['ums'], old_826)
                old_827 = s['uns']
                s['uns'] = F(jnp.minimum(s['uns'], _arith(F(1.2), s['dum'], 'mul')))
                s['uns'] = jnp.where(s["_pc"] == 0, s['uns'], old_827)
                old_828 = s['umr']
                s['umr'] = F(jnp.minimum(s['umr'], _arith(F(9.1), s['dum'], 'mul')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_828)
                old_829 = s['unr']
                s['unr'] = F(jnp.minimum(s['unr'], _arith(F(9.1), s['dum'], 'mul')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_829)
                old_830 = s['pracs']
                s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(_arith(s['cons41'], _arith(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(F(1.2), s['umr'], 'mul'), _arith(F(0.95), s['ums'], 'mul'), 'sub'), 2, 'pow'), _arith(_arith(F(0.08), s['ums'], 'mul'), s['umr'], 'mul'), 'add'), F(0.5), 'pow'), s['rho'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], 3, 'pow')), _arith(_arith(_div(F(5.0), _arith(_arith(s['lamr'][s['k'] - 1], 3, 'pow'), s['lams'][s['k'] - 1], 'mul')), _div(F(2.0), _arith(_arith(s['lamr'][s['k'] - 1], 2, 'pow'), _arith(s['lams'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(0.5), _arith(s['lamr'][s['k'] - 1], _arith(s['lams'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul'), 'mul')))
                s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_830)
                old_831 = s['npracs']
                s['npracs'] = s['npracs'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_arith(_arith(s['cons32'], s['rho'][s['k'] - 1], 'mul'), _arith(_arith(_arith(F(1.7), _arith(_arith(s['unr'], s['uns'], 'sub'), 2, 'pow'), 'mul'), _arith(_arith(F(0.3), s['unr'], 'mul'), s['uns'], 'mul'), 'add'), F(0.5), 'pow'), 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(_arith(_div(F(1.0), _arith(_arith(s['lamr'][s['k'] - 1], 3, 'pow'), s['lams'][s['k'] - 1], 'mul')), _div(F(1.0), _arith(_arith(s['lamr'][s['k'] - 1], 2, 'pow'), _arith(s['lams'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(1.0), _arith(s['lamr'][s['k'] - 1], _arith(s['lams'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul')))
                s['npracs'] = jnp.where(s["_pc"] == 0, s['npracs'], old_831)
                old_832 = s['pracs']
                s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(jnp.minimum(s['pracs'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']))))
                s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_832)
                def yes_833(s):
                    s = dict(s)
                    old_834 = s['psacr']
                    s['psacr'] = s['psacr'].at[s['k'] - 1].set(F(_arith(s['cons31'], _arith(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(F(1.2), s['umr'], 'mul'), _arith(F(0.95), s['ums'], 'mul'), 'sub'), 2, 'pow'), _arith(_arith(F(0.08), s['ums'], 'mul'), s['umr'], 'mul'), 'add'), F(0.5), 'pow'), s['rho'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(s['lams'][s['k'] - 1], 3, 'pow')), _arith(_arith(_div(F(5.0), _arith(_arith(s['lams'][s['k'] - 1], 3, 'pow'), s['lamr'][s['k'] - 1], 'mul')), _div(F(2.0), _arith(_arith(s['lams'][s['k'] - 1], 2, 'pow'), _arith(s['lamr'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(0.5), _arith(s['lams'][s['k'] - 1], _arith(s['lamr'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul'), 'mul')))
                    s['psacr'] = jnp.where(s["_pc"] == 0, s['psacr'], old_834)
                    return s
                def no_833(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((s['qni3d'][s['k'] - 1] >= F(0.0001))) & ((s['qr3d'][s['k'] - 1] >= F(0.0001))))), yes_833, no_833, s)
                return s
            def no_820(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qr3d'][s['k'] - 1] >= F(1e-08))) & ((s['qni3d'][s['k'] - 1] >= F(1e-08))))), yes_820, no_820, s)
            def yes_835(s):
                s = dict(s)
                old_836 = s['umg']
                s['umg'] = F(_div(_arith(s['agn'][s['k'] - 1], s['cons7'], 'mul'), _arith(s['lamg'][s['k'] - 1], s['bg'], 'pow')))
                s['umg'] = jnp.where(s["_pc"] == 0, s['umg'], old_836)
                old_837 = s['umr']
                s['umr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons4'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_837)
                old_838 = s['ung']
                s['ung'] = F(_div(_arith(s['agn'][s['k'] - 1], s['cons8'], 'mul'), _arith(s['lamg'][s['k'] - 1], s['bg'], 'pow')))
                s['ung'] = jnp.where(s["_pc"] == 0, s['ung'], old_838)
                old_839 = s['unr']
                s['unr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons6'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['br'], 'pow')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_839)
                old_840 = s['dum']
                s['dum'] = F(_arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.54), 'pow'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_840)
                old_841 = s['umg']
                s['umg'] = F(jnp.minimum(s['umg'], _arith(F(20.0), s['dum'], 'mul')))
                s['umg'] = jnp.where(s["_pc"] == 0, s['umg'], old_841)
                old_842 = s['ung']
                s['ung'] = F(jnp.minimum(s['ung'], _arith(F(20.0), s['dum'], 'mul')))
                s['ung'] = jnp.where(s["_pc"] == 0, s['ung'], old_842)
                old_843 = s['umr']
                s['umr'] = F(jnp.minimum(s['umr'], _arith(F(9.1), s['dum'], 'mul')))
                s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_843)
                old_844 = s['unr']
                s['unr'] = F(jnp.minimum(s['unr'], _arith(F(9.1), s['dum'], 'mul')))
                s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_844)
                old_845 = s['pracg']
                s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(_arith(s['cons41'], _arith(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(F(1.2), s['umr'], 'mul'), _arith(F(0.95), s['umg'], 'mul'), 'sub'), 2, 'pow'), _arith(_arith(F(0.08), s['umg'], 'mul'), s['umr'], 'mul'), 'add'), F(0.5), 'pow'), s['rho'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0g'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], 3, 'pow')), _arith(_arith(_div(F(5.0), _arith(_arith(s['lamr'][s['k'] - 1], 3, 'pow'), s['lamg'][s['k'] - 1], 'mul')), _div(F(2.0), _arith(_arith(s['lamr'][s['k'] - 1], 2, 'pow'), _arith(s['lamg'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(0.5), _arith(s['lamr'][s['k'] - 1], _arith(s['lamg'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul'), 'mul')))
                s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_845)
                old_846 = s['npracg']
                s['npracg'] = s['npracg'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_arith(_arith(s['cons32'], s['rho'][s['k'] - 1], 'mul'), _arith(_arith(_arith(F(1.7), _arith(_arith(s['unr'], s['ung'], 'sub'), 2, 'pow'), 'mul'), _arith(_arith(F(0.3), s['unr'], 'mul'), s['ung'], 'mul'), 'add'), F(0.5), 'pow'), 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['n0g'][s['k'] - 1], 'mul'), _arith(_arith(_div(F(1.0), _arith(_arith(s['lamr'][s['k'] - 1], 3, 'pow'), s['lamg'][s['k'] - 1], 'mul')), _div(F(1.0), _arith(_arith(s['lamr'][s['k'] - 1], 2, 'pow'), _arith(s['lamg'][s['k'] - 1], 2, 'pow'), 'mul')), 'add'), _div(F(1.0), _arith(s['lamr'][s['k'] - 1], _arith(s['lamg'][s['k'] - 1], 3, 'pow'), 'mul')), 'add'), 'mul')))
                s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_846)
                old_847 = s['pracg']
                s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(jnp.minimum(s['pracg'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']))))
                s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_847)
                return s
            def no_835(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qr3d'][s['k'] - 1] >= F(1e-08))) & ((s['qg3d'][s['k'] - 1] >= F(1e-08))))), yes_835, no_835, s)
            def yes_848(s):
                s = dict(s)
                def yes_849(s):
                    s = dict(s)
                    def yes_850(s):
                        s = dict(s)
                        def yes_851(s):
                            s = dict(s)
                            def yes_852(s):
                                s = dict(s)
                                old_853 = s['fmult']
                                s['fmult'] = F(F(0.0))
                                s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_853)
                                return s
                            def no_852(s):
                                s = dict(s)
                                def yes_854(s):
                                    s = dict(s)
                                    old_855 = s['fmult']
                                    s['fmult'] = F(_div(_arith(F(270.16), s['t3d'][s['k'] - 1], 'sub'), F(2.0)))
                                    s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_855)
                                    return s
                                def no_854(s):
                                    s = dict(s)
                                    def yes_856(s):
                                        s = dict(s)
                                        old_857 = s['fmult']
                                        s['fmult'] = F(_div(_arith(s['t3d'][s['k'] - 1], F(265.16), 'sub'), F(3.0)))
                                        s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_857)
                                        return s
                                    def no_856(s):
                                        s = dict(s)
                                        def yes_858(s):
                                            s = dict(s)
                                            old_859 = s['fmult']
                                            s['fmult'] = F(F(0.0))
                                            s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_859)
                                            return s
                                        def no_858(s):
                                            s = dict(s)
                                            return s
                                        s = lax.cond((s["_pc"] == 0) & ((s['t3d'][s['k'] - 1] < F(265.16))), yes_858, no_858, s)
                                        return s
                                    s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] >= F(265.16))) & ((s['t3d'][s['k'] - 1] <= F(268.16))))), yes_856, no_856, s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] <= F(270.16))) & ((s['t3d'][s['k'] - 1] > F(268.16))))), yes_854, no_854, s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['t3d'][s['k'] - 1] > F(270.16))), yes_852, no_852, s)
                            def yes_860(s):
                                s = dict(s)
                                old_861 = s['nmults']
                                s['nmults'] = s['nmults'].at[s['k'] - 1].set(F(_arith(_arith(_arith(F(350000.0), s['psacws'][s['k'] - 1], 'mul'), s['fmult'], 'mul'), F(1000.0), 'mul')))
                                s['nmults'] = jnp.where(s["_pc"] == 0, s['nmults'], old_861)
                                old_862 = s['qmults']
                                s['qmults'] = s['qmults'].at[s['k'] - 1].set(F(_arith(s['nmults'][s['k'] - 1], s['mmult'], 'mul')))
                                s['qmults'] = jnp.where(s["_pc"] == 0, s['qmults'], old_862)
                                old_863 = s['qmults']
                                s['qmults'] = s['qmults'].at[s['k'] - 1].set(F(jnp.minimum(s['qmults'][s['k'] - 1], s['psacws'][s['k'] - 1])))
                                s['qmults'] = jnp.where(s["_pc"] == 0, s['qmults'], old_863)
                                old_864 = s['psacws']
                                s['psacws'] = s['psacws'].at[s['k'] - 1].set(F(_arith(s['psacws'][s['k'] - 1], s['qmults'][s['k'] - 1], 'sub')))
                                s['psacws'] = jnp.where(s["_pc"] == 0, s['psacws'], old_864)
                                return s
                            def no_860(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['psacws'][s['k'] - 1] > F(0.0))), yes_860, no_860, s)
                            def yes_865(s):
                                s = dict(s)
                                old_866 = s['nmultr']
                                s['nmultr'] = s['nmultr'].at[s['k'] - 1].set(F(_arith(_arith(_arith(F(350000.0), s['pracs'][s['k'] - 1], 'mul'), s['fmult'], 'mul'), F(1000.0), 'mul')))
                                s['nmultr'] = jnp.where(s["_pc"] == 0, s['nmultr'], old_866)
                                old_867 = s['qmultr']
                                s['qmultr'] = s['qmultr'].at[s['k'] - 1].set(F(_arith(s['nmultr'][s['k'] - 1], s['mmult'], 'mul')))
                                s['qmultr'] = jnp.where(s["_pc"] == 0, s['qmultr'], old_867)
                                old_868 = s['qmultr']
                                s['qmultr'] = s['qmultr'].at[s['k'] - 1].set(F(jnp.minimum(s['qmultr'][s['k'] - 1], s['pracs'][s['k'] - 1])))
                                s['qmultr'] = jnp.where(s["_pc"] == 0, s['qmultr'], old_868)
                                old_869 = s['pracs']
                                s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(_arith(s['pracs'][s['k'] - 1], s['qmultr'][s['k'] - 1], 'sub')))
                                s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_869)
                                return s
                            def no_865(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['pracs'][s['k'] - 1] > F(0.0))), yes_865, no_865, s)
                            return s
                        def no_851(s):
                            s = dict(s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] < F(270.16))) & ((s['t3d'][s['k'] - 1] > F(265.16))))), yes_851, no_851, s)
                        return s
                    def no_850(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((((s['psacws'][s['k'] - 1] > F(0.0))) | ((s['pracs'][s['k'] - 1] > F(0.0))))), yes_850, no_850, s)
                    return s
                def no_849(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] >= F(0.0005))) | ((s['qr3d'][s['k'] - 1] >= F(0.0001))))), yes_849, no_849, s)
                return s
            def no_848(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= F(0.0001))), yes_848, no_848, s)
            def yes_870(s):
                s = dict(s)
                def yes_871(s):
                    s = dict(s)
                    def yes_872(s):
                        s = dict(s)
                        def yes_873(s):
                            s = dict(s)
                            def yes_874(s):
                                s = dict(s)
                                old_875 = s['fmult']
                                s['fmult'] = F(F(0.0))
                                s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_875)
                                return s
                            def no_874(s):
                                s = dict(s)
                                def yes_876(s):
                                    s = dict(s)
                                    old_877 = s['fmult']
                                    s['fmult'] = F(_div(_arith(F(270.16), s['t3d'][s['k'] - 1], 'sub'), F(2.0)))
                                    s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_877)
                                    return s
                                def no_876(s):
                                    s = dict(s)
                                    def yes_878(s):
                                        s = dict(s)
                                        old_879 = s['fmult']
                                        s['fmult'] = F(_div(_arith(s['t3d'][s['k'] - 1], F(265.16), 'sub'), F(3.0)))
                                        s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_879)
                                        return s
                                    def no_878(s):
                                        s = dict(s)
                                        def yes_880(s):
                                            s = dict(s)
                                            old_881 = s['fmult']
                                            s['fmult'] = F(F(0.0))
                                            s['fmult'] = jnp.where(s["_pc"] == 0, s['fmult'], old_881)
                                            return s
                                        def no_880(s):
                                            s = dict(s)
                                            return s
                                        s = lax.cond((s["_pc"] == 0) & ((s['t3d'][s['k'] - 1] < F(265.16))), yes_880, no_880, s)
                                        return s
                                    s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] >= F(265.16))) & ((s['t3d'][s['k'] - 1] <= F(268.16))))), yes_878, no_878, s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] <= F(270.16))) & ((s['t3d'][s['k'] - 1] > F(268.16))))), yes_876, no_876, s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['t3d'][s['k'] - 1] > F(270.16))), yes_874, no_874, s)
                            def yes_882(s):
                                s = dict(s)
                                old_883 = s['nmultg']
                                s['nmultg'] = s['nmultg'].at[s['k'] - 1].set(F(_arith(_arith(_arith(F(350000.0), s['psacwg'][s['k'] - 1], 'mul'), s['fmult'], 'mul'), F(1000.0), 'mul')))
                                s['nmultg'] = jnp.where(s["_pc"] == 0, s['nmultg'], old_883)
                                old_884 = s['qmultg']
                                s['qmultg'] = s['qmultg'].at[s['k'] - 1].set(F(_arith(s['nmultg'][s['k'] - 1], s['mmult'], 'mul')))
                                s['qmultg'] = jnp.where(s["_pc"] == 0, s['qmultg'], old_884)
                                old_885 = s['qmultg']
                                s['qmultg'] = s['qmultg'].at[s['k'] - 1].set(F(jnp.minimum(s['qmultg'][s['k'] - 1], s['psacwg'][s['k'] - 1])))
                                s['qmultg'] = jnp.where(s["_pc"] == 0, s['qmultg'], old_885)
                                old_886 = s['psacwg']
                                s['psacwg'] = s['psacwg'].at[s['k'] - 1].set(F(_arith(s['psacwg'][s['k'] - 1], s['qmultg'][s['k'] - 1], 'sub')))
                                s['psacwg'] = jnp.where(s["_pc"] == 0, s['psacwg'], old_886)
                                return s
                            def no_882(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['psacwg'][s['k'] - 1] > F(0.0))), yes_882, no_882, s)
                            def yes_887(s):
                                s = dict(s)
                                old_888 = s['nmultrg']
                                s['nmultrg'] = s['nmultrg'].at[s['k'] - 1].set(F(_arith(_arith(_arith(F(350000.0), s['pracg'][s['k'] - 1], 'mul'), s['fmult'], 'mul'), F(1000.0), 'mul')))
                                s['nmultrg'] = jnp.where(s["_pc"] == 0, s['nmultrg'], old_888)
                                old_889 = s['qmultrg']
                                s['qmultrg'] = s['qmultrg'].at[s['k'] - 1].set(F(_arith(s['nmultrg'][s['k'] - 1], s['mmult'], 'mul')))
                                s['qmultrg'] = jnp.where(s["_pc"] == 0, s['qmultrg'], old_889)
                                old_890 = s['qmultrg']
                                s['qmultrg'] = s['qmultrg'].at[s['k'] - 1].set(F(jnp.minimum(s['qmultrg'][s['k'] - 1], s['pracg'][s['k'] - 1])))
                                s['qmultrg'] = jnp.where(s["_pc"] == 0, s['qmultrg'], old_890)
                                old_891 = s['pracg']
                                s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(_arith(s['pracg'][s['k'] - 1], s['qmultrg'][s['k'] - 1], 'sub')))
                                s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_891)
                                return s
                            def no_887(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['pracg'][s['k'] - 1] > F(0.0))), yes_887, no_887, s)
                            return s
                        def no_873(s):
                            s = dict(s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] < F(270.16))) & ((s['t3d'][s['k'] - 1] > F(265.16))))), yes_873, no_873, s)
                        return s
                    def no_872(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((((s['psacwg'][s['k'] - 1] > F(0.0))) | ((s['pracg'][s['k'] - 1] > F(0.0))))), yes_872, no_872, s)
                    return s
                def no_871(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] >= F(0.0005))) | ((s['qr3d'][s['k'] - 1] >= F(0.0001))))), yes_871, no_871, s)
                return s
            def no_870(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] >= F(0.0001))), yes_870, no_870, s)
            def yes_892(s):
                s = dict(s)
                def yes_893(s):
                    s = dict(s)
                    old_894 = s['pgsacw']
                    s['pgsacw'] = s['pgsacw'].at[s['k'] - 1].set(F(jnp.minimum(s['psacws'][s['k'] - 1], _div(_arith(_arith(_arith(_arith(_arith(_arith(s['cons17'], s['dt'], 'mul'), s['n0s'][s['k'] - 1], 'mul'), s['qc3d'][s['k'] - 1], 'mul'), s['qc3d'][s['k'] - 1], 'mul'), s['asn'][s['k'] - 1], 'mul'), s['asn'][s['k'] - 1], 'mul'), _arith(s['rho'][s['k'] - 1], _arith(s['lams'][s['k'] - 1], _arith(_arith(F(2.0), s['bs'], 'mul'), F(2.0), 'add'), 'pow'), 'mul')))))
                    s['pgsacw'] = jnp.where(s["_pc"] == 0, s['pgsacw'], old_894)
                    old_895 = s['dum']
                    s['dum'] = F(jnp.maximum(_arith(_div(s['rhosn'], _arith(s['rhog'], s['rhosn'], 'sub')), s['pgsacw'][s['k'] - 1], 'mul'), F(0.0)))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_895)
                    old_896 = s['nscng']
                    s['nscng'] = s['nscng'].at[s['k'] - 1].set(F(_arith(_div(s['dum'], s['mg0']), s['rho'][s['k'] - 1], 'mul')))
                    s['nscng'] = jnp.where(s["_pc"] == 0, s['nscng'], old_896)
                    old_897 = s['nscng']
                    s['nscng'] = s['nscng'].at[s['k'] - 1].set(F(jnp.minimum(s['nscng'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']))))
                    s['nscng'] = jnp.where(s["_pc"] == 0, s['nscng'], old_897)
                    old_898 = s['psacws']
                    s['psacws'] = s['psacws'].at[s['k'] - 1].set(F(_arith(s['psacws'][s['k'] - 1], s['pgsacw'][s['k'] - 1], 'sub')))
                    s['psacws'] = jnp.where(s["_pc"] == 0, s['psacws'], old_898)
                    return s
                def no_893(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((s['qni3d'][s['k'] - 1] >= F(0.0001))) & ((s['qc3d'][s['k'] - 1] >= F(0.0005))))), yes_893, no_893, s)
                return s
            def no_892(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['psacws'][s['k'] - 1] > F(0.0))), yes_892, no_892, s)
            def yes_899(s):
                s = dict(s)
                def yes_900(s):
                    s = dict(s)
                    old_901 = s['dum']
                    s['dum'] = F(_div(_arith(_arith(s['cons18'], _arith(_div(F(4.0), s['lams'][s['k'] - 1]), 3, 'pow'), 'mul'), _arith(_div(F(4.0), s['lams'][s['k'] - 1]), 3, 'pow'), 'mul'), _arith(_arith(_arith(s['cons18'], _arith(_div(F(4.0), s['lams'][s['k'] - 1]), 3, 'pow'), 'mul'), _arith(_div(F(4.0), s['lams'][s['k'] - 1]), 3, 'pow'), 'mul'), _arith(_arith(s['cons19'], _arith(_div(F(4.0), s['lamr'][s['k'] - 1]), 3, 'pow'), 'mul'), _arith(_div(F(4.0), s['lamr'][s['k'] - 1]), 3, 'pow'), 'mul'), 'add')))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_901)
                    old_902 = s['dum']
                    s['dum'] = F(jnp.minimum(s['dum'], F(1.0)))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_902)
                    old_903 = s['dum']
                    s['dum'] = F(jnp.maximum(s['dum'], F(0.0)))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_903)
                    old_904 = s['pgracs']
                    s['pgracs'] = s['pgracs'].at[s['k'] - 1].set(F(_arith(_arith(F(1.0), s['dum'], 'sub'), s['pracs'][s['k'] - 1], 'mul')))
                    s['pgracs'] = jnp.where(s["_pc"] == 0, s['pgracs'], old_904)
                    old_905 = s['ngracs']
                    s['ngracs'] = s['ngracs'].at[s['k'] - 1].set(F(_arith(_arith(F(1.0), s['dum'], 'sub'), s['npracs'][s['k'] - 1], 'mul')))
                    s['ngracs'] = jnp.where(s["_pc"] == 0, s['ngracs'], old_905)
                    old_906 = s['ngracs']
                    s['ngracs'] = s['ngracs'].at[s['k'] - 1].set(F(jnp.minimum(s['ngracs'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']))))
                    s['ngracs'] = jnp.where(s["_pc"] == 0, s['ngracs'], old_906)
                    old_907 = s['ngracs']
                    s['ngracs'] = s['ngracs'].at[s['k'] - 1].set(F(jnp.minimum(s['ngracs'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']))))
                    s['ngracs'] = jnp.where(s["_pc"] == 0, s['ngracs'], old_907)
                    old_908 = s['pracs']
                    s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(_arith(s['pracs'][s['k'] - 1], s['pgracs'][s['k'] - 1], 'sub')))
                    s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_908)
                    old_909 = s['npracs']
                    s['npracs'] = s['npracs'].at[s['k'] - 1].set(F(_arith(s['npracs'][s['k'] - 1], s['ngracs'][s['k'] - 1], 'sub')))
                    s['npracs'] = jnp.where(s["_pc"] == 0, s['npracs'], old_909)
                    old_910 = s['psacr']
                    s['psacr'] = s['psacr'].at[s['k'] - 1].set(F(_arith(s['psacr'][s['k'] - 1], _arith(F(1.0), s['dum'], 'sub'), 'mul')))
                    s['psacr'] = jnp.where(s["_pc"] == 0, s['psacr'], old_910)
                    return s
                def no_900(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((s['qni3d'][s['k'] - 1] >= F(0.0001))) & ((s['qr3d'][s['k'] - 1] >= F(0.0001))))), yes_900, no_900, s)
                return s
            def no_899(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['pracs'][s['k'] - 1] > F(0.0))), yes_899, no_899, s)
            def yes_911(s):
                s = dict(s)
                old_912 = s['mnuccr']
                s['mnuccr'] = s['mnuccr'].at[s['k'] - 1].set(F(_div(_div(_arith(_arith(s['cons20'], s['nr3d'][s['k'] - 1], 'mul'), _intrinsic('exp', _arith(s['aimm'], _arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), 'mul')), 'mul'), _arith(s['lamr'][s['k'] - 1], 3, 'pow')), _arith(s['lamr'][s['k'] - 1], 3, 'pow'))))
                s['mnuccr'] = jnp.where(s["_pc"] == 0, s['mnuccr'], old_912)
                old_913 = s['nnuccr']
                s['nnuccr'] = s['nnuccr'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(s['pi'], s['nr3d'][s['k'] - 1], 'mul'), s['bimm'], 'mul'), _intrinsic('exp', _arith(s['aimm'], _arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), 'mul')), 'mul'), _arith(s['lamr'][s['k'] - 1], 3, 'pow'))))
                s['nnuccr'] = jnp.where(s["_pc"] == 0, s['nnuccr'], old_913)
                old_914 = s['nnuccr']
                s['nnuccr'] = s['nnuccr'].at[s['k'] - 1].set(F(jnp.minimum(s['nnuccr'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']))))
                s['nnuccr'] = jnp.where(s["_pc"] == 0, s['nnuccr'], old_914)
                return s
            def no_911(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] < F(269.15))) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])))), yes_911, no_911, s)
            def yes_915(s):
                s = dict(s)
                def yes_916(s):
                    s = dict(s)
                    old_917 = s['dum']
                    s['dum'] = F(_arith(s['qc3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'mul'))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_917)
                    old_918 = s['pra']
                    s['pra'] = s['pra'].at[s['k'] - 1].set(F(_arith(F(67.0), _arith(s['dum'], F(1.15), 'pow'), 'mul')))
                    s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_918)
                    old_919 = s['npra']
                    s['npra'] = s['npra'].at[s['k'] - 1].set(F(_div(s['pra'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['nc3d'][s['k'] - 1]))))
                    s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_919)
                    return s
                def no_916(s):
                    s = dict(s)
                    def yes_920(s):
                        s = dict(s)
                        old_921 = s['dum']
                        s['dum'] = F(_arith(F(1.0), _div(s['qc3d'][s['k'] - 1], _arith(s['qc3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')), 'sub'))
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_921)
                        old_922 = s['dum1']
                        s['dum1'] = F(_arith(_div(s['dum'], _arith(s['dum'], F(0.0005), 'add')), 4, 'pow'))
                        s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_922)
                        old_923 = s['pra']
                        s['pra'] = s['pra'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_div(_arith(F(5780.0), s['rho'][s['k'] - 1], 'mul'), F(1000.0)), s['qc3d'][s['k'] - 1], 'mul'), s['qr3d'][s['k'] - 1], 'mul'), s['dum1'], 'mul')))
                        s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_923)
                        old_924 = s['npra']
                        s['npra'] = s['npra'].at[s['k'] - 1].set(F(_div(_arith(_div(_arith(_div(_arith(s['pra'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), F(1000.0)), _div(_arith(s['nc3d'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), F(1000000.0)), 'mul'), _div(_arith(s['qc3d'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), F(1000.0))), F(1000000.0), 'mul'), s['rho'][s['k'] - 1])))
                        s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_924)
                        return s
                    def no_920(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 1)), yes_920, no_920, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['irain'] == 0)), yes_916, no_916, s)
                return s
            def no_915(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qr3d'][s['k'] - 1] >= F(1e-08))) & ((s['qc3d'][s['k'] - 1] >= F(1e-08))))), yes_915, no_915, s)
            def yes_925(s):
                s = dict(s)
                old_926 = s['dum1']
                s['dum1'] = F(F(0.0003))
                s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_926)
                def yes_927(s):
                    s = dict(s)
                    old_928 = s['dum']
                    s['dum'] = F(F(1.0))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_928)
                    return s
                def no_927(s):
                    s = dict(s)
                    def yes_929(s):
                        s = dict(s)
                        old_930 = s['dum']
                        s['dum'] = F(_arith(F(2.0), _intrinsic('exp', _arith(F(2300.0), _arith(_div(F(1.0), s['lamr'][s['k'] - 1]), s['dum1'], 'sub'), 'mul')), 'sub'))
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_930)
                        return s
                    def no_929(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((_div(F(1.0), s['lamr'][s['k'] - 1]) >= s['dum1'])), yes_929, no_929, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((_div(F(1.0), s['lamr'][s['k'] - 1]) < s['dum1'])), yes_927, no_927, s)
                old_931 = s['nragg']
                s['nragg'] = s['nragg'].at[s['k'] - 1].set(F(_arith(_arith(_arith(_arith((-F(5.78)), s['dum'], 'mul'), s['nr3d'][s['k'] - 1], 'mul'), s['qr3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul')))
                s['nragg'] = jnp.where(s["_pc"] == 0, s['nragg'], old_931)
                return s
            def no_925(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= F(1e-08))), yes_925, no_925, s)
            def yes_932(s):
                s = dict(s)
                old_933 = s['nprci']
                s['nprci'] = s['nprci'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(_arith(s['cons21'], _arith(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0i'][s['k'] - 1], 'mul'), _intrinsic('exp', _arith((-s['lami'][s['k'] - 1]), s['dcs'], 'mul')), 'mul'), s['dv'][s['k'] - 1], 'mul'), s['abi'][s['k'] - 1])))
                s['nprci'] = jnp.where(s["_pc"] == 0, s['nprci'], old_933)
                old_934 = s['prci']
                s['prci'] = s['prci'].at[s['k'] - 1].set(F(_arith(s['cons22'], s['nprci'][s['k'] - 1], 'mul')))
                s['prci'] = jnp.where(s["_pc"] == 0, s['prci'], old_934)
                old_935 = s['nprci']
                s['nprci'] = s['nprci'].at[s['k'] - 1].set(F(jnp.minimum(s['nprci'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']))))
                s['nprci'] = jnp.where(s["_pc"] == 0, s['nprci'], old_935)
                return s
            def no_932(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qi3d'][s['k'] - 1] >= F(1e-08))) & ((s['qvqvsi'][s['k'] - 1] >= F(1.0))))), yes_932, no_932, s)
            def yes_936(s):
                s = dict(s)
                old_937 = s['prai']
                s['prai'] = s['prai'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons23'], s['asn'][s['k'] - 1], 'mul'), s['qi3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(s['lams'][s['k'] - 1], _arith(s['bs'], F(3.0), 'add'), 'pow'))))
                s['prai'] = jnp.where(s["_pc"] == 0, s['prai'], old_937)
                old_938 = s['nprai']
                s['nprai'] = s['nprai'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(_arith(s['cons23'], s['asn'][s['k'] - 1], 'mul'), s['ni3d'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['n0s'][s['k'] - 1], 'mul'), _arith(s['lams'][s['k'] - 1], _arith(s['bs'], F(3.0), 'add'), 'pow'))))
                s['nprai'] = jnp.where(s["_pc"] == 0, s['nprai'], old_938)
                old_939 = s['nprai']
                s['nprai'] = s['nprai'].at[s['k'] - 1].set(F(jnp.minimum(s['nprai'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']))))
                s['nprai'] = jnp.where(s["_pc"] == 0, s['nprai'], old_939)
                return s
            def no_936(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qni3d'][s['k'] - 1] >= F(1e-08))) & ((s['qi3d'][s['k'] - 1] >= s['qsmall'])))), yes_936, no_936, s)
            def yes_940(s):
                s = dict(s)
                def yes_941(s):
                    s = dict(s)
                    old_942 = s['niacr']
                    s['niacr'] = s['niacr'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(_arith(s['cons24'], s['ni3d'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['arn'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], _arith(s['br'], F(3.0), 'add'), 'pow')), s['rho'][s['k'] - 1], 'mul')))
                    s['niacr'] = jnp.where(s["_pc"] == 0, s['niacr'], old_942)
                    old_943 = s['piacr']
                    s['piacr'] = s['piacr'].at[s['k'] - 1].set(F(_arith(_div(_div(_arith(_arith(_arith(s['cons25'], s['ni3d'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['arn'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], _arith(s['br'], F(3.0), 'add'), 'pow')), _arith(s['lamr'][s['k'] - 1], 3, 'pow')), s['rho'][s['k'] - 1], 'mul')))
                    s['piacr'] = jnp.where(s["_pc"] == 0, s['piacr'], old_943)
                    old_944 = s['praci']
                    s['praci'] = s['praci'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(_arith(s['cons24'], s['qi3d'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['arn'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], _arith(s['br'], F(3.0), 'add'), 'pow')), s['rho'][s['k'] - 1], 'mul')))
                    s['praci'] = jnp.where(s["_pc"] == 0, s['praci'], old_944)
                    old_945 = s['niacr']
                    s['niacr'] = s['niacr'].at[s['k'] - 1].set(F(jnp.minimum(s['niacr'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']))))
                    s['niacr'] = jnp.where(s["_pc"] == 0, s['niacr'], old_945)
                    old_946 = s['niacr']
                    s['niacr'] = s['niacr'].at[s['k'] - 1].set(F(jnp.minimum(s['niacr'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']))))
                    s['niacr'] = jnp.where(s["_pc"] == 0, s['niacr'], old_946)
                    return s
                def no_941(s):
                    s = dict(s)
                    old_947 = s['niacrs']
                    s['niacrs'] = s['niacrs'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(_arith(s['cons24'], s['ni3d'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['arn'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], _arith(s['br'], F(3.0), 'add'), 'pow')), s['rho'][s['k'] - 1], 'mul')))
                    s['niacrs'] = jnp.where(s["_pc"] == 0, s['niacrs'], old_947)
                    old_948 = s['piacrs']
                    s['piacrs'] = s['piacrs'].at[s['k'] - 1].set(F(_arith(_div(_div(_arith(_arith(_arith(s['cons25'], s['ni3d'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['arn'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], _arith(s['br'], F(3.0), 'add'), 'pow')), _arith(s['lamr'][s['k'] - 1], 3, 'pow')), s['rho'][s['k'] - 1], 'mul')))
                    s['piacrs'] = jnp.where(s["_pc"] == 0, s['piacrs'], old_948)
                    old_949 = s['pracis']
                    s['pracis'] = s['pracis'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(_arith(s['cons24'], s['qi3d'][s['k'] - 1], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['arn'][s['k'] - 1], 'mul'), _arith(s['lamr'][s['k'] - 1], _arith(s['br'], F(3.0), 'add'), 'pow')), s['rho'][s['k'] - 1], 'mul')))
                    s['pracis'] = jnp.where(s["_pc"] == 0, s['pracis'], old_949)
                    old_950 = s['niacrs']
                    s['niacrs'] = s['niacrs'].at[s['k'] - 1].set(F(jnp.minimum(s['niacrs'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']))))
                    s['niacrs'] = jnp.where(s["_pc"] == 0, s['niacrs'], old_950)
                    old_951 = s['niacrs']
                    s['niacrs'] = s['niacrs'].at[s['k'] - 1].set(F(jnp.minimum(s['niacrs'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']))))
                    s['niacrs'] = jnp.where(s["_pc"] == 0, s['niacrs'], old_951)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= F(0.0001))), yes_941, no_941, s)
                return s
            def no_940(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['qr3d'][s['k'] - 1] >= F(1e-08))) & ((s['qi3d'][s['k'] - 1] >= F(1e-08))) & ((s['t3d'][s['k'] - 1] <= s['tmelt'])))), yes_940, no_940, s)
            def yes_952(s):
                s = dict(s)
                def yes_953(s):
                    s = dict(s)
                    old_954 = s['kc2']
                    s['kc2'] = F(_arith(_arith(F(0.005), _intrinsic('exp', _arith(F(0.304), _arith(s['tmelt'], s['t3d'][s['k'] - 1], 'sub'), 'mul')), 'mul'), F(1000.0), 'mul'))
                    s['kc2'] = jnp.where(s["_pc"] == 0, s['kc2'], old_954)
                    old_955 = s['kc2']
                    s['kc2'] = F(jnp.minimum(s['kc2'], F(500000.0)))
                    s['kc2'] = jnp.where(s["_pc"] == 0, s['kc2'], old_955)
                    old_956 = s['kc2']
                    s['kc2'] = F(jnp.maximum(_div(s['kc2'], s['rho'][s['k'] - 1]), F(0.0)))
                    s['kc2'] = jnp.where(s["_pc"] == 0, s['kc2'], old_956)
                    def yes_957(s):
                        s = dict(s)
                        old_958 = s['nnuccd']
                        s['nnuccd'] = s['nnuccd'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(s['kc2'], s['ni3d'][s['k'] - 1], 'sub'), s['ns3d'][s['k'] - 1], 'sub'), s['ng3d'][s['k'] - 1], 'sub'), s['dt'])))
                        s['nnuccd'] = jnp.where(s["_pc"] == 0, s['nnuccd'], old_958)
                        old_959 = s['mnuccd']
                        s['mnuccd'] = s['mnuccd'].at[s['k'] - 1].set(F(_arith(s['nnuccd'][s['k'] - 1], s['mi0'], 'mul')))
                        s['mnuccd'] = jnp.where(s["_pc"] == 0, s['mnuccd'], old_959)
                        return s
                    def no_957(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['kc2'] > _arith(_arith(s['ni3d'][s['k'] - 1], s['ns3d'][s['k'] - 1], 'add'), s['ng3d'][s['k'] - 1], 'add'))), yes_957, no_957, s)
                    return s
                def no_953(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((((s['qvqvs'][s['k'] - 1] >= F(0.999))) & ((s['t3d'][s['k'] - 1] <= F(265.15))))) | ((s['qvqvsi'][s['k'] - 1] >= F(1.08))))), yes_953, no_953, s)
                return s
            def no_952(s):
                s = dict(s)
                def yes_960(s):
                    s = dict(s)
                    def yes_961(s):
                        s = dict(s)
                        old_962 = s['kc2']
                        s['kc2'] = F(_div(_arith(F(0.16), F(1000.0), 'mul'), s['rho'][s['k'] - 1]))
                        s['kc2'] = jnp.where(s["_pc"] == 0, s['kc2'], old_962)
                        def yes_963(s):
                            s = dict(s)
                            old_964 = s['nnuccd']
                            s['nnuccd'] = s['nnuccd'].at[s['k'] - 1].set(F(_div(_arith(_arith(_arith(s['kc2'], s['ni3d'][s['k'] - 1], 'sub'), s['ns3d'][s['k'] - 1], 'sub'), s['ng3d'][s['k'] - 1], 'sub'), s['dt'])))
                            s['nnuccd'] = jnp.where(s["_pc"] == 0, s['nnuccd'], old_964)
                            old_965 = s['mnuccd']
                            s['mnuccd'] = s['mnuccd'].at[s['k'] - 1].set(F(_arith(s['nnuccd'][s['k'] - 1], s['mi0'], 'mul')))
                            s['mnuccd'] = jnp.where(s["_pc"] == 0, s['mnuccd'], old_965)
                            return s
                        def no_963(s):
                            s = dict(s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['kc2'] > _arith(_arith(s['ni3d'][s['k'] - 1], s['ns3d'][s['k'] - 1], 'add'), s['ng3d'][s['k'] - 1], 'add'))), yes_963, no_963, s)
                        return s
                    def no_961(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] < s['tmelt'])) & ((s['qvqvsi'][s['k'] - 1] > F(1.0))))), yes_961, no_961, s)
                    return s
                def no_960(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['inuc'] == 1)), yes_960, no_960, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['inuc'] == 0)), yes_952, no_952, s)
            old_966 = s['nnuccd']
            s['nnuccd'] = s['nnuccd'].at[s['k'] - 1].set(F(_arith(s['nnuccd'][s['k'] - 1], s['nnuccd_reduce_coef'], 'mul')))
            s['nnuccd'] = jnp.where(s["_pc"] == 0, s['nnuccd'], old_966)
            s["_pc"] = jnp.where(s["_pc"] == 101, I(0), s["_pc"])
            def yes_967(s):
                s = dict(s)
                old_968 = s['epsi']
                s['epsi'] = F(_div(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0i'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), _arith(s['lami'][s['k'] - 1], s['lami'][s['k'] - 1], 'mul')))
                s['epsi'] = jnp.where(s["_pc"] == 0, s['epsi'], old_968)
                return s
            def no_967(s):
                s = dict(s)
                old_969 = s['epsi']
                s['epsi'] = F(F(0.0))
                s['epsi'] = jnp.where(s["_pc"] == 0, s['epsi'], old_969)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] >= s['qsmall'])), yes_967, no_967, s)
            def yes_970(s):
                s = dict(s)
                old_971 = s['epss']
                s['epss'] = F(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0s'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), _arith(_div(s['f1s'], _arith(s['lams'][s['k'] - 1], s['lams'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2s'], _arith(_div(_arith(s['asn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons10'], 'mul'), _arith(s['lams'][s['k'] - 1], s['cons35'], 'pow')), 'add'), 'mul'))
                s['epss'] = jnp.where(s["_pc"] == 0, s['epss'], old_971)
                return s
            def no_970(s):
                s = dict(s)
                old_972 = s['epss']
                s['epss'] = F(F(0.0))
                s['epss'] = jnp.where(s["_pc"] == 0, s['epss'], old_972)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])), yes_970, no_970, s)
            def yes_973(s):
                s = dict(s)
                old_974 = s['epsg']
                s['epsg'] = F(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0g'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), _arith(_div(s['f1s'], _arith(s['lamg'][s['k'] - 1], s['lamg'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2s'], _arith(_div(_arith(s['agn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons11'], 'mul'), _arith(s['lamg'][s['k'] - 1], s['cons36'], 'pow')), 'add'), 'mul'))
                s['epsg'] = jnp.where(s["_pc"] == 0, s['epsg'], old_974)
                return s
            def no_973(s):
                s = dict(s)
                old_975 = s['epsg']
                s['epsg'] = F(F(0.0))
                s['epsg'] = jnp.where(s["_pc"] == 0, s['epsg'], old_975)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] >= s['qsmall'])), yes_973, no_973, s)
            def yes_976(s):
                s = dict(s)
                old_977 = s['epsr']
                s['epsr'] = F(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['n0rr'][s['k'] - 1], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), _arith(_div(s['f1r'], _arith(s['lamr'][s['k'] - 1], s['lamr'][s['k'] - 1], 'mul')), _div(_arith(_arith(_arith(s['f2r'], _arith(_div(_arith(s['arn'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), s['mu'][s['k'] - 1]), F(0.5), 'pow'), 'mul'), _arith(s['sc'][s['k'] - 1], _div(F(1.0), F(3.0)), 'pow'), 'mul'), s['cons9'], 'mul'), _arith(s['lamr'][s['k'] - 1], s['cons34'], 'pow')), 'add'), 'mul'))
                s['epsr'] = jnp.where(s["_pc"] == 0, s['epsr'], old_977)
                return s
            def no_976(s):
                s = dict(s)
                old_978 = s['epsr']
                s['epsr'] = F(F(0.0))
                s['epsr'] = jnp.where(s["_pc"] == 0, s['epsr'], old_978)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])), yes_976, no_976, s)
            def yes_979(s):
                s = dict(s)
                old_980 = s['dum']
                s['dum'] = F(_arith(F(1.0), _arith(_intrinsic('exp', _arith((-s['lami'][s['k'] - 1]), s['dcs'], 'mul')), _arith(F(1.0), _arith(s['lami'][s['k'] - 1], s['dcs'], 'mul'), 'add'), 'mul'), 'sub'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_980)
                old_981 = s['prd']
                s['prd'] = s['prd'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['epsi'], _arith(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), 'mul'), s['abi'][s['k'] - 1]), s['dum'], 'mul')))
                s['prd'] = jnp.where(s["_pc"] == 0, s['prd'], old_981)
                return s
            def no_979(s):
                s = dict(s)
                old_982 = s['dum']
                s['dum'] = F(F(0.0))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_982)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] >= s['qsmall'])), yes_979, no_979, s)
            def yes_983(s):
                s = dict(s)
                old_984 = s['prds']
                s['prds'] = s['prds'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['epss'], _arith(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), 'mul'), s['abi'][s['k'] - 1]), _arith(_div(_arith(s['epsi'], _arith(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), 'mul'), s['abi'][s['k'] - 1]), _arith(F(1.0), s['dum'], 'sub'), 'mul'), 'add')))
                s['prds'] = jnp.where(s["_pc"] == 0, s['prds'], old_984)
                return s
            def no_983(s):
                s = dict(s)
                old_985 = s['prd']
                s['prd'] = s['prd'].at[s['k'] - 1].set(F(_arith(s['prd'][s['k'] - 1], _arith(_div(_arith(s['epsi'], _arith(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), 'mul'), s['abi'][s['k'] - 1]), _arith(F(1.0), s['dum'], 'sub'), 'mul'), 'add')))
                s['prd'] = jnp.where(s["_pc"] == 0, s['prd'], old_985)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])), yes_983, no_983, s)
            old_986 = s['prdg']
            s['prdg'] = s['prdg'].at[s['k'] - 1].set(F(_div(_arith(s['epsg'], _arith(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), 'mul'), s['abi'][s['k'] - 1])))
            s['prdg'] = jnp.where(s["_pc"] == 0, s['prdg'], old_986)
            def yes_987(s):
                s = dict(s)
                old_988 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(_div(_arith(s['epsr'], _arith(s['qv3d'][s['k'] - 1], s['qvs'][s['k'] - 1], 'sub'), 'mul'), s['ab'][s['k'] - 1])))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_988)
                old_989 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(jnp.minimum(s['pre'][s['k'] - 1], F(0.0))))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_989)
                return s
            def no_987(s):
                s = dict(s)
                old_990 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(F(0.0)))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_990)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qv3d'][s['k'] - 1] < s['qvs'][s['k'] - 1])), yes_987, no_987, s)
            old_991 = s['dum']
            s['dum'] = F(_div(_arith(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), s['dt']))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_991)
            old_992 = s['fudgef']
            s['fudgef'] = F(F(0.9999))
            s['fudgef'] = jnp.where(s["_pc"] == 0, s['fudgef'], old_992)
            old_993 = s['sum_dep']
            s['sum_dep'] = F(_arith(_arith(_arith(s['prd'][s['k'] - 1], s['prds'][s['k'] - 1], 'add'), s['mnuccd'][s['k'] - 1], 'add'), s['prdg'][s['k'] - 1], 'add'))
            s['sum_dep'] = jnp.where(s["_pc"] == 0, s['sum_dep'], old_993)
            def yes_994(s):
                s = dict(s)
                old_995 = s['mnuccd']
                s['mnuccd'] = s['mnuccd'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['fudgef'], s['mnuccd'][s['k'] - 1], 'mul'), s['dum'], 'mul'), s['sum_dep'])))
                s['mnuccd'] = jnp.where(s["_pc"] == 0, s['mnuccd'], old_995)
                old_996 = s['prd']
                s['prd'] = s['prd'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['fudgef'], s['prd'][s['k'] - 1], 'mul'), s['dum'], 'mul'), s['sum_dep'])))
                s['prd'] = jnp.where(s["_pc"] == 0, s['prd'], old_996)
                old_997 = s['prds']
                s['prds'] = s['prds'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['fudgef'], s['prds'][s['k'] - 1], 'mul'), s['dum'], 'mul'), s['sum_dep'])))
                s['prds'] = jnp.where(s["_pc"] == 0, s['prds'], old_997)
                old_998 = s['prdg']
                s['prdg'] = s['prdg'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['fudgef'], s['prdg'][s['k'] - 1], 'mul'), s['dum'], 'mul'), s['sum_dep'])))
                s['prdg'] = jnp.where(s["_pc"] == 0, s['prdg'], old_998)
                return s
            def no_994(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((((s['dum'] > F(0.0))) & ((s['sum_dep'] > _arith(s['dum'], s['fudgef'], 'mul'))))) | ((((s['dum'] < F(0.0))) & ((s['sum_dep'] < _arith(s['dum'], s['fudgef'], 'mul'))))))), yes_994, no_994, s)
            def yes_999(s):
                s = dict(s)
                old_1000 = s['eprd']
                s['eprd'] = s['eprd'].at[s['k'] - 1].set(F(s['prd'][s['k'] - 1]))
                s['eprd'] = jnp.where(s["_pc"] == 0, s['eprd'], old_1000)
                old_1001 = s['prd']
                s['prd'] = s['prd'].at[s['k'] - 1].set(F(F(0.0)))
                s['prd'] = jnp.where(s["_pc"] == 0, s['prd'], old_1001)
                return s
            def no_999(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['prd'][s['k'] - 1] < F(0.0))), yes_999, no_999, s)
            def yes_1002(s):
                s = dict(s)
                old_1003 = s['eprds']
                s['eprds'] = s['eprds'].at[s['k'] - 1].set(F(s['prds'][s['k'] - 1]))
                s['eprds'] = jnp.where(s["_pc"] == 0, s['eprds'], old_1003)
                old_1004 = s['prds']
                s['prds'] = s['prds'].at[s['k'] - 1].set(F(F(0.0)))
                s['prds'] = jnp.where(s["_pc"] == 0, s['prds'], old_1004)
                return s
            def no_1002(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['prds'][s['k'] - 1] < F(0.0))), yes_1002, no_1002, s)
            def yes_1005(s):
                s = dict(s)
                old_1006 = s['eprdg']
                s['eprdg'] = s['eprdg'].at[s['k'] - 1].set(F(s['prdg'][s['k'] - 1]))
                s['eprdg'] = jnp.where(s["_pc"] == 0, s['eprdg'], old_1006)
                old_1007 = s['prdg']
                s['prdg'] = s['prdg'].at[s['k'] - 1].set(F(F(0.0)))
                s['prdg'] = jnp.where(s["_pc"] == 0, s['prdg'], old_1007)
                return s
            def no_1005(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['prdg'][s['k'] - 1] < F(0.0))), yes_1005, no_1005, s)
            def yes_1008(s):
                s = dict(s)
                old_1009 = s['mnuccc']
                s['mnuccc'] = s['mnuccc'].at[s['k'] - 1].set(F(F(0.0)))
                s['mnuccc'] = jnp.where(s["_pc"] == 0, s['mnuccc'], old_1009)
                old_1010 = s['nnuccc']
                s['nnuccc'] = s['nnuccc'].at[s['k'] - 1].set(F(F(0.0)))
                s['nnuccc'] = jnp.where(s["_pc"] == 0, s['nnuccc'], old_1010)
                old_1011 = s['mnuccr']
                s['mnuccr'] = s['mnuccr'].at[s['k'] - 1].set(F(F(0.0)))
                s['mnuccr'] = jnp.where(s["_pc"] == 0, s['mnuccr'], old_1011)
                old_1012 = s['nnuccr']
                s['nnuccr'] = s['nnuccr'].at[s['k'] - 1].set(F(F(0.0)))
                s['nnuccr'] = jnp.where(s["_pc"] == 0, s['nnuccr'], old_1012)
                old_1013 = s['mnuccd']
                s['mnuccd'] = s['mnuccd'].at[s['k'] - 1].set(F(F(0.0)))
                s['mnuccd'] = jnp.where(s["_pc"] == 0, s['mnuccd'], old_1013)
                old_1014 = s['nnuccd']
                s['nnuccd'] = s['nnuccd'].at[s['k'] - 1].set(F(F(0.0)))
                s['nnuccd'] = jnp.where(s["_pc"] == 0, s['nnuccd'], old_1014)
                return s
            def no_1008(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['iliq'] == 1)), yes_1008, no_1008, s)
            def yes_1015(s):
                s = dict(s)
                old_1016 = s['pracg']
                s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(F(0.0)))
                s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_1016)
                old_1017 = s['psacr']
                s['psacr'] = s['psacr'].at[s['k'] - 1].set(F(F(0.0)))
                s['psacr'] = jnp.where(s["_pc"] == 0, s['psacr'], old_1017)
                old_1018 = s['psacwg']
                s['psacwg'] = s['psacwg'].at[s['k'] - 1].set(F(F(0.0)))
                s['psacwg'] = jnp.where(s["_pc"] == 0, s['psacwg'], old_1018)
                old_1019 = s['pgsacw']
                s['pgsacw'] = s['pgsacw'].at[s['k'] - 1].set(F(F(0.0)))
                s['pgsacw'] = jnp.where(s["_pc"] == 0, s['pgsacw'], old_1019)
                old_1020 = s['pgracs']
                s['pgracs'] = s['pgracs'].at[s['k'] - 1].set(F(F(0.0)))
                s['pgracs'] = jnp.where(s["_pc"] == 0, s['pgracs'], old_1020)
                old_1021 = s['prdg']
                s['prdg'] = s['prdg'].at[s['k'] - 1].set(F(F(0.0)))
                s['prdg'] = jnp.where(s["_pc"] == 0, s['prdg'], old_1021)
                old_1022 = s['eprdg']
                s['eprdg'] = s['eprdg'].at[s['k'] - 1].set(F(F(0.0)))
                s['eprdg'] = jnp.where(s["_pc"] == 0, s['eprdg'], old_1022)
                old_1023 = s['evpmg']
                s['evpmg'] = s['evpmg'].at[s['k'] - 1].set(F(F(0.0)))
                s['evpmg'] = jnp.where(s["_pc"] == 0, s['evpmg'], old_1023)
                old_1024 = s['pgmlt']
                s['pgmlt'] = s['pgmlt'].at[s['k'] - 1].set(F(F(0.0)))
                s['pgmlt'] = jnp.where(s["_pc"] == 0, s['pgmlt'], old_1024)
                old_1025 = s['npracg']
                s['npracg'] = s['npracg'].at[s['k'] - 1].set(F(F(0.0)))
                s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_1025)
                old_1026 = s['npsacwg']
                s['npsacwg'] = s['npsacwg'].at[s['k'] - 1].set(F(F(0.0)))
                s['npsacwg'] = jnp.where(s["_pc"] == 0, s['npsacwg'], old_1026)
                old_1027 = s['nscng']
                s['nscng'] = s['nscng'].at[s['k'] - 1].set(F(F(0.0)))
                s['nscng'] = jnp.where(s["_pc"] == 0, s['nscng'], old_1027)
                old_1028 = s['ngracs']
                s['ngracs'] = s['ngracs'].at[s['k'] - 1].set(F(F(0.0)))
                s['ngracs'] = jnp.where(s["_pc"] == 0, s['ngracs'], old_1028)
                old_1029 = s['nsubg']
                s['nsubg'] = s['nsubg'].at[s['k'] - 1].set(F(F(0.0)))
                s['nsubg'] = jnp.where(s["_pc"] == 0, s['nsubg'], old_1029)
                old_1030 = s['ngmltg']
                s['ngmltg'] = s['ngmltg'].at[s['k'] - 1].set(F(F(0.0)))
                s['ngmltg'] = jnp.where(s["_pc"] == 0, s['ngmltg'], old_1030)
                old_1031 = s['ngmltr']
                s['ngmltr'] = s['ngmltr'].at[s['k'] - 1].set(F(F(0.0)))
                s['ngmltr'] = jnp.where(s["_pc"] == 0, s['ngmltr'], old_1031)
                old_1032 = s['piacrs']
                s['piacrs'] = s['piacrs'].at[s['k'] - 1].set(F(_arith(s['piacrs'][s['k'] - 1], s['piacr'][s['k'] - 1], 'add')))
                s['piacrs'] = jnp.where(s["_pc"] == 0, s['piacrs'], old_1032)
                old_1033 = s['piacr']
                s['piacr'] = s['piacr'].at[s['k'] - 1].set(F(F(0.0)))
                s['piacr'] = jnp.where(s["_pc"] == 0, s['piacr'], old_1033)
                return s
            def no_1015(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 1)), yes_1015, no_1015, s)
            old_1034 = s['dum']
            s['dum'] = F(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['prc'][s['k'] - 1], s['pra'][s['k'] - 1], 'add'), s['mnuccc'][s['k'] - 1], 'add'), s['psacws'][s['k'] - 1], 'add'), s['psacwi'][s['k'] - 1], 'add'), s['qmults'][s['k'] - 1], 'add'), s['psacwg'][s['k'] - 1], 'add'), s['pgsacw'][s['k'] - 1], 'add'), s['qmultg'][s['k'] - 1], 'add'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1034)
            def yes_1035(s):
                s = dict(s)
                old_1036 = s['ratio']
                s['ratio'] = F(_div(s['qc3d'][s['k'] - 1], s['dum']))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_1036)
                old_1037 = s['prc']
                s['prc'] = s['prc'].at[s['k'] - 1].set(F(_arith(s['prc'][s['k'] - 1], s['ratio'], 'mul')))
                s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_1037)
                old_1038 = s['pra']
                s['pra'] = s['pra'].at[s['k'] - 1].set(F(_arith(s['pra'][s['k'] - 1], s['ratio'], 'mul')))
                s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_1038)
                old_1039 = s['mnuccc']
                s['mnuccc'] = s['mnuccc'].at[s['k'] - 1].set(F(_arith(s['mnuccc'][s['k'] - 1], s['ratio'], 'mul')))
                s['mnuccc'] = jnp.where(s["_pc"] == 0, s['mnuccc'], old_1039)
                old_1040 = s['psacws']
                s['psacws'] = s['psacws'].at[s['k'] - 1].set(F(_arith(s['psacws'][s['k'] - 1], s['ratio'], 'mul')))
                s['psacws'] = jnp.where(s["_pc"] == 0, s['psacws'], old_1040)
                old_1041 = s['psacwi']
                s['psacwi'] = s['psacwi'].at[s['k'] - 1].set(F(_arith(s['psacwi'][s['k'] - 1], s['ratio'], 'mul')))
                s['psacwi'] = jnp.where(s["_pc"] == 0, s['psacwi'], old_1041)
                old_1042 = s['qmults']
                s['qmults'] = s['qmults'].at[s['k'] - 1].set(F(_arith(s['qmults'][s['k'] - 1], s['ratio'], 'mul')))
                s['qmults'] = jnp.where(s["_pc"] == 0, s['qmults'], old_1042)
                old_1043 = s['qmultg']
                s['qmultg'] = s['qmultg'].at[s['k'] - 1].set(F(_arith(s['qmultg'][s['k'] - 1], s['ratio'], 'mul')))
                s['qmultg'] = jnp.where(s["_pc"] == 0, s['qmultg'], old_1043)
                old_1044 = s['psacwg']
                s['psacwg'] = s['psacwg'].at[s['k'] - 1].set(F(_arith(s['psacwg'][s['k'] - 1], s['ratio'], 'mul')))
                s['psacwg'] = jnp.where(s["_pc"] == 0, s['psacwg'], old_1044)
                old_1045 = s['pgsacw']
                s['pgsacw'] = s['pgsacw'].at[s['k'] - 1].set(F(_arith(s['pgsacw'][s['k'] - 1], s['ratio'], 'mul')))
                s['pgsacw'] = jnp.where(s["_pc"] == 0, s['pgsacw'], old_1045)
                return s
            def no_1035(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qc3d'][s['k'] - 1])) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])))), yes_1035, no_1035, s)
            old_1046 = s['dum']
            s['dum'] = F(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['prd'][s['k'] - 1]), s['mnuccc'][s['k'] - 1], 'sub'), s['prci'][s['k'] - 1], 'add'), s['prai'][s['k'] - 1], 'add'), s['qmults'][s['k'] - 1], 'sub'), s['qmultg'][s['k'] - 1], 'sub'), s['qmultr'][s['k'] - 1], 'sub'), s['qmultrg'][s['k'] - 1], 'sub'), s['mnuccd'][s['k'] - 1], 'sub'), s['praci'][s['k'] - 1], 'add'), s['pracis'][s['k'] - 1], 'add'), s['eprd'][s['k'] - 1], 'sub'), s['psacwi'][s['k'] - 1], 'sub'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1046)
            def yes_1047(s):
                s = dict(s)
                old_1048 = s['ratio']
                s['ratio'] = F(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_div(s['qi3d'][s['k'] - 1], s['dt']), s['prd'][s['k'] - 1], 'add'), s['mnuccc'][s['k'] - 1], 'add'), s['qmults'][s['k'] - 1], 'add'), s['qmultg'][s['k'] - 1], 'add'), s['qmultr'][s['k'] - 1], 'add'), s['qmultrg'][s['k'] - 1], 'add'), s['mnuccd'][s['k'] - 1], 'add'), s['psacwi'][s['k'] - 1], 'add'), _arith(_arith(_arith(_arith(s['prci'][s['k'] - 1], s['prai'][s['k'] - 1], 'add'), s['praci'][s['k'] - 1], 'add'), s['pracis'][s['k'] - 1], 'add'), s['eprd'][s['k'] - 1], 'sub')))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_1048)
                old_1049 = s['prci']
                s['prci'] = s['prci'].at[s['k'] - 1].set(F(_arith(s['prci'][s['k'] - 1], s['ratio'], 'mul')))
                s['prci'] = jnp.where(s["_pc"] == 0, s['prci'], old_1049)
                old_1050 = s['prai']
                s['prai'] = s['prai'].at[s['k'] - 1].set(F(_arith(s['prai'][s['k'] - 1], s['ratio'], 'mul')))
                s['prai'] = jnp.where(s["_pc"] == 0, s['prai'], old_1050)
                old_1051 = s['praci']
                s['praci'] = s['praci'].at[s['k'] - 1].set(F(_arith(s['praci'][s['k'] - 1], s['ratio'], 'mul')))
                s['praci'] = jnp.where(s["_pc"] == 0, s['praci'], old_1051)
                old_1052 = s['pracis']
                s['pracis'] = s['pracis'].at[s['k'] - 1].set(F(_arith(s['pracis'][s['k'] - 1], s['ratio'], 'mul')))
                s['pracis'] = jnp.where(s["_pc"] == 0, s['pracis'], old_1052)
                old_1053 = s['eprd']
                s['eprd'] = s['eprd'].at[s['k'] - 1].set(F(_arith(s['eprd'][s['k'] - 1], s['ratio'], 'mul')))
                s['eprd'] = jnp.where(s["_pc"] == 0, s['eprd'], old_1053)
                return s
            def no_1047(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qi3d'][s['k'] - 1])) & ((s['qi3d'][s['k'] - 1] >= s['qsmall'])))), yes_1047, no_1047, s)
            old_1054 = s['dum']
            s['dum'] = F(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['pracs'][s['k'] - 1], s['pre'][s['k'] - 1], 'sub'), _arith(_arith(s['qmultr'][s['k'] - 1], s['qmultrg'][s['k'] - 1], 'add'), s['prc'][s['k'] - 1], 'sub'), 'add'), _arith(s['mnuccr'][s['k'] - 1], s['pra'][s['k'] - 1], 'sub'), 'add'), s['piacr'][s['k'] - 1], 'add'), s['piacrs'][s['k'] - 1], 'add'), s['pgracs'][s['k'] - 1], 'add'), s['pracg'][s['k'] - 1], 'add'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1054)
            def yes_1055(s):
                s = dict(s)
                old_1056 = s['ratio']
                s['ratio'] = F(_div(_arith(_arith(_div(s['qr3d'][s['k'] - 1], s['dt']), s['prc'][s['k'] - 1], 'add'), s['pra'][s['k'] - 1], 'add'), _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['pre'][s['k'] - 1]), s['qmultr'][s['k'] - 1], 'add'), s['qmultrg'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'add'), s['mnuccr'][s['k'] - 1], 'add'), s['piacr'][s['k'] - 1], 'add'), s['piacrs'][s['k'] - 1], 'add'), s['pgracs'][s['k'] - 1], 'add'), s['pracg'][s['k'] - 1], 'add')))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_1056)
                old_1057 = s['pre']
                s['pre'] = s['pre'].at[s['k'] - 1].set(F(_arith(s['pre'][s['k'] - 1], s['ratio'], 'mul')))
                s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_1057)
                old_1058 = s['pracs']
                s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(_arith(s['pracs'][s['k'] - 1], s['ratio'], 'mul')))
                s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_1058)
                old_1059 = s['qmultr']
                s['qmultr'] = s['qmultr'].at[s['k'] - 1].set(F(_arith(s['qmultr'][s['k'] - 1], s['ratio'], 'mul')))
                s['qmultr'] = jnp.where(s["_pc"] == 0, s['qmultr'], old_1059)
                old_1060 = s['qmultrg']
                s['qmultrg'] = s['qmultrg'].at[s['k'] - 1].set(F(_arith(s['qmultrg'][s['k'] - 1], s['ratio'], 'mul')))
                s['qmultrg'] = jnp.where(s["_pc"] == 0, s['qmultrg'], old_1060)
                old_1061 = s['mnuccr']
                s['mnuccr'] = s['mnuccr'].at[s['k'] - 1].set(F(_arith(s['mnuccr'][s['k'] - 1], s['ratio'], 'mul')))
                s['mnuccr'] = jnp.where(s["_pc"] == 0, s['mnuccr'], old_1061)
                old_1062 = s['piacr']
                s['piacr'] = s['piacr'].at[s['k'] - 1].set(F(_arith(s['piacr'][s['k'] - 1], s['ratio'], 'mul')))
                s['piacr'] = jnp.where(s["_pc"] == 0, s['piacr'], old_1062)
                old_1063 = s['piacrs']
                s['piacrs'] = s['piacrs'].at[s['k'] - 1].set(F(_arith(s['piacrs'][s['k'] - 1], s['ratio'], 'mul')))
                s['piacrs'] = jnp.where(s["_pc"] == 0, s['piacrs'], old_1063)
                old_1064 = s['pgracs']
                s['pgracs'] = s['pgracs'].at[s['k'] - 1].set(F(_arith(s['pgracs'][s['k'] - 1], s['ratio'], 'mul')))
                s['pgracs'] = jnp.where(s["_pc"] == 0, s['pgracs'], old_1064)
                old_1065 = s['pracg']
                s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(_arith(s['pracg'][s['k'] - 1], s['ratio'], 'mul')))
                s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_1065)
                return s
            def no_1055(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qr3d'][s['k'] - 1])) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])))), yes_1055, no_1055, s)
            def yes_1066(s):
                s = dict(s)
                old_1067 = s['dum']
                s['dum'] = F(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['prds'][s['k'] - 1]), s['psacws'][s['k'] - 1], 'sub'), s['prai'][s['k'] - 1], 'sub'), s['prci'][s['k'] - 1], 'sub'), s['pracs'][s['k'] - 1], 'sub'), s['eprds'][s['k'] - 1], 'sub'), s['psacr'][s['k'] - 1], 'add'), s['piacrs'][s['k'] - 1], 'sub'), s['pracis'][s['k'] - 1], 'sub'), s['dt'], 'mul'))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1067)
                def yes_1068(s):
                    s = dict(s)
                    old_1069 = s['ratio']
                    s['ratio'] = F(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_div(s['qni3d'][s['k'] - 1], s['dt']), s['prds'][s['k'] - 1], 'add'), s['psacws'][s['k'] - 1], 'add'), s['prai'][s['k'] - 1], 'add'), s['prci'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'add'), s['piacrs'][s['k'] - 1], 'add'), s['pracis'][s['k'] - 1], 'add'), _arith((-s['eprds'][s['k'] - 1]), s['psacr'][s['k'] - 1], 'add')))
                    s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_1069)
                    old_1070 = s['eprds']
                    s['eprds'] = s['eprds'].at[s['k'] - 1].set(F(_arith(s['eprds'][s['k'] - 1], s['ratio'], 'mul')))
                    s['eprds'] = jnp.where(s["_pc"] == 0, s['eprds'], old_1070)
                    old_1071 = s['psacr']
                    s['psacr'] = s['psacr'].at[s['k'] - 1].set(F(_arith(s['psacr'][s['k'] - 1], s['ratio'], 'mul')))
                    s['psacr'] = jnp.where(s["_pc"] == 0, s['psacr'], old_1071)
                    return s
                def no_1068(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qni3d'][s['k'] - 1])) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])))), yes_1068, no_1068, s)
                return s
            def no_1066(s):
                s = dict(s)
                def yes_1072(s):
                    s = dict(s)
                    old_1073 = s['dum']
                    s['dum'] = F(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['prds'][s['k'] - 1]), s['psacws'][s['k'] - 1], 'sub'), s['prai'][s['k'] - 1], 'sub'), s['prci'][s['k'] - 1], 'sub'), s['pracs'][s['k'] - 1], 'sub'), s['eprds'][s['k'] - 1], 'sub'), s['psacr'][s['k'] - 1], 'add'), s['piacrs'][s['k'] - 1], 'sub'), s['pracis'][s['k'] - 1], 'sub'), s['mnuccr'][s['k'] - 1], 'sub'), s['dt'], 'mul'))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1073)
                    def yes_1074(s):
                        s = dict(s)
                        old_1075 = s['ratio']
                        s['ratio'] = F(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_div(s['qni3d'][s['k'] - 1], s['dt']), s['prds'][s['k'] - 1], 'add'), s['psacws'][s['k'] - 1], 'add'), s['prai'][s['k'] - 1], 'add'), s['prci'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'add'), s['piacrs'][s['k'] - 1], 'add'), s['pracis'][s['k'] - 1], 'add'), s['mnuccr'][s['k'] - 1], 'add'), _arith((-s['eprds'][s['k'] - 1]), s['psacr'][s['k'] - 1], 'add')))
                        s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_1075)
                        old_1076 = s['eprds']
                        s['eprds'] = s['eprds'].at[s['k'] - 1].set(F(_arith(s['eprds'][s['k'] - 1], s['ratio'], 'mul')))
                        s['eprds'] = jnp.where(s["_pc"] == 0, s['eprds'], old_1076)
                        old_1077 = s['psacr']
                        s['psacr'] = s['psacr'].at[s['k'] - 1].set(F(_arith(s['psacr'][s['k'] - 1], s['ratio'], 'mul')))
                        s['psacr'] = jnp.where(s["_pc"] == 0, s['psacr'], old_1077)
                        return s
                    def no_1074(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qni3d'][s['k'] - 1])) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])))), yes_1074, no_1074, s)
                    return s
                def no_1072(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 1)), yes_1072, no_1072, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 0)), yes_1066, no_1066, s)
            old_1078 = s['dum']
            s['dum'] = F(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['psacwg'][s['k'] - 1]), s['pracg'][s['k'] - 1], 'sub'), s['pgsacw'][s['k'] - 1], 'sub'), s['pgracs'][s['k'] - 1], 'sub'), s['prdg'][s['k'] - 1], 'sub'), s['mnuccr'][s['k'] - 1], 'sub'), s['eprdg'][s['k'] - 1], 'sub'), s['piacr'][s['k'] - 1], 'sub'), s['praci'][s['k'] - 1], 'sub'), s['psacr'][s['k'] - 1], 'sub'), s['dt'], 'mul'))
            s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1078)
            def yes_1079(s):
                s = dict(s)
                old_1080 = s['ratio']
                s['ratio'] = F(_div(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_div(s['qg3d'][s['k'] - 1], s['dt']), s['psacwg'][s['k'] - 1], 'add'), s['pracg'][s['k'] - 1], 'add'), s['pgsacw'][s['k'] - 1], 'add'), s['pgracs'][s['k'] - 1], 'add'), s['prdg'][s['k'] - 1], 'add'), s['mnuccr'][s['k'] - 1], 'add'), s['psacr'][s['k'] - 1], 'add'), s['piacr'][s['k'] - 1], 'add'), s['praci'][s['k'] - 1], 'add'), (-s['eprdg'][s['k'] - 1])))
                s['ratio'] = jnp.where(s["_pc"] == 0, s['ratio'], old_1080)
                old_1081 = s['eprdg']
                s['eprdg'] = s['eprdg'].at[s['k'] - 1].set(F(_arith(s['eprdg'][s['k'] - 1], s['ratio'], 'mul')))
                s['eprdg'] = jnp.where(s["_pc"] == 0, s['eprdg'], old_1081)
                return s
            def no_1079(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['dum'] > s['qg3d'][s['k'] - 1])) & ((s['qg3d'][s['k'] - 1] >= s['qsmall'])))), yes_1079, no_1079, s)
            old_1082 = s['qv3dten']
            s['qv3dten'] = s['qv3dten'].at[s['k'] - 1].set(F(_arith(s['qv3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['pre'][s['k'] - 1]), s['prd'][s['k'] - 1], 'sub'), s['prds'][s['k'] - 1], 'sub'), s['mnuccd'][s['k'] - 1], 'sub'), s['eprd'][s['k'] - 1], 'sub'), s['eprds'][s['k'] - 1], 'sub'), s['prdg'][s['k'] - 1], 'sub'), s['eprdg'][s['k'] - 1], 'sub'), 'add')))
            s['qv3dten'] = jnp.where(s["_pc"] == 0, s['qv3dten'], old_1082)
            old_1083 = s['t3dten']
            s['t3dten'] = s['t3dten'].at[s['k'] - 1].set(F(_arith(s['t3dten'][s['k'] - 1], _div(_arith(_arith(_arith(s['pre'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(_arith(_arith(_arith(_arith(_arith(s['prd'][s['k'] - 1], s['prds'][s['k'] - 1], 'add'), s['mnuccd'][s['k'] - 1], 'add'), s['eprd'][s['k'] - 1], 'add'), s['eprds'][s['k'] - 1], 'add'), s['prdg'][s['k'] - 1], 'add'), s['eprdg'][s['k'] - 1], 'add'), s['xxls'][s['k'] - 1], 'mul'), 'add'), _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['psacws'][s['k'] - 1], s['psacwi'][s['k'] - 1], 'add'), s['mnuccc'][s['k'] - 1], 'add'), s['mnuccr'][s['k'] - 1], 'add'), s['qmults'][s['k'] - 1], 'add'), s['qmultg'][s['k'] - 1], 'add'), s['qmultr'][s['k'] - 1], 'add'), s['qmultrg'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'add'), s['psacwg'][s['k'] - 1], 'add'), s['pracg'][s['k'] - 1], 'add'), s['pgsacw'][s['k'] - 1], 'add'), s['pgracs'][s['k'] - 1], 'add'), s['piacr'][s['k'] - 1], 'add'), s['piacrs'][s['k'] - 1], 'add'), s['xlf'][s['k'] - 1], 'mul'), 'add'), s['cpm'][s['k'] - 1]), 'add')))
            s['t3dten'] = jnp.where(s["_pc"] == 0, s['t3dten'], old_1083)
            def yes_1084(s):
                s = dict(s)
                lax.cond(s["_pc"] == 0, lambda _: jax.debug.callback(_fatal, s["_pc"], ordered=True), lambda _: None, None)
                return s
            def no_1084(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['pcc'][s['k'] - 1] != F(0.0))), yes_1084, no_1084, s)
            old_1085 = s['qc3dten']
            s['qc3dten'] = s['qc3dten'].at[s['k'] - 1].set(F(_arith(s['qc3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith((-s['pra'][s['k'] - 1]), s['prc'][s['k'] - 1], 'sub'), s['mnuccc'][s['k'] - 1], 'sub'), s['pcc'][s['k'] - 1], 'add'), s['psacws'][s['k'] - 1], 'sub'), s['psacwi'][s['k'] - 1], 'sub'), s['qmults'][s['k'] - 1], 'sub'), s['qmultg'][s['k'] - 1], 'sub'), s['psacwg'][s['k'] - 1], 'sub'), s['pgsacw'][s['k'] - 1], 'sub'), 'add')))
            s['qc3dten'] = jnp.where(s["_pc"] == 0, s['qc3dten'], old_1085)
            old_1086 = s['qi3dten']
            s['qi3dten'] = s['qi3dten'].at[s['k'] - 1].set(F(_arith(s['qi3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['prd'][s['k'] - 1], s['eprd'][s['k'] - 1], 'add'), s['psacwi'][s['k'] - 1], 'add'), s['mnuccc'][s['k'] - 1], 'add'), s['prci'][s['k'] - 1], 'sub'), s['prai'][s['k'] - 1], 'sub'), s['qmults'][s['k'] - 1], 'add'), s['qmultg'][s['k'] - 1], 'add'), s['qmultr'][s['k'] - 1], 'add'), s['qmultrg'][s['k'] - 1], 'add'), s['mnuccd'][s['k'] - 1], 'add'), s['praci'][s['k'] - 1], 'sub'), s['pracis'][s['k'] - 1], 'sub'), 'add')))
            s['qi3dten'] = jnp.where(s["_pc"] == 0, s['qi3dten'], old_1086)
            old_1087 = s['qr3dten']
            s['qr3dten'] = s['qr3dten'].at[s['k'] - 1].set(F(_arith(s['qr3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['pre'][s['k'] - 1], s['pra'][s['k'] - 1], 'add'), s['prc'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'sub'), s['mnuccr'][s['k'] - 1], 'sub'), s['qmultr'][s['k'] - 1], 'sub'), s['qmultrg'][s['k'] - 1], 'sub'), s['piacr'][s['k'] - 1], 'sub'), s['piacrs'][s['k'] - 1], 'sub'), s['pracg'][s['k'] - 1], 'sub'), s['pgracs'][s['k'] - 1], 'sub'), 'add')))
            s['qr3dten'] = jnp.where(s["_pc"] == 0, s['qr3dten'], old_1087)
            def yes_1088(s):
                s = dict(s)
                old_1089 = s['qni3dten']
                s['qni3dten'] = s['qni3dten'].at[s['k'] - 1].set(F(_arith(s['qni3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['prai'][s['k'] - 1], s['psacws'][s['k'] - 1], 'add'), s['prds'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'add'), s['prci'][s['k'] - 1], 'add'), s['eprds'][s['k'] - 1], 'add'), s['psacr'][s['k'] - 1], 'sub'), s['piacrs'][s['k'] - 1], 'add'), s['pracis'][s['k'] - 1], 'add'), 'add')))
                s['qni3dten'] = jnp.where(s["_pc"] == 0, s['qni3dten'], old_1089)
                old_1090 = s['ns3dten']
                s['ns3dten'] = s['ns3dten'].at[s['k'] - 1].set(F(_arith(s['ns3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(s['nsagg'][s['k'] - 1], s['nprci'][s['k'] - 1], 'add'), s['nscng'][s['k'] - 1], 'sub'), s['ngracs'][s['k'] - 1], 'sub'), s['niacrs'][s['k'] - 1], 'add'), 'add')))
                s['ns3dten'] = jnp.where(s["_pc"] == 0, s['ns3dten'], old_1090)
                old_1091 = s['qg3dten']
                s['qg3dten'] = s['qg3dten'].at[s['k'] - 1].set(F(_arith(s['qg3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['pracg'][s['k'] - 1], s['psacwg'][s['k'] - 1], 'add'), s['pgsacw'][s['k'] - 1], 'add'), s['pgracs'][s['k'] - 1], 'add'), s['prdg'][s['k'] - 1], 'add'), s['eprdg'][s['k'] - 1], 'add'), s['mnuccr'][s['k'] - 1], 'add'), s['piacr'][s['k'] - 1], 'add'), s['praci'][s['k'] - 1], 'add'), s['psacr'][s['k'] - 1], 'add'), 'add')))
                s['qg3dten'] = jnp.where(s["_pc"] == 0, s['qg3dten'], old_1091)
                old_1092 = s['ng3dten']
                s['ng3dten'] = s['ng3dten'].at[s['k'] - 1].set(F(_arith(s['ng3dten'][s['k'] - 1], _arith(_arith(_arith(s['nscng'][s['k'] - 1], s['ngracs'][s['k'] - 1], 'add'), s['nnuccr'][s['k'] - 1], 'add'), s['niacr'][s['k'] - 1], 'add'), 'add')))
                s['ng3dten'] = jnp.where(s["_pc"] == 0, s['ng3dten'], old_1092)
                return s
            def no_1088(s):
                s = dict(s)
                def yes_1093(s):
                    s = dict(s)
                    old_1094 = s['qni3dten']
                    s['qni3dten'] = s['qni3dten'].at[s['k'] - 1].set(F(_arith(s['qni3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['prai'][s['k'] - 1], s['psacws'][s['k'] - 1], 'add'), s['prds'][s['k'] - 1], 'add'), s['pracs'][s['k'] - 1], 'add'), s['prci'][s['k'] - 1], 'add'), s['eprds'][s['k'] - 1], 'add'), s['psacr'][s['k'] - 1], 'sub'), s['piacrs'][s['k'] - 1], 'add'), s['pracis'][s['k'] - 1], 'add'), s['mnuccr'][s['k'] - 1], 'add'), 'add')))
                    s['qni3dten'] = jnp.where(s["_pc"] == 0, s['qni3dten'], old_1094)
                    old_1095 = s['ns3dten']
                    s['ns3dten'] = s['ns3dten'].at[s['k'] - 1].set(F(_arith(s['ns3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(s['nsagg'][s['k'] - 1], s['nprci'][s['k'] - 1], 'add'), s['nscng'][s['k'] - 1], 'sub'), s['ngracs'][s['k'] - 1], 'sub'), s['niacrs'][s['k'] - 1], 'add'), s['nnuccr'][s['k'] - 1], 'add'), 'add')))
                    s['ns3dten'] = jnp.where(s["_pc"] == 0, s['ns3dten'], old_1095)
                    return s
                def no_1093(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 1)), yes_1093, no_1093, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 0)), yes_1088, no_1088, s)
            old_1096 = s['nc3dten']
            s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith((-s['nnuccc'][s['k'] - 1]), s['npsacws'][s['k'] - 1], 'sub'), s['npra'][s['k'] - 1], 'sub'), s['nprc'][s['k'] - 1], 'sub'), s['npsacwi'][s['k'] - 1], 'sub'), s['npsacwg'][s['k'] - 1], 'sub'), 'add')))
            s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1096)
            old_1097 = s['ni3dten']
            s['ni3dten'] = s['ni3dten'].at[s['k'] - 1].set(F(_arith(s['ni3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(_arith(s['nnuccc'][s['k'] - 1], s['nprci'][s['k'] - 1], 'sub'), s['nprai'][s['k'] - 1], 'sub'), s['nmults'][s['k'] - 1], 'add'), s['nmultg'][s['k'] - 1], 'add'), s['nmultr'][s['k'] - 1], 'add'), s['nmultrg'][s['k'] - 1], 'add'), s['nnuccd'][s['k'] - 1], 'add'), s['niacr'][s['k'] - 1], 'sub'), s['niacrs'][s['k'] - 1], 'sub'), 'add')))
            s['ni3dten'] = jnp.where(s["_pc"] == 0, s['ni3dten'], old_1097)
            old_1098 = s['nr3dten']
            s['nr3dten'] = s['nr3dten'].at[s['k'] - 1].set(F(_arith(s['nr3dten'][s['k'] - 1], _arith(_arith(_arith(_arith(_arith(_arith(_arith(s['nprc1'][s['k'] - 1], s['npracs'][s['k'] - 1], 'sub'), s['nnuccr'][s['k'] - 1], 'sub'), s['nragg'][s['k'] - 1], 'add'), s['niacr'][s['k'] - 1], 'sub'), s['niacrs'][s['k'] - 1], 'sub'), s['npracg'][s['k'] - 1], 'sub'), s['ngracs'][s['k'] - 1], 'sub'), 'add')))
            s['nr3dten'] = jnp.where(s["_pc"] == 0, s['nr3dten'], old_1098)
            def yes_1099(s):
                s = dict(s)
                old_1100 = s['dumt']
                s['dumt'] = F(_arith(s['t3d'][s['k'] - 1], _arith(s['dt'], s['t3dten'][s['k'] - 1], 'mul'), 'add'))
                s['dumt'] = jnp.where(s["_pc"] == 0, s['dumt'], old_1100)
                old_1101 = s['dumqv']
                s['dumqv'] = F(_arith(s['qv3d'][s['k'] - 1], _arith(s['dt'], s['qv3dten'][s['k'] - 1], 'mul'), 'add'))
                s['dumqv'] = jnp.where(s["_pc"] == 0, s['dumqv'], old_1101)
                old_1102 = s['dum']
                s['dum'] = F(jnp.minimum(_arith(F(0.99), s['pres'][s['k'] - 1], 'mul'), POLYSVP(s['dumt'], 0)))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1102)
                old_1103 = s['dumqss']
                s['dumqss'] = F(_div(_arith(s['ep_2'], s['dum'], 'mul'), _arith(s['pres'][s['k'] - 1], s['dum'], 'sub')))
                s['dumqss'] = jnp.where(s["_pc"] == 0, s['dumqss'], old_1103)
                old_1104 = s['dumqc']
                s['dumqc'] = F(_arith(s['qc3d'][s['k'] - 1], _arith(s['dt'], s['qc3dten'][s['k'] - 1], 'mul'), 'add'))
                s['dumqc'] = jnp.where(s["_pc"] == 0, s['dumqc'], old_1104)
                old_1105 = s['dumqc']
                s['dumqc'] = F(jnp.maximum(s['dumqc'], F(0.0)))
                s['dumqc'] = jnp.where(s["_pc"] == 0, s['dumqc'], old_1105)
                old_1106 = s['dums']
                s['dums'] = F(_arith(s['dumqv'], s['dumqss'], 'sub'))
                s['dums'] = jnp.where(s["_pc"] == 0, s['dums'], old_1106)
                old_1107 = s['pcc']
                s['pcc'] = s['pcc'].at[s['k'] - 1].set(F(_div(_div(s['dums'], _arith(F(1.0), _div(_arith(_arith(s['xxlv'][s['k'] - 1], 2, 'pow'), s['dumqss'], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rv'], 'mul'), _arith(s['dumt'], 2, 'pow'), 'mul')), 'add')), s['dt'])))
                s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_1107)
                def yes_1108(s):
                    s = dict(s)
                    old_1109 = s['pcc']
                    s['pcc'] = s['pcc'].at[s['k'] - 1].set(F(_div((-_arith(s['qc3d'][s['k'] - 1], _arith(s['dt'], s['qc3dten'][s['k'] - 1], 'mul'), 'add')), s['dt'])))
                    s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_1109)
                    return s
                def no_1108(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((_arith(_arith(_arith(s['pcc'][s['k'] - 1], s['dt'], 'mul'), s['qc3d'][s['k'] - 1], 'add'), _arith(s['dt'], s['qc3dten'][s['k'] - 1], 'mul'), 'add') < F(0.0))), yes_1108, no_1108, s)
                old_1110 = s['qv3dten']
                s['qv3dten'] = s['qv3dten'].at[s['k'] - 1].set(F(_arith(s['qv3dten'][s['k'] - 1], s['pcc'][s['k'] - 1], 'sub')))
                s['qv3dten'] = jnp.where(s["_pc"] == 0, s['qv3dten'], old_1110)
                old_1111 = s['t3dten']
                s['t3dten'] = s['t3dten'].at[s['k'] - 1].set(F(_arith(s['t3dten'][s['k'] - 1], _div(_arith(s['pcc'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'add')))
                s['t3dten'] = jnp.where(s["_pc"] == 0, s['t3dten'], old_1111)
                old_1112 = s['qc3dten']
                s['qc3dten'] = s['qc3dten'].at[s['k'] - 1].set(F(_arith(s['qc3dten'][s['k'] - 1], s['pcc'][s['k'] - 1], 'add')))
                s['qc3dten'] = jnp.where(s["_pc"] == 0, s['qc3dten'], old_1112)
                return s
            def no_1099(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['isatadj'] == 0)), yes_1099, no_1099, s)
            def yes_1113(s):
                s = dict(s)
                def yes_1114(s):
                    s = dict(s)
                    old_1115 = s['dum']
                    s['dum'] = F(_arith(s['w3d'][s['k'] - 1], s['wvar'][s['k'] - 1], 'add'))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1115)
                    old_1116 = s['dum']
                    s['dum'] = F(jnp.maximum(s['dum'], F(0.1)))
                    s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1116)
                    return s
                def no_1114(s):
                    s = dict(s)
                    def yes_1117(s):
                        s = dict(s)
                        old_1118 = s['dum']
                        s['dum'] = F(s['w3d'][s['k'] - 1])
                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1118)
                        return s
                    def no_1117(s):
                        s = dict(s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['isub'] == 1)), yes_1117, no_1117, s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['isub'] == 0)), yes_1114, no_1114, s)
                def yes_1119(s):
                    s = dict(s)
                    def yes_1120(s):
                        s = dict(s)
                        old_1121 = s['idrop']
                        s['idrop'] = I(0)
                        s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_1121)
                        def yes_1122(s):
                            s = dict(s)
                            old_1123 = s['idrop']
                            s['idrop'] = I(1)
                            s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_1123)
                            return s
                        def no_1122(s):
                            s = dict(s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] <= _div(F(5e-05), s['rho'][s['k'] - 1]))), yes_1122, no_1122, s)
                        def yes_1124(s):
                            s = dict(s)
                            old_1125 = s['idrop']
                            s['idrop'] = I(1)
                            s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_1125)
                            return s
                        def no_1124(s):
                            s = dict(s)
                            def yes_1126(s):
                                s = dict(s)
                                def yes_1127(s):
                                    s = dict(s)
                                    old_1128 = s['idrop']
                                    s['idrop'] = I(1)
                                    s['idrop'] = jnp.where(s["_pc"] == 0, s['idrop'], old_1128)
                                    return s
                                def no_1127(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] > _div(F(5e-05), s['rho'][s['k'] - 1]))) & ((s['qc3d'][_arith(s['k'], 1, 'sub') - 1] <= _div(F(5e-05), s['rho'][_arith(s['k'], 1, 'sub') - 1]))))), yes_1127, no_1127, s)
                                return s
                            def no_1126(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['k'] >= 2)), yes_1126, no_1126, s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['k'] == 1)), yes_1124, no_1124, s)
                        def yes_1129(s):
                            s = dict(s)
                            def yes_1130(s):
                                s = dict(s)
                                old_1131 = s['dum']
                                s['dum'] = F(_arith(s['dum'], F(100.0), 'mul'))
                                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1131)
                                old_1132 = s['dum2']
                                s['dum2'] = F(_arith(_arith(F(0.88), _arith(s['c1'], _div(F(2.0), _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'), _arith(_arith(F(0.07), _arith(s['dum'], F(1.5), 'pow'), 'mul'), _div(s['k1'], _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1132)
                                old_1133 = s['dum2']
                                s['dum2'] = F(_arith(s['dum2'], F(1000000.0), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1133)
                                old_1134 = s['dum2']
                                s['dum2'] = F(_div(s['dum2'], s['rho'][s['k'] - 1]))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1134)
                                old_1135 = s['dum2']
                                s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1135)
                                old_1136 = s['dum2']
                                s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1136)
                                old_1137 = s['nc3dten']
                                s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1137)
                                old_1138 = s['nact']
                                s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1138)
                                return s
                            def no_1130(s):
                                s = dict(s)
                                def yes_1139(s):
                                    s = dict(s)
                                    old_1140 = s['sigvl']
                                    s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                    s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_1140)
                                    old_1141 = s['aact']
                                    s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                    s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_1141)
                                    old_1142 = s['alpha']
                                    s['alpha'] = F(_arith(_div(_arith(_arith(s['g'], s['mw'], 'mul'), s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rr'], 'mul'), _arith(s['t3d'][s['k'] - 1], 2, 'pow'), 'mul')), _div(_arith(s['g'], s['ma'], 'mul'), _arith(s['rr'], s['t3d'][s['k'] - 1], 'mul')), 'sub'))
                                    s['alpha'] = jnp.where(s["_pc"] == 0, s['alpha'], old_1142)
                                    old_1143 = s['gamm']
                                    s['gamm'] = F(_arith(_div(_arith(s['rr'], s['t3d'][s['k'] - 1], 'mul'), _arith(s['evs'][s['k'] - 1], s['mw'], 'mul')), _div(_arith(s['mw'], _arith(s['xxlv'][s['k'] - 1], 2, 'pow'), 'mul'), _arith(_arith(_arith(s['cpm'][s['k'] - 1], s['pres'][s['k'] - 1], 'mul'), s['ma'], 'mul'), s['t3d'][s['k'] - 1], 'mul')), 'add'))
                                    s['gamm'] = jnp.where(s["_pc"] == 0, s['gamm'], old_1143)
                                    old_1144 = s['gg']
                                    s['gg'] = F(_div(F(1.0), _arith(_div(_arith(_arith(s['rhow'], s['rr'], 'mul'), s['t3d'][s['k'] - 1], 'mul'), _arith(_arith(s['evs'][s['k'] - 1], s['dv'][s['k'] - 1], 'mul'), s['mw'], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['rhow'], 'mul'), _arith(s['kap'][s['k'] - 1], s['t3d'][s['k'] - 1], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['mw'], 'mul'), _arith(s['t3d'][s['k'] - 1], s['rr'], 'mul')), F(1.0), 'sub'), 'mul'), 'add')))
                                    s['gg'] = jnp.where(s["_pc"] == 0, s['gg'], old_1144)
                                    old_1145 = s['psi']
                                    s['psi'] = F(_arith(_arith(_div(F(2.0), F(3.0)), _arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(0.5), 'pow'), 'mul'), s['aact'], 'mul'))
                                    s['psi'] = jnp.where(s["_pc"] == 0, s['psi'], old_1145)
                                    old_1146 = s['eta1']
                                    s['eta1'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew1'], 'mul')))
                                    s['eta1'] = jnp.where(s["_pc"] == 0, s['eta1'], old_1146)
                                    old_1147 = s['eta2']
                                    s['eta2'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew2'], 'mul')))
                                    s['eta2'] = jnp.where(s["_pc"] == 0, s['eta2'], old_1147)
                                    old_1148 = s['sm1']
                                    s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_1148)
                                    old_1149 = s['sm2']
                                    s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_1149)
                                    old_1150 = s['dum1']
                                    s['dum1'] = F(_arith(_div(F(1.0), _arith(s['sm1'], 2, 'pow')), _arith(_arith(s['f11'], _arith(_div(s['psi'], s['eta1']), F(1.5), 'pow'), 'mul'), _arith(s['f21'], _arith(_div(_arith(s['sm1'], 2, 'pow'), _arith(s['eta1'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_1150)
                                    old_1151 = s['dum2']
                                    s['dum2'] = F(_arith(_div(F(1.0), _arith(s['sm2'], 2, 'pow')), _arith(_arith(s['f12'], _arith(_div(s['psi'], s['eta2']), F(1.5), 'pow'), 'mul'), _arith(s['f22'], _arith(_div(_arith(s['sm2'], 2, 'pow'), _arith(s['eta2'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1151)
                                    old_1152 = s['smax']
                                    s['smax'] = F(_div(F(1.0), _arith(_arith(s['dum1'], s['dum2'], 'add'), F(0.5), 'pow')))
                                    s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_1152)
                                    old_1153 = s['uu1']
                                    s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                    s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_1153)
                                    old_1154 = s['uu2']
                                    s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                    s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_1154)
                                    old_1155 = s['dum1']
                                    s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_1155)
                                    old_1156 = s['dum2']
                                    s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1156)
                                    old_1157 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1157)
                                    old_1158 = s['dum2']
                                    s['dum2'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1158)
                                    old_1159 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1159)
                                    old_1160 = s['dum2']
                                    s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1160)
                                    old_1161 = s['nc3dten']
                                    s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1161)
                                    old_1162 = s['nact']
                                    s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1162)
                                    return s
                                def no_1139(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 2)), yes_1139, no_1139, s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 1)), yes_1130, no_1130, s)
                            return s
                        def no_1129(s):
                            s = dict(s)
                            def yes_1163(s):
                                s = dict(s)
                                old_1164 = s['tauc']
                                s['tauc'] = F(_div(F(1.0), _div(_arith(_arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rho'][s['k'] - 1], 'mul'), s['dv'][s['k'] - 1], 'mul'), s['nc3d'][s['k'] - 1], 'mul'), _arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), 'mul'), s['lamc'][s['k'] - 1])))
                                s['tauc'] = jnp.where(s["_pc"] == 0, s['tauc'], old_1164)
                                def yes_1165(s):
                                    s = dict(s)
                                    old_1166 = s['taur']
                                    s['taur'] = F(_div(F(1.0), s['epsr']))
                                    s['taur'] = jnp.where(s["_pc"] == 0, s['taur'], old_1166)
                                    return s
                                def no_1165(s):
                                    s = dict(s)
                                    old_1167 = s['taur']
                                    s['taur'] = F(F(100000000.0))
                                    s['taur'] = jnp.where(s["_pc"] == 0, s['taur'], old_1167)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['epsr'] > F(1e-08))), yes_1165, no_1165, s)
                                def yes_1168(s):
                                    s = dict(s)
                                    old_1169 = s['taui']
                                    s['taui'] = F(_div(F(1.0), s['epsi']))
                                    s['taui'] = jnp.where(s["_pc"] == 0, s['taui'], old_1169)
                                    return s
                                def no_1168(s):
                                    s = dict(s)
                                    old_1170 = s['taui']
                                    s['taui'] = F(F(100000000.0))
                                    s['taui'] = jnp.where(s["_pc"] == 0, s['taui'], old_1170)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['epsi'] > F(1e-08))), yes_1168, no_1168, s)
                                def yes_1171(s):
                                    s = dict(s)
                                    old_1172 = s['taus']
                                    s['taus'] = F(_div(F(1.0), s['epss']))
                                    s['taus'] = jnp.where(s["_pc"] == 0, s['taus'], old_1172)
                                    return s
                                def no_1171(s):
                                    s = dict(s)
                                    old_1173 = s['taus']
                                    s['taus'] = F(F(100000000.0))
                                    s['taus'] = jnp.where(s["_pc"] == 0, s['taus'], old_1173)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['epss'] > F(1e-08))), yes_1171, no_1171, s)
                                def yes_1174(s):
                                    s = dict(s)
                                    old_1175 = s['taug']
                                    s['taug'] = F(_div(F(1.0), s['epsg']))
                                    s['taug'] = jnp.where(s["_pc"] == 0, s['taug'], old_1175)
                                    return s
                                def no_1174(s):
                                    s = dict(s)
                                    old_1176 = s['taug']
                                    s['taug'] = F(F(100000000.0))
                                    s['taug'] = jnp.where(s["_pc"] == 0, s['taug'], old_1176)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['epsg'] > F(1e-08))), yes_1174, no_1174, s)
                                old_1177 = s['dum3']
                                s['dum3'] = F(_arith(_arith(_arith(_div(_arith(s['qvs'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul'), _arith(s['pres'][s['k'] - 1], s['evs'][s['k'] - 1], 'sub')), _div(s['dqsdt'], s['cp']), 'add'), s['g'], 'mul'), s['dum'], 'mul'))
                                s['dum3'] = jnp.where(s["_pc"] == 0, s['dum3'], old_1177)
                                old_1178 = s['dum3']
                                s['dum3'] = F(_div(_arith(_arith(_arith(_arith(_arith(_arith(s['dum3'], s['tauc'], 'mul'), s['taur'], 'mul'), s['taui'], 'mul'), s['taus'], 'mul'), s['taug'], 'mul'), _arith(_arith(s['qvs'][s['k'] - 1], s['qvi'][s['k'] - 1], 'sub'), _arith(_arith(_arith(_arith(_arith(s['tauc'], s['taur'], 'mul'), s['taui'], 'mul'), s['taug'], 'mul'), _arith(_arith(_arith(s['tauc'], s['taur'], 'mul'), s['taus'], 'mul'), s['taug'], 'mul'), 'add'), _arith(_arith(_arith(s['tauc'], s['taur'], 'mul'), s['taui'], 'mul'), s['taus'], 'mul'), 'add'), 'mul'), 'sub'), _arith(_arith(_arith(_arith(_arith(_arith(_arith(s['tauc'], s['taur'], 'mul'), s['taui'], 'mul'), s['taug'], 'mul'), _arith(_arith(_arith(s['tauc'], s['taur'], 'mul'), s['taus'], 'mul'), s['taug'], 'mul'), 'add'), _arith(_arith(_arith(s['tauc'], s['taur'], 'mul'), s['taui'], 'mul'), s['taus'], 'mul'), 'add'), _arith(_arith(_arith(s['taur'], s['taui'], 'mul'), s['taus'], 'mul'), s['taug'], 'mul'), 'add'), _arith(_arith(_arith(s['tauc'], s['taui'], 'mul'), s['taus'], 'mul'), s['taug'], 'mul'), 'add')))
                                s['dum3'] = jnp.where(s["_pc"] == 0, s['dum3'], old_1178)
                                def yes_1179(s):
                                    s = dict(s)
                                    def yes_1180(s):
                                        s = dict(s)
                                        old_1181 = s['dum']
                                        s['dum'] = F(_arith(s['dum'], F(100.0), 'mul'))
                                        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1181)
                                        old_1182 = s['dumact']
                                        s['dumact'] = F(_arith(_arith(F(0.88), _arith(s['c1'], _div(F(2.0), _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'), _arith(_arith(F(0.07), _arith(s['dum'], F(1.5), 'pow'), 'mul'), _div(s['k1'], _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'))
                                        s['dumact'] = jnp.where(s["_pc"] == 0, s['dumact'], old_1182)
                                        old_1183 = s['dum3']
                                        s['dum3'] = F(_arith(_div(s['dum3'], s['qvs'][s['k'] - 1]), F(100.0), 'mul'))
                                        s['dum3'] = jnp.where(s["_pc"] == 0, s['dum3'], old_1183)
                                        old_1184 = s['dum2']
                                        s['dum2'] = F(_arith(s['c1'], _arith(s['dum3'], s['k1'], 'pow'), 'mul'))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1184)
                                        old_1185 = s['dum2']
                                        s['dum2'] = F(jnp.minimum(s['dum2'], s['dumact']))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1185)
                                        old_1186 = s['dum2']
                                        s['dum2'] = F(_arith(s['dum2'], F(1000000.0), 'mul'))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1186)
                                        old_1187 = s['dum2']
                                        s['dum2'] = F(_div(s['dum2'], s['rho'][s['k'] - 1]))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1187)
                                        old_1188 = s['dum2']
                                        s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1188)
                                        old_1189 = s['dum2']
                                        s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                        s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1189)
                                        old_1190 = s['nc3dten']
                                        s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                        s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1190)
                                        old_1191 = s['nact']
                                        s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                        s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1191)
                                        return s
                                    def no_1180(s):
                                        s = dict(s)
                                        def yes_1192(s):
                                            s = dict(s)
                                            old_1193 = s['sigvl']
                                            s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                            s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_1193)
                                            old_1194 = s['aact']
                                            s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                            s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_1194)
                                            old_1195 = s['alpha']
                                            s['alpha'] = F(_arith(_div(_arith(_arith(s['g'], s['mw'], 'mul'), s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rr'], 'mul'), _arith(s['t3d'][s['k'] - 1], 2, 'pow'), 'mul')), _div(_arith(s['g'], s['ma'], 'mul'), _arith(s['rr'], s['t3d'][s['k'] - 1], 'mul')), 'sub'))
                                            s['alpha'] = jnp.where(s["_pc"] == 0, s['alpha'], old_1195)
                                            old_1196 = s['gamm']
                                            s['gamm'] = F(_arith(_div(_arith(s['rr'], s['t3d'][s['k'] - 1], 'mul'), _arith(s['evs'][s['k'] - 1], s['mw'], 'mul')), _div(_arith(s['mw'], _arith(s['xxlv'][s['k'] - 1], 2, 'pow'), 'mul'), _arith(_arith(_arith(s['cpm'][s['k'] - 1], s['pres'][s['k'] - 1], 'mul'), s['ma'], 'mul'), s['t3d'][s['k'] - 1], 'mul')), 'add'))
                                            s['gamm'] = jnp.where(s["_pc"] == 0, s['gamm'], old_1196)
                                            old_1197 = s['gg']
                                            s['gg'] = F(_div(F(1.0), _arith(_div(_arith(_arith(s['rhow'], s['rr'], 'mul'), s['t3d'][s['k'] - 1], 'mul'), _arith(_arith(s['evs'][s['k'] - 1], s['dv'][s['k'] - 1], 'mul'), s['mw'], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['rhow'], 'mul'), _arith(s['kap'][s['k'] - 1], s['t3d'][s['k'] - 1], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['mw'], 'mul'), _arith(s['t3d'][s['k'] - 1], s['rr'], 'mul')), F(1.0), 'sub'), 'mul'), 'add')))
                                            s['gg'] = jnp.where(s["_pc"] == 0, s['gg'], old_1197)
                                            old_1198 = s['psi']
                                            s['psi'] = F(_arith(_arith(_div(F(2.0), F(3.0)), _arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(0.5), 'pow'), 'mul'), s['aact'], 'mul'))
                                            s['psi'] = jnp.where(s["_pc"] == 0, s['psi'], old_1198)
                                            old_1199 = s['eta1']
                                            s['eta1'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew1'], 'mul')))
                                            s['eta1'] = jnp.where(s["_pc"] == 0, s['eta1'], old_1199)
                                            old_1200 = s['eta2']
                                            s['eta2'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew2'], 'mul')))
                                            s['eta2'] = jnp.where(s["_pc"] == 0, s['eta2'], old_1200)
                                            old_1201 = s['sm1']
                                            s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_1201)
                                            old_1202 = s['sm2']
                                            s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_1202)
                                            old_1203 = s['dum1']
                                            s['dum1'] = F(_arith(_div(F(1.0), _arith(s['sm1'], 2, 'pow')), _arith(_arith(s['f11'], _arith(_div(s['psi'], s['eta1']), F(1.5), 'pow'), 'mul'), _arith(s['f21'], _arith(_div(_arith(s['sm1'], 2, 'pow'), _arith(s['eta1'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                            s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_1203)
                                            old_1204 = s['dum2']
                                            s['dum2'] = F(_arith(_div(F(1.0), _arith(s['sm2'], 2, 'pow')), _arith(_arith(s['f12'], _arith(_div(s['psi'], s['eta2']), F(1.5), 'pow'), 'mul'), _arith(s['f22'], _arith(_div(_arith(s['sm2'], 2, 'pow'), _arith(s['eta2'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1204)
                                            old_1205 = s['smax']
                                            s['smax'] = F(_div(F(1.0), _arith(_arith(s['dum1'], s['dum2'], 'add'), F(0.5), 'pow')))
                                            s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_1205)
                                            old_1206 = s['uu1']
                                            s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                            s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_1206)
                                            old_1207 = s['uu2']
                                            s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                            s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_1207)
                                            old_1208 = s['dum1']
                                            s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                            s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_1208)
                                            old_1209 = s['dum2']
                                            s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1209)
                                            old_1210 = s['dum2']
                                            s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1210)
                                            old_1211 = s['dumact']
                                            s['dumact'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                            s['dumact'] = jnp.where(s["_pc"] == 0, s['dumact'], old_1211)
                                            old_1212 = s['sigvl']
                                            s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                            s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_1212)
                                            old_1213 = s['aact']
                                            s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                            s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_1213)
                                            old_1214 = s['sm1']
                                            s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_1214)
                                            old_1215 = s['sm2']
                                            s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                            s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_1215)
                                            old_1216 = s['smax']
                                            s['smax'] = F(_div(s['dum3'], s['qvs'][s['k'] - 1]))
                                            s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_1216)
                                            old_1217 = s['uu1']
                                            s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                            s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_1217)
                                            old_1218 = s['uu2']
                                            s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                            s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_1218)
                                            old_1219 = s['dum1']
                                            s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                            s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_1219)
                                            old_1220 = s['dum2']
                                            s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1220)
                                            old_1221 = s['dum2']
                                            s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1221)
                                            old_1222 = s['dum2']
                                            s['dum2'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1222)
                                            old_1223 = s['dum2']
                                            s['dum2'] = F(jnp.minimum(s['dum2'], s['dumact']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1223)
                                            old_1224 = s['dum2']
                                            s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1224)
                                            old_1225 = s['dum2']
                                            s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                            s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1225)
                                            old_1226 = s['nc3dten']
                                            s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                            s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1226)
                                            old_1227 = s['nact']
                                            s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                            s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1227)
                                            return s
                                        def no_1192(s):
                                            s = dict(s)
                                            return s
                                        s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 2)), yes_1192, no_1192, s)
                                        return s
                                    s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 1)), yes_1180, no_1180, s)
                                    return s
                                def no_1179(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((_div(s['dum3'], s['qvs'][s['k'] - 1]) >= F(1e-06))), yes_1179, no_1179, s)
                                return s
                            def no_1163(s):
                                s = dict(s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['idrop'] == 0)), yes_1163, no_1163, s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['idrop'] == 1)), yes_1129, no_1129, s)
                        return s
                    def no_1120(s):
                        s = dict(s)
                        def yes_1228(s):
                            s = dict(s)
                            def yes_1229(s):
                                s = dict(s)
                                old_1230 = s['dum']
                                s['dum'] = F(_arith(s['dum'], F(100.0), 'mul'))
                                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1230)
                                old_1231 = s['dum2']
                                s['dum2'] = F(_arith(_arith(F(0.88), _arith(s['c1'], _div(F(2.0), _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'), _arith(_arith(F(0.07), _arith(s['dum'], F(1.5), 'pow'), 'mul'), _div(s['k1'], _arith(s['k1'], F(2.0), 'add')), 'pow'), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1231)
                                old_1232 = s['dum2']
                                s['dum2'] = F(_arith(s['dum2'], F(1000000.0), 'mul'))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1232)
                                old_1233 = s['dum2']
                                s['dum2'] = F(_div(s['dum2'], s['rho'][s['k'] - 1]))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1233)
                                old_1234 = s['dum2']
                                s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1234)
                                old_1235 = s['dum2']
                                s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1235)
                                old_1236 = s['nc3dten']
                                s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1236)
                                old_1237 = s['nact']
                                s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1237)
                                return s
                            def no_1229(s):
                                s = dict(s)
                                def yes_1238(s):
                                    s = dict(s)
                                    old_1239 = s['sigvl']
                                    s['sigvl'] = F(_arith(F(0.0761), _arith(F(0.000155), _arith(s['t3d'][s['k'] - 1], s['tmelt'], 'sub'), 'mul'), 'sub'))
                                    s['sigvl'] = jnp.where(s["_pc"] == 0, s['sigvl'], old_1239)
                                    old_1240 = s['aact']
                                    s['aact'] = F(_div(_arith(_div(_arith(F(2.0), s['mw'], 'mul'), _arith(s['rhow'], s['rr'], 'mul')), s['sigvl'], 'mul'), s['t3d'][s['k'] - 1]))
                                    s['aact'] = jnp.where(s["_pc"] == 0, s['aact'], old_1240)
                                    old_1241 = s['alpha']
                                    s['alpha'] = F(_arith(_div(_arith(_arith(s['g'], s['mw'], 'mul'), s['xxlv'][s['k'] - 1], 'mul'), _arith(_arith(s['cpm'][s['k'] - 1], s['rr'], 'mul'), _arith(s['t3d'][s['k'] - 1], 2, 'pow'), 'mul')), _div(_arith(s['g'], s['ma'], 'mul'), _arith(s['rr'], s['t3d'][s['k'] - 1], 'mul')), 'sub'))
                                    s['alpha'] = jnp.where(s["_pc"] == 0, s['alpha'], old_1241)
                                    old_1242 = s['gamm']
                                    s['gamm'] = F(_arith(_div(_arith(s['rr'], s['t3d'][s['k'] - 1], 'mul'), _arith(s['evs'][s['k'] - 1], s['mw'], 'mul')), _div(_arith(s['mw'], _arith(s['xxlv'][s['k'] - 1], 2, 'pow'), 'mul'), _arith(_arith(_arith(s['cpm'][s['k'] - 1], s['pres'][s['k'] - 1], 'mul'), s['ma'], 'mul'), s['t3d'][s['k'] - 1], 'mul')), 'add'))
                                    s['gamm'] = jnp.where(s["_pc"] == 0, s['gamm'], old_1242)
                                    old_1243 = s['gg']
                                    s['gg'] = F(_div(F(1.0), _arith(_div(_arith(_arith(s['rhow'], s['rr'], 'mul'), s['t3d'][s['k'] - 1], 'mul'), _arith(_arith(s['evs'][s['k'] - 1], s['dv'][s['k'] - 1], 'mul'), s['mw'], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['rhow'], 'mul'), _arith(s['kap'][s['k'] - 1], s['t3d'][s['k'] - 1], 'mul')), _arith(_div(_arith(s['xxlv'][s['k'] - 1], s['mw'], 'mul'), _arith(s['t3d'][s['k'] - 1], s['rr'], 'mul')), F(1.0), 'sub'), 'mul'), 'add')))
                                    s['gg'] = jnp.where(s["_pc"] == 0, s['gg'], old_1243)
                                    old_1244 = s['psi']
                                    s['psi'] = F(_arith(_arith(_div(F(2.0), F(3.0)), _arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(0.5), 'pow'), 'mul'), s['aact'], 'mul'))
                                    s['psi'] = jnp.where(s["_pc"] == 0, s['psi'], old_1244)
                                    old_1245 = s['eta1']
                                    s['eta1'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew1'], 'mul')))
                                    s['eta1'] = jnp.where(s["_pc"] == 0, s['eta1'], old_1245)
                                    old_1246 = s['eta2']
                                    s['eta2'] = F(_div(_arith(_div(_arith(s['alpha'], s['dum'], 'mul'), s['gg']), F(1.5), 'pow'), _arith(_arith(_arith(_arith(F(2.0), s['pi'], 'mul'), s['rhow'], 'mul'), s['gamm'], 'mul'), s['nanew2'], 'mul')))
                                    s['eta2'] = jnp.where(s["_pc"] == 0, s['eta2'], old_1246)
                                    old_1247 = s['sm1']
                                    s['sm1'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm1'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm1'] = jnp.where(s["_pc"] == 0, s['sm1'], old_1247)
                                    old_1248 = s['sm2']
                                    s['sm2'] = F(_arith(_div(F(2.0), _arith(s['bact'], F(0.5), 'pow')), _arith(_div(s['aact'], _arith(F(3.0), s['rm2'], 'mul')), F(1.5), 'pow'), 'mul'))
                                    s['sm2'] = jnp.where(s["_pc"] == 0, s['sm2'], old_1248)
                                    old_1249 = s['dum1']
                                    s['dum1'] = F(_arith(_div(F(1.0), _arith(s['sm1'], 2, 'pow')), _arith(_arith(s['f11'], _arith(_div(s['psi'], s['eta1']), F(1.5), 'pow'), 'mul'), _arith(s['f21'], _arith(_div(_arith(s['sm1'], 2, 'pow'), _arith(s['eta1'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_1249)
                                    old_1250 = s['dum2']
                                    s['dum2'] = F(_arith(_div(F(1.0), _arith(s['sm2'], 2, 'pow')), _arith(_arith(s['f12'], _arith(_div(s['psi'], s['eta2']), F(1.5), 'pow'), 'mul'), _arith(s['f22'], _arith(_div(_arith(s['sm2'], 2, 'pow'), _arith(s['eta2'], _arith(F(3.0), s['psi'], 'mul'), 'add')), F(0.75), 'pow'), 'mul'), 'add'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1250)
                                    old_1251 = s['smax']
                                    s['smax'] = F(_div(F(1.0), _arith(_arith(s['dum1'], s['dum2'], 'add'), F(0.5), 'pow')))
                                    s['smax'] = jnp.where(s["_pc"] == 0, s['smax'], old_1251)
                                    old_1252 = s['uu1']
                                    s['uu1'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm1'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig1']), 'mul')))
                                    s['uu1'] = jnp.where(s["_pc"] == 0, s['uu1'], old_1252)
                                    old_1253 = s['uu2']
                                    s['uu2'] = F(_div(_arith(F(2.0), _intrinsic('log', _div(s['sm2'], s['smax'])), 'mul'), _arith(F(4.242), _intrinsic('log', s['sig2']), 'mul')))
                                    s['uu2'] = jnp.where(s["_pc"] == 0, s['uu2'], old_1253)
                                    old_1254 = s['dum1']
                                    s['dum1'] = F(_arith(_div(s['nanew1'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu1'])), 'sub'), 'mul'))
                                    s['dum1'] = jnp.where(s["_pc"] == 0, s['dum1'], old_1254)
                                    old_1255 = s['dum2']
                                    s['dum2'] = F(_arith(_div(s['nanew2'], F(2.0)), _arith(F(1.0), F(DERF1(s['uu2'])), 'sub'), 'mul'))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1255)
                                    old_1256 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum1'], s['dum2'], 'add'), s['rho'][s['k'] - 1]))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1256)
                                    old_1257 = s['dum2']
                                    s['dum2'] = F(jnp.minimum(_div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1257)
                                    old_1258 = s['dum2']
                                    s['dum2'] = F(_div(_arith(s['dum2'], s['nc3d'][s['k'] - 1], 'sub'), s['dt']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1258)
                                    old_1259 = s['dum2']
                                    s['dum2'] = F(jnp.maximum(F(0.0), s['dum2']))
                                    s['dum2'] = jnp.where(s["_pc"] == 0, s['dum2'], old_1259)
                                    old_1260 = s['nc3dten']
                                    s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1260)
                                    old_1261 = s['nact']
                                    s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['dum2'], 'add')))
                                    s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1261)
                                    return s
                                def no_1238(s):
                                    s = dict(s)
                                    return s
                                s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 2)), yes_1238, no_1238, s)
                                return s
                            s = lax.cond((s["_pc"] == 0) & ((s['iact'] == 1)), yes_1229, no_1229, s)
                            return s
                        def no_1228(s):
                            s = dict(s)
                            return s
                        s = lax.cond((s["_pc"] == 0) & ((s['ibase'] == 2)), yes_1228, no_1228, s)
                        return s
                    s = lax.cond((s["_pc"] == 0) & ((s['ibase'] == 1)), yes_1120, no_1120, s)
                    return s
                def no_1119(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['dum'] >= F(0.001))), yes_1119, no_1119, s)
                return s
            def no_1113(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((_arith(s['qc3d'][s['k'] - 1], _arith(s['qc3dten'][s['k'] - 1], s['dt'], 'mul'), 'add') >= s['qsmall'])) & ((s['inum'] == 0)))), yes_1113, no_1113, s)
            def yes_1262(s):
                s = dict(s)
                old_1263 = s['dum']
                s['dum'] = F(_div(_arith(s['eprd'][s['k'] - 1], s['dt'], 'mul'), s['qi3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1263)
                old_1264 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1264)
                old_1265 = s['nsubi']
                s['nsubi'] = s['nsubi'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['ni3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['nsubi'] = jnp.where(s["_pc"] == 0, s['nsubi'], old_1265)
                return s
            def no_1262(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['eprd'][s['k'] - 1] < F(0.0))), yes_1262, no_1262, s)
            def yes_1266(s):
                s = dict(s)
                old_1267 = s['dum']
                s['dum'] = F(_div(_arith(s['eprds'][s['k'] - 1], s['dt'], 'mul'), s['qni3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1267)
                old_1268 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1268)
                old_1269 = s['nsubs']
                s['nsubs'] = s['nsubs'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['ns3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['nsubs'] = jnp.where(s["_pc"] == 0, s['nsubs'], old_1269)
                return s
            def no_1266(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['eprds'][s['k'] - 1] < F(0.0))), yes_1266, no_1266, s)
            def yes_1270(s):
                s = dict(s)
                old_1271 = s['dum']
                s['dum'] = F(_div(_arith(s['pre'][s['k'] - 1], s['dt'], 'mul'), s['qr3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1271)
                old_1272 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1272)
                old_1273 = s['nsubr']
                s['nsubr'] = s['nsubr'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['nr3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['nsubr'] = jnp.where(s["_pc"] == 0, s['nsubr'], old_1273)
                return s
            def no_1270(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['pre'][s['k'] - 1] < F(0.0))), yes_1270, no_1270, s)
            def yes_1274(s):
                s = dict(s)
                old_1275 = s['dum']
                s['dum'] = F(_div(_arith(s['eprdg'][s['k'] - 1], s['dt'], 'mul'), s['qg3d'][s['k'] - 1]))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1275)
                old_1276 = s['dum']
                s['dum'] = F(jnp.maximum((-F(1.0)), s['dum']))
                s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1276)
                old_1277 = s['nsubg']
                s['nsubg'] = s['nsubg'].at[s['k'] - 1].set(F(_div(_arith(s['dum'], s['ng3d'][s['k'] - 1], 'mul'), s['dt'])))
                s['nsubg'] = jnp.where(s["_pc"] == 0, s['nsubg'], old_1277)
                return s
            def no_1274(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['eprdg'][s['k'] - 1] < F(0.0))), yes_1274, no_1274, s)
            old_1278 = s['ni3dten']
            s['ni3dten'] = s['ni3dten'].at[s['k'] - 1].set(F(_arith(s['ni3dten'][s['k'] - 1], s['nsubi'][s['k'] - 1], 'add')))
            s['ni3dten'] = jnp.where(s["_pc"] == 0, s['ni3dten'], old_1278)
            old_1279 = s['ns3dten']
            s['ns3dten'] = s['ns3dten'].at[s['k'] - 1].set(F(_arith(s['ns3dten'][s['k'] - 1], s['nsubs'][s['k'] - 1], 'add')))
            s['ns3dten'] = jnp.where(s["_pc"] == 0, s['ns3dten'], old_1279)
            old_1280 = s['ng3dten']
            s['ng3dten'] = s['ng3dten'].at[s['k'] - 1].set(F(_arith(s['ng3dten'][s['k'] - 1], s['nsubg'][s['k'] - 1], 'add')))
            s['ng3dten'] = jnp.where(s["_pc"] == 0, s['ng3dten'], old_1280)
            old_1281 = s['nr3dten']
            s['nr3dten'] = s['nr3dten'].at[s['k'] - 1].set(F(_arith(s['nr3dten'][s['k'] - 1], s['nsubr'][s['k'] - 1], 'add')))
            s['nr3dten'] = jnp.where(s["_pc"] == 0, s['nr3dten'], old_1281)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['t3d'][s['k'] - 1] >= s['tmelt'])), yes_203, no_203, s)
        old_1282 = s['ltrue']
        s['ltrue'] = I(1)
        s['ltrue'] = jnp.where(s["_pc"] == 0, s['ltrue'], old_1282)
        s["_pc"] = jnp.where(s["_pc"] == 200, I(0), s["_pc"])
        def yes_1283(s):
            s = dict(s)
            old_1284 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d_init'], _arith(_arith(s['t3d'][s['k'] - 1], s['t3d_init'], 'sub'), s['cf3d'][s['k'] - 1], 'mul'), 'add')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1284)
            old_1285 = s['t3dten']
            s['t3dten'] = s['t3dten'].at[s['k'] - 1].set(F(_arith(s['t3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['t3dten'] = jnp.where(s["_pc"] == 0, s['t3dten'], old_1285)
            old_1286 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv_init'], _arith(_arith(s['qv3d'][s['k'] - 1], s['qsat_init'], 'sub'), s['cf3d'][s['k'] - 1], 'mul'), 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1286)
            old_1287 = s['qv3dten']
            s['qv3dten'] = s['qv3dten'].at[s['k'] - 1].set(F(_arith(s['qv3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qv3dten'] = jnp.where(s["_pc"] == 0, s['qv3dten'], old_1287)
            old_1288 = s['qc3d']
            s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(_arith(s['qc3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_1288)
            old_1289 = s['qc3dten']
            s['qc3dten'] = s['qc3dten'].at[s['k'] - 1].set(F(_arith(s['qc3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qc3dten'] = jnp.where(s["_pc"] == 0, s['qc3dten'], old_1289)
            def yes_1290(s):
                s = dict(s)
                old_1291 = s['nc3d']
                s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_arith(s['nc3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1291)
                old_1292 = s['nc3dten']
                s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1292)
                return s
            def no_1290(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['inum'] == 0)), yes_1290, no_1290, s)
            old_1293 = s['qr3d']
            s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(_arith(s['qr3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_1293)
            old_1294 = s['qr3dten']
            s['qr3dten'] = s['qr3dten'].at[s['k'] - 1].set(F(_arith(s['qr3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qr3dten'] = jnp.where(s["_pc"] == 0, s['qr3dten'], old_1294)
            old_1295 = s['nr3d']
            s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1295)
            old_1296 = s['nr3dten']
            s['nr3dten'] = s['nr3dten'].at[s['k'] - 1].set(F(_arith(s['nr3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nr3dten'] = jnp.where(s["_pc"] == 0, s['nr3dten'], old_1296)
            def yes_1297(s):
                s = dict(s)
                old_1298 = s['qi3d']
                s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(_arith(s['qi3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_1298)
                old_1299 = s['qi3dten']
                s['qi3dten'] = s['qi3dten'].at[s['k'] - 1].set(F(_arith(s['qi3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['qi3dten'] = jnp.where(s["_pc"] == 0, s['qi3dten'], old_1299)
                old_1300 = s['ni3d']
                s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_arith(s['ni3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1300)
                old_1301 = s['ni3dten']
                s['ni3dten'] = s['ni3dten'].at[s['k'] - 1].set(F(_arith(s['ni3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['ni3dten'] = jnp.where(s["_pc"] == 0, s['ni3dten'], old_1301)
                old_1302 = s['qni3d']
                s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(_arith(s['qni3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_1302)
                old_1303 = s['qni3dten']
                s['qni3dten'] = s['qni3dten'].at[s['k'] - 1].set(F(_arith(s['qni3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['qni3dten'] = jnp.where(s["_pc"] == 0, s['qni3dten'], old_1303)
                old_1304 = s['ns3d']
                s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_arith(s['ns3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_1304)
                old_1305 = s['ns3dten']
                s['ns3dten'] = s['ns3dten'].at[s['k'] - 1].set(F(_arith(s['ns3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['ns3dten'] = jnp.where(s["_pc"] == 0, s['ns3dten'], old_1305)
                return s
            def no_1297(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['iliq'] == 0)), yes_1297, no_1297, s)
            def yes_1306(s):
                s = dict(s)
                old_1307 = s['qg3d']
                s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(_arith(s['qg3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_1307)
                old_1308 = s['qg3dten']
                s['qg3dten'] = s['qg3dten'].at[s['k'] - 1].set(F(_arith(s['qg3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['qg3dten'] = jnp.where(s["_pc"] == 0, s['qg3dten'], old_1308)
                old_1309 = s['ng3d']
                s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_arith(s['ng3d'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_1309)
                old_1310 = s['ng3dten']
                s['ng3dten'] = s['ng3dten'].at[s['k'] - 1].set(F(_arith(s['ng3dten'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
                s['ng3dten'] = jnp.where(s["_pc"] == 0, s['ng3dten'], old_1310)
                return s
            def no_1306(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 0)), yes_1306, no_1306, s)
            old_1311 = s['prc']
            s['prc'] = s['prc'].at[s['k'] - 1].set(F(_arith(s['prc'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['prc'] = jnp.where(s["_pc"] == 0, s['prc'], old_1311)
            old_1312 = s['pra']
            s['pra'] = s['pra'].at[s['k'] - 1].set(F(_arith(s['pra'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pra'] = jnp.where(s["_pc"] == 0, s['pra'], old_1312)
            old_1313 = s['psmlt']
            s['psmlt'] = s['psmlt'].at[s['k'] - 1].set(F(_arith(s['psmlt'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['psmlt'] = jnp.where(s["_pc"] == 0, s['psmlt'], old_1313)
            old_1314 = s['evpms']
            s['evpms'] = s['evpms'].at[s['k'] - 1].set(F(_arith(s['evpms'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['evpms'] = jnp.where(s["_pc"] == 0, s['evpms'], old_1314)
            old_1315 = s['pracs']
            s['pracs'] = s['pracs'].at[s['k'] - 1].set(F(_arith(s['pracs'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pracs'] = jnp.where(s["_pc"] == 0, s['pracs'], old_1315)
            old_1316 = s['evpmg']
            s['evpmg'] = s['evpmg'].at[s['k'] - 1].set(F(_arith(s['evpmg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['evpmg'] = jnp.where(s["_pc"] == 0, s['evpmg'], old_1316)
            old_1317 = s['pracg']
            s['pracg'] = s['pracg'].at[s['k'] - 1].set(F(_arith(s['pracg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pracg'] = jnp.where(s["_pc"] == 0, s['pracg'], old_1317)
            old_1318 = s['pre']
            s['pre'] = s['pre'].at[s['k'] - 1].set(F(_arith(s['pre'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pre'] = jnp.where(s["_pc"] == 0, s['pre'], old_1318)
            old_1319 = s['pgmlt']
            s['pgmlt'] = s['pgmlt'].at[s['k'] - 1].set(F(_arith(s['pgmlt'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pgmlt'] = jnp.where(s["_pc"] == 0, s['pgmlt'], old_1319)
            old_1320 = s['mnuccc']
            s['mnuccc'] = s['mnuccc'].at[s['k'] - 1].set(F(_arith(s['mnuccc'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['mnuccc'] = jnp.where(s["_pc"] == 0, s['mnuccc'], old_1320)
            old_1321 = s['psacws']
            s['psacws'] = s['psacws'].at[s['k'] - 1].set(F(_arith(s['psacws'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['psacws'] = jnp.where(s["_pc"] == 0, s['psacws'], old_1321)
            old_1322 = s['psacwi']
            s['psacwi'] = s['psacwi'].at[s['k'] - 1].set(F(_arith(s['psacwi'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['psacwi'] = jnp.where(s["_pc"] == 0, s['psacwi'], old_1322)
            old_1323 = s['qmults']
            s['qmults'] = s['qmults'].at[s['k'] - 1].set(F(_arith(s['qmults'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qmults'] = jnp.where(s["_pc"] == 0, s['qmults'], old_1323)
            old_1324 = s['qmultg']
            s['qmultg'] = s['qmultg'].at[s['k'] - 1].set(F(_arith(s['qmultg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qmultg'] = jnp.where(s["_pc"] == 0, s['qmultg'], old_1324)
            old_1325 = s['psacwg']
            s['psacwg'] = s['psacwg'].at[s['k'] - 1].set(F(_arith(s['psacwg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['psacwg'] = jnp.where(s["_pc"] == 0, s['psacwg'], old_1325)
            old_1326 = s['pgsacw']
            s['pgsacw'] = s['pgsacw'].at[s['k'] - 1].set(F(_arith(s['pgsacw'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pgsacw'] = jnp.where(s["_pc"] == 0, s['pgsacw'], old_1326)
            old_1327 = s['prd']
            s['prd'] = s['prd'].at[s['k'] - 1].set(F(_arith(s['prd'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['prd'] = jnp.where(s["_pc"] == 0, s['prd'], old_1327)
            old_1328 = s['prci']
            s['prci'] = s['prci'].at[s['k'] - 1].set(F(_arith(s['prci'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['prci'] = jnp.where(s["_pc"] == 0, s['prci'], old_1328)
            old_1329 = s['prai']
            s['prai'] = s['prai'].at[s['k'] - 1].set(F(_arith(s['prai'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['prai'] = jnp.where(s["_pc"] == 0, s['prai'], old_1329)
            old_1330 = s['qmultr']
            s['qmultr'] = s['qmultr'].at[s['k'] - 1].set(F(_arith(s['qmultr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qmultr'] = jnp.where(s["_pc"] == 0, s['qmultr'], old_1330)
            old_1331 = s['qmultrg']
            s['qmultrg'] = s['qmultrg'].at[s['k'] - 1].set(F(_arith(s['qmultrg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qmultrg'] = jnp.where(s["_pc"] == 0, s['qmultrg'], old_1331)
            old_1332 = s['mnuccd']
            s['mnuccd'] = s['mnuccd'].at[s['k'] - 1].set(F(_arith(s['mnuccd'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['mnuccd'] = jnp.where(s["_pc"] == 0, s['mnuccd'], old_1332)
            old_1333 = s['praci']
            s['praci'] = s['praci'].at[s['k'] - 1].set(F(_arith(s['praci'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['praci'] = jnp.where(s["_pc"] == 0, s['praci'], old_1333)
            old_1334 = s['pracis']
            s['pracis'] = s['pracis'].at[s['k'] - 1].set(F(_arith(s['pracis'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pracis'] = jnp.where(s["_pc"] == 0, s['pracis'], old_1334)
            old_1335 = s['eprd']
            s['eprd'] = s['eprd'].at[s['k'] - 1].set(F(_arith(s['eprd'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['eprd'] = jnp.where(s["_pc"] == 0, s['eprd'], old_1335)
            old_1336 = s['mnuccr']
            s['mnuccr'] = s['mnuccr'].at[s['k'] - 1].set(F(_arith(s['mnuccr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['mnuccr'] = jnp.where(s["_pc"] == 0, s['mnuccr'], old_1336)
            old_1337 = s['piacr']
            s['piacr'] = s['piacr'].at[s['k'] - 1].set(F(_arith(s['piacr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['piacr'] = jnp.where(s["_pc"] == 0, s['piacr'], old_1337)
            old_1338 = s['piacrs']
            s['piacrs'] = s['piacrs'].at[s['k'] - 1].set(F(_arith(s['piacrs'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['piacrs'] = jnp.where(s["_pc"] == 0, s['piacrs'], old_1338)
            old_1339 = s['pgracs']
            s['pgracs'] = s['pgracs'].at[s['k'] - 1].set(F(_arith(s['pgracs'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pgracs'] = jnp.where(s["_pc"] == 0, s['pgracs'], old_1339)
            old_1340 = s['prds']
            s['prds'] = s['prds'].at[s['k'] - 1].set(F(_arith(s['prds'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['prds'] = jnp.where(s["_pc"] == 0, s['prds'], old_1340)
            old_1341 = s['eprds']
            s['eprds'] = s['eprds'].at[s['k'] - 1].set(F(_arith(s['eprds'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['eprds'] = jnp.where(s["_pc"] == 0, s['eprds'], old_1341)
            old_1342 = s['psacr']
            s['psacr'] = s['psacr'].at[s['k'] - 1].set(F(_arith(s['psacr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['psacr'] = jnp.where(s["_pc"] == 0, s['psacr'], old_1342)
            old_1343 = s['prdg']
            s['prdg'] = s['prdg'].at[s['k'] - 1].set(F(_arith(s['prdg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['prdg'] = jnp.where(s["_pc"] == 0, s['prdg'], old_1343)
            old_1344 = s['eprdg']
            s['eprdg'] = s['eprdg'].at[s['k'] - 1].set(F(_arith(s['eprdg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['eprdg'] = jnp.where(s["_pc"] == 0, s['eprdg'], old_1344)
            old_1345 = s['nprc1']
            s['nprc1'] = s['nprc1'].at[s['k'] - 1].set(F(_arith(s['nprc1'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nprc1'] = jnp.where(s["_pc"] == 0, s['nprc1'], old_1345)
            old_1346 = s['nragg']
            s['nragg'] = s['nragg'].at[s['k'] - 1].set(F(_arith(s['nragg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nragg'] = jnp.where(s["_pc"] == 0, s['nragg'], old_1346)
            old_1347 = s['npracg']
            s['npracg'] = s['npracg'].at[s['k'] - 1].set(F(_arith(s['npracg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['npracg'] = jnp.where(s["_pc"] == 0, s['npracg'], old_1347)
            old_1348 = s['nsubr']
            s['nsubr'] = s['nsubr'].at[s['k'] - 1].set(F(_arith(s['nsubr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nsubr'] = jnp.where(s["_pc"] == 0, s['nsubr'], old_1348)
            old_1349 = s['nsmltr']
            s['nsmltr'] = s['nsmltr'].at[s['k'] - 1].set(F(_arith(s['nsmltr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nsmltr'] = jnp.where(s["_pc"] == 0, s['nsmltr'], old_1349)
            old_1350 = s['ngmltr']
            s['ngmltr'] = s['ngmltr'].at[s['k'] - 1].set(F(_arith(s['ngmltr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['ngmltr'] = jnp.where(s["_pc"] == 0, s['ngmltr'], old_1350)
            old_1351 = s['npracs']
            s['npracs'] = s['npracs'].at[s['k'] - 1].set(F(_arith(s['npracs'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['npracs'] = jnp.where(s["_pc"] == 0, s['npracs'], old_1351)
            old_1352 = s['nnuccr']
            s['nnuccr'] = s['nnuccr'].at[s['k'] - 1].set(F(_arith(s['nnuccr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nnuccr'] = jnp.where(s["_pc"] == 0, s['nnuccr'], old_1352)
            old_1353 = s['niacr']
            s['niacr'] = s['niacr'].at[s['k'] - 1].set(F(_arith(s['niacr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['niacr'] = jnp.where(s["_pc"] == 0, s['niacr'], old_1353)
            old_1354 = s['niacrs']
            s['niacrs'] = s['niacrs'].at[s['k'] - 1].set(F(_arith(s['niacrs'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['niacrs'] = jnp.where(s["_pc"] == 0, s['niacrs'], old_1354)
            old_1355 = s['ngracs']
            s['ngracs'] = s['ngracs'].at[s['k'] - 1].set(F(_arith(s['ngracs'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['ngracs'] = jnp.where(s["_pc"] == 0, s['ngracs'], old_1355)
            old_1356 = s['nsmlts']
            s['nsmlts'] = s['nsmlts'].at[s['k'] - 1].set(F(_arith(s['nsmlts'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nsmlts'] = jnp.where(s["_pc"] == 0, s['nsmlts'], old_1356)
            old_1357 = s['nsagg']
            s['nsagg'] = s['nsagg'].at[s['k'] - 1].set(F(_arith(s['nsagg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nsagg'] = jnp.where(s["_pc"] == 0, s['nsagg'], old_1357)
            old_1358 = s['nprci']
            s['nprci'] = s['nprci'].at[s['k'] - 1].set(F(_arith(s['nprci'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nprci'] = jnp.where(s["_pc"] == 0, s['nprci'], old_1358)
            old_1359 = s['nscng']
            s['nscng'] = s['nscng'].at[s['k'] - 1].set(F(_arith(s['nscng'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nscng'] = jnp.where(s["_pc"] == 0, s['nscng'], old_1359)
            old_1360 = s['nsubs']
            s['nsubs'] = s['nsubs'].at[s['k'] - 1].set(F(_arith(s['nsubs'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nsubs'] = jnp.where(s["_pc"] == 0, s['nsubs'], old_1360)
            old_1361 = s['pcc']
            s['pcc'] = s['pcc'].at[s['k'] - 1].set(F(_arith(s['pcc'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['pcc'] = jnp.where(s["_pc"] == 0, s['pcc'], old_1361)
            old_1362 = s['nnuccc']
            s['nnuccc'] = s['nnuccc'].at[s['k'] - 1].set(F(_arith(s['nnuccc'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nnuccc'] = jnp.where(s["_pc"] == 0, s['nnuccc'], old_1362)
            old_1363 = s['npsacws']
            s['npsacws'] = s['npsacws'].at[s['k'] - 1].set(F(_arith(s['npsacws'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['npsacws'] = jnp.where(s["_pc"] == 0, s['npsacws'], old_1363)
            old_1364 = s['npra']
            s['npra'] = s['npra'].at[s['k'] - 1].set(F(_arith(s['npra'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['npra'] = jnp.where(s["_pc"] == 0, s['npra'], old_1364)
            old_1365 = s['nprc']
            s['nprc'] = s['nprc'].at[s['k'] - 1].set(F(_arith(s['nprc'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nprc'] = jnp.where(s["_pc"] == 0, s['nprc'], old_1365)
            old_1366 = s['npsacwi']
            s['npsacwi'] = s['npsacwi'].at[s['k'] - 1].set(F(_arith(s['npsacwi'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['npsacwi'] = jnp.where(s["_pc"] == 0, s['npsacwi'], old_1366)
            old_1367 = s['npsacwg']
            s['npsacwg'] = s['npsacwg'].at[s['k'] - 1].set(F(_arith(s['npsacwg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['npsacwg'] = jnp.where(s["_pc"] == 0, s['npsacwg'], old_1367)
            old_1368 = s['nprai']
            s['nprai'] = s['nprai'].at[s['k'] - 1].set(F(_arith(s['nprai'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nprai'] = jnp.where(s["_pc"] == 0, s['nprai'], old_1368)
            old_1369 = s['nmults']
            s['nmults'] = s['nmults'].at[s['k'] - 1].set(F(_arith(s['nmults'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nmults'] = jnp.where(s["_pc"] == 0, s['nmults'], old_1369)
            old_1370 = s['nmultg']
            s['nmultg'] = s['nmultg'].at[s['k'] - 1].set(F(_arith(s['nmultg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nmultg'] = jnp.where(s["_pc"] == 0, s['nmultg'], old_1370)
            old_1371 = s['nmultr']
            s['nmultr'] = s['nmultr'].at[s['k'] - 1].set(F(_arith(s['nmultr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nmultr'] = jnp.where(s["_pc"] == 0, s['nmultr'], old_1371)
            old_1372 = s['nmultrg']
            s['nmultrg'] = s['nmultrg'].at[s['k'] - 1].set(F(_arith(s['nmultrg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nmultrg'] = jnp.where(s["_pc"] == 0, s['nmultrg'], old_1372)
            old_1373 = s['nnuccd']
            s['nnuccd'] = s['nnuccd'].at[s['k'] - 1].set(F(_arith(s['nnuccd'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nnuccd'] = jnp.where(s["_pc"] == 0, s['nnuccd'], old_1373)
            old_1374 = s['nsubi']
            s['nsubi'] = s['nsubi'].at[s['k'] - 1].set(F(_arith(s['nsubi'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nsubi'] = jnp.where(s["_pc"] == 0, s['nsubi'], old_1374)
            old_1375 = s['ngmltg']
            s['ngmltg'] = s['ngmltg'].at[s['k'] - 1].set(F(_arith(s['ngmltg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['ngmltg'] = jnp.where(s["_pc"] == 0, s['ngmltg'], old_1375)
            old_1376 = s['nsubg']
            s['nsubg'] = s['nsubg'].at[s['k'] - 1].set(F(_arith(s['nsubg'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nsubg'] = jnp.where(s["_pc"] == 0, s['nsubg'], old_1376)
            old_1377 = s['nact']
            s['nact'] = s['nact'].at[s['k'] - 1].set(F(_arith(s['nact'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1377)
            old_1378 = s['negfix_ni']
            s['negfix_ni'] = s['negfix_ni'].at[s['k'] - 1].set(F(_arith(s['negfix_ni'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['negfix_ni'] = jnp.where(s["_pc"] == 0, s['negfix_ni'], old_1378)
            old_1379 = s['negfix_ns']
            s['negfix_ns'] = s['negfix_ns'].at[s['k'] - 1].set(F(_arith(s['negfix_ns'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['negfix_ns'] = jnp.where(s["_pc"] == 0, s['negfix_ns'], old_1379)
            old_1380 = s['negfix_nc']
            s['negfix_nc'] = s['negfix_nc'].at[s['k'] - 1].set(F(_arith(s['negfix_nc'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['negfix_nc'] = jnp.where(s["_pc"] == 0, s['negfix_nc'], old_1380)
            old_1381 = s['negfix_nr']
            s['negfix_nr'] = s['negfix_nr'].at[s['k'] - 1].set(F(_arith(s['negfix_nr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['negfix_nr'] = jnp.where(s["_pc"] == 0, s['negfix_nr'], old_1381)
            old_1382 = s['negfix_ng']
            s['negfix_ng'] = s['negfix_ng'].at[s['k'] - 1].set(F(_arith(s['negfix_ng'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['negfix_ng'] = jnp.where(s["_pc"] == 0, s['negfix_ng'], old_1382)
            old_1383 = s['sizefix_nr']
            s['sizefix_nr'] = s['sizefix_nr'].at[s['k'] - 1].set(F(_arith(s['sizefix_nr'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_1383)
            old_1384 = s['sizefix_nc']
            s['sizefix_nc'] = s['sizefix_nc'].at[s['k'] - 1].set(F(_arith(s['sizefix_nc'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_1384)
            old_1385 = s['sizefix_ni']
            s['sizefix_ni'] = s['sizefix_ni'].at[s['k'] - 1].set(F(_arith(s['sizefix_ni'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['sizefix_ni'] = jnp.where(s["_pc"] == 0, s['sizefix_ni'], old_1385)
            old_1386 = s['sizefix_ns']
            s['sizefix_ns'] = s['sizefix_ns'].at[s['k'] - 1].set(F(_arith(s['sizefix_ns'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_1386)
            old_1387 = s['sizefix_ng']
            s['sizefix_ng'] = s['sizefix_ng'].at[s['k'] - 1].set(F(_arith(s['sizefix_ng'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_1387)
            old_1388 = s['qc_inst']
            s['qc_inst'] = s['qc_inst'].at[s['k'] - 1].set(F(_arith(s['qc_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qc_inst'] = jnp.where(s["_pc"] == 0, s['qc_inst'], old_1388)
            old_1389 = s['qr_inst']
            s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_1389)
            old_1390 = s['qi_inst']
            s['qi_inst'] = s['qi_inst'].at[s['k'] - 1].set(F(_arith(s['qi_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_1390)
            old_1391 = s['qs_inst']
            s['qs_inst'] = s['qs_inst'].at[s['k'] - 1].set(F(_arith(s['qs_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_1391)
            old_1392 = s['qg_inst']
            s['qg_inst'] = s['qg_inst'].at[s['k'] - 1].set(F(_arith(s['qg_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_1392)
            old_1393 = s['nc_inst']
            s['nc_inst'] = s['nc_inst'].at[s['k'] - 1].set(F(_arith(s['nc_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nc_inst'] = jnp.where(s["_pc"] == 0, s['nc_inst'], old_1393)
            old_1394 = s['nr_inst']
            s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_1394)
            old_1395 = s['ni_inst']
            s['ni_inst'] = s['ni_inst'].at[s['k'] - 1].set(F(_arith(s['ni_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['ni_inst'] = jnp.where(s["_pc"] == 0, s['ni_inst'], old_1395)
            old_1396 = s['ns_inst']
            s['ns_inst'] = s['ns_inst'].at[s['k'] - 1].set(F(_arith(s['ns_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['ns_inst'] = jnp.where(s["_pc"] == 0, s['ns_inst'], old_1396)
            old_1397 = s['ng_inst']
            s['ng_inst'] = s['ng_inst'].at[s['k'] - 1].set(F(_arith(s['ng_inst'][s['k'] - 1], s['cf3d'][s['k'] - 1], 'mul')))
            s['ng_inst'] = jnp.where(s["_pc"] == 0, s['ng_inst'], old_1397)
            return s
        def no_1283(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['cf3d'][s['k'] - 1] > s['cloud_frac_thresh'])), yes_1283, no_1283, s)
        old_1398 = s['nact']
        s['nact'] = s['nact'].at[s['k'] - 1].set(F(_div(s['nact'][s['k'] - 1], s['dt'])))
        s['nact'] = jnp.where(s["_pc"] == 0, s['nact'], old_1398)
        return s
    s = lax.fori_loop(0, I((s['kte'] - s['kts']) // (1) + 1), loop_76, s)
    old_1399 = s['precrt']
    s['precrt'] = F(F(0.0))
    s['precrt'] = jnp.where(s["_pc"] == 0, s['precrt'], old_1399)
    old_1400 = s['snowrt']
    s['snowrt'] = F(F(0.0))
    s['snowrt'] = jnp.where(s["_pc"] == 0, s['snowrt'], old_1400)
    def yes_1401(s):
        s = dict(s)
        s["_pc"] = jnp.where(s["_pc"] == 0, I(400), s["_pc"])
        return s
    def no_1401(s):
        s = dict(s)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['ltrue'] == 0)), yes_1401, no_1401, s)
    old_1402 = s['nstep']
    s['nstep'] = I(1)
    s['nstep'] = jnp.where(s["_pc"] == 0, s['nstep'], old_1402)
    def loop_1403(iteration, s):
        s = dict(s)
        s['k'] = I(s['kte'] + iteration * ((-1)))
        old_1404 = s['dumi']
        s['dumi'] = s['dumi'].at[s['k'] - 1].set(F(_arith(s['qi3d'][s['k'] - 1], _arith(s['qi3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumi'] = jnp.where(s["_pc"] == 0, s['dumi'], old_1404)
        old_1405 = s['dumqs']
        s['dumqs'] = s['dumqs'].at[s['k'] - 1].set(F(_arith(s['qni3d'][s['k'] - 1], _arith(s['qni3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumqs'] = jnp.where(s["_pc"] == 0, s['dumqs'], old_1405)
        old_1406 = s['dumr']
        s['dumr'] = s['dumr'].at[s['k'] - 1].set(F(_arith(s['qr3d'][s['k'] - 1], _arith(s['qr3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumr'] = jnp.where(s["_pc"] == 0, s['dumr'], old_1406)
        old_1407 = s['dumfni']
        s['dumfni'] = s['dumfni'].at[s['k'] - 1].set(F(_arith(s['ni3d'][s['k'] - 1], _arith(s['ni3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumfni'] = jnp.where(s["_pc"] == 0, s['dumfni'], old_1407)
        old_1408 = s['dumfns']
        s['dumfns'] = s['dumfns'].at[s['k'] - 1].set(F(_arith(s['ns3d'][s['k'] - 1], _arith(s['ns3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumfns'] = jnp.where(s["_pc"] == 0, s['dumfns'], old_1408)
        old_1409 = s['dumfnr']
        s['dumfnr'] = s['dumfnr'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], _arith(s['nr3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumfnr'] = jnp.where(s["_pc"] == 0, s['dumfnr'], old_1409)
        old_1410 = s['dumc']
        s['dumc'] = s['dumc'].at[s['k'] - 1].set(F(_arith(s['qc3d'][s['k'] - 1], _arith(s['qc3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumc'] = jnp.where(s["_pc"] == 0, s['dumc'], old_1410)
        old_1411 = s['dumfnc']
        s['dumfnc'] = s['dumfnc'].at[s['k'] - 1].set(F(_arith(s['nc3d'][s['k'] - 1], _arith(s['nc3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumfnc'] = jnp.where(s["_pc"] == 0, s['dumfnc'], old_1411)
        old_1412 = s['dumg']
        s['dumg'] = s['dumg'].at[s['k'] - 1].set(F(_arith(s['qg3d'][s['k'] - 1], _arith(s['qg3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumg'] = jnp.where(s["_pc"] == 0, s['dumg'], old_1412)
        old_1413 = s['dumfng']
        s['dumfng'] = s['dumfng'].at[s['k'] - 1].set(F(_arith(s['ng3d'][s['k'] - 1], _arith(s['ng3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['dumfng'] = jnp.where(s["_pc"] == 0, s['dumfng'], old_1413)
        def yes_1414(s):
            s = dict(s)
            old_1415 = s['dumfnc']
            s['dumfnc'] = s['dumfnc'].at[s['k'] - 1].set(F(s['nc3d'][s['k'] - 1]))
            s['dumfnc'] = jnp.where(s["_pc"] == 0, s['dumfnc'], old_1415)
            return s
        def no_1414(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['inum'] == 1)), yes_1414, no_1414, s)
        old_1416 = s['dumfni']
        s['dumfni'] = s['dumfni'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['dumfni'][s['k'] - 1])))
        s['dumfni'] = jnp.where(s["_pc"] == 0, s['dumfni'], old_1416)
        old_1417 = s['dumfns']
        s['dumfns'] = s['dumfns'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['dumfns'][s['k'] - 1])))
        s['dumfns'] = jnp.where(s["_pc"] == 0, s['dumfns'], old_1417)
        old_1418 = s['dumfnc']
        s['dumfnc'] = s['dumfnc'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['dumfnc'][s['k'] - 1])))
        s['dumfnc'] = jnp.where(s["_pc"] == 0, s['dumfnc'], old_1418)
        old_1419 = s['dumfnr']
        s['dumfnr'] = s['dumfnr'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['dumfnr'][s['k'] - 1])))
        s['dumfnr'] = jnp.where(s["_pc"] == 0, s['dumfnr'], old_1419)
        old_1420 = s['dumfng']
        s['dumfng'] = s['dumfng'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['dumfng'][s['k'] - 1])))
        s['dumfng'] = jnp.where(s["_pc"] == 0, s['dumfng'], old_1420)
        def yes_1421(s):
            s = dict(s)
            old_1422 = s['dlami']
            s['dlami'] = F(_arith(_div(_arith(s['cons12'], s['dumfni'][s['k'] - 1], 'mul'), s['dumi'][s['k'] - 1]), _div(F(1.0), s['di']), 'pow'))
            s['dlami'] = jnp.where(s["_pc"] == 0, s['dlami'], old_1422)
            old_1423 = s['dlami']
            s['dlami'] = F(jnp.maximum(s['dlami'], s['lammini']))
            s['dlami'] = jnp.where(s["_pc"] == 0, s['dlami'], old_1423)
            old_1424 = s['dlami']
            s['dlami'] = F(jnp.minimum(s['dlami'], s['lammaxi']))
            s['dlami'] = jnp.where(s["_pc"] == 0, s['dlami'], old_1424)
            return s
        def no_1421(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumi'][s['k'] - 1] >= s['qsmall'])), yes_1421, no_1421, s)
        def yes_1425(s):
            s = dict(s)
            old_1426 = s['dlamr']
            s['dlamr'] = F(_arith(_div(_arith(_arith(s['pi'], s['rhow'], 'mul'), s['dumfnr'][s['k'] - 1], 'mul'), s['dumr'][s['k'] - 1]), _div(F(1.0), F(3.0)), 'pow'))
            s['dlamr'] = jnp.where(s["_pc"] == 0, s['dlamr'], old_1426)
            old_1427 = s['dlamr']
            s['dlamr'] = F(jnp.maximum(s['dlamr'], s['lamminr']))
            s['dlamr'] = jnp.where(s["_pc"] == 0, s['dlamr'], old_1427)
            old_1428 = s['dlamr']
            s['dlamr'] = F(jnp.minimum(s['dlamr'], s['lammaxr']))
            s['dlamr'] = jnp.where(s["_pc"] == 0, s['dlamr'], old_1428)
            return s
        def no_1425(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumr'][s['k'] - 1] >= s['qsmall'])), yes_1425, no_1425, s)
        def yes_1429(s):
            s = dict(s)
            def yes_1430(s):
                s = dict(s)
                old_1431 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(s['pgam_fixed']))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1431)
                return s
            def no_1430(s):
                s = dict(s)
                old_1432 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_arith(F(0.0005714), _arith(_div(s['nc3d'][s['k'] - 1], F(1000000.0)), s['rho'][s['k'] - 1], 'mul'), 'mul'), F(0.2714), 'add')))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1432)
                old_1433 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_div(F(1.0), _arith(s['pgam'][s['k'] - 1], 2, 'pow')), F(1.0), 'sub')))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1433)
                old_1434 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.maximum(s['pgam'][s['k'] - 1], F(2.0))))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1434)
                old_1435 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.minimum(s['pgam'][s['k'] - 1], F(10.0))))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1435)
                return s
            s = lax.cond((s["_pc"] == 0) & (s['dofix_pgam']), yes_1430, no_1430, s)
            old_1436 = s['dlamc']
            s['dlamc'] = F(_arith(_div(_arith(_arith(s['cons26'], s['dumfnc'][s['k'] - 1], 'mul'), GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add')), 'mul'), _arith(s['dumc'][s['k'] - 1], GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add')), 'mul')), _div(F(1.0), F(3.0)), 'pow'))
            s['dlamc'] = jnp.where(s["_pc"] == 0, s['dlamc'], old_1436)
            old_1437 = s['lammin']
            s['lammin'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(6e-05)))
            s['lammin'] = jnp.where(s["_pc"] == 0, s['lammin'], old_1437)
            old_1438 = s['lammax']
            s['lammax'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(1e-06)))
            s['lammax'] = jnp.where(s["_pc"] == 0, s['lammax'], old_1438)
            old_1439 = s['dlamc']
            s['dlamc'] = F(jnp.maximum(s['dlamc'], s['lammin']))
            s['dlamc'] = jnp.where(s["_pc"] == 0, s['dlamc'], old_1439)
            old_1440 = s['dlamc']
            s['dlamc'] = F(jnp.minimum(s['dlamc'], s['lammax']))
            s['dlamc'] = jnp.where(s["_pc"] == 0, s['dlamc'], old_1440)
            return s
        def no_1429(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumc'][s['k'] - 1] >= s['qsmall'])), yes_1429, no_1429, s)
        def yes_1441(s):
            s = dict(s)
            old_1442 = s['dlams']
            s['dlams'] = F(_arith(_div(_arith(s['cons1'], s['dumfns'][s['k'] - 1], 'mul'), s['dumqs'][s['k'] - 1]), _div(F(1.0), s['ds']), 'pow'))
            s['dlams'] = jnp.where(s["_pc"] == 0, s['dlams'], old_1442)
            old_1443 = s['dlams']
            s['dlams'] = F(jnp.maximum(s['dlams'], s['lammins']))
            s['dlams'] = jnp.where(s["_pc"] == 0, s['dlams'], old_1443)
            old_1444 = s['dlams']
            s['dlams'] = F(jnp.minimum(s['dlams'], s['lammaxs']))
            s['dlams'] = jnp.where(s["_pc"] == 0, s['dlams'], old_1444)
            return s
        def no_1441(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumqs'][s['k'] - 1] >= s['qsmall'])), yes_1441, no_1441, s)
        def yes_1445(s):
            s = dict(s)
            old_1446 = s['dlamg']
            s['dlamg'] = F(_arith(_div(_arith(s['cons2'], s['dumfng'][s['k'] - 1], 'mul'), s['dumg'][s['k'] - 1]), _div(F(1.0), s['dg']), 'pow'))
            s['dlamg'] = jnp.where(s["_pc"] == 0, s['dlamg'], old_1446)
            old_1447 = s['dlamg']
            s['dlamg'] = F(jnp.maximum(s['dlamg'], s['lamming']))
            s['dlamg'] = jnp.where(s["_pc"] == 0, s['dlamg'], old_1447)
            old_1448 = s['dlamg']
            s['dlamg'] = F(jnp.minimum(s['dlamg'], s['lammaxg']))
            s['dlamg'] = jnp.where(s["_pc"] == 0, s['dlamg'], old_1448)
            return s
        def no_1445(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumg'][s['k'] - 1] >= s['qsmall'])), yes_1445, no_1445, s)
        def yes_1449(s):
            s = dict(s)
            old_1450 = s['unc']
            s['unc'] = F(_div(_arith(s['acn'][s['k'] - 1], GAMMA(_arith(_arith(F(1.0), s['bc'], 'add'), s['pgam'][s['k'] - 1], 'add')), 'mul'), _arith(_arith(s['dlamc'], s['bc'], 'pow'), GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add')), 'mul')))
            s['unc'] = jnp.where(s["_pc"] == 0, s['unc'], old_1450)
            old_1451 = s['umc']
            s['umc'] = F(_div(_arith(s['acn'][s['k'] - 1], GAMMA(_arith(_arith(F(4.0), s['bc'], 'add'), s['pgam'][s['k'] - 1], 'add')), 'mul'), _arith(_arith(s['dlamc'], s['bc'], 'pow'), GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add')), 'mul')))
            s['umc'] = jnp.where(s["_pc"] == 0, s['umc'], old_1451)
            return s
        def no_1449(s):
            s = dict(s)
            old_1452 = s['umc']
            s['umc'] = F(F(0.0))
            s['umc'] = jnp.where(s["_pc"] == 0, s['umc'], old_1452)
            old_1453 = s['unc']
            s['unc'] = F(F(0.0))
            s['unc'] = jnp.where(s["_pc"] == 0, s['unc'], old_1453)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumc'][s['k'] - 1] >= s['qsmall'])), yes_1449, no_1449, s)
        def yes_1454(s):
            s = dict(s)
            old_1455 = s['uni']
            s['uni'] = F(_div(_arith(s['ain'][s['k'] - 1], s['cons27'], 'mul'), _arith(s['dlami'], s['bi'], 'pow')))
            s['uni'] = jnp.where(s["_pc"] == 0, s['uni'], old_1455)
            old_1456 = s['umi']
            s['umi'] = F(_div(_arith(s['ain'][s['k'] - 1], s['cons28'], 'mul'), _arith(s['dlami'], s['bi'], 'pow')))
            s['umi'] = jnp.where(s["_pc"] == 0, s['umi'], old_1456)
            return s
        def no_1454(s):
            s = dict(s)
            old_1457 = s['umi']
            s['umi'] = F(F(0.0))
            s['umi'] = jnp.where(s["_pc"] == 0, s['umi'], old_1457)
            old_1458 = s['uni']
            s['uni'] = F(F(0.0))
            s['uni'] = jnp.where(s["_pc"] == 0, s['uni'], old_1458)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumi'][s['k'] - 1] >= s['qsmall'])), yes_1454, no_1454, s)
        def yes_1459(s):
            s = dict(s)
            old_1460 = s['unr']
            s['unr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons6'], 'mul'), _arith(s['dlamr'], s['br'], 'pow')))
            s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_1460)
            old_1461 = s['umr']
            s['umr'] = F(_div(_arith(s['arn'][s['k'] - 1], s['cons4'], 'mul'), _arith(s['dlamr'], s['br'], 'pow')))
            s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_1461)
            return s
        def no_1459(s):
            s = dict(s)
            old_1462 = s['umr']
            s['umr'] = F(F(0.0))
            s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_1462)
            old_1463 = s['unr']
            s['unr'] = F(F(0.0))
            s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_1463)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumr'][s['k'] - 1] >= s['qsmall'])), yes_1459, no_1459, s)
        def yes_1464(s):
            s = dict(s)
            old_1465 = s['ums']
            s['ums'] = F(_div(_arith(s['asn'][s['k'] - 1], s['cons3'], 'mul'), _arith(s['dlams'], s['bs'], 'pow')))
            s['ums'] = jnp.where(s["_pc"] == 0, s['ums'], old_1465)
            old_1466 = s['uns']
            s['uns'] = F(_div(_arith(s['asn'][s['k'] - 1], s['cons5'], 'mul'), _arith(s['dlams'], s['bs'], 'pow')))
            s['uns'] = jnp.where(s["_pc"] == 0, s['uns'], old_1466)
            return s
        def no_1464(s):
            s = dict(s)
            old_1467 = s['ums']
            s['ums'] = F(F(0.0))
            s['ums'] = jnp.where(s["_pc"] == 0, s['ums'], old_1467)
            old_1468 = s['uns']
            s['uns'] = F(F(0.0))
            s['uns'] = jnp.where(s["_pc"] == 0, s['uns'], old_1468)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumqs'][s['k'] - 1] >= s['qsmall'])), yes_1464, no_1464, s)
        def yes_1469(s):
            s = dict(s)
            old_1470 = s['umg']
            s['umg'] = F(_div(_arith(s['agn'][s['k'] - 1], s['cons7'], 'mul'), _arith(s['dlamg'], s['bg'], 'pow')))
            s['umg'] = jnp.where(s["_pc"] == 0, s['umg'], old_1470)
            old_1471 = s['ung']
            s['ung'] = F(_div(_arith(s['agn'][s['k'] - 1], s['cons8'], 'mul'), _arith(s['dlamg'], s['bg'], 'pow')))
            s['ung'] = jnp.where(s["_pc"] == 0, s['ung'], old_1471)
            return s
        def no_1469(s):
            s = dict(s)
            old_1472 = s['umg']
            s['umg'] = F(F(0.0))
            s['umg'] = jnp.where(s["_pc"] == 0, s['umg'], old_1472)
            old_1473 = s['ung']
            s['ung'] = F(F(0.0))
            s['ung'] = jnp.where(s["_pc"] == 0, s['ung'], old_1473)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['dumg'][s['k'] - 1] >= s['qsmall'])), yes_1469, no_1469, s)
        old_1474 = s['dum']
        s['dum'] = F(_arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.54), 'pow'))
        s['dum'] = jnp.where(s["_pc"] == 0, s['dum'], old_1474)
        old_1475 = s['ums']
        s['ums'] = F(jnp.minimum(s['ums'], _arith(F(1.2), s['dum'], 'mul')))
        s['ums'] = jnp.where(s["_pc"] == 0, s['ums'], old_1475)
        old_1476 = s['uns']
        s['uns'] = F(jnp.minimum(s['uns'], _arith(F(1.2), s['dum'], 'mul')))
        s['uns'] = jnp.where(s["_pc"] == 0, s['uns'], old_1476)
        old_1477 = s['umi']
        s['umi'] = F(jnp.minimum(s['umi'], _arith(F(1.2), _arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.35), 'pow'), 'mul')))
        s['umi'] = jnp.where(s["_pc"] == 0, s['umi'], old_1477)
        old_1478 = s['uni']
        s['uni'] = F(jnp.minimum(s['uni'], _arith(F(1.2), _arith(_div(s['rhosu'], s['rho'][s['k'] - 1]), F(0.35), 'pow'), 'mul')))
        s['uni'] = jnp.where(s["_pc"] == 0, s['uni'], old_1478)
        old_1479 = s['umr']
        s['umr'] = F(jnp.minimum(s['umr'], _arith(F(9.1), s['dum'], 'mul')))
        s['umr'] = jnp.where(s["_pc"] == 0, s['umr'], old_1479)
        old_1480 = s['unr']
        s['unr'] = F(jnp.minimum(s['unr'], _arith(F(9.1), s['dum'], 'mul')))
        s['unr'] = jnp.where(s["_pc"] == 0, s['unr'], old_1480)
        old_1481 = s['umg']
        s['umg'] = F(jnp.minimum(s['umg'], _arith(F(20.0), s['dum'], 'mul')))
        s['umg'] = jnp.where(s["_pc"] == 0, s['umg'], old_1481)
        old_1482 = s['ung']
        s['ung'] = F(jnp.minimum(s['ung'], _arith(F(20.0), s['dum'], 'mul')))
        s['ung'] = jnp.where(s["_pc"] == 0, s['ung'], old_1482)
        old_1483 = s['fr']
        s['fr'] = s['fr'].at[s['k'] - 1].set(F(s['umr']))
        s['fr'] = jnp.where(s["_pc"] == 0, s['fr'], old_1483)
        old_1484 = s['fi']
        s['fi'] = s['fi'].at[s['k'] - 1].set(F(s['umi']))
        s['fi'] = jnp.where(s["_pc"] == 0, s['fi'], old_1484)
        old_1485 = s['fni']
        s['fni'] = s['fni'].at[s['k'] - 1].set(F(s['uni']))
        s['fni'] = jnp.where(s["_pc"] == 0, s['fni'], old_1485)
        old_1486 = s['fs']
        s['fs'] = s['fs'].at[s['k'] - 1].set(F(s['ums']))
        s['fs'] = jnp.where(s["_pc"] == 0, s['fs'], old_1486)
        old_1487 = s['fns']
        s['fns'] = s['fns'].at[s['k'] - 1].set(F(s['uns']))
        s['fns'] = jnp.where(s["_pc"] == 0, s['fns'], old_1487)
        old_1488 = s['fnr']
        s['fnr'] = s['fnr'].at[s['k'] - 1].set(F(s['unr']))
        s['fnr'] = jnp.where(s["_pc"] == 0, s['fnr'], old_1488)
        old_1489 = s['fc']
        s['fc'] = s['fc'].at[s['k'] - 1].set(F(s['umc']))
        s['fc'] = jnp.where(s["_pc"] == 0, s['fc'], old_1489)
        old_1490 = s['fnc']
        s['fnc'] = s['fnc'].at[s['k'] - 1].set(F(s['unc']))
        s['fnc'] = jnp.where(s["_pc"] == 0, s['fnc'], old_1490)
        old_1491 = s['fg']
        s['fg'] = s['fg'].at[s['k'] - 1].set(F(s['umg']))
        s['fg'] = jnp.where(s["_pc"] == 0, s['fg'], old_1491)
        old_1492 = s['fng']
        s['fng'] = s['fng'].at[s['k'] - 1].set(F(s['ung']))
        s['fng'] = jnp.where(s["_pc"] == 0, s['fng'], old_1492)
        def yes_1493(s):
            s = dict(s)
            def yes_1494(s):
                s = dict(s)
                old_1495 = s['fr']
                s['fr'] = s['fr'].at[s['k'] - 1].set(F(s['fr'][_arith(s['k'], 1, 'add') - 1]))
                s['fr'] = jnp.where(s["_pc"] == 0, s['fr'], old_1495)
                return s
            def no_1494(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fr'][s['k'] - 1] < F(1e-10))), yes_1494, no_1494, s)
            def yes_1496(s):
                s = dict(s)
                old_1497 = s['fi']
                s['fi'] = s['fi'].at[s['k'] - 1].set(F(s['fi'][_arith(s['k'], 1, 'add') - 1]))
                s['fi'] = jnp.where(s["_pc"] == 0, s['fi'], old_1497)
                return s
            def no_1496(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fi'][s['k'] - 1] < F(1e-10))), yes_1496, no_1496, s)
            def yes_1498(s):
                s = dict(s)
                old_1499 = s['fni']
                s['fni'] = s['fni'].at[s['k'] - 1].set(F(s['fni'][_arith(s['k'], 1, 'add') - 1]))
                s['fni'] = jnp.where(s["_pc"] == 0, s['fni'], old_1499)
                return s
            def no_1498(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fni'][s['k'] - 1] < F(1e-10))), yes_1498, no_1498, s)
            def yes_1500(s):
                s = dict(s)
                old_1501 = s['fs']
                s['fs'] = s['fs'].at[s['k'] - 1].set(F(s['fs'][_arith(s['k'], 1, 'add') - 1]))
                s['fs'] = jnp.where(s["_pc"] == 0, s['fs'], old_1501)
                return s
            def no_1500(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fs'][s['k'] - 1] < F(1e-10))), yes_1500, no_1500, s)
            def yes_1502(s):
                s = dict(s)
                old_1503 = s['fns']
                s['fns'] = s['fns'].at[s['k'] - 1].set(F(s['fns'][_arith(s['k'], 1, 'add') - 1]))
                s['fns'] = jnp.where(s["_pc"] == 0, s['fns'], old_1503)
                return s
            def no_1502(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fns'][s['k'] - 1] < F(1e-10))), yes_1502, no_1502, s)
            def yes_1504(s):
                s = dict(s)
                old_1505 = s['fnr']
                s['fnr'] = s['fnr'].at[s['k'] - 1].set(F(s['fnr'][_arith(s['k'], 1, 'add') - 1]))
                s['fnr'] = jnp.where(s["_pc"] == 0, s['fnr'], old_1505)
                return s
            def no_1504(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fnr'][s['k'] - 1] < F(1e-10))), yes_1504, no_1504, s)
            def yes_1506(s):
                s = dict(s)
                old_1507 = s['fc']
                s['fc'] = s['fc'].at[s['k'] - 1].set(F(s['fc'][_arith(s['k'], 1, 'add') - 1]))
                s['fc'] = jnp.where(s["_pc"] == 0, s['fc'], old_1507)
                return s
            def no_1506(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fc'][s['k'] - 1] < F(1e-10))), yes_1506, no_1506, s)
            def yes_1508(s):
                s = dict(s)
                old_1509 = s['fnc']
                s['fnc'] = s['fnc'].at[s['k'] - 1].set(F(s['fnc'][_arith(s['k'], 1, 'add') - 1]))
                s['fnc'] = jnp.where(s["_pc"] == 0, s['fnc'], old_1509)
                return s
            def no_1508(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fnc'][s['k'] - 1] < F(1e-10))), yes_1508, no_1508, s)
            def yes_1510(s):
                s = dict(s)
                old_1511 = s['fg']
                s['fg'] = s['fg'].at[s['k'] - 1].set(F(s['fg'][_arith(s['k'], 1, 'add') - 1]))
                s['fg'] = jnp.where(s["_pc"] == 0, s['fg'], old_1511)
                return s
            def no_1510(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fg'][s['k'] - 1] < F(1e-10))), yes_1510, no_1510, s)
            def yes_1512(s):
                s = dict(s)
                old_1513 = s['fng']
                s['fng'] = s['fng'].at[s['k'] - 1].set(F(s['fng'][_arith(s['k'], 1, 'add') - 1]))
                s['fng'] = jnp.where(s["_pc"] == 0, s['fng'], old_1513)
                return s
            def no_1512(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['fng'][s['k'] - 1] < F(1e-10))), yes_1512, no_1512, s)
            return s
        def no_1493(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['k'] <= _arith(s['kte'], 1, 'sub'))), yes_1493, no_1493, s)
        old_1514 = s['rgvm']
        s['rgvm'] = F(jnp.maximum(jnp.maximum(jnp.maximum(jnp.maximum(jnp.maximum(jnp.maximum(jnp.maximum(jnp.maximum(jnp.maximum(s['fr'][s['k'] - 1], s['fi'][s['k'] - 1]), s['fs'][s['k'] - 1]), s['fc'][s['k'] - 1]), s['fni'][s['k'] - 1]), s['fnr'][s['k'] - 1]), s['fns'][s['k'] - 1]), s['fnc'][s['k'] - 1]), s['fg'][s['k'] - 1]), s['fng'][s['k'] - 1]))
        s['rgvm'] = jnp.where(s["_pc"] == 0, s['rgvm'], old_1514)
        old_1515 = s['nstep']
        s['nstep'] = I(jnp.maximum(I(_arith(_div(_arith(s['rgvm'], s['dt'], 'mul'), s['dzq'][s['k'] - 1]), F(1.0), 'add')), s['nstep']))
        s['nstep'] = jnp.where(s["_pc"] == 0, s['nstep'], old_1515)
        old_1516 = s['dumr']
        s['dumr'] = s['dumr'].at[s['k'] - 1].set(F(_arith(s['dumr'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumr'] = jnp.where(s["_pc"] == 0, s['dumr'], old_1516)
        old_1517 = s['dumi']
        s['dumi'] = s['dumi'].at[s['k'] - 1].set(F(_arith(s['dumi'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumi'] = jnp.where(s["_pc"] == 0, s['dumi'], old_1517)
        old_1518 = s['dumfni']
        s['dumfni'] = s['dumfni'].at[s['k'] - 1].set(F(_arith(s['dumfni'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumfni'] = jnp.where(s["_pc"] == 0, s['dumfni'], old_1518)
        old_1519 = s['dumqs']
        s['dumqs'] = s['dumqs'].at[s['k'] - 1].set(F(_arith(s['dumqs'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumqs'] = jnp.where(s["_pc"] == 0, s['dumqs'], old_1519)
        old_1520 = s['dumfns']
        s['dumfns'] = s['dumfns'].at[s['k'] - 1].set(F(_arith(s['dumfns'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumfns'] = jnp.where(s["_pc"] == 0, s['dumfns'], old_1520)
        old_1521 = s['dumfnr']
        s['dumfnr'] = s['dumfnr'].at[s['k'] - 1].set(F(_arith(s['dumfnr'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumfnr'] = jnp.where(s["_pc"] == 0, s['dumfnr'], old_1521)
        old_1522 = s['dumc']
        s['dumc'] = s['dumc'].at[s['k'] - 1].set(F(_arith(s['dumc'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumc'] = jnp.where(s["_pc"] == 0, s['dumc'], old_1522)
        old_1523 = s['dumfnc']
        s['dumfnc'] = s['dumfnc'].at[s['k'] - 1].set(F(_arith(s['dumfnc'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumfnc'] = jnp.where(s["_pc"] == 0, s['dumfnc'], old_1523)
        old_1524 = s['dumg']
        s['dumg'] = s['dumg'].at[s['k'] - 1].set(F(_arith(s['dumg'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumg'] = jnp.where(s["_pc"] == 0, s['dumg'], old_1524)
        old_1525 = s['dumfng']
        s['dumfng'] = s['dumfng'].at[s['k'] - 1].set(F(_arith(s['dumfng'][s['k'] - 1], s['rho'][s['k'] - 1], 'mul')))
        s['dumfng'] = jnp.where(s["_pc"] == 0, s['dumfng'], old_1525)
        return s
    s = lax.fori_loop(0, I((s['kts'] - s['kte']) // ((-1)) + 1), loop_1403, s)
    def loop_1526(iteration, s):
        s = dict(s)
        s['n'] = I(1 + iteration * (1))
        def loop_1527(iteration, s):
            s = dict(s)
            s['k'] = I(s['kts'] + iteration * (1))
            old_1528 = s['faloutr']
            s['faloutr'] = s['faloutr'].at[s['k'] - 1].set(F(_arith(s['fr'][s['k'] - 1], s['dumr'][s['k'] - 1], 'mul')))
            s['faloutr'] = jnp.where(s["_pc"] == 0, s['faloutr'], old_1528)
            old_1529 = s['falouti']
            s['falouti'] = s['falouti'].at[s['k'] - 1].set(F(_arith(s['fi'][s['k'] - 1], s['dumi'][s['k'] - 1], 'mul')))
            s['falouti'] = jnp.where(s["_pc"] == 0, s['falouti'], old_1529)
            old_1530 = s['faloutni']
            s['faloutni'] = s['faloutni'].at[s['k'] - 1].set(F(_arith(s['fni'][s['k'] - 1], s['dumfni'][s['k'] - 1], 'mul')))
            s['faloutni'] = jnp.where(s["_pc"] == 0, s['faloutni'], old_1530)
            old_1531 = s['falouts']
            s['falouts'] = s['falouts'].at[s['k'] - 1].set(F(_arith(s['fs'][s['k'] - 1], s['dumqs'][s['k'] - 1], 'mul')))
            s['falouts'] = jnp.where(s["_pc"] == 0, s['falouts'], old_1531)
            old_1532 = s['faloutns']
            s['faloutns'] = s['faloutns'].at[s['k'] - 1].set(F(_arith(s['fns'][s['k'] - 1], s['dumfns'][s['k'] - 1], 'mul')))
            s['faloutns'] = jnp.where(s["_pc"] == 0, s['faloutns'], old_1532)
            old_1533 = s['faloutnr']
            s['faloutnr'] = s['faloutnr'].at[s['k'] - 1].set(F(_arith(s['fnr'][s['k'] - 1], s['dumfnr'][s['k'] - 1], 'mul')))
            s['faloutnr'] = jnp.where(s["_pc"] == 0, s['faloutnr'], old_1533)
            old_1534 = s['faloutc']
            s['faloutc'] = s['faloutc'].at[s['k'] - 1].set(F(_arith(s['fc'][s['k'] - 1], s['dumc'][s['k'] - 1], 'mul')))
            s['faloutc'] = jnp.where(s["_pc"] == 0, s['faloutc'], old_1534)
            old_1535 = s['faloutnc']
            s['faloutnc'] = s['faloutnc'].at[s['k'] - 1].set(F(_arith(s['fnc'][s['k'] - 1], s['dumfnc'][s['k'] - 1], 'mul')))
            s['faloutnc'] = jnp.where(s["_pc"] == 0, s['faloutnc'], old_1535)
            old_1536 = s['faloutg']
            s['faloutg'] = s['faloutg'].at[s['k'] - 1].set(F(_arith(s['fg'][s['k'] - 1], s['dumg'][s['k'] - 1], 'mul')))
            s['faloutg'] = jnp.where(s["_pc"] == 0, s['faloutg'], old_1536)
            old_1537 = s['faloutng']
            s['faloutng'] = s['faloutng'].at[s['k'] - 1].set(F(_arith(s['fng'][s['k'] - 1], s['dumfng'][s['k'] - 1], 'mul')))
            s['faloutng'] = jnp.where(s["_pc"] == 0, s['faloutng'], old_1537)
            return s
        s = lax.fori_loop(0, I((s['kte'] - s['kts']) // (1) + 1), loop_1527, s)
        old_1538 = s['k']
        s['k'] = I(s['kte'])
        s['k'] = jnp.where(s["_pc"] == 0, s['k'], old_1538)
        old_1539 = s['faltndr']
        s['faltndr'] = F(_div(s['faloutr'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndr'] = jnp.where(s["_pc"] == 0, s['faltndr'], old_1539)
        old_1540 = s['faltndi']
        s['faltndi'] = F(_div(s['falouti'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndi'] = jnp.where(s["_pc"] == 0, s['faltndi'], old_1540)
        old_1541 = s['faltndni']
        s['faltndni'] = F(_div(s['faloutni'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndni'] = jnp.where(s["_pc"] == 0, s['faltndni'], old_1541)
        old_1542 = s['faltnds']
        s['faltnds'] = F(_div(s['falouts'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltnds'] = jnp.where(s["_pc"] == 0, s['faltnds'], old_1542)
        old_1543 = s['faltndns']
        s['faltndns'] = F(_div(s['faloutns'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndns'] = jnp.where(s["_pc"] == 0, s['faltndns'], old_1543)
        old_1544 = s['faltndnr']
        s['faltndnr'] = F(_div(s['faloutnr'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndnr'] = jnp.where(s["_pc"] == 0, s['faltndnr'], old_1544)
        old_1545 = s['faltndc']
        s['faltndc'] = F(_div(s['faloutc'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndc'] = jnp.where(s["_pc"] == 0, s['faltndc'], old_1545)
        old_1546 = s['faltndnc']
        s['faltndnc'] = F(_div(s['faloutnc'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndnc'] = jnp.where(s["_pc"] == 0, s['faltndnc'], old_1546)
        old_1547 = s['faltndg']
        s['faltndg'] = F(_div(s['faloutg'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndg'] = jnp.where(s["_pc"] == 0, s['faltndg'], old_1547)
        old_1548 = s['faltndng']
        s['faltndng'] = F(_div(s['faloutng'][s['k'] - 1], s['dzq'][s['k'] - 1]))
        s['faltndng'] = jnp.where(s["_pc"] == 0, s['faltndng'], old_1548)
        old_1549 = s['qrsten']
        s['qrsten'] = s['qrsten'].at[s['k'] - 1].set(F(_arith(s['qrsten'][s['k'] - 1], _div(_div(s['faltndr'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['qrsten'] = jnp.where(s["_pc"] == 0, s['qrsten'], old_1549)
        old_1550 = s['qisten']
        s['qisten'] = s['qisten'].at[s['k'] - 1].set(F(_arith(s['qisten'][s['k'] - 1], _div(_div(s['faltndi'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['qisten'] = jnp.where(s["_pc"] == 0, s['qisten'], old_1550)
        old_1551 = s['ni3dten']
        s['ni3dten'] = s['ni3dten'].at[s['k'] - 1].set(F(_arith(s['ni3dten'][s['k'] - 1], _div(_div(s['faltndni'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['ni3dten'] = jnp.where(s["_pc"] == 0, s['ni3dten'], old_1551)
        old_1552 = s['qnisten']
        s['qnisten'] = s['qnisten'].at[s['k'] - 1].set(F(_arith(s['qnisten'][s['k'] - 1], _div(_div(s['faltnds'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['qnisten'] = jnp.where(s["_pc"] == 0, s['qnisten'], old_1552)
        old_1553 = s['ns3dten']
        s['ns3dten'] = s['ns3dten'].at[s['k'] - 1].set(F(_arith(s['ns3dten'][s['k'] - 1], _div(_div(s['faltndns'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['ns3dten'] = jnp.where(s["_pc"] == 0, s['ns3dten'], old_1553)
        old_1554 = s['nr3dten']
        s['nr3dten'] = s['nr3dten'].at[s['k'] - 1].set(F(_arith(s['nr3dten'][s['k'] - 1], _div(_div(s['faltndnr'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['nr3dten'] = jnp.where(s["_pc"] == 0, s['nr3dten'], old_1554)
        old_1555 = s['qcsten']
        s['qcsten'] = s['qcsten'].at[s['k'] - 1].set(F(_arith(s['qcsten'][s['k'] - 1], _div(_div(s['faltndc'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['qcsten'] = jnp.where(s["_pc"] == 0, s['qcsten'], old_1555)
        old_1556 = s['nc3dten']
        s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], _div(_div(s['faltndnc'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1556)
        old_1557 = s['qgsten']
        s['qgsten'] = s['qgsten'].at[s['k'] - 1].set(F(_arith(s['qgsten'][s['k'] - 1], _div(_div(s['faltndg'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['qgsten'] = jnp.where(s["_pc"] == 0, s['qgsten'], old_1557)
        old_1558 = s['ng3dten']
        s['ng3dten'] = s['ng3dten'].at[s['k'] - 1].set(F(_arith(s['ng3dten'][s['k'] - 1], _div(_div(s['faltndng'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['ng3dten'] = jnp.where(s["_pc"] == 0, s['ng3dten'], old_1558)
        old_1559 = s['nisten']
        s['nisten'] = s['nisten'].at[s['k'] - 1].set(F(_arith(s['nisten'][s['k'] - 1], _div(_div(s['faltndni'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['nisten'] = jnp.where(s["_pc"] == 0, s['nisten'], old_1559)
        old_1560 = s['nssten']
        s['nssten'] = s['nssten'].at[s['k'] - 1].set(F(_arith(s['nssten'][s['k'] - 1], _div(_div(s['faltndns'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['nssten'] = jnp.where(s["_pc"] == 0, s['nssten'], old_1560)
        old_1561 = s['nrsten']
        s['nrsten'] = s['nrsten'].at[s['k'] - 1].set(F(_arith(s['nrsten'][s['k'] - 1], _div(_div(s['faltndnr'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['nrsten'] = jnp.where(s["_pc"] == 0, s['nrsten'], old_1561)
        old_1562 = s['ncsten']
        s['ncsten'] = s['ncsten'].at[s['k'] - 1].set(F(_arith(s['ncsten'][s['k'] - 1], _div(_div(s['faltndnc'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['ncsten'] = jnp.where(s["_pc"] == 0, s['ncsten'], old_1562)
        old_1563 = s['ngsten']
        s['ngsten'] = s['ngsten'].at[s['k'] - 1].set(F(_arith(s['ngsten'][s['k'] - 1], _div(_div(s['faltndng'], s['nstep']), s['rho'][s['k'] - 1]), 'sub')))
        s['ngsten'] = jnp.where(s["_pc"] == 0, s['ngsten'], old_1563)
        old_1564 = s['dumr']
        s['dumr'] = s['dumr'].at[s['k'] - 1].set(F(_arith(s['dumr'][s['k'] - 1], _div(_arith(s['faltndr'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumr'] = jnp.where(s["_pc"] == 0, s['dumr'], old_1564)
        old_1565 = s['dumi']
        s['dumi'] = s['dumi'].at[s['k'] - 1].set(F(_arith(s['dumi'][s['k'] - 1], _div(_arith(s['faltndi'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumi'] = jnp.where(s["_pc"] == 0, s['dumi'], old_1565)
        old_1566 = s['dumfni']
        s['dumfni'] = s['dumfni'].at[s['k'] - 1].set(F(_arith(s['dumfni'][s['k'] - 1], _div(_arith(s['faltndni'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumfni'] = jnp.where(s["_pc"] == 0, s['dumfni'], old_1566)
        old_1567 = s['dumqs']
        s['dumqs'] = s['dumqs'].at[s['k'] - 1].set(F(_arith(s['dumqs'][s['k'] - 1], _div(_arith(s['faltnds'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumqs'] = jnp.where(s["_pc"] == 0, s['dumqs'], old_1567)
        old_1568 = s['dumfns']
        s['dumfns'] = s['dumfns'].at[s['k'] - 1].set(F(_arith(s['dumfns'][s['k'] - 1], _div(_arith(s['faltndns'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumfns'] = jnp.where(s["_pc"] == 0, s['dumfns'], old_1568)
        old_1569 = s['dumfnr']
        s['dumfnr'] = s['dumfnr'].at[s['k'] - 1].set(F(_arith(s['dumfnr'][s['k'] - 1], _div(_arith(s['faltndnr'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumfnr'] = jnp.where(s["_pc"] == 0, s['dumfnr'], old_1569)
        old_1570 = s['dumc']
        s['dumc'] = s['dumc'].at[s['k'] - 1].set(F(_arith(s['dumc'][s['k'] - 1], _div(_arith(s['faltndc'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumc'] = jnp.where(s["_pc"] == 0, s['dumc'], old_1570)
        old_1571 = s['dumfnc']
        s['dumfnc'] = s['dumfnc'].at[s['k'] - 1].set(F(_arith(s['dumfnc'][s['k'] - 1], _div(_arith(s['faltndnc'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumfnc'] = jnp.where(s["_pc"] == 0, s['dumfnc'], old_1571)
        old_1572 = s['dumg']
        s['dumg'] = s['dumg'].at[s['k'] - 1].set(F(_arith(s['dumg'][s['k'] - 1], _div(_arith(s['faltndg'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumg'] = jnp.where(s["_pc"] == 0, s['dumg'], old_1572)
        old_1573 = s['dumfng']
        s['dumfng'] = s['dumfng'].at[s['k'] - 1].set(F(_arith(s['dumfng'][s['k'] - 1], _div(_arith(s['faltndng'], s['dt'], 'mul'), s['nstep']), 'sub')))
        s['dumfng'] = jnp.where(s["_pc"] == 0, s['dumfng'], old_1573)
        def loop_1574(iteration, s):
            s = dict(s)
            s['k'] = I(_arith(s['kte'], 1, 'sub') + iteration * ((-1)))
            old_1575 = s['faltndr']
            s['faltndr'] = F(_div(_arith(s['faloutr'][_arith(s['k'], 1, 'add') - 1], s['faloutr'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndr'] = jnp.where(s["_pc"] == 0, s['faltndr'], old_1575)
            old_1576 = s['faltndi']
            s['faltndi'] = F(_div(_arith(s['falouti'][_arith(s['k'], 1, 'add') - 1], s['falouti'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndi'] = jnp.where(s["_pc"] == 0, s['faltndi'], old_1576)
            old_1577 = s['faltndni']
            s['faltndni'] = F(_div(_arith(s['faloutni'][_arith(s['k'], 1, 'add') - 1], s['faloutni'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndni'] = jnp.where(s["_pc"] == 0, s['faltndni'], old_1577)
            old_1578 = s['faltnds']
            s['faltnds'] = F(_div(_arith(s['falouts'][_arith(s['k'], 1, 'add') - 1], s['falouts'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltnds'] = jnp.where(s["_pc"] == 0, s['faltnds'], old_1578)
            old_1579 = s['faltndns']
            s['faltndns'] = F(_div(_arith(s['faloutns'][_arith(s['k'], 1, 'add') - 1], s['faloutns'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndns'] = jnp.where(s["_pc"] == 0, s['faltndns'], old_1579)
            old_1580 = s['faltndnr']
            s['faltndnr'] = F(_div(_arith(s['faloutnr'][_arith(s['k'], 1, 'add') - 1], s['faloutnr'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndnr'] = jnp.where(s["_pc"] == 0, s['faltndnr'], old_1580)
            old_1581 = s['faltndc']
            s['faltndc'] = F(_div(_arith(s['faloutc'][_arith(s['k'], 1, 'add') - 1], s['faloutc'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndc'] = jnp.where(s["_pc"] == 0, s['faltndc'], old_1581)
            old_1582 = s['faltndnc']
            s['faltndnc'] = F(_div(_arith(s['faloutnc'][_arith(s['k'], 1, 'add') - 1], s['faloutnc'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndnc'] = jnp.where(s["_pc"] == 0, s['faltndnc'], old_1582)
            old_1583 = s['faltndg']
            s['faltndg'] = F(_div(_arith(s['faloutg'][_arith(s['k'], 1, 'add') - 1], s['faloutg'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndg'] = jnp.where(s["_pc"] == 0, s['faltndg'], old_1583)
            old_1584 = s['faltndng']
            s['faltndng'] = F(_div(_arith(s['faloutng'][_arith(s['k'], 1, 'add') - 1], s['faloutng'][s['k'] - 1], 'sub'), s['dzq'][s['k'] - 1]))
            s['faltndng'] = jnp.where(s["_pc"] == 0, s['faltndng'], old_1584)
            old_1585 = s['qrsten']
            s['qrsten'] = s['qrsten'].at[s['k'] - 1].set(F(_arith(s['qrsten'][s['k'] - 1], _div(_div(s['faltndr'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['qrsten'] = jnp.where(s["_pc"] == 0, s['qrsten'], old_1585)
            old_1586 = s['qisten']
            s['qisten'] = s['qisten'].at[s['k'] - 1].set(F(_arith(s['qisten'][s['k'] - 1], _div(_div(s['faltndi'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['qisten'] = jnp.where(s["_pc"] == 0, s['qisten'], old_1586)
            old_1587 = s['ni3dten']
            s['ni3dten'] = s['ni3dten'].at[s['k'] - 1].set(F(_arith(s['ni3dten'][s['k'] - 1], _div(_div(s['faltndni'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['ni3dten'] = jnp.where(s["_pc"] == 0, s['ni3dten'], old_1587)
            old_1588 = s['qnisten']
            s['qnisten'] = s['qnisten'].at[s['k'] - 1].set(F(_arith(s['qnisten'][s['k'] - 1], _div(_div(s['faltnds'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['qnisten'] = jnp.where(s["_pc"] == 0, s['qnisten'], old_1588)
            old_1589 = s['ns3dten']
            s['ns3dten'] = s['ns3dten'].at[s['k'] - 1].set(F(_arith(s['ns3dten'][s['k'] - 1], _div(_div(s['faltndns'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['ns3dten'] = jnp.where(s["_pc"] == 0, s['ns3dten'], old_1589)
            old_1590 = s['nr3dten']
            s['nr3dten'] = s['nr3dten'].at[s['k'] - 1].set(F(_arith(s['nr3dten'][s['k'] - 1], _div(_div(s['faltndnr'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['nr3dten'] = jnp.where(s["_pc"] == 0, s['nr3dten'], old_1590)
            old_1591 = s['qcsten']
            s['qcsten'] = s['qcsten'].at[s['k'] - 1].set(F(_arith(s['qcsten'][s['k'] - 1], _div(_div(s['faltndc'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['qcsten'] = jnp.where(s["_pc"] == 0, s['qcsten'], old_1591)
            old_1592 = s['nc3dten']
            s['nc3dten'] = s['nc3dten'].at[s['k'] - 1].set(F(_arith(s['nc3dten'][s['k'] - 1], _div(_div(s['faltndnc'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['nc3dten'] = jnp.where(s["_pc"] == 0, s['nc3dten'], old_1592)
            old_1593 = s['qgsten']
            s['qgsten'] = s['qgsten'].at[s['k'] - 1].set(F(_arith(s['qgsten'][s['k'] - 1], _div(_div(s['faltndg'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['qgsten'] = jnp.where(s["_pc"] == 0, s['qgsten'], old_1593)
            old_1594 = s['ng3dten']
            s['ng3dten'] = s['ng3dten'].at[s['k'] - 1].set(F(_arith(s['ng3dten'][s['k'] - 1], _div(_div(s['faltndng'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['ng3dten'] = jnp.where(s["_pc"] == 0, s['ng3dten'], old_1594)
            old_1595 = s['nisten']
            s['nisten'] = s['nisten'].at[s['k'] - 1].set(F(_arith(s['nisten'][s['k'] - 1], _div(_div(s['faltndni'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['nisten'] = jnp.where(s["_pc"] == 0, s['nisten'], old_1595)
            old_1596 = s['nssten']
            s['nssten'] = s['nssten'].at[s['k'] - 1].set(F(_arith(s['nssten'][s['k'] - 1], _div(_div(s['faltndns'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['nssten'] = jnp.where(s["_pc"] == 0, s['nssten'], old_1596)
            old_1597 = s['nrsten']
            s['nrsten'] = s['nrsten'].at[s['k'] - 1].set(F(_arith(s['nrsten'][s['k'] - 1], _div(_div(s['faltndnr'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['nrsten'] = jnp.where(s["_pc"] == 0, s['nrsten'], old_1597)
            old_1598 = s['ncsten']
            s['ncsten'] = s['ncsten'].at[s['k'] - 1].set(F(_arith(s['ncsten'][s['k'] - 1], _div(_div(s['faltndnc'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['ncsten'] = jnp.where(s["_pc"] == 0, s['ncsten'], old_1598)
            old_1599 = s['ngsten']
            s['ngsten'] = s['ngsten'].at[s['k'] - 1].set(F(_arith(s['ngsten'][s['k'] - 1], _div(_div(s['faltndng'], s['nstep']), s['rho'][s['k'] - 1]), 'add')))
            s['ngsten'] = jnp.where(s["_pc"] == 0, s['ngsten'], old_1599)
            old_1600 = s['dumr']
            s['dumr'] = s['dumr'].at[s['k'] - 1].set(F(_arith(s['dumr'][s['k'] - 1], _div(_arith(s['faltndr'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumr'] = jnp.where(s["_pc"] == 0, s['dumr'], old_1600)
            old_1601 = s['dumi']
            s['dumi'] = s['dumi'].at[s['k'] - 1].set(F(_arith(s['dumi'][s['k'] - 1], _div(_arith(s['faltndi'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumi'] = jnp.where(s["_pc"] == 0, s['dumi'], old_1601)
            old_1602 = s['dumfni']
            s['dumfni'] = s['dumfni'].at[s['k'] - 1].set(F(_arith(s['dumfni'][s['k'] - 1], _div(_arith(s['faltndni'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumfni'] = jnp.where(s["_pc"] == 0, s['dumfni'], old_1602)
            old_1603 = s['dumqs']
            s['dumqs'] = s['dumqs'].at[s['k'] - 1].set(F(_arith(s['dumqs'][s['k'] - 1], _div(_arith(s['faltnds'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumqs'] = jnp.where(s["_pc"] == 0, s['dumqs'], old_1603)
            old_1604 = s['dumfns']
            s['dumfns'] = s['dumfns'].at[s['k'] - 1].set(F(_arith(s['dumfns'][s['k'] - 1], _div(_arith(s['faltndns'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumfns'] = jnp.where(s["_pc"] == 0, s['dumfns'], old_1604)
            old_1605 = s['dumfnr']
            s['dumfnr'] = s['dumfnr'].at[s['k'] - 1].set(F(_arith(s['dumfnr'][s['k'] - 1], _div(_arith(s['faltndnr'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumfnr'] = jnp.where(s["_pc"] == 0, s['dumfnr'], old_1605)
            old_1606 = s['dumc']
            s['dumc'] = s['dumc'].at[s['k'] - 1].set(F(_arith(s['dumc'][s['k'] - 1], _div(_arith(s['faltndc'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumc'] = jnp.where(s["_pc"] == 0, s['dumc'], old_1606)
            old_1607 = s['dumfnc']
            s['dumfnc'] = s['dumfnc'].at[s['k'] - 1].set(F(_arith(s['dumfnc'][s['k'] - 1], _div(_arith(s['faltndnc'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumfnc'] = jnp.where(s["_pc"] == 0, s['dumfnc'], old_1607)
            old_1608 = s['dumg']
            s['dumg'] = s['dumg'].at[s['k'] - 1].set(F(_arith(s['dumg'][s['k'] - 1], _div(_arith(s['faltndg'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumg'] = jnp.where(s["_pc"] == 0, s['dumg'], old_1608)
            old_1609 = s['dumfng']
            s['dumfng'] = s['dumfng'].at[s['k'] - 1].set(F(_arith(s['dumfng'][s['k'] - 1], _div(_arith(s['faltndng'], s['dt'], 'mul'), s['nstep']), 'add')))
            s['dumfng'] = jnp.where(s["_pc"] == 0, s['dumfng'], old_1609)
            return s
        s = lax.fori_loop(0, I((s['kts'] - _arith(s['kte'], 1, 'sub')) // ((-1)) + 1), loop_1574, s)
        old_1610 = s['precrt']
        s['precrt'] = F(_arith(s['precrt'], _div(_arith(_arith(_arith(_arith(_arith(s['faloutr'][s['kts'] - 1], s['faloutc'][s['kts'] - 1], 'add'), s['falouts'][s['kts'] - 1], 'add'), s['falouti'][s['kts'] - 1], 'add'), s['faloutg'][s['kts'] - 1], 'add'), s['dt'], 'mul'), s['nstep']), 'add'))
        s['precrt'] = jnp.where(s["_pc"] == 0, s['precrt'], old_1610)
        old_1611 = s['snowrt']
        s['snowrt'] = F(_arith(s['snowrt'], _div(_arith(_arith(_arith(s['falouts'][s['kts'] - 1], s['falouti'][s['kts'] - 1], 'add'), s['faloutg'][s['kts'] - 1], 'add'), s['dt'], 'mul'), s['nstep']), 'add'))
        s['snowrt'] = jnp.where(s["_pc"] == 0, s['snowrt'], old_1611)
        return s
    s = lax.fori_loop(0, I((s['nstep'] - 1) // (1) + 1), loop_1526, s)
    def loop_1612(iteration, s):
        s = dict(s)
        s['k'] = I(s['kts'] + iteration * (1))
        old_1613 = s['qr3dten']
        s['qr3dten'] = s['qr3dten'].at[s['k'] - 1].set(F(_arith(s['qr3dten'][s['k'] - 1], s['qrsten'][s['k'] - 1], 'add')))
        s['qr3dten'] = jnp.where(s["_pc"] == 0, s['qr3dten'], old_1613)
        old_1614 = s['qi3dten']
        s['qi3dten'] = s['qi3dten'].at[s['k'] - 1].set(F(_arith(s['qi3dten'][s['k'] - 1], s['qisten'][s['k'] - 1], 'add')))
        s['qi3dten'] = jnp.where(s["_pc"] == 0, s['qi3dten'], old_1614)
        old_1615 = s['qc3dten']
        s['qc3dten'] = s['qc3dten'].at[s['k'] - 1].set(F(_arith(s['qc3dten'][s['k'] - 1], s['qcsten'][s['k'] - 1], 'add')))
        s['qc3dten'] = jnp.where(s["_pc"] == 0, s['qc3dten'], old_1615)
        old_1616 = s['qg3dten']
        s['qg3dten'] = s['qg3dten'].at[s['k'] - 1].set(F(_arith(s['qg3dten'][s['k'] - 1], s['qgsten'][s['k'] - 1], 'add')))
        s['qg3dten'] = jnp.where(s["_pc"] == 0, s['qg3dten'], old_1616)
        old_1617 = s['qni3dten']
        s['qni3dten'] = s['qni3dten'].at[s['k'] - 1].set(F(_arith(s['qni3dten'][s['k'] - 1], s['qnisten'][s['k'] - 1], 'add')))
        s['qni3dten'] = jnp.where(s["_pc"] == 0, s['qni3dten'], old_1617)
        old_1618 = s['qc3d']
        s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(_arith(s['qc3d'][s['k'] - 1], _arith(s['qc3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_1618)
        old_1619 = s['qi3d']
        s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(_arith(s['qi3d'][s['k'] - 1], _arith(s['qi3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_1619)
        old_1620 = s['qni3d']
        s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(_arith(s['qni3d'][s['k'] - 1], _arith(s['qni3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_1620)
        old_1621 = s['qr3d']
        s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(_arith(s['qr3d'][s['k'] - 1], _arith(s['qr3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_1621)
        old_1622 = s['nc3d']
        s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_arith(s['nc3d'][s['k'] - 1], _arith(s['nc3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1622)
        old_1623 = s['ni3d']
        s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_arith(s['ni3d'][s['k'] - 1], _arith(s['ni3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1623)
        old_1624 = s['ns3d']
        s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_arith(s['ns3d'][s['k'] - 1], _arith(s['ns3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_1624)
        old_1625 = s['nr3d']
        s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], _arith(s['nr3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1625)
        def yes_1626(s):
            s = dict(s)
            old_1627 = s['qg3d']
            s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(_arith(s['qg3d'][s['k'] - 1], _arith(s['qg3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
            s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_1627)
            old_1628 = s['ng3d']
            s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_arith(s['ng3d'][s['k'] - 1], _arith(s['ng3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
            s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_1628)
            return s
        def no_1626(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 0)), yes_1626, no_1626, s)
        old_1629 = s['t3d']
        s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _arith(s['t3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1629)
        old_1630 = s['qv3d']
        s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], _arith(s['qv3dten'][s['k'] - 1], s['dt'], 'mul'), 'add')))
        s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1630)
        old_1631 = s['evs']
        s['evs'] = s['evs'].at[s['k'] - 1].set(F(jnp.minimum(_arith(F(0.99), s['pres'][s['k'] - 1], 'mul'), POLYSVP(s['t3d'][s['k'] - 1], 0))))
        s['evs'] = jnp.where(s["_pc"] == 0, s['evs'], old_1631)
        old_1632 = s['eis']
        s['eis'] = s['eis'].at[s['k'] - 1].set(F(jnp.minimum(_arith(F(0.99), s['pres'][s['k'] - 1], 'mul'), POLYSVP(s['t3d'][s['k'] - 1], 1))))
        s['eis'] = jnp.where(s["_pc"] == 0, s['eis'], old_1632)
        def yes_1633(s):
            s = dict(s)
            old_1634 = s['eis']
            s['eis'] = s['eis'].at[s['k'] - 1].set(F(s['evs'][s['k'] - 1]))
            s['eis'] = jnp.where(s["_pc"] == 0, s['eis'], old_1634)
            return s
        def no_1633(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['eis'][s['k'] - 1] > s['evs'][s['k'] - 1])), yes_1633, no_1633, s)
        old_1635 = s['qvs']
        s['qvs'] = s['qvs'].at[s['k'] - 1].set(F(_div(_arith(s['ep_2'], s['evs'][s['k'] - 1], 'mul'), _arith(s['pres'][s['k'] - 1], s['evs'][s['k'] - 1], 'sub'))))
        s['qvs'] = jnp.where(s["_pc"] == 0, s['qvs'], old_1635)
        old_1636 = s['qvi']
        s['qvi'] = s['qvi'].at[s['k'] - 1].set(F(_div(_arith(s['ep_2'], s['eis'][s['k'] - 1], 'mul'), _arith(s['pres'][s['k'] - 1], s['eis'][s['k'] - 1], 'sub'))))
        s['qvi'] = jnp.where(s["_pc"] == 0, s['qvi'], old_1636)
        old_1637 = s['qvqvs']
        s['qvqvs'] = s['qvqvs'].at[s['k'] - 1].set(F(_div(s['qv3d'][s['k'] - 1], s['qvs'][s['k'] - 1])))
        s['qvqvs'] = jnp.where(s["_pc"] == 0, s['qvqvs'], old_1637)
        old_1638 = s['qvqvsi']
        s['qvqvsi'] = s['qvqvsi'].at[s['k'] - 1].set(F(_div(s['qv3d'][s['k'] - 1], s['qvi'][s['k'] - 1])))
        s['qvqvsi'] = jnp.where(s["_pc"] == 0, s['qvqvsi'], old_1638)
        def yes_1639(s):
            s = dict(s)
            def yes_1640(s):
                s = dict(s)
                old_1641 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1641)
                old_1642 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qr3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1642)
                old_1643 = s['qr_inst']
                s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_1643)
                old_1644 = s['qr3d']
                s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_1644)
                return s
            def no_1640(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] < F(1e-08))), yes_1640, no_1640, s)
            def yes_1645(s):
                s = dict(s)
                old_1646 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qc3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1646)
                old_1647 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qc3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1647)
                old_1648 = s['qc_inst']
                s['qc_inst'] = s['qc_inst'].at[s['k'] - 1].set(F(_arith(s['qc_inst'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qc_inst'] = jnp.where(s["_pc"] == 0, s['qc_inst'], old_1648)
                old_1649 = s['qc3d']
                s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_1649)
                return s
            def no_1645(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] < F(1e-08))), yes_1645, no_1645, s)
            return s
        def no_1639(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qvqvs'][s['k'] - 1] < F(0.9))), yes_1639, no_1639, s)
        def yes_1650(s):
            s = dict(s)
            def yes_1651(s):
                s = dict(s)
                old_1652 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qi3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1652)
                old_1653 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qi3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1653)
                old_1654 = s['qi_inst']
                s['qi_inst'] = s['qi_inst'].at[s['k'] - 1].set(F(_arith(s['qi_inst'][s['k'] - 1], _div(s['qi3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_1654)
                old_1655 = s['qi3d']
                s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_1655)
                return s
            def no_1651(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] < F(1e-08))), yes_1651, no_1651, s)
            def yes_1656(s):
                s = dict(s)
                old_1657 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qni3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1657)
                old_1658 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qni3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1658)
                old_1659 = s['qs_inst']
                s['qs_inst'] = s['qs_inst'].at[s['k'] - 1].set(F(_arith(s['qs_inst'][s['k'] - 1], _div(s['qni3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_1659)
                old_1660 = s['qni3d']
                s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_1660)
                return s
            def no_1656(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] < F(1e-08))), yes_1656, no_1656, s)
            def yes_1661(s):
                s = dict(s)
                old_1662 = s['qv3d']
                s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qg3d'][s['k'] - 1], 'add')))
                s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1662)
                old_1663 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qg3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1663)
                old_1664 = s['qg_inst']
                s['qg_inst'] = s['qg_inst'].at[s['k'] - 1].set(F(_arith(s['qg_inst'][s['k'] - 1], _div(s['qg3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_1664)
                old_1665 = s['qg3d']
                s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_1665)
                return s
            def no_1661(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] < F(1e-08))), yes_1661, no_1661, s)
            return s
        def no_1650(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qvqvsi'][s['k'] - 1] < F(0.9))), yes_1650, no_1650, s)
        def yes_1666(s):
            s = dict(s)
            old_1667 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qc3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1667)
            old_1668 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qc3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1668)
            old_1669 = s['nc_inst']
            s['nc_inst'] = s['nc_inst'].at[s['k'] - 1].set(F(_arith(s['nc_inst'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']), 'sub')))
            s['nc_inst'] = jnp.where(s["_pc"] == 0, s['nc_inst'], old_1669)
            old_1670 = s['qc_inst']
            s['qc_inst'] = s['qc_inst'].at[s['k'] - 1].set(F(_arith(s['qc_inst'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qc_inst'] = jnp.where(s["_pc"] == 0, s['qc_inst'], old_1670)
            old_1671 = s['qc3d']
            s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_1671)
            old_1672 = s['nc3d']
            s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1672)
            old_1673 = s['effc']
            s['effc'] = s['effc'].at[s['k'] - 1].set(F(F(0.0)))
            s['effc'] = jnp.where(s["_pc"] == 0, s['effc'], old_1673)
            return s
        def no_1666(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] < s['qsmall'])), yes_1666, no_1666, s)
        def yes_1674(s):
            s = dict(s)
            old_1675 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1675)
            old_1676 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qr3d'][s['k'] - 1], s['xxlv'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1676)
            old_1677 = s['nr_inst']
            s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'sub')))
            s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_1677)
            old_1678 = s['qr_inst']
            s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_1678)
            old_1679 = s['qr3d']
            s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_1679)
            old_1680 = s['nr3d']
            s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1680)
            old_1681 = s['effr']
            s['effr'] = s['effr'].at[s['k'] - 1].set(F(F(0.0)))
            s['effr'] = jnp.where(s["_pc"] == 0, s['effr'], old_1681)
            return s
        def no_1674(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] < s['qsmall'])), yes_1674, no_1674, s)
        def yes_1682(s):
            s = dict(s)
            old_1683 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qi3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1683)
            old_1684 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qi3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1684)
            old_1685 = s['ni_inst']
            s['ni_inst'] = s['ni_inst'].at[s['k'] - 1].set(F(_arith(s['ni_inst'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']), 'sub')))
            s['ni_inst'] = jnp.where(s["_pc"] == 0, s['ni_inst'], old_1685)
            old_1686 = s['qi_inst']
            s['qi_inst'] = s['qi_inst'].at[s['k'] - 1].set(F(_arith(s['qi_inst'][s['k'] - 1], _div(s['qi3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_1686)
            old_1687 = s['qi3d']
            s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_1687)
            old_1688 = s['ni3d']
            s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1688)
            old_1689 = s['effi']
            s['effi'] = s['effi'].at[s['k'] - 1].set(F(F(0.0)))
            s['effi'] = jnp.where(s["_pc"] == 0, s['effi'], old_1689)
            return s
        def no_1682(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] < s['qsmall'])), yes_1682, no_1682, s)
        def yes_1690(s):
            s = dict(s)
            old_1691 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qni3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1691)
            old_1692 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qni3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1692)
            old_1693 = s['ns_inst']
            s['ns_inst'] = s['ns_inst'].at[s['k'] - 1].set(F(_arith(s['ns_inst'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']), 'sub')))
            s['ns_inst'] = jnp.where(s["_pc"] == 0, s['ns_inst'], old_1693)
            old_1694 = s['qs_inst']
            s['qs_inst'] = s['qs_inst'].at[s['k'] - 1].set(F(_arith(s['qs_inst'][s['k'] - 1], _div(s['qni3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_1694)
            old_1695 = s['qni3d']
            s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_1695)
            old_1696 = s['ns3d']
            s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_1696)
            old_1697 = s['effs']
            s['effs'] = s['effs'].at[s['k'] - 1].set(F(F(0.0)))
            s['effs'] = jnp.where(s["_pc"] == 0, s['effs'], old_1697)
            return s
        def no_1690(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] < s['qsmall'])), yes_1690, no_1690, s)
        def yes_1698(s):
            s = dict(s)
            old_1699 = s['qv3d']
            s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(F(_arith(s['qv3d'][s['k'] - 1], s['qg3d'][s['k'] - 1], 'add')))
            s['qv3d'] = jnp.where(s["_pc"] == 0, s['qv3d'], old_1699)
            old_1700 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qg3d'][s['k'] - 1], s['xxls'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1700)
            old_1701 = s['ng_inst']
            s['ng_inst'] = s['ng_inst'].at[s['k'] - 1].set(F(_arith(s['ng_inst'][s['k'] - 1], _div(s['ng3d'][s['k'] - 1], s['dt']), 'sub')))
            s['ng_inst'] = jnp.where(s["_pc"] == 0, s['ng_inst'], old_1701)
            old_1702 = s['qg_inst']
            s['qg_inst'] = s['qg_inst'].at[s['k'] - 1].set(F(_arith(s['qg_inst'][s['k'] - 1], _div(s['qg3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_1702)
            old_1703 = s['qg3d']
            s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_1703)
            old_1704 = s['ng3d']
            s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_1704)
            old_1705 = s['effg']
            s['effg'] = s['effg'].at[s['k'] - 1].set(F(F(0.0)))
            s['effg'] = jnp.where(s["_pc"] == 0, s['effg'], old_1705)
            return s
        def no_1698(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] < s['qsmall'])), yes_1698, no_1698, s)
        def yes_1706(s):
            s = dict(s)
            s["_pc"] = jnp.where(s["_pc"] == 0, I(500), s["_pc"])
            return s
        def no_1706(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((((s['qc3d'][s['k'] - 1] < s['qsmall'])) & ((s['qi3d'][s['k'] - 1] < s['qsmall'])) & ((s['qni3d'][s['k'] - 1] < s['qsmall'])) & ((s['qr3d'][s['k'] - 1] < s['qsmall'])) & ((s['qg3d'][s['k'] - 1] < s['qsmall'])))), yes_1706, no_1706, s)
        def yes_1707(s):
            s = dict(s)
            old_1708 = s['qr_inst']
            s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qi3d'][s['k'] - 1], s['dt']), 'add')))
            s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_1708)
            old_1709 = s['qi_inst']
            s['qi_inst'] = s['qi_inst'].at[s['k'] - 1].set(F(_arith(s['qi_inst'][s['k'] - 1], _div(s['qi3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_1709)
            old_1710 = s['nr_inst']
            s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']), 'add')))
            s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_1710)
            old_1711 = s['ni_inst']
            s['ni_inst'] = s['ni_inst'].at[s['k'] - 1].set(F(_arith(s['ni_inst'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']), 'sub')))
            s['ni_inst'] = jnp.where(s["_pc"] == 0, s['ni_inst'], old_1711)
            old_1712 = s['qr3d']
            s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(_arith(s['qr3d'][s['k'] - 1], s['qi3d'][s['k'] - 1], 'add')))
            s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_1712)
            old_1713 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qi3d'][s['k'] - 1], s['xlf'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'sub')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1713)
            old_1714 = s['qi3d']
            s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_1714)
            old_1715 = s['nr3d']
            s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_arith(s['nr3d'][s['k'] - 1], s['ni3d'][s['k'] - 1], 'add')))
            s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1715)
            old_1716 = s['ni3d']
            s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1716)
            return s
        def no_1707(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((((s['qi3d'][s['k'] - 1] >= s['qsmall'])) & ((s['t3d'][s['k'] - 1] >= s['tmelt'])))), yes_1707, no_1707, s)
        def yes_1717(s):
            s = dict(s)
            s["_pc"] = jnp.where(s["_pc"] == 0, I(778), s["_pc"])
            return s
        def no_1717(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['iliq'] == 1)), yes_1717, no_1717, s)
        def yes_1718(s):
            s = dict(s)
            old_1719 = s['qi_inst']
            s['qi_inst'] = s['qi_inst'].at[s['k'] - 1].set(F(_arith(s['qi_inst'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['dt']), 'add')))
            s['qi_inst'] = jnp.where(s["_pc"] == 0, s['qi_inst'], old_1719)
            old_1720 = s['qc_inst']
            s['qc_inst'] = s['qc_inst'].at[s['k'] - 1].set(F(_arith(s['qc_inst'][s['k'] - 1], _div(s['qc3d'][s['k'] - 1], s['dt']), 'sub')))
            s['qc_inst'] = jnp.where(s["_pc"] == 0, s['qc_inst'], old_1720)
            old_1721 = s['ni_inst']
            s['ni_inst'] = s['ni_inst'].at[s['k'] - 1].set(F(_arith(s['ni_inst'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']), 'add')))
            s['ni_inst'] = jnp.where(s["_pc"] == 0, s['ni_inst'], old_1721)
            old_1722 = s['nc_inst']
            s['nc_inst'] = s['nc_inst'].at[s['k'] - 1].set(F(_arith(s['nc_inst'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']), 'sub')))
            s['nc_inst'] = jnp.where(s["_pc"] == 0, s['nc_inst'], old_1722)
            old_1723 = s['qi3d']
            s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(F(_arith(s['qi3d'][s['k'] - 1], s['qc3d'][s['k'] - 1], 'add')))
            s['qi3d'] = jnp.where(s["_pc"] == 0, s['qi3d'], old_1723)
            old_1724 = s['t3d']
            s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qc3d'][s['k'] - 1], s['xlf'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'add')))
            s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1724)
            old_1725 = s['qc3d']
            s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['qc3d'] = jnp.where(s["_pc"] == 0, s['qc3d'], old_1725)
            old_1726 = s['ni3d']
            s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_arith(s['ni3d'][s['k'] - 1], s['nc3d'][s['k'] - 1], 'add')))
            s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1726)
            old_1727 = s['nc3d']
            s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(F(0.0)))
            s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1727)
            return s
        def no_1718(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] <= F(233.15))) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])))), yes_1718, no_1718, s)
        def yes_1728(s):
            s = dict(s)
            def yes_1729(s):
                s = dict(s)
                old_1730 = s['qg_inst']
                s['qg_inst'] = s['qg_inst'].at[s['k'] - 1].set(F(_arith(s['qg_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'add')))
                s['qg_inst'] = jnp.where(s["_pc"] == 0, s['qg_inst'], old_1730)
                old_1731 = s['qr_inst']
                s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'sub')))
                s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_1731)
                old_1732 = s['ng_inst']
                s['ng_inst'] = s['ng_inst'].at[s['k'] - 1].set(F(_arith(s['ng_inst'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'add')))
                s['ng_inst'] = jnp.where(s["_pc"] == 0, s['ng_inst'], old_1732)
                old_1733 = s['nr_inst']
                s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'sub')))
                s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_1733)
                old_1734 = s['qg3d']
                s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(F(_arith(s['qg3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')))
                s['qg3d'] = jnp.where(s["_pc"] == 0, s['qg3d'], old_1734)
                old_1735 = s['t3d']
                s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qr3d'][s['k'] - 1], s['xlf'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'add')))
                s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1735)
                old_1736 = s['qr3d']
                s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_1736)
                old_1737 = s['ng3d']
                s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_arith(s['ng3d'][s['k'] - 1], s['nr3d'][s['k'] - 1], 'add')))
                s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_1737)
                old_1738 = s['nr3d']
                s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(F(0.0)))
                s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1738)
                return s
            def no_1729(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] <= F(233.15))) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])))), yes_1729, no_1729, s)
            return s
        def no_1728(s):
            s = dict(s)
            def yes_1739(s):
                s = dict(s)
                def yes_1740(s):
                    s = dict(s)
                    old_1741 = s['qs_inst']
                    s['qs_inst'] = s['qs_inst'].at[s['k'] - 1].set(F(_arith(s['qs_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'add')))
                    s['qs_inst'] = jnp.where(s["_pc"] == 0, s['qs_inst'], old_1741)
                    old_1742 = s['qr_inst']
                    s['qr_inst'] = s['qr_inst'].at[s['k'] - 1].set(F(_arith(s['qr_inst'][s['k'] - 1], _div(s['qr3d'][s['k'] - 1], s['dt']), 'sub')))
                    s['qr_inst'] = jnp.where(s["_pc"] == 0, s['qr_inst'], old_1742)
                    old_1743 = s['ns_inst']
                    s['ns_inst'] = s['ns_inst'].at[s['k'] - 1].set(F(_arith(s['ns_inst'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'add')))
                    s['ns_inst'] = jnp.where(s["_pc"] == 0, s['ns_inst'], old_1743)
                    old_1744 = s['nr_inst']
                    s['nr_inst'] = s['nr_inst'].at[s['k'] - 1].set(F(_arith(s['nr_inst'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'sub')))
                    s['nr_inst'] = jnp.where(s["_pc"] == 0, s['nr_inst'], old_1744)
                    old_1745 = s['qni3d']
                    s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(F(_arith(s['qni3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], 'add')))
                    s['qni3d'] = jnp.where(s["_pc"] == 0, s['qni3d'], old_1745)
                    old_1746 = s['t3d']
                    s['t3d'] = s['t3d'].at[s['k'] - 1].set(F(_arith(s['t3d'][s['k'] - 1], _div(_arith(s['qr3d'][s['k'] - 1], s['xlf'][s['k'] - 1], 'mul'), s['cpm'][s['k'] - 1]), 'add')))
                    s['t3d'] = jnp.where(s["_pc"] == 0, s['t3d'], old_1746)
                    old_1747 = s['qr3d']
                    s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(F(F(0.0)))
                    s['qr3d'] = jnp.where(s["_pc"] == 0, s['qr3d'], old_1747)
                    old_1748 = s['ns3d']
                    s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_arith(s['ns3d'][s['k'] - 1], s['nr3d'][s['k'] - 1], 'add')))
                    s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_1748)
                    old_1749 = s['nr3d']
                    s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(F(0.0)))
                    s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1749)
                    return s
                def no_1740(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((((s['t3d'][s['k'] - 1] <= F(233.15))) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])))), yes_1740, no_1740, s)
                return s
            def no_1739(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 1)), yes_1739, no_1739, s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['igraup'] == 0)), yes_1728, no_1728, s)
        s["_pc"] = jnp.where(s["_pc"] == 778, I(0), s["_pc"])
        def yes_1750(s):
            s = dict(s)
            old_1751 = s['negfix_ni']
            s['negfix_ni'] = s['negfix_ni'].at[s['k'] - 1].set(F(_arith(s['negfix_ni'][s['k'] - 1], _div(s['ni3d'][s['k'] - 1], s['dt']), 'add')))
            s['negfix_ni'] = jnp.where(s["_pc"] == 0, s['negfix_ni'], old_1751)
            return s
        def no_1750(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['ni3d'][s['k'] - 1] < F(0.0))), yes_1750, no_1750, s)
        def yes_1752(s):
            s = dict(s)
            old_1753 = s['negfix_ns']
            s['negfix_ns'] = s['negfix_ns'].at[s['k'] - 1].set(F(_arith(s['negfix_ns'][s['k'] - 1], _div(s['ns3d'][s['k'] - 1], s['dt']), 'add')))
            s['negfix_ns'] = jnp.where(s["_pc"] == 0, s['negfix_ns'], old_1753)
            return s
        def no_1752(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['ns3d'][s['k'] - 1] < F(0.0))), yes_1752, no_1752, s)
        def yes_1754(s):
            s = dict(s)
            old_1755 = s['negfix_nc']
            s['negfix_nc'] = s['negfix_nc'].at[s['k'] - 1].set(F(_arith(s['negfix_nc'][s['k'] - 1], _div(s['nc3d'][s['k'] - 1], s['dt']), 'add')))
            s['negfix_nc'] = jnp.where(s["_pc"] == 0, s['negfix_nc'], old_1755)
            return s
        def no_1754(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['nc3d'][s['k'] - 1] < F(0.0))), yes_1754, no_1754, s)
        def yes_1756(s):
            s = dict(s)
            old_1757 = s['negfix_nr']
            s['negfix_nr'] = s['negfix_nr'].at[s['k'] - 1].set(F(_arith(s['negfix_nr'][s['k'] - 1], _div(s['nr3d'][s['k'] - 1], s['dt']), 'add')))
            s['negfix_nr'] = jnp.where(s["_pc"] == 0, s['negfix_nr'], old_1757)
            return s
        def no_1756(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['nr3d'][s['k'] - 1] < F(0.0))), yes_1756, no_1756, s)
        def yes_1758(s):
            s = dict(s)
            old_1759 = s['negfix_ng']
            s['negfix_ng'] = s['negfix_ng'].at[s['k'] - 1].set(F(_arith(s['negfix_ng'][s['k'] - 1], _div(s['ng3d'][s['k'] - 1], s['dt']), 'add')))
            s['negfix_ng'] = jnp.where(s["_pc"] == 0, s['negfix_ng'], old_1759)
            return s
        def no_1758(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['ng3d'][s['k'] - 1] < F(0.0))), yes_1758, no_1758, s)
        old_1760 = s['ni3d']
        s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ni3d'][s['k'] - 1])))
        s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1760)
        old_1761 = s['ns3d']
        s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ns3d'][s['k'] - 1])))
        s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_1761)
        old_1762 = s['nc3d']
        s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['nc3d'][s['k'] - 1])))
        s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1762)
        old_1763 = s['nr3d']
        s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['nr3d'][s['k'] - 1])))
        s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1763)
        old_1764 = s['ng3d']
        s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(jnp.maximum(F(0.0), s['ng3d'][s['k'] - 1])))
        s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_1764)
        def yes_1765(s):
            s = dict(s)
            old_1766 = s['lami']
            s['lami'] = s['lami'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons12'], s['ni3d'][s['k'] - 1], 'mul'), s['qi3d'][s['k'] - 1]), _div(F(1.0), s['di']), 'pow')))
            s['lami'] = jnp.where(s["_pc"] == 0, s['lami'], old_1766)
            old_1767 = s['tmpnum']
            s['tmpnum'] = F(F(0.0))
            s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1767)
            def yes_1768(s):
                s = dict(s)
                old_1769 = s['lami']
                s['lami'] = s['lami'].at[s['k'] - 1].set(F(s['lammini']))
                s['lami'] = jnp.where(s["_pc"] == 0, s['lami'], old_1769)
                old_1770 = s['n0i']
                s['n0i'] = s['n0i'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lami'][s['k'] - 1], _arith(s['di'], F(1.0), 'add'), 'pow'), s['qi3d'][s['k'] - 1], 'mul'), s['cons12'])))
                s['n0i'] = jnp.where(s["_pc"] == 0, s['n0i'], old_1770)
                old_1771 = s['tmpnum']
                s['tmpnum'] = F(s['ni3d'][s['k'] - 1])
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1771)
                old_1772 = s['ni3d']
                s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_div(s['n0i'][s['k'] - 1], s['lami'][s['k'] - 1])))
                s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1772)
                old_1773 = s['sizefix_ni']
                s['sizefix_ni'] = s['sizefix_ni'].at[s['k'] - 1].set(F(_arith(s['sizefix_ni'][s['k'] - 1], _div(_arith(s['ni3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                s['sizefix_ni'] = jnp.where(s["_pc"] == 0, s['sizefix_ni'], old_1773)
                return s
            def no_1768(s):
                s = dict(s)
                def yes_1774(s):
                    s = dict(s)
                    old_1775 = s['lami']
                    s['lami'] = s['lami'].at[s['k'] - 1].set(F(s['lammaxi']))
                    s['lami'] = jnp.where(s["_pc"] == 0, s['lami'], old_1775)
                    old_1776 = s['n0i']
                    s['n0i'] = s['n0i'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lami'][s['k'] - 1], _arith(s['di'], F(1.0), 'add'), 'pow'), s['qi3d'][s['k'] - 1], 'mul'), s['cons12'])))
                    s['n0i'] = jnp.where(s["_pc"] == 0, s['n0i'], old_1776)
                    old_1777 = s['tmpnum']
                    s['tmpnum'] = F(s['ni3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1777)
                    old_1778 = s['ni3d']
                    s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(_div(s['n0i'][s['k'] - 1], s['lami'][s['k'] - 1])))
                    s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1778)
                    old_1779 = s['sizefix_ni']
                    s['sizefix_ni'] = s['sizefix_ni'].at[s['k'] - 1].set(F(_arith(s['sizefix_ni'][s['k'] - 1], _div(_arith(s['ni3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ni'] = jnp.where(s["_pc"] == 0, s['sizefix_ni'], old_1779)
                    return s
                def no_1774(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lami'][s['k'] - 1] > s['lammaxi'])), yes_1774, no_1774, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['lami'][s['k'] - 1] < s['lammini'])), yes_1768, no_1768, s)
            return s
        def no_1765(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] >= s['qsmall'])), yes_1765, no_1765, s)
        def yes_1780(s):
            s = dict(s)
            old_1781 = s['lamr']
            s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(s['pi'], s['rhow'], 'mul'), s['nr3d'][s['k'] - 1], 'mul'), s['qr3d'][s['k'] - 1]), _div(F(1.0), F(3.0)), 'pow')))
            s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_1781)
            old_1782 = s['tmpnum']
            s['tmpnum'] = F(F(0.0))
            s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1782)
            def yes_1783(s):
                s = dict(s)
                old_1784 = s['tmpnum']
                s['tmpnum'] = F(s['nr3d'][s['k'] - 1])
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1784)
                old_1785 = s['lamr']
                s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(s['lamminr']))
                s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_1785)
                old_1786 = s['n0rr']
                s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamr'][s['k'] - 1], 4, 'pow'), s['qr3d'][s['k'] - 1], 'mul'), _arith(s['pi'], s['rhow'], 'mul'))))
                s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_1786)
                old_1787 = s['nr3d']
                s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_div(s['n0rr'][s['k'] - 1], s['lamr'][s['k'] - 1])))
                s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1787)
                old_1788 = s['sizefix_nr']
                s['sizefix_nr'] = s['sizefix_nr'].at[s['k'] - 1].set(F(_arith(s['sizefix_nr'][s['k'] - 1], _div(_arith(s['nr3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_1788)
                return s
            def no_1783(s):
                s = dict(s)
                def yes_1789(s):
                    s = dict(s)
                    old_1790 = s['tmpnum']
                    s['tmpnum'] = F(s['nr3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1790)
                    old_1791 = s['lamr']
                    s['lamr'] = s['lamr'].at[s['k'] - 1].set(F(s['lammaxr']))
                    s['lamr'] = jnp.where(s["_pc"] == 0, s['lamr'], old_1791)
                    old_1792 = s['n0rr']
                    s['n0rr'] = s['n0rr'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamr'][s['k'] - 1], 4, 'pow'), s['qr3d'][s['k'] - 1], 'mul'), _arith(s['pi'], s['rhow'], 'mul'))))
                    s['n0rr'] = jnp.where(s["_pc"] == 0, s['n0rr'], old_1792)
                    old_1793 = s['nr3d']
                    s['nr3d'] = s['nr3d'].at[s['k'] - 1].set(F(_div(s['n0rr'][s['k'] - 1], s['lamr'][s['k'] - 1])))
                    s['nr3d'] = jnp.where(s["_pc"] == 0, s['nr3d'], old_1793)
                    old_1794 = s['sizefix_nr']
                    s['sizefix_nr'] = s['sizefix_nr'].at[s['k'] - 1].set(F(_arith(s['sizefix_nr'][s['k'] - 1], _div(_arith(s['nr3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_nr'] = jnp.where(s["_pc"] == 0, s['sizefix_nr'], old_1794)
                    return s
                def no_1789(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamr'][s['k'] - 1] > s['lammaxr'])), yes_1789, no_1789, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['lamr'][s['k'] - 1] < s['lamminr'])), yes_1783, no_1783, s)
            return s
        def no_1780(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= s['qsmall'])), yes_1780, no_1780, s)
        def yes_1795(s):
            s = dict(s)
            def yes_1796(s):
                s = dict(s)
                old_1797 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(s['pgam_fixed']))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1797)
                return s
            def no_1796(s):
                s = dict(s)
                old_1798 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_arith(F(0.0005714), _arith(_div(s['nc3d'][s['k'] - 1], F(1000000.0)), s['rho'][s['k'] - 1], 'mul'), 'mul'), F(0.2714), 'add')))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1798)
                old_1799 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(_arith(_div(F(1.0), _arith(s['pgam'][s['k'] - 1], 2, 'pow')), F(1.0), 'sub')))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1799)
                old_1800 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.maximum(s['pgam'][s['k'] - 1], F(2.0))))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1800)
                old_1801 = s['pgam']
                s['pgam'] = s['pgam'].at[s['k'] - 1].set(F(jnp.minimum(s['pgam'][s['k'] - 1], F(10.0))))
                s['pgam'] = jnp.where(s["_pc"] == 0, s['pgam'], old_1801)
                return s
            s = lax.cond((s["_pc"] == 0) & (s['dofix_pgam']), yes_1796, no_1796, s)
            old_1802 = s['lamc']
            s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(_arith(_div(_arith(_arith(s['cons26'], s['nc3d'][s['k'] - 1], 'mul'), GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add')), 'mul'), _arith(s['qc3d'][s['k'] - 1], GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add')), 'mul')), _div(F(1.0), F(3.0)), 'pow')))
            s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_1802)
            old_1803 = s['lammin']
            s['lammin'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(6e-05)))
            s['lammin'] = jnp.where(s["_pc"] == 0, s['lammin'], old_1803)
            old_1804 = s['lammax']
            s['lammax'] = F(_div(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'), F(1e-06)))
            s['lammax'] = jnp.where(s["_pc"] == 0, s['lammax'], old_1804)
            old_1805 = s['tmpnum']
            s['tmpnum'] = F(F(0.0))
            s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1805)
            def yes_1806(s):
                s = dict(s)
                old_1807 = s['lamc']
                s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(s['lammin']))
                s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_1807)
                old_1808 = s['tmpnum']
                s['tmpnum'] = F(s['nc3d'][s['k'] - 1])
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1808)
                old_1809 = s['nc3d']
                s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_intrinsic('exp', _arith(_arith(_arith(_arith(F(3.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), _intrinsic('log', s['qc3d'][s['k'] - 1]), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'))), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add'))), 'sub')), s['cons26'])))
                s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1809)
                old_1810 = s['sizefix_nc']
                s['sizefix_nc'] = s['sizefix_nc'].at[s['k'] - 1].set(F(_arith(s['sizefix_nc'][s['k'] - 1], _div(_arith(s['nc3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_1810)
                return s
            def no_1806(s):
                s = dict(s)
                def yes_1811(s):
                    s = dict(s)
                    old_1812 = s['lamc']
                    s['lamc'] = s['lamc'].at[s['k'] - 1].set(F(s['lammax']))
                    s['lamc'] = jnp.where(s["_pc"] == 0, s['lamc'], old_1812)
                    old_1813 = s['tmpnum']
                    s['tmpnum'] = F(s['nc3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1813)
                    old_1814 = s['nc3d']
                    s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_intrinsic('exp', _arith(_arith(_arith(_arith(F(3.0), _intrinsic('log', s['lamc'][s['k'] - 1]), 'mul'), _intrinsic('log', s['qc3d'][s['k'] - 1]), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(1.0), 'add'))), 'add'), _intrinsic('log', GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add'))), 'sub')), s['cons26'])))
                    s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1814)
                    old_1815 = s['sizefix_nc']
                    s['sizefix_nc'] = s['sizefix_nc'].at[s['k'] - 1].set(F(_arith(s['sizefix_nc'][s['k'] - 1], _div(_arith(s['nc3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_nc'] = jnp.where(s["_pc"] == 0, s['sizefix_nc'], old_1815)
                    return s
                def no_1811(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamc'][s['k'] - 1] > s['lammax'])), yes_1811, no_1811, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['lamc'][s['k'] - 1] < s['lammin'])), yes_1806, no_1806, s)
            return s
        def no_1795(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] >= s['qsmall'])), yes_1795, no_1795, s)
        def yes_1816(s):
            s = dict(s)
            old_1817 = s['lams']
            s['lams'] = s['lams'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons1'], s['ns3d'][s['k'] - 1], 'mul'), s['qni3d'][s['k'] - 1]), _div(F(1.0), s['ds']), 'pow')))
            s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_1817)
            old_1818 = s['tmpnum']
            s['tmpnum'] = F(F(0.0))
            s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1818)
            def yes_1819(s):
                s = dict(s)
                old_1820 = s['lams']
                s['lams'] = s['lams'].at[s['k'] - 1].set(F(s['lammins']))
                s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_1820)
                old_1821 = s['n0s']
                s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lams'][s['k'] - 1], _arith(s['ds'], F(1.0), 'add'), 'pow'), s['qni3d'][s['k'] - 1], 'mul'), s['cons1'])))
                s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_1821)
                old_1822 = s['tmpnum']
                s['tmpnum'] = F(s['ns3d'][s['k'] - 1])
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1822)
                old_1823 = s['ns3d']
                s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_div(s['n0s'][s['k'] - 1], s['lams'][s['k'] - 1])))
                s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_1823)
                old_1824 = s['sizefix_ns']
                s['sizefix_ns'] = s['sizefix_ns'].at[s['k'] - 1].set(F(_arith(s['sizefix_ns'][s['k'] - 1], _div(_arith(s['ns3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_1824)
                return s
            def no_1819(s):
                s = dict(s)
                def yes_1825(s):
                    s = dict(s)
                    old_1826 = s['lams']
                    s['lams'] = s['lams'].at[s['k'] - 1].set(F(s['lammaxs']))
                    s['lams'] = jnp.where(s["_pc"] == 0, s['lams'], old_1826)
                    old_1827 = s['n0s']
                    s['n0s'] = s['n0s'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lams'][s['k'] - 1], _arith(s['ds'], F(1.0), 'add'), 'pow'), s['qni3d'][s['k'] - 1], 'mul'), s['cons1'])))
                    s['n0s'] = jnp.where(s["_pc"] == 0, s['n0s'], old_1827)
                    old_1828 = s['tmpnum']
                    s['tmpnum'] = F(s['ns3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1828)
                    old_1829 = s['ns3d']
                    s['ns3d'] = s['ns3d'].at[s['k'] - 1].set(F(_div(s['n0s'][s['k'] - 1], s['lams'][s['k'] - 1])))
                    s['ns3d'] = jnp.where(s["_pc"] == 0, s['ns3d'], old_1829)
                    old_1830 = s['sizefix_ns']
                    s['sizefix_ns'] = s['sizefix_ns'].at[s['k'] - 1].set(F(_arith(s['sizefix_ns'][s['k'] - 1], _div(_arith(s['ns3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ns'] = jnp.where(s["_pc"] == 0, s['sizefix_ns'], old_1830)
                    return s
                def no_1825(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lams'][s['k'] - 1] > s['lammaxs'])), yes_1825, no_1825, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['lams'][s['k'] - 1] < s['lammins'])), yes_1819, no_1819, s)
            return s
        def no_1816(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= s['qsmall'])), yes_1816, no_1816, s)
        def yes_1831(s):
            s = dict(s)
            old_1832 = s['lamg']
            s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(_arith(_div(_arith(s['cons2'], s['ng3d'][s['k'] - 1], 'mul'), s['qg3d'][s['k'] - 1]), _div(F(1.0), s['dg']), 'pow')))
            s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_1832)
            old_1833 = s['tmpnum']
            s['tmpnum'] = F(F(0.0))
            s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1833)
            def yes_1834(s):
                s = dict(s)
                old_1835 = s['lamg']
                s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(s['lamming']))
                s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_1835)
                old_1836 = s['n0g']
                s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamg'][s['k'] - 1], _arith(s['dg'], F(1.0), 'add'), 'pow'), s['qg3d'][s['k'] - 1], 'mul'), s['cons2'])))
                s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_1836)
                old_1837 = s['tmpnum']
                s['tmpnum'] = F(s['ng3d'][s['k'] - 1])
                s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1837)
                old_1838 = s['ng3d']
                s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_div(s['n0g'][s['k'] - 1], s['lamg'][s['k'] - 1])))
                s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_1838)
                old_1839 = s['sizefix_ng']
                s['sizefix_ng'] = s['sizefix_ng'].at[s['k'] - 1].set(F(_arith(s['sizefix_ng'][s['k'] - 1], _div(_arith(s['ng3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_1839)
                return s
            def no_1834(s):
                s = dict(s)
                def yes_1840(s):
                    s = dict(s)
                    old_1841 = s['lamg']
                    s['lamg'] = s['lamg'].at[s['k'] - 1].set(F(s['lammaxg']))
                    s['lamg'] = jnp.where(s["_pc"] == 0, s['lamg'], old_1841)
                    old_1842 = s['n0g']
                    s['n0g'] = s['n0g'].at[s['k'] - 1].set(F(_div(_arith(_arith(s['lamg'][s['k'] - 1], _arith(s['dg'], F(1.0), 'add'), 'pow'), s['qg3d'][s['k'] - 1], 'mul'), s['cons2'])))
                    s['n0g'] = jnp.where(s["_pc"] == 0, s['n0g'], old_1842)
                    old_1843 = s['tmpnum']
                    s['tmpnum'] = F(s['ng3d'][s['k'] - 1])
                    s['tmpnum'] = jnp.where(s["_pc"] == 0, s['tmpnum'], old_1843)
                    old_1844 = s['ng3d']
                    s['ng3d'] = s['ng3d'].at[s['k'] - 1].set(F(_div(s['n0g'][s['k'] - 1], s['lamg'][s['k'] - 1])))
                    s['ng3d'] = jnp.where(s["_pc"] == 0, s['ng3d'], old_1844)
                    old_1845 = s['sizefix_ng']
                    s['sizefix_ng'] = s['sizefix_ng'].at[s['k'] - 1].set(F(_arith(s['sizefix_ng'][s['k'] - 1], _div(_arith(s['ng3d'][s['k'] - 1], s['tmpnum'], 'sub'), s['dt']), 'add')))
                    s['sizefix_ng'] = jnp.where(s["_pc"] == 0, s['sizefix_ng'], old_1845)
                    return s
                def no_1840(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['lamg'][s['k'] - 1] > s['lammaxg'])), yes_1840, no_1840, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['lamg'][s['k'] - 1] < s['lamming'])), yes_1834, no_1834, s)
            return s
        def no_1831(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] >= s['qsmall'])), yes_1831, no_1831, s)
        s["_pc"] = jnp.where(s["_pc"] == 500, I(0), s["_pc"])
        def yes_1846(s):
            s = dict(s)
            old_1847 = s['tmpqsmall']
            s['tmpqsmall'] = F(_div(s['qsmall'], s['cf3d'][s['k'] - 1]))
            s['tmpqsmall'] = jnp.where(s["_pc"] == 0, s['tmpqsmall'], old_1847)
            return s
        def no_1846(s):
            s = dict(s)
            old_1848 = s['tmpqsmall']
            s['tmpqsmall'] = F(s['qsmall'])
            s['tmpqsmall'] = jnp.where(s["_pc"] == 0, s['tmpqsmall'], old_1848)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['cf3d'][s['k'] - 1] > s['cloud_frac_thresh'])), yes_1846, no_1846, s)
        def yes_1849(s):
            s = dict(s)
            old_1850 = s['effi']
            s['effi'] = s['effi'].at[s['k'] - 1].set(F(_arith(_div(_div(F(3.0), s['lami'][s['k'] - 1]), F(2.0)), F(1000000.0), 'mul')))
            s['effi'] = jnp.where(s["_pc"] == 0, s['effi'], old_1850)
            return s
        def no_1849(s):
            s = dict(s)
            old_1851 = s['effi']
            s['effi'] = s['effi'].at[s['k'] - 1].set(F(F(25.0)))
            s['effi'] = jnp.where(s["_pc"] == 0, s['effi'], old_1851)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qi3d'][s['k'] - 1] >= s['tmpqsmall'])), yes_1849, no_1849, s)
        def yes_1852(s):
            s = dict(s)
            old_1853 = s['effs']
            s['effs'] = s['effs'].at[s['k'] - 1].set(F(_arith(_div(_div(F(3.0), s['lams'][s['k'] - 1]), F(2.0)), F(1000000.0), 'mul')))
            s['effs'] = jnp.where(s["_pc"] == 0, s['effs'], old_1853)
            return s
        def no_1852(s):
            s = dict(s)
            old_1854 = s['effs']
            s['effs'] = s['effs'].at[s['k'] - 1].set(F(F(25.0)))
            s['effs'] = jnp.where(s["_pc"] == 0, s['effs'], old_1854)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qni3d'][s['k'] - 1] >= s['tmpqsmall'])), yes_1852, no_1852, s)
        def yes_1855(s):
            s = dict(s)
            old_1856 = s['effr']
            s['effr'] = s['effr'].at[s['k'] - 1].set(F(_arith(_div(_div(F(3.0), s['lamr'][s['k'] - 1]), F(2.0)), F(1000000.0), 'mul')))
            s['effr'] = jnp.where(s["_pc"] == 0, s['effr'], old_1856)
            return s
        def no_1855(s):
            s = dict(s)
            old_1857 = s['effr']
            s['effr'] = s['effr'].at[s['k'] - 1].set(F(F(25.0)))
            s['effr'] = jnp.where(s["_pc"] == 0, s['effr'], old_1857)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qr3d'][s['k'] - 1] >= s['tmpqsmall'])), yes_1855, no_1855, s)
        def yes_1858(s):
            s = dict(s)
            old_1859 = s['effc']
            s['effc'] = s['effc'].at[s['k'] - 1].set(F(_arith(_div(_div(_div(GAMMA(_arith(s['pgam'][s['k'] - 1], F(4.0), 'add')), GAMMA(_arith(s['pgam'][s['k'] - 1], F(3.0), 'add'))), s['lamc'][s['k'] - 1]), F(2.0)), F(1000000.0), 'mul')))
            s['effc'] = jnp.where(s["_pc"] == 0, s['effc'], old_1859)
            return s
        def no_1858(s):
            s = dict(s)
            old_1860 = s['effc']
            s['effc'] = s['effc'].at[s['k'] - 1].set(F(F(25.0)))
            s['effc'] = jnp.where(s["_pc"] == 0, s['effc'], old_1860)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qc3d'][s['k'] - 1] >= s['tmpqsmall'])), yes_1858, no_1858, s)
        def yes_1861(s):
            s = dict(s)
            old_1862 = s['effg']
            s['effg'] = s['effg'].at[s['k'] - 1].set(F(_arith(_div(_div(F(3.0), s['lamg'][s['k'] - 1]), F(2.0)), F(1000000.0), 'mul')))
            s['effg'] = jnp.where(s["_pc"] == 0, s['effg'], old_1862)
            return s
        def no_1861(s):
            s = dict(s)
            old_1863 = s['effg']
            s['effg'] = s['effg'].at[s['k'] - 1].set(F(F(25.0)))
            s['effg'] = jnp.where(s["_pc"] == 0, s['effg'], old_1863)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['qg3d'][s['k'] - 1] >= s['tmpqsmall'])), yes_1861, no_1861, s)
        def yes_1864(s):
            s = dict(s)
            old_1865 = s['nim_morr_cl']
            s['nim_morr_cl'] = s['nim_morr_cl'].at[s['k'] - 1].set(F(_div(_arith(s['ni3d'][s['k'] - 1], _div(F(10000000.0), s['rho'][s['k'] - 1]), 'sub'), s['dt'])))
            s['nim_morr_cl'] = jnp.where(s["_pc"] == 0, s['nim_morr_cl'], old_1865)
            return s
        def no_1864(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['ni3d'][s['k'] - 1] > _div(F(10000000.0), s['rho'][s['k'] - 1]))), yes_1864, no_1864, s)
        old_1866 = s['ni3d']
        s['ni3d'] = s['ni3d'].at[s['k'] - 1].set(F(jnp.minimum(s['ni3d'][s['k'] - 1], _div(F(10000000.0), s['rho'][s['k'] - 1]))))
        s['ni3d'] = jnp.where(s["_pc"] == 0, s['ni3d'], old_1866)
        def yes_1867(s):
            s = dict(s)
            old_1868 = s['nc3d']
            s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(jnp.minimum(s['nc3d'][s['k'] - 1], _div(_arith(s['nanew1'], s['nanew2'], 'add'), s['rho'][s['k'] - 1]))))
            s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1868)
            return s
        def no_1867(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((((s['inum'] == 0)) & ((s['iact'] == 2)))), yes_1867, no_1867, s)
        def yes_1869(s):
            s = dict(s)
            old_1870 = s['nc3d']
            s['nc3d'] = s['nc3d'].at[s['k'] - 1].set(F(_div(_arith(s['ndcnst'], F(1000000.0), 'mul'), s['rho'][s['k'] - 1])))
            s['nc3d'] = jnp.where(s["_pc"] == 0, s['nc3d'], old_1870)
            return s
        def no_1869(s):
            s = dict(s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['inum'] == 1)), yes_1869, no_1869, s)
        adj = _positive_qv_adj(s['qv3d'][s['k'] - 1], s['qc3d'][s['k'] - 1], s['qr3d'][s['k'] - 1], s['qi3d'][s['k'] - 1], s['qni3d'][s['k'] - 1], s['qg3d'][s['k'] - 1], s['t3d'][s['k'] - 1], s["iliq"], s["igraup"])
        s['qv3d'] = s['qv3d'].at[s['k'] - 1].set(jnp.where(s["_pc"] == 0, adj[0], s['qv3d'][s['k'] - 1]))
        s['qc3d'] = s['qc3d'].at[s['k'] - 1].set(jnp.where(s["_pc"] == 0, adj[1], s['qc3d'][s['k'] - 1]))
        s['qr3d'] = s['qr3d'].at[s['k'] - 1].set(jnp.where(s["_pc"] == 0, adj[2], s['qr3d'][s['k'] - 1]))
        s['qi3d'] = s['qi3d'].at[s['k'] - 1].set(jnp.where(s["_pc"] == 0, adj[3], s['qi3d'][s['k'] - 1]))
        s['qni3d'] = s['qni3d'].at[s['k'] - 1].set(jnp.where(s["_pc"] == 0, adj[4], s['qni3d'][s['k'] - 1]))
        s['qg3d'] = s['qg3d'].at[s['k'] - 1].set(jnp.where(s["_pc"] == 0, adj[5], s['qg3d'][s['k'] - 1]))
        s['t3d'] = s['t3d'].at[s['k'] - 1].set(jnp.where(s["_pc"] == 0, adj[6], s['t3d'][s['k'] - 1]))
        return s
    s = lax.fori_loop(0, I((s['kte'] - s['kts']) // (1) + 1), loop_1612, s)
    s["_pc"] = jnp.where(s["_pc"] == 400, I(0), s["_pc"])
    return {'qc3dten':s['qc3dten'], 'qi3dten':s['qi3dten'], 'qni3dten':s['qni3dten'], 'qr3dten':s['qr3dten'], 'nc3dten':s['nc3dten'], 'ni3dten':s['ni3dten'], 'ns3dten':s['ns3dten'], 'nr3dten':s['nr3dten'], 'qc3d':s['qc3d'], 'qi3d':s['qi3d'], 'qni3d':s['qni3d'], 'qr3d':s['qr3d'], 'nc3d':s['nc3d'], 'ni3d':s['ni3d'], 'ns3d':s['ns3d'], 'nr3d':s['nr3d'], 't3dten':s['t3dten'], 'qv3dten':s['qv3dten'], 't3d':s['t3d'], 'qv3d':s['qv3d'], 'pres':s['pres'], 'rho':s['rho'], 'dzq':s['dzq'], 'w3d':s['w3d'], 'wvar':s['wvar'], 'fr':s['fr'], 'effc':s['effc'], 'effi':s['effi'], 'effs':s['effs'], 'effr':s['effr'], 'qg3dten':s['qg3dten'], 'ng3dten':s['ng3dten'], 'qg3d':s['qg3d'], 'ng3d':s['ng3d'], 'effg':s['effg'], 'qgsten':s['qgsten'], 'qrsten':s['qrsten'], 'qisten':s['qisten'], 'qnisten':s['qnisten'], 'qcsten':s['qcsten'], 'ngsten':s['ngsten'], 'nrsten':s['nrsten'], 'nisten':s['nisten'], 'nssten':s['nssten'], 'ncsten':s['ncsten'], 'cf3d':s['cf3d'], 'prc':s['prc'], 'pra':s['pra'], 'psmlt':s['psmlt'], 'evpms':s['evpms'], 'pracs':s['pracs'], 'evpmg':s['evpmg'], 'pracg':s['pracg'], 'pre':s['pre'], 'pgmlt':s['pgmlt'], 'mnuccc':s['mnuccc'], 'psacws':s['psacws'], 'psacwi':s['psacwi'], 'qmults':s['qmults'], 'qmultg':s['qmultg'], 'psacwg':s['psacwg'], 'pgsacw':s['pgsacw'], 'prd':s['prd'], 'prci':s['prci'], 'prai':s['prai'], 'qmultr':s['qmultr'], 'qmultrg':s['qmultrg'], 'mnuccd':s['mnuccd'], 'praci':s['praci'], 'pracis':s['pracis'], 'eprd':s['eprd'], 'mnuccr':s['mnuccr'], 'piacr':s['piacr'], 'piacrs':s['piacrs'], 'pgracs':s['pgracs'], 'prds':s['prds'], 'eprds':s['eprds'], 'psacr':s['psacr'], 'prdg':s['prdg'], 'eprdg':s['eprdg'], 'nprc1':s['nprc1'], 'nragg':s['nragg'], 'npracg':s['npracg'], 'nsubr':s['nsubr'], 'nsmltr':s['nsmltr'], 'ngmltr':s['ngmltr'], 'npracs':s['npracs'], 'nnuccr':s['nnuccr'], 'niacr':s['niacr'], 'niacrs':s['niacrs'], 'ngracs':s['ngracs'], 'nsmlts':s['nsmlts'], 'nsagg':s['nsagg'], 'nprci':s['nprci'], 'nscng':s['nscng'], 'nsubs':s['nsubs'], 'pcc':s['pcc'], 'nnuccc':s['nnuccc'], 'npsacws':s['npsacws'], 'npra':s['npra'], 'nprc':s['nprc'], 'npsacwi':s['npsacwi'], 'npsacwg':s['npsacwg'], 'nprai':s['nprai'], 'nmults':s['nmults'], 'nmultg':s['nmultg'], 'nmultr':s['nmultr'], 'nmultrg':s['nmultrg'], 'nnuccd':s['nnuccd'], 'nsubi':s['nsubi'], 'ngmltg':s['ngmltg'], 'nsubg':s['nsubg'], 'nact':s['nact'], 'sizefix_nr':s['sizefix_nr'], 'sizefix_nc':s['sizefix_nc'], 'sizefix_ni':s['sizefix_ni'], 'sizefix_ns':s['sizefix_ns'], 'sizefix_ng':s['sizefix_ng'], 'negfix_nr':s['negfix_nr'], 'negfix_nc':s['negfix_nc'], 'negfix_ni':s['negfix_ni'], 'negfix_ns':s['negfix_ns'], 'negfix_ng':s['negfix_ng'], 'nim_morr_cl':s['nim_morr_cl'], 'qc_inst':s['qc_inst'], 'qr_inst':s['qr_inst'], 'qi_inst':s['qi_inst'], 'qs_inst':s['qs_inst'], 'qg_inst':s['qg_inst'], 'nc_inst':s['nc_inst'], 'nr_inst':s['nr_inst'], 'ni_inst':s['ni_inst'], 'ns_inst':s['ns_inst'], 'ng_inst':s['ng_inst'], 'precrt':s['precrt'], 'snowrt':s['snowrt']}


def POLYSVP(t, type):
    s = {}
    s['dum'] = jnp.zeros((), dtype=jnp.float32)
    s['a0i'] = jnp.zeros((), dtype=jnp.float32)
    s['a1i'] = jnp.zeros((), dtype=jnp.float32)
    s['a2i'] = jnp.zeros((), dtype=jnp.float32)
    s['a3i'] = jnp.zeros((), dtype=jnp.float32)
    s['a4i'] = jnp.zeros((), dtype=jnp.float32)
    s['a5i'] = jnp.zeros((), dtype=jnp.float32)
    s['a6i'] = jnp.zeros((), dtype=jnp.float32)
    s['a7i'] = jnp.zeros((), dtype=jnp.float32)
    s['a8i'] = jnp.zeros((), dtype=jnp.float32)
    s['a0'] = jnp.zeros((), dtype=jnp.float32)
    s['a1'] = jnp.zeros((), dtype=jnp.float32)
    s['a2'] = jnp.zeros((), dtype=jnp.float32)
    s['a3'] = jnp.zeros((), dtype=jnp.float32)
    s['a4'] = jnp.zeros((), dtype=jnp.float32)
    s['a5'] = jnp.zeros((), dtype=jnp.float32)
    s['a6'] = jnp.zeros((), dtype=jnp.float32)
    s['a7'] = jnp.zeros((), dtype=jnp.float32)
    s['a8'] = jnp.zeros((), dtype=jnp.float32)
    s['dt'] = jnp.zeros((), dtype=jnp.float32)
    s['_pc'] = jnp.zeros((), dtype=jnp.int32)
    s['polysvp'] = jnp.zeros((), dtype=jnp.float32)
    s['a0i'] = F(F(6.11147274))
    s['a1i'] = F(F(0.50316082))
    s['a2i'] = F(F(0.0188439774))
    s['a3i'] = F(F(0.000420895665))
    s['a4i'] = F(F(6.15021634e-06))
    s['a5i'] = F(F(6.02588177e-08))
    s['a6i'] = F(F(3.85852041e-10))
    s['a7i'] = F(F(1.46898966e-12))
    s['a8i'] = F(F(2.52751365e-15))
    s['a0'] = F(F(6.11239921))
    s['a1'] = F(F(0.443987641))
    s['a2'] = F(F(0.0142986287))
    s['a3'] = F(F(0.00026484743))
    s['a4'] = F(F(3.02950461e-06))
    s['a5'] = F(F(2.06739458e-08))
    s['a6'] = F(F(6.40689451e-11))
    s['a7'] = F((-F(9.52447341e-14)))
    s['a8'] = F((-F(9.76195544e-16)))
    s['t'] = F(t)
    s['type'] = I(type)
    def yes_1(s):
        s = dict(s)
        old_2 = s['dt']
        s['dt'] = F(jnp.maximum((-F(80.0)), _arith(s['t'], F(273.16), 'sub')))
        s['dt'] = jnp.where(s["_pc"] == 0, s['dt'], old_2)
        old_3 = s['polysvp']
        s['polysvp'] = F(_arith(s['a0i'], _arith(s['dt'], _arith(s['a1i'], _arith(s['dt'], _arith(s['a2i'], _arith(s['dt'], _arith(s['a3i'], _arith(s['dt'], _arith(s['a4i'], _arith(s['dt'], _arith(s['a5i'], _arith(s['dt'], _arith(s['a6i'], _arith(s['dt'], _arith(s['a7i'], _arith(s['a8i'], s['dt'], 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'))
        s['polysvp'] = jnp.where(s["_pc"] == 0, s['polysvp'], old_3)
        old_4 = s['polysvp']
        s['polysvp'] = F(_arith(s['polysvp'], F(100.0), 'mul'))
        s['polysvp'] = jnp.where(s["_pc"] == 0, s['polysvp'], old_4)
        return s
    def no_1(s):
        s = dict(s)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['type'] == 1)), yes_1, no_1, s)
    def yes_5(s):
        s = dict(s)
        old_6 = s['dt']
        s['dt'] = F(jnp.maximum((-F(80.0)), _arith(s['t'], F(273.16), 'sub')))
        s['dt'] = jnp.where(s["_pc"] == 0, s['dt'], old_6)
        old_7 = s['polysvp']
        s['polysvp'] = F(_arith(s['a0'], _arith(s['dt'], _arith(s['a1'], _arith(s['dt'], _arith(s['a2'], _arith(s['dt'], _arith(s['a3'], _arith(s['dt'], _arith(s['a4'], _arith(s['dt'], _arith(s['a5'], _arith(s['dt'], _arith(s['a6'], _arith(s['dt'], _arith(s['a7'], _arith(s['a8'], s['dt'], 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'), 'mul'), 'add'))
        s['polysvp'] = jnp.where(s["_pc"] == 0, s['polysvp'], old_7)
        old_8 = s['polysvp']
        s['polysvp'] = F(_arith(s['polysvp'], F(100.0), 'mul'))
        s['polysvp'] = jnp.where(s["_pc"] == 0, s['polysvp'], old_8)
        return s
    def no_5(s):
        s = dict(s)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['type'] == 0)), yes_5, no_5, s)
    return s['polysvp']


def GAMMA(x):
    s = {}
    s['pi'] = jnp.zeros((), dtype=jnp.float32)
    s['sqrtpi'] = jnp.zeros((), dtype=jnp.float32)
    s['i'] = jnp.zeros((), dtype=jnp.int32)
    s['n'] = jnp.zeros((), dtype=jnp.int32)
    s['parity'] = jnp.zeros((), dtype=jnp.bool_)
    s['conv'] = jnp.zeros((), dtype=jnp.float32)
    s['eps'] = jnp.zeros((), dtype=jnp.float32)
    s['fact'] = jnp.zeros((), dtype=jnp.float32)
    s['half'] = jnp.zeros((), dtype=jnp.float32)
    s['one'] = jnp.zeros((), dtype=jnp.float32)
    s['res'] = jnp.zeros((), dtype=jnp.float32)
    s['sum'] = jnp.zeros((), dtype=jnp.float32)
    s['twelve'] = jnp.zeros((), dtype=jnp.float32)
    s['two'] = jnp.zeros((), dtype=jnp.float32)
    s['xbig'] = jnp.zeros((), dtype=jnp.float32)
    s['xden'] = jnp.zeros((), dtype=jnp.float32)
    s['xinf'] = jnp.zeros((), dtype=jnp.float32)
    s['xminin'] = jnp.zeros((), dtype=jnp.float32)
    s['xnum'] = jnp.zeros((), dtype=jnp.float32)
    s['y'] = jnp.zeros((), dtype=jnp.float32)
    s['y1'] = jnp.zeros((), dtype=jnp.float32)
    s['ysq'] = jnp.zeros((), dtype=jnp.float32)
    s['z'] = jnp.zeros((), dtype=jnp.float32)
    s['zero'] = jnp.zeros((), dtype=jnp.float32)
    s['c'] = jnp.zeros((7,), dtype=jnp.float32)
    s['p'] = jnp.zeros((8,), dtype=jnp.float32)
    s['q'] = jnp.zeros((8,), dtype=jnp.float32)
    s['_pc'] = jnp.zeros((), dtype=jnp.int32)
    s['gamma'] = jnp.zeros((), dtype=jnp.float32)
    s['one'] = F(F(1.0))
    s['half'] = F(F(0.5))
    s['twelve'] = F(F(12.0))
    s['two'] = F(F(2.0))
    s['zero'] = F(F(0.0))
    s['xbig'] = F(F(35.04))
    s['xminin'] = F(F(1.18e-38))
    s['eps'] = F(F(1.19e-07))
    s['xinf'] = F(F(3.4e+38))
    s['p'] = s['p'].at[1 - 1].set(F((-F(1.716185138865495))))
    s['p'] = s['p'].at[2 - 1].set(F(F(24.76565080557592)))
    s['p'] = s['p'].at[3 - 1].set(F((-F(379.80425647094563))))
    s['p'] = s['p'].at[4 - 1].set(F(F(629.3311553128184)))
    s['p'] = s['p'].at[5 - 1].set(F(F(866.9662027904133)))
    s['p'] = s['p'].at[6 - 1].set(F((-F(31451.272968848367))))
    s['p'] = s['p'].at[7 - 1].set(F((-F(36144.413418691176))))
    s['p'] = s['p'].at[8 - 1].set(F(F(66456.14382024054)))
    s['q'] = s['q'].at[1 - 1].set(F((-F(30.840230011973897))))
    s['q'] = s['q'].at[2 - 1].set(F(F(315.35062697960416)))
    s['q'] = s['q'].at[3 - 1].set(F((-F(1015.1563674902192))))
    s['q'] = s['q'].at[4 - 1].set(F((-F(3107.771671572311))))
    s['q'] = s['q'].at[5 - 1].set(F(F(22538.11842098015)))
    s['q'] = s['q'].at[6 - 1].set(F(F(4755.846277527881)))
    s['q'] = s['q'].at[7 - 1].set(F((-F(134659.9598649693))))
    s['q'] = s['q'].at[8 - 1].set(F((-F(115132.25967555349))))
    s['c'] = s['c'].at[1 - 1].set(F((-F(0.001910444077728))))
    s['c'] = s['c'].at[2 - 1].set(F(F(0.00084171387781295)))
    s['c'] = s['c'].at[3 - 1].set(F((-F(0.0005952379913043012))))
    s['c'] = s['c'].at[4 - 1].set(F(F(0.0007936507935003503)))
    s['c'] = s['c'].at[5 - 1].set(F((-F(0.0027777777777776816))))
    s['c'] = s['c'].at[6 - 1].set(F(F(0.08333333333333333)))
    s['c'] = s['c'].at[7 - 1].set(F(F(0.0057083835261)))
    s['x'] = F(x)
    s["pi"] = F(3.1415926535897932384626434)
    s["sqrtpi"] = F(0.9189385332046727417803297)
    old_1 = s['parity']
    s['parity'] = B(False)
    s['parity'] = jnp.where(s["_pc"] == 0, s['parity'], old_1)
    old_2 = s['fact']
    s['fact'] = F(s['one'])
    s['fact'] = jnp.where(s["_pc"] == 0, s['fact'], old_2)
    old_3 = s['n']
    s['n'] = I(0)
    s['n'] = jnp.where(s["_pc"] == 0, s['n'], old_3)
    old_4 = s['y']
    s['y'] = F(s['x'])
    s['y'] = jnp.where(s["_pc"] == 0, s['y'], old_4)
    def yes_5(s):
        s = dict(s)
        old_6 = s['y']
        s['y'] = F((-s['x']))
        s['y'] = jnp.where(s["_pc"] == 0, s['y'], old_6)
        old_7 = s['y1']
        s['y1'] = F(jnp.trunc(s['y']))
        s['y1'] = jnp.where(s["_pc"] == 0, s['y1'], old_7)
        old_8 = s['res']
        s['res'] = F(_arith(s['y'], s['y1'], 'sub'))
        s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_8)
        def yes_9(s):
            s = dict(s)
            def yes_10(s):
                s = dict(s)
                old_11 = s['parity']
                s['parity'] = B(True)
                s['parity'] = jnp.where(s["_pc"] == 0, s['parity'], old_11)
                return s
            def no_10(s):
                s = dict(s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['y1'] != _arith(jnp.trunc(_arith(s['y1'], s['half'], 'mul')), s['two'], 'mul'))), yes_10, no_10, s)
            old_12 = s['fact']
            s['fact'] = F(_div((-s['pi']), _intrinsic('sin', _arith(s['pi'], s['res'], 'mul'))))
            s['fact'] = jnp.where(s["_pc"] == 0, s['fact'], old_12)
            old_13 = s['y']
            s['y'] = F(_arith(s['y'], s['one'], 'add'))
            s['y'] = jnp.where(s["_pc"] == 0, s['y'], old_13)
            return s
        def no_9(s):
            s = dict(s)
            old_14 = s['res']
            s['res'] = F(s['xinf'])
            s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_14)
            s["_pc"] = jnp.where(s["_pc"] == 0, I(900), s["_pc"])
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['res'] != s['zero'])), yes_9, no_9, s)
        return s
    def no_5(s):
        s = dict(s)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['y'] <= s['zero'])), yes_5, no_5, s)
    def yes_15(s):
        s = dict(s)
        def yes_16(s):
            s = dict(s)
            old_17 = s['res']
            s['res'] = F(_div(s['one'], s['y']))
            s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_17)
            return s
        def no_16(s):
            s = dict(s)
            old_18 = s['res']
            s['res'] = F(s['xinf'])
            s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_18)
            s["_pc"] = jnp.where(s["_pc"] == 0, I(900), s["_pc"])
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['y'] >= s['xminin'])), yes_16, no_16, s)
        return s
    def no_15(s):
        s = dict(s)
        def yes_19(s):
            s = dict(s)
            old_20 = s['y1']
            s['y1'] = F(s['y'])
            s['y1'] = jnp.where(s["_pc"] == 0, s['y1'], old_20)
            def yes_21(s):
                s = dict(s)
                old_22 = s['z']
                s['z'] = F(s['y'])
                s['z'] = jnp.where(s["_pc"] == 0, s['z'], old_22)
                old_23 = s['y']
                s['y'] = F(_arith(s['y'], s['one'], 'add'))
                s['y'] = jnp.where(s["_pc"] == 0, s['y'], old_23)
                return s
            def no_21(s):
                s = dict(s)
                old_24 = s['n']
                s['n'] = I(_arith(I(s['y']), 1, 'sub'))
                s['n'] = jnp.where(s["_pc"] == 0, s['n'], old_24)
                old_25 = s['y']
                s['y'] = F(_arith(s['y'], F(s['n']), 'sub'))
                s['y'] = jnp.where(s["_pc"] == 0, s['y'], old_25)
                old_26 = s['z']
                s['z'] = F(_arith(s['y'], s['one'], 'sub'))
                s['z'] = jnp.where(s["_pc"] == 0, s['z'], old_26)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['y'] < s['one'])), yes_21, no_21, s)
            old_27 = s['xnum']
            s['xnum'] = F(s['zero'])
            s['xnum'] = jnp.where(s["_pc"] == 0, s['xnum'], old_27)
            old_28 = s['xden']
            s['xden'] = F(s['one'])
            s['xden'] = jnp.where(s["_pc"] == 0, s['xden'], old_28)
            def loop_29(iteration, s):
                s = dict(s)
                s['i'] = I(1 + iteration * (1))
                old_30 = s['xnum']
                s['xnum'] = F(_arith(_arith(s['xnum'], s['p'][s['i'] - 1], 'add'), s['z'], 'mul'))
                s['xnum'] = jnp.where(s["_pc"] == 0, s['xnum'], old_30)
                old_31 = s['xden']
                s['xden'] = F(_arith(_arith(s['xden'], s['z'], 'mul'), s['q'][s['i'] - 1], 'add'))
                s['xden'] = jnp.where(s["_pc"] == 0, s['xden'], old_31)
                return s
            s = lax.fori_loop(0, I((8 - 1) // (1) + 1), loop_29, s)
            old_32 = s['res']
            s['res'] = F(_arith(_div(s['xnum'], s['xden']), s['one'], 'add'))
            s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_32)
            def yes_33(s):
                s = dict(s)
                old_34 = s['res']
                s['res'] = F(_div(s['res'], s['y1']))
                s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_34)
                return s
            def no_33(s):
                s = dict(s)
                def yes_35(s):
                    s = dict(s)
                    def loop_36(iteration, s):
                        s = dict(s)
                        s['i'] = I(1 + iteration * (1))
                        old_37 = s['res']
                        s['res'] = F(_arith(s['res'], s['y'], 'mul'))
                        s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_37)
                        old_38 = s['y']
                        s['y'] = F(_arith(s['y'], s['one'], 'add'))
                        s['y'] = jnp.where(s["_pc"] == 0, s['y'], old_38)
                        return s
                    s = lax.fori_loop(0, I((s['n'] - 1) // (1) + 1), loop_36, s)
                    return s
                def no_35(s):
                    s = dict(s)
                    return s
                s = lax.cond((s["_pc"] == 0) & ((s['y1'] > s['y'])), yes_35, no_35, s)
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['y1'] < s['y'])), yes_33, no_33, s)
            return s
        def no_19(s):
            s = dict(s)
            def yes_39(s):
                s = dict(s)
                old_40 = s['ysq']
                s['ysq'] = F(_arith(s['y'], s['y'], 'mul'))
                s['ysq'] = jnp.where(s["_pc"] == 0, s['ysq'], old_40)
                old_41 = s['sum']
                s['sum'] = F(s['c'][7 - 1])
                s['sum'] = jnp.where(s["_pc"] == 0, s['sum'], old_41)
                def loop_42(iteration, s):
                    s = dict(s)
                    s['i'] = I(1 + iteration * (1))
                    old_43 = s['sum']
                    s['sum'] = F(_arith(_div(s['sum'], s['ysq']), s['c'][s['i'] - 1], 'add'))
                    s['sum'] = jnp.where(s["_pc"] == 0, s['sum'], old_43)
                    return s
                s = lax.fori_loop(0, I((6 - 1) // (1) + 1), loop_42, s)
                old_44 = s['sum']
                s['sum'] = F(_arith(_arith(_div(s['sum'], s['y']), s['y'], 'sub'), s['sqrtpi'], 'add'))
                s['sum'] = jnp.where(s["_pc"] == 0, s['sum'], old_44)
                old_45 = s['sum']
                s['sum'] = F(_arith(s['sum'], _arith(_arith(s['y'], s['half'], 'sub'), _intrinsic('log', s['y']), 'mul'), 'add'))
                s['sum'] = jnp.where(s["_pc"] == 0, s['sum'], old_45)
                old_46 = s['res']
                s['res'] = F(_intrinsic('exp', s['sum']))
                s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_46)
                return s
            def no_39(s):
                s = dict(s)
                old_47 = s['res']
                s['res'] = F(s['xinf'])
                s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_47)
                s["_pc"] = jnp.where(s["_pc"] == 0, I(900), s["_pc"])
                return s
            s = lax.cond((s["_pc"] == 0) & ((s['y'] <= s['xbig'])), yes_39, no_39, s)
            return s
        s = lax.cond((s["_pc"] == 0) & ((s['y'] < s['twelve'])), yes_19, no_19, s)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['y'] < s['eps'])), yes_15, no_15, s)
    def yes_48(s):
        s = dict(s)
        old_49 = s['res']
        s['res'] = F((-s['res']))
        s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_49)
        return s
    def no_48(s):
        s = dict(s)
        return s
    s = lax.cond((s["_pc"] == 0) & (s['parity']), yes_48, no_48, s)
    def yes_50(s):
        s = dict(s)
        old_51 = s['res']
        s['res'] = F(_div(s['fact'], s['res']))
        s['res'] = jnp.where(s["_pc"] == 0, s['res'], old_51)
        return s
    def no_50(s):
        s = dict(s)
        return s
    s = lax.cond((s["_pc"] == 0) & ((s['fact'] != s['one'])), yes_50, no_50, s)
    s["_pc"] = jnp.where(s["_pc"] == 900, I(0), s["_pc"])
    old_52 = s['gamma']
    s['gamma'] = F(s['res'])
    s['gamma'] = jnp.where(s["_pc"] == 0, s['gamma'], old_52)
    return s['gamma']


def DERF1(x):
    """Error function erf(x), faithful to module_mp_graupel.F90:DERF1 (Ooura approximation).
    Array-capable; the coefficient block is gathered per element by the argument range.
    """
    x = jnp.asarray(x, dtype=jnp.float64)
    w = jnp.abs(x)

    # Branch 1: w < 2.2 — k=int(w^2) in 0..4, polynomial in t=frac(w^2), ×w
    t1 = w * w
    k1 = jnp.clip(jnp.floor(t1).astype(jnp.int32), 0, 4)
    t1f = t1 - k1
    y1 = _horner13(_DERF1_A_BLOCKS[k1], t1f) * w

    # Branch 2: 2.2 <= w < 6.9 — k=int(w) in 2..6, block k-2, polynomial in t=w-k, then 1-poly^16
    k2 = jnp.clip(jnp.floor(w).astype(jnp.int32), 2, 6)
    t2 = w - k2
    y2 = _horner13(_DERF1_B_BLOCKS[k2 - 2], t2)
    y2 = y2 * y2
    y2 = y2 * y2
    y2 = y2 * y2
    y2 = 1.0 - y2 * y2

    y = jnp.where(w < 2.2, y1, jnp.where(w < 6.9, y2, 1.0))
    return jnp.where(x < 0.0, -y, y)


def _positive_qv_adj(qv,qc,qr,qi,qs,qg,t,iliq,igraup):
    initial=(qv,qc,qr,qi,qs,qg,t)
    def adjust(_):
        liq=F(qc+qr)
        ice=F(jnp.where(iliq==0,F(qs+qi),F(0)))
        ice=F(jnp.where(igraup==0,F(ice+qg),ice))
        total=F(F(qv+liq)+ice)
        def enough(_):
            delta=F(F(1.e-12)-qv)
            factor=F(F(1)-F(delta/F(liq+ice)))
            c,r=F(qc*factor),F(qr*factor)
            i,s=jnp.where(iliq==0,F(qi*factor),qi),jnp.where(iliq==0,F(qs*factor),qs)
            g=jnp.where(igraup==0,F(qg*factor),qg)
            dl=F(liq-F(c+r))
            di=jnp.where(iliq==0,F(ice-F(i+s)),F(0))
            di=jnp.where(igraup==0,F(ice-F(F(i+s)+g)),di)
            temp=F(D(t)-D(C.Lv/C.Cp)*D(dl)-D(C.Ls/C.Cp)*D(di))
            return F(1.e-12),c,r,i,s,g,temp
        return lax.cond(total>=F(2)*F(1.e-12),enough,lambda _:initial,None)
    return lax.cond(qv<F(1.e-12),adjust,lambda _:initial,None)


# Internal arithmetic and polynomial helpers.
def _arith(a,b,op):
    if op=='add':return a+b
    if op=='sub':return a-b
    if op=='mul':return a*b
    return a**b

def _intrinsic(name,x):
    return getattr(jnp,name)(x)

def _div(a,b):
    return a/b

def _fatal(_):
    raise FloatingPointError("Morrison core conservation assertion failed")


def _horner13(coeffs, t):
    """Horner eval of 13 coeffs (coeffs[..., 13]) at t: c0·t^12 + ... + c12 (Fortran order)."""
    y = coeffs[..., 0]
    for i in range(1, 13):
        y = y * t + coeffs[..., i]
    return y
