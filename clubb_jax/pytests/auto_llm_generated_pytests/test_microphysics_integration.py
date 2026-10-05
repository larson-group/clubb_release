"""Focused interface contracts; Morrison checks do not establish case accuracy."""
from types import SimpleNamespace
import numpy as np
import pytest
import jax
import jax.numpy as jnp

from clubb_jax.src.Microphys import parameters_microphys as parameters
from clubb_jax.src.Microphys.microphys_init_cleanup import init_microphys, cleanup_microphys
from clubb_jax.src.Microphys.morrison_microphys_module import morrison_microphys_driver
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.pytests.microphysics_test_inputs import initialize




@pytest.mark.parametrize('scheme,ice,graupel,dimension', [
    ('none', False, False, 0), ('khairoutdinov_kogan', False, False, 2),
    ('morrison', False, False, 2), ('morrison', True, False, 6),
    ('morrison', True, True, 8)])
def test_species_lifecycle(scheme, ice, graupel, dimension):
    size, pdf_size, metadata, *rest = initialize(microphys_scheme=scheme,
        l_ice_microphys=ice, l_graupel=graupel)
    assert size == dimension
    assert len(parameters.l_hydromet_sed) == dimension
    if dimension:
        assert (metadata.iirr, metadata.iiNr) == (0, 1)
        assert metadata.iiri == (2 if ice else -1)
        assert metadata.iirg == (6 if graupel else -1)
    cleanup_microphys()
    assert parameters.l_hydromet_sed == ()


@pytest.mark.parametrize(
    "config,reason",
    [
        ({"lh_microphys_type": "interactive", "microphys_scheme": "none"}, "SILHS"),
        ({"microphys_scheme": "coamps"}, "Unsupported"),
        ({"microphys_scheme": "simplified_ice"}, "Unsupported"),
        ({"l_gfdl_activation": True}, "GFDL"),
        ({"microphys_scheme": "morrison", "l_cloud_sed": True}, "sedimentation"),
        ({"microphys_scheme": "khairoutdinov_kogan", "l_predict_Nc": True}, "l_predict_Nc"),
        ({"l_morr_xp2_mc": True}, "l_morr_xp2_mc"),
    ],
)
def test_disabled_dependencies_rejected(config, reason):
    with pytest.raises(ValueError, match=reason):
        initialize(**config)


@pytest.mark.parametrize('dimension,ncol', [(2,1),(6,2),(8,2)])
def test_morrison_species_interface_eager_and_jit(dimension, ncol):
    _, _, metadata, *_ = initialize(microphys_scheme='morrison',
        l_ice_microphys=dimension>2, l_graupel=dimension>6)
    nzt=5
    gr = SimpleNamespace(zt=jnp.broadcast_to(jnp.arange(nzt)*100.,(ncol,nzt)),
                         dzt=jnp.full((ncol,nzt),100.))
    one=jnp.ones((ncol,nzt)); zero=jnp.zeros_like(one)
    hydromet=jnp.zeros((ncol,nzt,dimension))
    for i in range(0,dimension,2):
        hydromet=hydromet.at[...,i].set(1.e-5)
        hydromet=hydromet.at[...,i+1].set(2.e4)
    # Columns differ to catch accidental flattening or cross-column transport.
    hydromet=hydromet*jnp.arange(1,ncol+1)[:,None,None]
    stats=JaxStats.empty(l_sample=True, names=("rrm_auto", "precip_rate_sfc"),
        grids=("zt", "sfc"), ncol=ncol, max_nlev=nzt)
    def run(hm):
        return morrison_microphys_driver(
            gr, ncol, 10., nzt,                                                # In
            dimension, metadata,                                               # In
            False, jnp.linspace(258.,285.,nzt)[None,:]*one, zero, 90000.*one,  # In
            one, one, 0.5*one, 0.2*one,                                        # In
            100.*one, 1.e-4*one, 1.e8*one, zero, 0.007*one, hm,                # In
            1,                                                                 # In
            one,                                                               # In
            stats,                                                             # InOut
        )
    eager=run(hydromet)
    compiled=jax.jit(run)(hydromet)
    assert len(eager)==12
    for a,b in zip(eager[1:],compiled[1:]):
        assert np.all(np.isfinite(a))
        np.testing.assert_allclose(a,b,rtol=1.e-6,atol=1.e-9)
    assert eager[1].shape==hydromet.shape
    assert np.all(np.asarray(hydromet+10.*eager[1])>=-1.e-12)
    assert np.all(np.asarray(eager[2][...,1:])==0.)
    assert np.all(np.asarray(eager[2][...,0])<=0.)
    for a,b in zip(eager[0].buffers,compiled[0].buffers):
        np.testing.assert_allclose(a,b,rtol=1.e-6,atol=1.e-9)
    np.testing.assert_allclose(eager[0].buffers[0][0],eager[7])
    assert np.all(np.asarray(eager[0].nsamples[0]) == 1)
    assert np.all(np.asarray(eager[0].nsamples[2]) == 1)


def test_previous_tendencies_feed_next_core_step(monkeypatch):
    from clubb_jax.src import advance_clubb_to_end as driver
    zero=jnp.zeros((2,3)); one=jnp.ones_like(zero)
    state=dict(dt_main=60.,dt_rad=60.,time_initial=0.,ifinal=2,l_stats=False,
        ngrdcol=2,nzt=3,nzm=4,thlm=290.*one,rtm=.01*one,rcm=zero,
        exner=one,thv_ds_zt=290.*one,l_calc_thlp2_rad=False,radht=zero,
        err_info=SimpleNamespace(is_fatal=lambda:False))
    moments=('wprtp','wpthlp','rtp2','thlp2','rtpthlp')
    for name in ('rcm','rvm','thlm')+moments:
        state[name+'_mc']=zero
    observed=[];order=[]
    def forcing(state,time):
        for name in ('rtm','thlm')+moments:
            state[name+'_forcing']=zero
    def core(state):
        order.append('core')
        observed.append([np.asarray(state[name+'_forcing']) for name in ('rtm','thlm')+moments])
    def micro(state,*args):
        order.append('micro')
        for i,name in enumerate(('rcm','rvm','thlm')+moments,1):
            state[name+'_mc']=i*one
    def radiation(**kwargs):
        order.append('radiation')
    monkeypatch.setattr(driver,'_prescribe_forcings',forcing)
    monkeypatch.setattr(driver,'_advance_clubb_core',core)
    monkeypatch.setattr(driver,'_advance_microphysics',micro)
    monkeypatch.setattr(driver,'_advance_radiation',radiation)
    driver.advance_clubb_to_end(state,l_stdout=False)
    assert order==['core','micro','radiation']*2
    for field in observed[0]: np.testing.assert_array_equal(field,zero)
    for field,value in zip(observed[1],(3.,3.,4.,5.,6.,7.,8.)):
        np.testing.assert_array_equal(field,value*one)


@pytest.mark.parametrize('ncol',[1,3])
def test_cloud_sedimentation_closed_boundaries_conserve_column(ncol):
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
    from clubb_jax.src.Microphys.cloud_sed_module import cloud_drop_sed
    from clubb_jax.src.CLUBB_core.constants_clubb import Lv,Cp
    gr=setup_grid(ncol,100.,0.,600.)
    one=jnp.ones((ncol,gr.nzt));zero=jnp.zeros_like(one)
    rcm=one*jnp.linspace(0.,2.e-4,gr.nzt)[None,:]
    stats=JaxStats.empty(l_sample=True,names=('Fcsed','sed_rcm'),grids=('zm','zt'),
        ncol=ncol,max_nlev=gr.nzm,grid_nlev=(gr.nzt,gr.nzm,1,gr.nzt,1,gr.nzt,gr.nzt))
    def run(rc):
        return cloud_drop_sed(gr, gr.ngrdcol, rc, 1.e8*one, jnp.ones((ncol,gr.nzm)),
                   one, one, 1.5, stats, zero,
                   zero)
    eager=run(rcm);compiled=jax.jit(run)(rcm)
    for a,b in zip(eager[1:],compiled[1:]): np.testing.assert_allclose(a,b,rtol=1.e-13,atol=1.e-20)
    np.testing.assert_allclose(jnp.sum(eager[1]*gr.dzt,axis=-1),0.,atol=1.e-18)
    np.testing.assert_allclose(eager[2],-Lv/Cp*eager[1],rtol=1.e-14,atol=1.e-20)
    assert np.any(np.asarray(eager[1])!=0.)
    assert np.all(np.asarray(run(zero)[1])==0.)


def test_kk_adjustment_limits_depletion_and_conserves_water():
    from clubb_jax.src.Microphys.KK_microphys_module import KK_microphys_adjust
    rcm=jnp.array([[1.e-4,1.e-5]]);rrm=jnp.array([[1.e-5,1.e-8]])
    nr=jnp.array([[1.e4,1.e3]]);dt=60.
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import KK_Nrm_auto_mean
    def run(rc):
        return KK_microphys_adjust(
            dt, jnp.ones_like(rc), rc, rrm, nr,  # In
            -rrm, rc,                            # In
            rc, -nr,                             # In
            KK_Nrm_auto_mean(rc), True,          # In
            True,                                # In
        )
    out=jax.jit(run)(rcm)
    rr_t,nr_t,rv_t,rc_t,th_t,_=out
    assert np.all(np.asarray(rcm+dt*rc_t)>=-1.e-18)
    assert np.all(np.asarray(rrm+dt*rr_t)>=-1.e-18)
    assert np.all(np.asarray(nr+dt*nr_t)>=-1.e-10)
    np.testing.assert_allclose(rr_t+rv_t+rc_t,0.,atol=1.e-20)
    for a,b in zip(run(rcm)[:5],out[:5]):np.testing.assert_allclose(a,b,rtol=1.e-14,atol=1.e-20)


def test_scheme_startup_returns_zero_tendencies_before_start():
    import inspect
    from clubb_jax.src.Microphys.microphys_driver import calc_microphys_scheme_tendcies
    initialize(microphys_scheme='khairoutdinov_kogan',microphys_start_time=120.)
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.parameter_indices import nparams
    gr=setup_grid(2,100.,0.,600.)
    args={name:None for name in inspect.signature(calc_microphys_scheme_tendcies).parameters}
    args.update(gr=gr,ngrdcol=2,time_current=60.,rcm=jnp.ones((2,gr.nzt)),
        hydromet=jnp.ones((2,gr.nzt,2)), wp2=jnp.ones((2,gr.nzm)),
        wp3=jnp.ones((2,gr.nzt)), clubb_params=jnp.ones((2,nparams)))
    # Skewness is prepared even before startup, but PDF/core inputs stay unused.
    result=calc_microphys_scheme_tendcies(**args)
    for tendency in result[2:-2]:assert np.all(np.asarray(tendency)==0.)
    assert np.all(np.asarray(result[-2]) > 0.)
    assert not np.any(np.asarray(result[-1]))


@pytest.mark.parametrize('upwind',[False,True])
def test_hydrometeor_transport_surface_loss_and_column_independence(upwind):
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid,zt2zm
    from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    from clubb_jax.src.Microphys.advance_microphys_module import microphys_lhs,microphys_rhs,microphys_solve
    initialize(microphys_scheme='khairoutdinov_kogan',l_upwind_diff_sed=upwind)
    gr=setup_grid(2,100.,0.,600.)
    zero=jnp.zeros((2,gr.nzt));zm=jnp.zeros((2,gr.nzm));one=jnp.ones_like(zero)
    q=one*jnp.array([1.e-4,2.e-4])[:,None];dt=10.
    vt=(-one).at[:,-1].set(0.)
    vm=zt2zm(gr.nzm,gr.nzt,2,gr,vt).at[:,-1].set(0.)
    stats=JaxStats.empty(l_sample=False,names=(),ncol=2,max_nlev=gr.nzm)
    def run(q):
        st,ta,ma,ts,sd,lhs=microphys_lhs(gr, gr.ngrdcol, 'rrm', True, dt,
            zm, jnp.zeros(2), zero, vm, vt,
            zero, jnp.ones_like(zm), one, one, False,
            stats)
        st,rhs=microphys_rhs(gr, gr.ngrdcol, 'rrm', dt, True,
                   q, zero, zm, jnp.zeros(2), one,
                   zero, jnp.ones_like(zm), one, one, st)
        return microphys_solve(gr, gr.ngrdcol, 'rrm', True, ta,
                   ma, ts, sd, one, 2,
                   st, lhs, rhs, q, ErrInfo.initialized(2))[3:]
    eager,err=run(q);compiled,jit_err=jax.jit(run)(q)
    assert not err.is_fatal() and not jit_err.is_fatal()
    np.testing.assert_allclose(eager,compiled,rtol=1.e-13,atol=1.e-20)
    np.testing.assert_allclose(eager[1],2*eager[0],rtol=1.e-13)
    loss=jnp.sum((q-eager)*gr.dzt,axis=-1)
    np.testing.assert_allclose(loss,-dt*vt[:,0]*eager[:,0],rtol=1.e-12,atol=1.e-17)
    assert np.all(np.asarray(eager)>=0.)


def test_parameter_override_does_not_reuse_previous_jit_constants():
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import KK_Nrm_auto_mean
    compiled=jax.jit(KK_Nrm_auto_mean)
    initialize(microphys_scheme='khairoutdinov_kogan',r_0=25.e-6)
    first=compiled(jnp.array(1.e-8))
    initialize(microphys_scheme='khairoutdinov_kogan',r_0=50.e-6)
    second=compiled(jnp.array(1.e-8))
    np.testing.assert_allclose(second,first/8.,rtol=1.e-14)
    initialize(microphys_scheme='khairoutdinov_kogan')
    assert parameters.l_hydromet_sed==(True,True)


@pytest.mark.parametrize('in_cloud',[True,False])
def test_predicted_cloud_number_transport_and_clipping(in_cloud):
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    from clubb_jax.src.Microphys.advance_microphys_module import advance_Ncm
    from clubb_jax.src.CLUBB_core.constants_clubb import Nc_in_cloud_min
    initialize(microphys_scheme='morrison', l_predict_Nc=True,
        specify_aerosol='morrison_no_aerosol',l_in_cloud_Nc_diff=in_cloud)
    gr=setup_grid(2,100.,0.,600.)
    zero=jnp.zeros((2,gr.nzt)); one=jnp.ones_like(zero)
    cloud=.5*one; nc=1.e8*one
    stats=JaxStats.empty(l_sample=False,names=(),ncol=2,max_nlev=gr.nzm)
    def run(source):
        return advance_Ncm(gr, gr.ngrdcol, 10., zero, cloud,
                   jnp.zeros((2,gr.nzm)), zero, jnp.ones((2,gr.nzm)), one, one,
                   source, SimpleNamespace(nu_hm=jnp.zeros(2)), False, 2, stats,
                   cloud*nc, nc, ErrInfo.initialized(2))
    unchanged=run(zero)
    np.testing.assert_allclose(unchanged[1],cloud*nc,rtol=1.e-14)
    assert not unchanged[3].is_fatal()
    clipped=jax.jit(run)(-nc)
    np.testing.assert_allclose(clipped[1],cloud*Nc_in_cloud_min,rtol=1.e-14)
    np.testing.assert_allclose(clipped[2],Nc_in_cloud_min,rtol=1.e-14)
    assert np.all(np.asarray(clipped[4][:,[0,-1]])==0.)


@pytest.mark.parametrize('ncol',[1,3])
def test_local_kk_direct_statistics_and_core_contract(ncol):
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.Microphys.KK_microphys_module import KK_local_microphys_driver, KK_local_microphys_core
    _,_,metadata,*_=initialize(microphys_scheme='khairoutdinov_kogan',l_local_kk=True)
    gr=setup_grid(ncol,100.,0.,600.)
    one=jnp.ones((ncol,gr.nzt));zero=jnp.zeros_like(one)
    rcm=one*jnp.arange(1,ncol+1)[:,None]*1.e-4
    hm=jnp.stack((rcm/10.,one*2.e4),axis=-1)
    stats=JaxStats.empty(l_sample=True,names=('rrm_auto','rrm_evap','Nrm_auto'),
        ncol=ncol,max_nlev=gr.nzt)
    def run(hydromet):
        return KK_local_microphys_driver(
            gr, ncol, 10., gr.nzt,                 # In
            2, metadata,                           # In
            False,                                 # In
            280.*one, zero, 90000.*one, one, one,  # In
            one, zero, 100.*one, rcm,              # In
            1.e8*one, rcm, .007*one, hydromet,     # In
            1, one,                                # In
            stats,                                 # InOut
        )
    eager=run(hm);compiled=jax.jit(run)(hm)
    core=KK_local_microphys_core(
        gr, ncol, 10., gr.nzt,           # In
        2, metadata,                     # In
        False,                           # In
        280.*one, 90000.*one, one, one,  # In
        rcm, 1.e8*one, rcm, hm,          # In
        1,                               # In
    )
    for driver_field,core_field in zip(eager[1:],core[:11]):
        np.testing.assert_array_equal(driver_field,core_field)
    for a,b in zip(eager[1:],compiled[1:]):
        np.testing.assert_allclose(a,b,rtol=1.e-12,atol=1.e-20)
    for slot,diag in enumerate((7,9,10)):
        np.testing.assert_array_equal(eager[0].buffers[0][slot],eager[diag])
    assert np.all(np.asarray(eager[0].nsamples[0])==1)
    np.testing.assert_allclose(jnp.sum(eager[1][...,0]+eager[4]+eager[5]),0.,atol=1.e-20)


@pytest.mark.parametrize("mode", ["interactive", "non-interactive"])
@pytest.mark.parametrize("compiled", [False, True])
def test_morrison_silhs_feedback_boundary_and_velocity_statistics(monkeypatch, mode, compiled):
    import inspect
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.parameter_indices import nparams
    from clubb_jax.src.Microphys import microphys_driver as driver
    from clubb_jax.src.Microphys import lh_microphys_driver_module as sampled
    from clubb_jax.src.Microphys import morrison_microphys_module as ordinary

    _, pdf_dim, metadata, *_ = initialize(microphys_scheme="morrison", lh_microphys_type=mode)
    gr = setup_grid(1, 100.0, 0.0, 500.0)
    one = jnp.ones((1, gr.nzt))
    zm = jnp.ones((1, gr.nzm))
    hm = jnp.ones((1, gr.nzt, 2))
    stats = JaxStats.empty(l_sample=True, names=("Vrr",), grids=("zm",), ncol=1, max_nlev=gr.nzm)

    def sampled_tendencies(*args):
        return (
            args[-2], 3 * hm, 4 * hm, one, one, one, one,
            *(zm,) * 5, *(one,) * 7, jnp.zeros(1, dtype=bool),
        )

    def mean_tendencies(*args):
        return (args[-1], jnp.zeros_like(hm), 2 * hm, *(jnp.zeros_like(one),) * 9)

    monkeypatch.setattr(sampled, "lh_microphys_driver", sampled_tendencies)
    monkeypatch.setattr(ordinary, "morrison_microphys_driver", mean_tendencies)
    args = {
        name: None for name in inspect.signature(driver.calc_microphys_scheme_tendcies).parameters
    }
    args.update(
        gr=gr,
        ngrdcol=1,
        dt=10.0,
        time_current=0.0,
        pdf_dim=pdf_dim,
        hydromet_dim=2,
        runtype="synthetic",
        thlm=280 * one,
        p_in_Pa=90000 * one,
        exner=one,
        rho=one,
        rho_zm=zm,
        rtm=0.01 * one,
        rcm=0.001 * one,
        cloud_frac=0.5 * one,
        wm_zt=one,
        wm_zm=zm,
        wp2=zm,
        wp3=one,
        clubb_params=jnp.ones((1, nparams)),
        hydromet=hm,
        Nc_in_cloud=1.0e8 * one,
        hm_metadata=metadata,
        pdf_params=SimpleNamespace(mixt_frac=0.5 * one, chi_1=one, chi_2=one),
        stats=stats,
        Nccnm=one,
    )

    def run(q):
        # Lightweight PDF doubles stay in the closure. Compile the body with
        # jax.jit(run) below; standalone tests use the decorated entry point.
        return driver.calc_microphys_scheme_tendcies.__wrapped__(
            **{**args, "hydromet": q}
        )

    result = (jax.jit(run) if compiled else run)(hm)
    np.testing.assert_array_equal(result[2], 3 * hm if mode == "interactive" else 0.0)
    for tendency in result[10:15]:
        np.testing.assert_array_equal(tendency, 1.0 if mode == "interactive" else 0.0)
    # Both modes update the ordinary velocity statistic once, using the active
    # sedimentation tendency. This catches a source update inside the wrong IF.
    np.testing.assert_array_equal(result[0].buffers[1], 4.0 if mode == "interactive" else 2.0)
    np.testing.assert_array_equal(result[0].nsamples[1], 1)


@pytest.mark.parametrize("nested", [False, True])
def test_silhs_prescribed_probabilities_and_flags_reset_between_cases(nested):
    from clubb_jax.src.SILHS import parameters_silhs

    config = (
        {"eight_cluster_presc_probs": {"cloud_precip_comp1": 0.35}}
        if nested
        else {"eight_cluster_presc_probs%cloud_precip_comp1": 0.35}
    )
    _, _, _, flags, decorr, _, _ = initialize(
        microphys_scheme="khairoutdinov_kogan",
        lh_microphys_type="interactive",
        l_local_kk=True,
        cluster_allocation_strategy=1,
        l_lh_deterministic_test=True,
        l_lh_importance_sampling=False,
        vert_decorr_coef=0.2,
        **config
    )
    assert flags.cluster_allocation_strategy == 1
    assert flags.l_lh_deterministic_test
    assert decorr == 0.2
    assert parameters_silhs.eight_cluster_presc_probs.cloud_precip_comp1 == 0.35
    _, _, _, defaults, decorr, _, _ = initialize(microphys_scheme="khairoutdinov_kogan")
    assert defaults.cluster_allocation_strategy == 3
    assert not defaults.l_lh_deterministic_test
    assert decorr == 0.03
    assert parameters_silhs.eight_cluster_presc_probs.cloud_precip_comp1 == 0.15


@pytest.mark.parametrize("random_option", ["l_lh_importance_sampling", "l_random_k_lh_start"])
def test_deterministic_silhs_rejects_independent_random_draws(random_option):
    config = {"l_lh_importance_sampling": False, "l_random_k_lh_start": False}
    config[random_option] = True
    with pytest.raises(ValueError, match="importance sampling and random starts disabled"):
        initialize(
            microphys_scheme="khairoutdinov_kogan",
            lh_microphys_type="interactive",
            l_local_kk=True,
            l_lh_deterministic_test=True,
            **config
        )


@pytest.mark.parametrize("debug_level", [0, 1])
@pytest.mark.parametrize("zeta", [0.0, 0.25])
def test_initialization_appends_source_configuration_and_correlations(
    monkeypatch, tmp_path, capsys, debug_level, zeta,
):
    from clubb_jax.src.CLUBB_core import error_code
    from clubb_jax.src.CLUBB_core.parameter_indices import nparams, iomicron, izeta_vrnce_rat

    monkeypatch.setattr(error_code, "_debug_level", debug_level)
    case_info = tmp_path / "case_setup.txt"
    case_info.write_text("Existing standalone settings\n")
    clubb_params = jnp.zeros((1, nparams)).at[0, iomicron].set(0.2)
    clubb_params = clubb_params.at[0, izeta_vrnce_rat].set(zeta)
    config = {
        "microphys_scheme": "khairoutdinov_kogan",
        "lh_microphys_type": "interactive", "l_local_kk": True,
        "lh_num_samples": 8, "l_lh_importance_sampling": False,
        "l_lh_deterministic_test": True, "c_evap": 0.9,
    }
    init_microphys(
        0, "rico_silhs", config, case_info,  # In
        1000.0, 1000.0,                      # In
        clubb_params,                        # In
        False,                               # In
        True,                                # InOut
        True,                                # InOut
    )
    text = case_info.read_text()
    output = capsys.readouterr().out
    assert text.startswith("Existing standalone settings\n")
    if debug_level == 0:
        assert text == "Existing standalone settings\n"
        assert output == ""
    else:
        assert text.index("&microphysics_setting") < text.index("&SILHS_setting")
        assert "lh_microphys_type = interactive" in text
        assert "lh_num_samples = 8" in text
        assert "C_evap = 0.9" in text
        assert "l_lh_deterministic_test = True" in text
        assert "l_lh_importance_sampling = False" in text
        assert text.count("hmp2_ip_on_hmm2_ip_slope%Ni =") == 2
        for location in ("in cloud", "below cloud"):
            label = f"Correlation array (approximate); {location}:"
            assert (label in text) == (zeta == 0.0)
            assert (label in output) == (zeta == 0.0)
        if zeta == 0.0:
            rows = text.split("Correlation array (approximate); in cloud:\n")[1].splitlines()[:6]
            matrix = np.array([[float(value) for value in row.split()] for row in rows])
            np.testing.assert_array_equal(np.diag(matrix), 1.0)
            np.testing.assert_array_equal(matrix, matrix.T)
