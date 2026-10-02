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


def initialize(**config):
    return init_microphys(0, 'synthetic', config, None, 1000., 1000.,
                          jnp.zeros((1, 1)), False, True, True)


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


@pytest.mark.parametrize('config,reason', [
    ({'lh_microphys_type':'interactive'}, 'SILHS'),
    ({'microphys_scheme':'coamps'}, 'Unsupported'),
    ({'microphys_scheme':'simplified_ice'}, 'Unsupported'),
    ({'l_gfdl_activation':True}, 'GFDL'),
    ({'microphys_scheme':'morrison','l_cloud_sed':True}, 'sedimentation'),
    ({'microphys_scheme':'khairoutdinov_kogan','l_predict_Nc':True}, 'l_predict_Nc'),
    ({'l_morr_xp2_mc':True}, 'l_morr_xp2_mc')])
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
        return morrison_microphys_driver(gr, ncol, 10., nzt, dimension, metadata, False,
            jnp.linspace(258.,285.,nzt)[None,:]*one, zero, 90000.*one, one, one, 0.5*one, 0.2*one, 100.*one,
            1.e-4*one, 1.e8*one, zero, 0.007*one, hm, 1, one, stats)
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
        return KK_microphys_adjust(dt,jnp.ones_like(rc),rc,rrm,nr,
            -rrm,rc,rc,-nr,KK_Nrm_auto_mean(rc),True,True)
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
    for tendency in result[2:-1]:assert np.all(np.asarray(tendency)==0.)
    assert np.all(np.asarray(result[-1]) > 0.)


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
        return KK_local_microphys_driver(gr,ncol,10.,gr.nzt,2,metadata,False,
            280.*one,zero,90000.*one,one,one,one,zero,100.*one,rcm,
            1.e8*one,rcm,.007*one,hydromet,1,one,stats)
    eager=run(hm);compiled=jax.jit(run)(hm)
    core=KK_local_microphys_core(gr,ncol,10.,gr.nzt,2,metadata,False,
        280.*one,90000.*one,one,one,rcm,1.e8*one,rcm,hm,1)
    for driver_field,core_field in zip(eager[1:],core[:11]):
        np.testing.assert_array_equal(driver_field,core_field)
    for a,b in zip(eager[1:],compiled[1:]):
        np.testing.assert_allclose(a,b,rtol=1.e-12,atol=1.e-20)
    for slot,diag in enumerate((7,9,10)):
        np.testing.assert_array_equal(eager[0].buffers[0][slot],eager[diag])
    assert np.all(np.asarray(eager[0].nsamples[0])==1)
    np.testing.assert_allclose(jnp.sum(eager[1][...,0]+eager[4]+eager[5]),0.,atol=1.e-20)
