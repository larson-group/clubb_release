"""GFDL is deliberately disabled; dummy interfaces must fail explicitly."""
import inspect
import pytest
from clubb_jax.src.Microphys import gfdl_activation


@pytest.mark.parametrize('name', ['aer_act_clubb_quadrature_Gauss', 'erff', 'Loading', 'get_unit'])
def test_disabled_activation_interfaces_cannot_simulate(name):
    routine=getattr(gfdl_activation,name)
    args={name:None for name in inspect.signature(routine).parameters}
    with pytest.raises(NotImplementedError,match='GFDL'):
        routine(**args)
