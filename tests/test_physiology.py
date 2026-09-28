import numpy as np
import pytest
from cardiotwin.physiology import Patient, simulate, fit_patient
from cardiotwin.benchmark import conformal_radius


def test_constant_flow_matches_analytic_solution():
    # HR near zero gives nearly constant sinusoidal inflow over this window.
    from cardiotwin.physiology import flow
    t=np.arange(0,5,.01); p=Patient(1.2,1.5,70,75)
    # Cancel the sinusoid so the integrator sees precisely constant inflow.
    mult=87.5/flow(t,75,70)
    y,_=simulate(t,p,mult)
    expected=110+(85-110)*np.exp(-t/(1.2*1.5))
    np.testing.assert_allclose(y,expected,atol=1e-10)


def test_noiseless_parameter_recovery():
    t=np.arange(0,8,.04); p=Patient(1.3,1.7,72,80)
    y,_=simulate(t,p); y[::3]=np.nan
    fitted=fit_patient(t,y,80,72)
    np.testing.assert_allclose(fitted,[1.3,1.7,85],rtol=1e-5)


def test_pressure_only_scale_symmetry():
    t=np.arange(0,8,.04)
    y,_=simulate(t,Patient(1.2,1.5,70,75))
    other,_=simulate(t,Patient(.6,3.,140,75))
    np.testing.assert_allclose(y,other,atol=1e-10)


def test_future_intervention_cannot_change_history():
    t=np.arange(0,16,.04); p=Patient(1.2,1.5,70,75,.02)
    baseline,_=simulate(t,p); changed,_=simulate(t,p,np.where(t>8,1.5,1))
    np.testing.assert_array_equal(baseline[t<=8],changed[t<=8])
    assert np.max(np.abs(baseline[t>8]-changed[t>8]))>10


def test_conformal_order_statistic():
    assert conformal_radius(np.arange(1,20),.1)==18
    assert np.isinf(conformal_radius([1,2],.1))


def test_invalid_inputs():
    with pytest.raises(ValueError): Patient(-1,1,70,75)
    with pytest.raises(ValueError): simulate([0,.1,.3],Patient(1,1,70,75))
    with pytest.raises(ValueError): fit_patient(np.arange(10),np.full(10,np.nan),75,70)


def test_timestep_refinement_converges():
    patient=Patient(1.2,1.5,70,75,.02)
    grids=[np.arange(0,4.001,dt) for dt in [.04,.02,.01]]
    curves=[simulate(t,patient)[0] for t in grids]
    coarse_error=np.mean(np.abs(curves[0]-curves[2][::4]))
    fine_error=np.mean(np.abs(curves[1]-curves[2][::2]))
    assert fine_error < coarse_error


def test_missing_history_yields_finite_features_and_forecast():
    from cardiotwin.benchmark import make_case
    case=make_case(np.random.default_rng(7),'test',missing=.85)
    assert np.isnan(case['obs']).any()
    assert np.isfinite(case['x']).all()
    assert np.isfinite(case['xh']).all()
    assert np.isfinite(case['physics']).all()
    assert case['x'].shape[0]==case['y'].shape[0]
