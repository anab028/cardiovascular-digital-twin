"""Reproducible patient-disjoint synthetic experiments."""
import argparse
import json
import platform
from pathlib import Path
import importlib.metadata
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesRegressor
from .physiology import Patient, simulate, fit_patient

T = np.arange(0, 16.001, .04)
HISTORY = T <= 8
FUTURE = ~HISTORY


def make_case(rng, pid, scenario="standard", missing=.5, noise=2.):
    r, c, sv, hr = rng.uniform(.8,1.5), rng.uniform(.9,2.2), rng.uniform(55,90), rng.uniform(55,100)
    z = rng.uniform(.035,.075) if scenario == "shift" else rng.uniform(.005,.025)
    patient = Patient(r,c,sv,hr,z)
    intervention = 1.5 if scenario == "intervention" else rng.uniform(.8,1.2)
    multiplier = np.where(HISTORY, 1., intervention)
    truth, q = simulate(T, patient, multiplier)
    obs = truth[HISTORY] + rng.normal(0,noise,HISTORY.sum())
    obs[rng.random(len(obs)) < missing] = np.nan
    # Preserve enough observations to make the task defined, without future data.
    obs[:8] = truth[:8] + rng.normal(0,noise,8)
    estimated = fit_patient(T[HISTORY], obs, hr, sv)
    fitted = Patient(estimated[0], estimated[1], sv, hr)
    mechanistic, _ = simulate(T, fitted, multiplier, initial_pressure=estimated[2])
    history_filled = np.interp(T[HISTORY], T[HISTORY][np.isfinite(obs)], obs[np.isfinite(obs)])
    history_summary = history_filled[np.linspace(0,len(obs)-1,24,dtype=int)]
    tf = T[FUTURE]
    phase = 2*np.pi*hr*tf/60
    # Features only from observed history, known inflow, and planned intervention.
    base = np.r_[history_summary, hr, sv, intervention]
    x = np.column_stack([np.tile(base,(len(tf),1)), tf-8, np.sin(phase), np.cos(phase)])
    xh = np.column_stack([x,mechanistic[FUTURE],np.tile(estimated[:2],(len(tf),1))])
    return dict(id=pid, x=x, xh=xh, y=truth[FUTURE], physics=mechanistic[FUTURE],
                obs=obs, truth=truth, estimated=estimated[:2], parameters=np.array([r,c]),
                last=obs[np.isfinite(obs)][-1], scenario=scenario)


def stack(cases,key):
    return np.concatenate([c[key] for c in cases])


def conformal_radius(scores, alpha=.1):
    """Finite-sample split-conformal order statistic over patient max errors."""
    scores = np.asarray(scores)
    if scores.ndim != 1 or len(scores)==0 or not np.all(np.isfinite(scores)) or not 0 < alpha < 1:
        raise ValueError("Need finite patient scores and alpha in (0,1)")
    k = int(np.ceil((len(scores)+1)*(1-alpha)))
    return float(np.sort(scores)[k-1]) if k <= len(scores) else float("inf")


def predictions(case, data_model, hybrid_model):
    return {"persistence":np.full_like(case["y"],case["last"]),
            "physiology":case["physics"],
            "data_driven":data_model.predict(case["x"]),
            "hybrid":case["physics"]+hybrid_model.predict(case["xh"])}


def run(output, seed=42, train_n=80, calibration_n=30, test_n=30):
    if min(train_n,calibration_n,test_n)<10:
        raise ValueError("Use at least ten patients per split")
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    train=[make_case(rng,f"train-{i}") for i in range(train_n)]
    calibration=[make_case(rng,f"cal-{i}") for i in range(calibration_n)]
    data=ExtraTreesRegressor(n_estimators=80,min_samples_leaf=6,random_state=seed,n_jobs=1)
    hybrid=ExtraTreesRegressor(n_estimators=80,min_samples_leaf=6,random_state=seed,n_jobs=1)
    data.fit(stack(train,"x"),stack(train,"y"))
    hybrid.fit(stack(train,"xh"),stack(train,"y")-stack(train,"physics"))
    names=["persistence","physiology","data_driven","hybrid"]
    radii={name:conformal_radius([np.max(np.abs(predictions(c,data,hybrid)[name]-c["y"])) for c in calibration]) for name in names}
    rows=[]; examples={}; patient_rows=[]; manifests={"train":[c["id"] for c in train],"calibration":[c["id"] for c in calibration]}
    for scenario,missing,noise in [("standard",.5,2),("sparse",.85,2),("noisy",.5,5),("intervention",.5,2),("shift",.5,2)]:
        cases=[make_case(rng,f"{scenario}-{i}",scenario,missing,noise) for i in range(test_n)]
        manifests[scenario]=[c["id"] for c in cases]
        preds=[predictions(c,data,hybrid) for c in cases]
        for name in names:
            errors=np.array([p[name]-c["y"] for p,c in zip(preds,cases)])
            maes=np.mean(np.abs(errors),axis=1)
            boot=np.mean(rng.choice(maes,size=(1000,len(maes)),replace=True),axis=1)
            rows.append(dict(scenario=scenario,model=name,mae=float(maes.mean()),
                mae_ci95=np.quantile(boot,[.025,.975]).tolist(),
                trajectory_coverage=float(np.mean(np.max(np.abs(errors),axis=1)<=radii[name])),
                point_coverage=float(np.mean(np.abs(errors)<=radii[name])),interval_width=2*radii[name]))
            patient_rows.extend(dict(patient=c["id"],scenario=scenario,model=name,mae=float(m)) for c,m in zip(cases,maes))
        examples[scenario]=(cases[0],preds[0])
    fig,axs=plt.subplots(1,2,figsize=(13,4.5),layout="constrained")
    for ax,scenario in zip(axs,["standard","shift"]):
        c,p=examples[scenario]
        ax.plot(T,c["truth"],color="#142b49",label="Synthetic truth")
        ax.scatter(T[HISTORY],c["obs"],s=7,alpha=.5,color="#667085",label="Observed history")
        ax.plot(T[FUTURE],p["physiology"],color="#ee9234",label="Fitted physiology")
        ax.plot(T[FUTURE],p["hybrid"],color="#008c95",label="Hybrid forecast")
        ax.fill_between(T[FUTURE],p["hybrid"]-radii["hybrid"],p["hybrid"]+radii["hybrid"],color="#008c95",alpha=.15,label="90% nominal band")
        ax.axvline(8,color="gray",ls="--"); ax.set(title=scenario.title(),xlabel="Time (s)",ylabel="Pressure (mmHg)")
    axs[0].legend(fontsize=8); fig.savefig(output/"forecast.png",dpi=160); plt.close(fig)
    # Exact scale symmetry: R->R/k, C->kC, SV->kSV yields identical pressure.
    base=Patient(1.2,1.5,70,75)
    fig,ax=plt.subplots(figsize=(8,4),layout="constrained")
    curves=[]
    for k in [.7,1.,1.4]:
        p,_=simulate(T,Patient(base.resistance/k,base.compliance*k,base.stroke_volume*k,base.heart_rate))
        curves.append(p); ax.plot(T,p,label=f"R={base.resistance/k:.2f}, C={base.compliance*k:.2f}, SV={base.stroke_volume*k:.0f}",ls={.7:"-",1.:"--",1.4:":"}[k])
    ax.set(title="Different parameters, identical pressure: structural ambiguity",xlabel="Time (s)",ylabel="Pressure (mmHg)"); ax.legend(); fig.savefig(output/"identifiability.png",dpi=160); plt.close(fig)
    # Matched-model noiseless recovery is separate from misspecified fitting.
    recovery=[]
    for i in range(20):
        r,c,sv,hr=rng.uniform(.8,1.5),rng.uniform(.9,2.2),rng.uniform(55,90),rng.uniform(55,100)
        p,_=simulate(T[HISTORY],Patient(r,c,sv,hr))
        obs=p+rng.normal(0,2,len(p)); obs[rng.random(len(obs))<.5]=np.nan
        est=fit_patient(T[HISTORY],obs,hr,sv)[:2]
        recovery.append(np.abs(est-[r,c])/[r,c])
    result={"seed":seed,"counts":{"train":train_n,"calibration":calibration_n,"test_per_scenario":test_n},
        "environment":{"python":platform.python_version(),**{x:importlib.metadata.version(x) for x in ["numpy","scipy","scikit-learn","matplotlib"]}},
        "metrics":rows,"matched_model_mean_relative_parameter_error":dict(zip(["R","C"],np.mean(recovery,axis=0).tolist())),
        "scale_symmetry_max_difference":float(np.max(np.abs(np.array(curves)-curves[1])))}
    (output/"metrics.json").write_text(json.dumps(result,indent=2)+"\n")
    (output/"patient_metrics.json").write_text(json.dumps(patient_rows,indent=2)+"\n")
    (output/"splits.json").write_text(json.dumps(manifests,indent=2)+"\n")
    lines=["# Benchmark results","","Synthetic, single-seed proof of concept. Lower MAE is better. Confidence intervals bootstrap patients, not time points.","","| Scenario | Model | MAE (mmHg), 95% CI | Trajectory coverage | Band width (mmHg) |","|---|---|---|---|---|"]
    for row in rows:
        lo,hi=row["mae_ci95"]
        lines.append(f'| {row["scenario"]} | {row["model"]} | {row["mae"]:.2f} [{lo:.2f}, {hi:.2f}] | {row["trajectory_coverage"]:.0%} | {row["interval_width"]:.2f} |')
    lines += ["","Coverage targets the complete discrete forecast trajectory at nominal 90% under exchangeability. Shifted scenarios have no coverage guarantee. Bands are calibrated against noiseless synthetic truth, not future noisy measurements.","",f"Matched two-element generator, known HR/SV: mean relative parameter errors {result['matched_model_mean_relative_parameter_error']}.","","The three-element generator creates intentional model mismatch. Its fitted two-element parameters are effective parameters, not necessarily the generating physiological values."]
    (output/"REPORT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({"output":str(output),"metrics":rows[:4]},indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",default="results")
    parser.add_argument("--seed",type=int,default=42)
    parser.add_argument("--train",type=int,default=80)
    parser.add_argument("--calibration",type=int,default=30)
    parser.add_argument("--test",type=int,default=30)
    args=parser.parse_args(); run(args.output,args.seed,args.train,args.calibration,args.test)

if __name__=="__main__":
    main()
