"""Direct forward witnesses and independent numerical release checks."""
import json, math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from region_geometry import lp_network,planar_network
from finite_budget import B,rate
DATA=Path(__file__).resolve().parent/"data"
s=json.loads((DATA/"s2_regions.json").read_text(encoding="utf-8"))
assert s["schema"]=="strict-regions-v2"
assert len(s["region_counts"])==9 and len(s["K_at_alpha1"])==8
witnesses=0; crosschecks=0
for row in s["region_counts"]:
 counts=[]
 for i,d in enumerate(row["draws"]):
  Ws=[np.asarray(W) for W in d["weights"]]
  assert d["strict_count"]==len(d["strict_masks"])==len(d["strict_witnesses"])
  assert d["ordinary_count"]==d["strict_count"]+d["forced_zero_count"]
  for mask,x in zip(d["strict_masks"],d["strict_witnesses"]):
   x=np.asarray(x);seen=[]
   for W in Ws:
    y=W@x
    assert np.all(y!=0)
    seen.extend((y>0).astype(int).tolist());x=np.maximum(y,0)
   assert seen==mask
   witnesses+=1
  if row["n"]==2 and row["L"]==2:
   assert d["strict_count"]==3+int(np.sum(Ws[1][:,0]*Ws[1][:,1]<0))
  if row["n"]==2 and i==0:
   b=lp_network(Ws)
   assert set(b["strict"])==set(map(tuple,d["strict_masks"]))
   crosschecks+=1
  counts.append(d["strict_count"])
 assert abs(sum(counts)/len(counts)-row["strict_mean"])<1e-12
 assert min(counts)==row["strict_min"] and max(counts)==row["strict_max"]
for r in s["K_at_alpha1"].values():
 k=np.asarray(r["K_values"])
 assert len(k)==r["trials"] and np.isfinite(k).all() and (k>=0).all()
 assert float(np.median(k))==r["median_K1"]
# A separate minimizer checks both sides of the corrected integer depth.
for L,positive in [(237380,False),(237381,True)]:
 opt=minimize_scalar(lambda t:B(L,1.95,t),bounds=(1e-10,.249999),method="bounded",options={"xatol":1e-14})
 assert (-opt.fun>0)==positive
 assert abs(rate(L,1.95)+opt.fun)<1e-7
print(json.dumps(dict(strict_forward_witnesses=witnesses,all_mask_LP_crosschecks=crosschecks,
    width_two_L2_analytic_trials=60,depth_195_adjacent_integer_check="passed"),indent=2))
