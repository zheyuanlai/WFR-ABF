"""Smoke test: validate eb_abffr_core helpers vs the notebook's numpy versions,
then run a short batched simulation to confirm it executes end-to-end."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import torch
import eb_abffr_core as eb

dev = eb.DEVICE
dt = eb.DTYPE
print("device:", dev, "dtype:", dt)

x_grid, dx, eval_mask, idx0 = eb.build_grid(dev, dt)
xg_np = x_grid.cpu().numpy()

# ---- reference numpy implementations from the notebook ----
beta, H, oout, oin, s = 8.0, 2.5, 1.0, 25.0, 0.25
def omega_np(x): return oout + (oin - oout) * np.exp(-x*x/(2*s*s))
def smooth_np(v, k, r):
    pad = min(r, len(v)-1)
    return np.convolve(np.pad(v, (pad,pad), mode="reflect"),
                       k[r-pad:len(k)-(r-pad)], mode="valid")
def gk_np(bw):
    r = max(1, int(round(4.0*bw/dx)))
    k = np.exp(-0.5*(np.arange(-r,r+1)*dx/bw)**2)
    return k/(k.sum()*dx), r

# 1) gaussian kernel matches
k_np, r_np = gk_np(0.07)
k_t, r_t = eb.gaussian_kernel(0.07, dx, dev, dt)
assert r_np == r_t
print("kernel max abs diff:", np.abs(k_np - k_t.cpu().numpy()).max())

# 2) smooth matches on a random vector
rng = np.random.default_rng(0)
v = rng.standard_normal(eb.N_GRID)
sm_np = smooth_np(v, k_np, r_np)
sm_t = eb.smooth(torch.as_tensor(v, device=dev, dtype=dt).unsqueeze(0), k_t, r_t, dx)[0].cpu().numpy()
print("smooth max abs diff:", np.abs(sm_np - sm_t).max())

# 3) cumtrapz matches
def cumtrapz_np(y): return np.concatenate([[0.0], np.cumsum(0.5*(y[1:]+y[:-1])*dx)])
ct_np = cumtrapz_np(v)
ct_t = eb.cumtrapz(torch.as_tensor(v, device=dev, dtype=dt).unsqueeze(0), dx)[0].cpu().numpy()
print("cumtrapz max abs diff:", np.abs(ct_np - ct_t).max())

# 4) binned_density vs notebook
def binned_density_np(X, k, r):
    idx = np.clip(np.round((X - eb.XMIN)/dx).astype(int), 0, eb.N_GRID-1)
    p = smooth_np(np.bincount(idx, minlength=eb.N_GRID).astype(float), k, r)/len(X)
    return np.maximum(p/max(np.trapezoid(p, dx=dx), eb.EPS), eb.EPS)
keta_np, reta_np = gk_np(0.10)
keta_t, reta_t = eb.gaussian_kernel(0.10, dx, dev, dt)
X = rng.uniform(-1.5, 1.5, 256)
bd_np = binned_density_np(X, keta_np, reta_np)
bd_t = eb.binned_density(torch.as_tensor(X, device=dev, dtype=dt).unsqueeze(0), keta_t, reta_t, dx)[0].cpu().numpy()
print("binned_density max abs diff:", np.abs(bd_np - bd_t).max())

# 5) interp1d vs np.interp
prof = rng.standard_normal(eb.N_GRID)
xq = rng.uniform(-1.8, 1.8, 100)
ip_np = np.interp(xq, xg_np, prof)
ip_t = eb.interp1d(torch.as_tensor(xq, device=dev, dtype=dt).unsqueeze(0),
                   torch.as_tensor(prof, device=dev, dtype=dt).unsqueeze(0), dx)[0].cpu().numpy()
print("interp1d max abs diff:", np.abs(ip_np - ip_t).max())

# 6) reference profiles centering
F_ref_b, Fp_ref_b = eb.reference_profiles(
    x_grid, eval_mask,
    torch.tensor([[beta]], device=dev, dtype=dt),
    torch.tensor([[H]], device=dev, dtype=dt),
    torch.tensor([[oout]], device=dev, dtype=dt),
    torch.tensor([[oin]], device=dev, dtype=dt),
    torch.tensor([[s]], device=dev, dtype=dt))
Fref_np = H*(xg_np**2-1)**2 + np.log(omega_np(xg_np))/beta
Fref_np = Fref_np - Fref_np[eval_mask.cpu().numpy()].mean()
print("F_ref max abs diff:", np.abs(Fref_np - F_ref_b[0].cpu().numpy()).max())

# 7) short batched run: 2 configs x 1 seed x [abf, fr_estimated], 2000 steps
print("\n--- short batched run (2000 steps) ---")
cfgA = eb.PhysConfig(n_steps=2000, save_every=200, N=256)
cfgB = eb.PhysConfig(n_steps=2000, save_every=200, N=256, omega_in=5.0)
spec = eb.BatchSpec(configs=[cfgA, cfgB], seeds=[0, 0],
                    methods=[eb.ABF, eb.FR_ESTIMATED], batch_seed=0)
import time
t0 = time.time()
recs = eb.simulate_batch(spec)
torch.cuda.synchronize() if torch.cuda.is_available() else None
print(f"ran {len(recs)} runs in {time.time()-t0:.1f}s")
for r in recs:
    print(f"  oin={r['config']['omega_in']:5.1f} {r['method']:13s} "
          f"L2(F)={r['final_l2_f']:.4f} L2(F')={r['final_l2_fp']:.4f} "
          f"ESS={r['final_ess']:.0f} repl_frac={r['repl_fraction']:.4f} "
          f"n_die={r['n_die']:.0f} n_clone={r['n_clone']:.0f}")
    assert np.isfinite(r['final_l2_f']), "NaN in L2(F)!"
print("\nsmoke OK")
