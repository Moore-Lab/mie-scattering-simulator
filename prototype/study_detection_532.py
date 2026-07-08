"""Detection physics at the REAL detection wavelength (532 nm imaging beam),
silica sphere in vacuum, radius 3-20 um. Trap is 1064 nm but detection is 532.
Adds backward collection (available on this apparatus) and a grid-convergence
check (angle-integrated g vs coefficient g)."""
import numpy as np, mie_core as mc, json
from information_pattern import (sphere_grid, scattered_intensity_map,
                                 info_density, cone_mask)

LAM = 0.532          # um, detection wavelength
N_SILICA = 1.46      # fused silica at 532 nm
NMED = 1.0           # vacuum

def wiscombe(x): return int(np.ceil(x + 4.05*x**(1/3) + 2))

def g_from_map(I, TH, W):
    # <cos theta> weighted by scattered intensity -> should match coefficient g
    num = np.sum(I*np.cos(TH)*W); den = np.sum(I*W)
    return num/den

def eta_curve(F, TH, W, direction, na_max=0.95, npts=80):
    Ftot = np.sum(F*W); NA = np.linspace(0.05, na_max, npts); eta=np.empty(npts)
    for i,na in enumerate(NA):
        a = np.arcsin(min(na/NMED,1.0)); m = cone_mask(TH,a,direction)
        eta[i] = np.sum(F[m]*W[m])/Ftot
    return NA, eta

def summarize(a_um, ntheta, nphi):
    x = 2*np.pi*a_um*NMED/LAM; m = N_SILICA+0j
    qext,qsca,qback,g = mc.efficiencies(m, x)
    theta,phi,TH,PH,W = sphere_grid(ntheta,nphi)
    I = scattered_intensity_map(m,x,theta,phi)
    g_map = g_from_map(I,TH,W)
    fwd = TH<=np.pi/2
    out = dict(a_um=a_um, x=round(float(x),2), n_max=wiscombe(x),
               g_coeff=round(float(g),4), g_map=round(float(g_map),4),
               g_relerr=abs(g_map-g)/abs(g),
               pwr_fwd=float(np.sum(I[fwd]*W[fwd])/np.sum(I*W)))
    for axis in ("x","z"):
        F = info_density(I,TH,PH,axis); Ftot=np.sum(F*W)
        out[f"info_{axis}_fwd"]=float(np.sum(F[fwd]*W[fwd])/Ftot)
        Fpol=np.sum(F*W,axis=1); out[f"{axis}_peak_deg"]=float(np.degrees(theta[np.argmax(Fpol)]))
        for d in ("forward","backward"):
            NA,eta=eta_curve(F,TH,W,d)
            for na in (0.5,0.8,0.95):
                out[f"eta_{axis}_{d[:3]}_NA{na}"]=round(float(np.interp(na,NA,eta)),3)
    return out

SIZES=[3.0,5.0,8.0,12.0,20.0]
def grid_for(a): return (1400 if a>=12 else (1000 if a>=8 else 700), 420)

def run():
    res=[]
    print(f"{'a(um)':>6}{'x':>8}{'n_max':>6}{'g_relerr':>10}  {'pwrF':>6} "
          f"{'iX_F':>6}{'Xpk':>6}{'etaXfwd.8':>10}{'etaXbwd.8':>10}  "
          f"{'iZ_F':>6}{'Zpk':>6}{'etaZfwd.8':>10}{'etaZbwd.8':>10}")
    for a in SIZES:
        nt,nph = grid_for(a)
        r=summarize(a,nt,nph); res.append(r)
        print(f"{a:>6}{r['x']:>8.1f}{r['n_max']:>6}{r['g_relerr']:>10.1e}  "
              f"{r['pwr_fwd']:>6.3f} {r['info_x_fwd']:>6.3f}{r['x_peak_deg']:>6.0f}"
              f"{r['eta_x_for_NA0.8']:>10.3f}{r['eta_x_bac_NA0.8']:>10.3f}  "
              f"{r['info_z_fwd']:>6.3f}{r['z_peak_deg']:>6.0f}"
              f"{r['eta_z_for_NA0.8']:>10.3f}{r['eta_z_bac_NA0.8']:>10.3f}")
    return res

if __name__ == "__main__":
    res = run()
    json.dump(res, open("information_pattern_results_532.json","w"), indent=2)
    print("\nwrote information_pattern_results_532.json  (g_relerr = grid convergence check)")
