"""
Information radiation pattern for position sensing of a levitated sphere.

Physics
-------
Sphere at position r_s in an incident plane wave (k along +z, x-polarized).
Displacing the sphere by r_s imprints a direction-dependent phase on the
far-field scattered wave:

    E_s(s_hat; r_s) = E_s(s_hat; 0) * exp[i k (k_hat - s_hat) . r_s]

so the derivative w.r.t. a displacement component j is

    d E_s / d r_{s,j} = i k (k_hat - s_hat)_j E_s

The (classical/quantum) Fisher information *density* for estimating r_{s,j} from
an ideal coherent (mode-matched homodyne) measurement of the scattered field is

    dF_j/dOmega  proportional to  |d E_s/d r_{s,j}|^2
                 = k^2 (k_hat - s_hat)_j^2 |E_s(s_hat)|^2

This is the "information radiation pattern". It is NOT the intensity pattern
|E_s|^2 -- it is weighted by (k_hat - s_hat)_j^2, which vanishes in the forward
direction for transverse motion and is largest to the sides/back.

The collection efficiency of a detector spanning solid angle Omega is

    eta_j(Omega) = F_j(Omega) / F_j(4pi)

which upper-bounds how close a real detector gets to the imprecision-noise SQL.
Measurement backaction (recoil heating) is set by total scattered power and is
independent of Omega, so eta_j is the central optimization figure of merit.

This module uses plane-wave Mie (the trap is really a focused beam -> GLMT in the
full project) but the information-pattern *structure* and the NA-collection
tradeoffs it reveals are the quantities the project must reproduce and refine.
"""
import json
import numpy as np
import mie_core as mc

LAMBDA = 1.064e-6  # vacuum wavelength, m


def sphere_grid(ntheta=400, nphi=240):
    """Gauss-Legendre in cos(theta), uniform in phi. Returns grids + weights."""
    mu, wmu = np.polynomial.legendre.leggauss(ntheta)   # mu in [-1,1]
    theta = np.arccos(mu)                                # (ntheta,)
    phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    dphi = 2 * np.pi / nphi
    TH, PH = np.meshgrid(theta, phi, indexing="ij")      # (ntheta, nphi)
    # solid-angle weight dOmega = dphi * w_mu   (w_mu already integrates dcos)
    W = np.outer(wmu, np.full(nphi, dphi))               # (ntheta, nphi)
    return theta, phi, TH, PH, W


def scattered_intensity_map(m, x, theta, phi):
    """
    |E_s|^2 on the (theta,phi) grid for x-polarized incident plane wave.
    E_theta ~ S2(theta) cos(phi), E_phi ~ -S1(theta) sin(phi)  (B&H amplitude
    scattering matrix). Returns array (ntheta, nphi), arbitrary common units.
    """
    S1, S2 = mc.scattering_amplitudes(m, x, theta)      # each (ntheta,)
    I1 = np.abs(S1) ** 2
    I2 = np.abs(S2) ** 2
    cph = np.cos(phi)[None, :]
    sph = np.sin(phi)[None, :]
    # |E_theta|^2 + |E_phi|^2
    return I2[:, None] * cph ** 2 + I1[:, None] * sph ** 2


def info_density(intensity, TH, PH, axis):
    """
    Fisher-information density (up to common k^2 and constants) for displacement
    along 'axis' in {'x','y','z'}. Weight factor (k_hat - s_hat)_j^2 with k_hat=z.
      s_hat = (sinT cosP, sinT sinP, cosT)
      x: (0 - sinT cosP)^2 ;  y: (0 - sinT sinP)^2 ;  z: (1 - cosT)^2
    """
    sinT = np.sin(TH)
    if axis == "x":
        wj = (sinT * np.cos(PH)) ** 2
    elif axis == "y":
        wj = (sinT * np.sin(PH)) ** 2
    elif axis == "z":
        wj = (1.0 - np.cos(TH)) ** 2
    else:
        raise ValueError(axis)
    return wj * intensity


def cone_mask(TH, half_angle, direction):
    """Boolean mask for a collection cone of given half-angle about +z (forward)
    or -z (backward)."""
    if direction == "forward":
        return TH <= half_angle
    elif direction == "backward":
        return TH >= (np.pi - half_angle)
    raise ValueError(direction)


def collection_efficiency(m, x, axis, direction, na_medium=1.0, npts=200):
    """
    eta_axis(NA) = fraction of total Fisher info for 'axis' displacement collected
    by a cone in 'direction' as a function of numerical aperture NA = n*sin(alpha).
    Returns (NA_array, eta_array).
    """
    theta, phi, TH, PH, W = sphere_grid(600, 360)
    I = scattered_intensity_map(m, x, theta, phi)
    F = info_density(I, TH, PH, axis)
    F_total = np.sum(F * W)
    NA = np.linspace(0.05, na_medium, npts)
    eta = np.empty_like(NA)
    for i, na in enumerate(NA):
        alpha = np.arcsin(min(na / na_medium, 1.0))
        mask = cone_mask(TH, alpha, direction)
        eta[i] = np.sum(F[mask] * W[mask]) / F_total
    return NA, eta


def summarize(m, x, label):
    theta, phi, TH, PH, W = sphere_grid(600, 360)
    I = scattered_intensity_map(m, x, theta, phi)
    out = {"label": label, "m": m.real, "x": float(x)}
    # forward vs backward power fraction
    fwd = TH <= np.pi / 2
    out["power_forward_frac"] = float(np.sum(I[fwd] * W[fwd]) / np.sum(I * W))
    for axis in ("x", "z"):
        F = info_density(I, TH, PH, axis)
        Ftot = np.sum(F * W)
        # where is the info? forward hemisphere fraction
        out[f"info_{axis}_forward_frac"] = float(np.sum(F[fwd] * W[fwd]) / Ftot)
        # peak polar angle of the info pattern (phi-averaged)
        Fpol = np.sum(F * W, axis=1)
        out[f"info_{axis}_peak_theta_deg"] = float(np.degrees(theta[np.argmax(Fpol)]))
        # efficiency at a few representative NAs, forward collection
        for na in (0.5, 0.8, 0.95):
            NA, eta = collection_efficiency(m, x, axis, "forward", 1.0, 60)
            out[f"eta_{axis}_fwd_NA{na}"] = float(np.interp(na, NA, eta))
    return out


if __name__ == "__main__":
    LAMBDA_um = LAMBDA * 1e6
    results = []
    for a_um in [0.1, 0.5, 2.5, 5.0]:
        x = 2 * np.pi * a_um / LAMBDA_um
        results.append(summarize(1.45 + 0j, x, f"silica a={a_um}um"))

    print(f"{'case':>18} {'x':>7} {'pwr_fwd':>8} "
          f"{'infoX_fwd':>9} {'Xpeak°':>7} {'etaX@0.8':>9} "
          f"{'infoZ_fwd':>9} {'Zpeak°':>7} {'etaZ@0.8':>9}")
    print("-" * 95)
    for r in results:
        print(f"{r['label']:>18} {r['x']:>7.2f} {r['power_forward_frac']:>8.3f} "
              f"{r['info_x_forward_frac']:>9.3f} {r['info_x_peak_theta_deg']:>7.1f} "
              f"{r['eta_x_fwd_NA0.8']:>9.3f} "
              f"{r['info_z_forward_frac']:>9.3f} {r['info_z_peak_theta_deg']:>7.1f} "
              f"{r['eta_z_fwd_NA0.8']:>9.3f}")

    with open("information_pattern_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nwrote information_pattern_results.json")
