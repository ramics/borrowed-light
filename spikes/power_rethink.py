"""Spike: is the transfer penalty measured on the right scale?

Three parts, all cheap (minutes on a laptop, no EC2):

A. Simulation. For one standardised marker, compare the precision and power of
   three estimands of "the correlate got weaker in the target":
     - AUC penalty (source AUC - target AUC), the current headline for ranking;
     - calibration slope of the source model in the target (beta_T / beta_S);
     - the marker-coefficient difference, delta = beta_T - beta_S (log-OR per
       SD), estimated by a likelihood fit with a marker x stratum interaction.
   Both sides share a 10% attack rate. "Attenuation" means beta_T = beta_S / 2.

B. Real data (HK H3N2, PNG P. vivax, HVTN 505). For each split, the
   interaction estimate of delta with a robust SE, plus a random-split
   (permutation) null for both attenuation and the AUC penalty, so each observed
   penalty is compared with the penalty you get from splitting the same people
   at random into groups of the same sizes.

C. Random-effects pooling of delta across independent datasets
   (DerSimonian-Laird), to show what pooling buys and what it cannot.

Usage:
    python spikes/power_rethink.py --hk DIR --png XLSX --h505 RDA [--png-loader PY]
        [--out DIR] [--reps N]

No data is written into the repository; --out should sit outside it.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


# --------------------------------------------------------------------------
# Small weighted logistic fit with sandwich (optionally cluster-robust) SEs.
# --------------------------------------------------------------------------
def logit_fit(X, y, w=None, cluster=None, iters=50):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    b = np.zeros(X.shape[1])
    b[0] = np.log((np.average(y, weights=w) + 1e-9) / (1 - np.average(y, weights=w) + 1e-9))
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        g = X.T @ (w * (y - p))
        H = (X * (w * p * (1 - p))[:, None]).T @ X
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            return None
        b += step
        if np.max(np.abs(step)) < 1e-9:
            break
    eta = np.clip(X @ b, -30, 30)
    p = 1 / (1 + np.exp(-eta))
    H = (X * (w * p * (1 - p))[:, None]).T @ X
    try:
        Hinv = np.linalg.inv(H)
    except np.linalg.LinAlgError:
        return None
    U = X * (w * (y - p))[:, None]
    if cluster is not None:
        U = pd.DataFrame(U).groupby(np.asarray(cluster)).sum().to_numpy()
    V = Hinv @ (U.T @ U) @ Hinv
    return b, V


def auc(score, y, w=None):
    """Weighted AUC (Mann-Whitney), ties counted half."""
    score, y = np.asarray(score, float), np.asarray(y, int)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    pos, neg = y == 1, y == 0
    if pos.sum() == 0 or neg.sum() == 0:
        return np.nan
    order = np.argsort(score)
    s, yy, ww = score[order], y[order], w[order]
    # cumulative negative weight strictly below each score, ties half
    uniq, inv = np.unique(s, return_inverse=True)
    negw = np.bincount(inv, weights=ww * (yy == 0), minlength=len(uniq))
    posw = np.bincount(inv, weights=ww * (yy == 1), minlength=len(uniq))
    below = np.concatenate([[0], np.cumsum(negw)[:-1]])
    num = np.sum(posw * (below + 0.5 * negw))
    return num / (posw.sum() * negw.sum())


# --------------------------------------------------------------------------
# Part A: simulation
# --------------------------------------------------------------------------
def beta_for_auc(target_auc, rng, n=400_000, prev=0.10):
    """Find the per-SD log-OR giving a population AUC for a N(0,1) marker."""
    x = rng.standard_normal(n)
    lo, hi = 0.0, 4.0
    for _ in range(40):
        b = (lo + hi) / 2
        a = intercept_for(b, prev, x)
        y = rng.random(n) < 1 / (1 + np.exp(-(a - b * x)))
        if auc(-x, y) < target_auc:
            lo = b
        else:
            hi = b
    return (lo + hi) / 2


def intercept_for(b, prev, x):
    lo, hi = -15.0, 5.0
    for _ in range(60):
        a = (lo + hi) / 2
        if np.mean(1 / (1 + np.exp(-(a - b * x)))) < prev:
            lo = a
        else:
            hi = a
    return (lo + hi) / 2


def sim_side(n, a, b, rng):
    x = rng.standard_normal(n)
    y = (rng.random(n) < 1 / (1 + np.exp(-(a - b * x)))).astype(int)
    return x, y


def one_rep(nS, nT, aS, bS, aT, bT, rng):
    xs, ys = sim_side(nS, aS, bS, rng)
    xt, yt = sim_side(nT, aT, bT, rng)
    fs = logit_fit(np.c_[np.ones(nS), xs], ys)
    ft = logit_fit(np.c_[np.ones(nT), xt], yt)
    if fs is None or ft is None:
        return None
    bs, vs = fs
    bt, vt = ft
    d = bt[1] - bs[1]
    se_d = np.sqrt(vs[1, 1] + vt[1, 1])
    # AUC: a single-marker model ranks by the marker, so source AUC needs no
    # CV for the point estimate; use the sign fitted in the source.
    sgn = np.sign(bs[1]) or 1.0
    auc_pen = auc(sgn * xs, ys) - auc(sgn * xt, yt)
    # Calibration slope of the source model in the target: refit y_T on LP_S.
    lp = bs[0] + bs[1] * xt
    fc = logit_fit(np.c_[np.ones(nT), lp], yt)
    slope = fc[0][1] if fc is not None else np.nan
    return d, se_d, auc_pen, slope


def part_a(rng, reps):
    prev = 0.10
    rows = []
    for home_auc in (0.60, 0.70, 0.80):
        b = beta_for_auc(home_auc, rng)
        xx = rng.standard_normal(200_000)
        a_full = intercept_for(b, prev, xx)
        a_half = intercept_for(b / 2, prev, xx)
        for events in (25, 50, 100, 200, 400):
            n = int(round(events / prev))
            null = [one_rep(n, n, a_full, b, a_full, b, rng) for _ in range(reps)]
            alt = [one_rep(n, n, a_full, b, a_half, b / 2, rng) for _ in range(reps)]
            null = np.array([r for r in null if r is not None])
            alt = np.array([r for r in alt if r is not None])
            sd_d, sd_auc = null[:, 0].std(), null[:, 2].std()
            slope_q = np.nanpercentile(null[:, 3], [2.5, 97.5])
            # Power at two-sided 5%: AUC penalty uses its null SD (as a
            # bootstrap would); delta uses its Wald SE from each fit.
            # Sign convention: attenuation makes |beta| smaller. The marker
            # is protective (beta enters as -b), so delta > 0 is attenuation.
            pow_d = np.mean(np.abs(alt[:, 0] / alt[:, 1]) > 1.96)
            pow_auc = np.mean(np.abs(alt[:, 2] / sd_auc) > 1.96)
            rows.append({
                "home_auc": home_auc, "beta_per_sd": round(b, 3), "events_per_side": events,
                "n_per_side": n,
                "sd_auc_penalty_null": sd_auc, "sd_delta_null": sd_d,
                "mean_wald_se_delta": null[:, 1].mean(),
                "cal_slope_null_2.5": slope_q[0], "cal_slope_null_97.5": slope_q[1],
                "true_auc_penalty_attenuated": float(np.mean(alt[:, 2])),
                "power_auc_penalty": pow_auc, "power_delta_wald": pow_d,
            })
            print(f"A: AUC {home_auc} E={events}: sd(AUCpen)={sd_auc:.3f} "
                  f"sd(delta)={sd_d:.3f} slope95=[{slope_q[0]:.2f},{slope_q[1]:.2f}] "
                  f"power AUC={pow_auc:.2f} delta={pow_d:.2f}", flush=True)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Part B: real-data splits
# --------------------------------------------------------------------------
def split_stats(M, y, s, w, cl):
    """Attenuation of the marker coefficient from source to target.

    M holds one or more candidate markers (columns). If more than one, the
    marker is chosen in the source alone (largest |z| in a single-marker fit),
    so a permutation null that calls this function repeats the selection.
    Each marker is standardised to its pooled SD. Returns
    (attenuation, se, auc_penalty, beta_source, se_source, chosen column), where
    attenuation = sign(beta_S) * (beta_S - beta_T): positive means the marker's
    log-OR per SD is weaker in the target, the same sign as the other penalties.
    """
    M = np.atleast_2d(np.asarray(M, float).T).T
    src = s == 0
    Z = (M - np.average(M, axis=0, weights=w)) / np.sqrt(np.average((M - np.average(M, axis=0, weights=w))**2, axis=0, weights=w))
    j = 0
    if Z.shape[1] > 1:
        best = -1.0
        for k in range(Z.shape[1]):
            f = logit_fit(np.c_[np.ones(src.sum()), Z[src, k]], y[src], w[src])
            if f is not None and abs(f[0][1] / np.sqrt(f[1][1, 1])) > best:
                best, j = abs(f[0][1] / np.sqrt(f[1][1, 1])), k
    z = Z[:, j]
    X = np.c_[np.ones(len(z)), z, s, z * s]
    f = logit_fit(X, y, w, cl)
    nan = (np.nan,) * 5 + (j,)
    if f is None or not np.all(np.isfinite(f[1])) or f[1][3, 3] <= 0:
        return nan
    b, V = f
    sg = np.sign(b[1]) or 1.0
    a_s = auc(sg * z[src], y[src], w[src])
    a_t = auc(sg * z[~src], y[~src], w[~src])
    return -sg * b[3], np.sqrt(V[3, 3]), a_s - a_t, b[1], np.sqrt(V[1, 1]), j


def robust_sd(v):
    """1.4826 x MAD: a few permutations with separation blow up the plain SD."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def permute_labels(s, cl, rng):
    if cl is None:
        return rng.permutation(s)
    return rng.permutation(s)  # row-level; HK overlap noted in the doc


def part_b(datasets, rng, reps):
    rows = []
    for name, (df, marker, src_label, cluster_col) in datasets.items():
        markers = [marker] if isinstance(marker, str) else list(marker)
        x = df[markers].to_numpy(float)
        y = df["endpoint"].to_numpy(int)
        w = df["sampling_weight"].to_numpy(float)
        s = (df["stratum"] != src_label).to_numpy(int)
        cl = df[cluster_col].to_numpy() if cluster_col else None
        d, se, auc_pen, bS, seS, j = split_stats(x, y, s, w, cl)
        nd, na = [], []
        for _ in range(reps):
            sp = permute_labels(s, cl, rng)
            r = split_stats(x, y, sp, w, cl)
            nd.append(r[0])
            na.append(r[2])
        nd, na = np.array(nd), np.array(na)
        rows.append({
            "split": name, "marker": markers[j], "selected_from": len(markers),
            "events_source": int(y[s == 0].sum()), "events_target": int(y[s == 1].sum()),
            "beta_source_per_sd": bS, "se_beta_source": seS,
            "attenuation": d, "se_attenuation": se, "p_wald": 2 * stats.norm.sf(abs(d / se)),
            "p_perm_attenuation": float(np.mean(np.abs(nd) >= abs(d))),
            "null_sd_attenuation": robust_sd(nd),
            "auc_penalty": auc_pen,
            "null_sd_auc_penalty": robust_sd(na),
            "null_mean_auc_penalty": float(np.nanmean(na)),
            "p_perm_auc": float(np.mean(np.abs(na - np.nanmean(na)) >= abs(auc_pen - np.nanmean(na)))),
        })
        print("B:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    return pd.DataFrame(rows)


def dersimonian_laird(est, se):
    est, se = np.asarray(est), np.asarray(se)
    w = 1 / se**2
    fixed = np.sum(w * est) / w.sum()
    Q = np.sum(w * (est - fixed) ** 2)
    k = len(est)
    tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - np.sum(w**2) / w.sum()))
    wr = 1 / (se**2 + tau2)
    mu = np.sum(wr * est) / wr.sum()
    se_mu = np.sqrt(1 / wr.sum())
    # Q-profile CI for tau^2
    def q_at(t2):
        ww = 1 / (se**2 + t2)
        m = np.sum(ww * est) / ww.sum()
        return np.sum(ww * (est - m) ** 2)
    def solve(target):
        lo, hi = 0.0, 50.0
        if q_at(lo) < target:
            return 0.0
        for _ in range(100):
            mid = (lo + hi) / 2
            if q_at(mid) > target:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2
    tau2_lo = solve(stats.chi2.ppf(0.975, k - 1))
    tau2_hi = solve(stats.chi2.ppf(0.025, k - 1))
    tq = stats.t.ppf(0.975, max(k - 2, 1))
    pi = (mu - tq * np.sqrt(tau2 + se_mu**2), mu + tq * np.sqrt(tau2 + se_mu**2))
    return {"k": k, "mu": mu, "se_mu": se_mu, "Q": Q, "tau": np.sqrt(tau2),
            "tau_ci": [np.sqrt(tau2_lo), np.sqrt(tau2_hi)], "prediction_interval": list(pi)}


def load_png(loader_py, path, split):
    spec = importlib.util.spec_from_file_location("png_loaders", loader_py)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.png_pv(path, split)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hk")
    ap.add_argument("--png")
    ap.add_argument("--png-loader", help="loaders.py that defines png_pv (PR 7) if main lacks it")
    ap.add_argument("--h505")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--sim-reps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261122)
    ap.add_argument("--skip-sim", action="store_true")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    if not a.skip_sim:
        part_a(rng, a.sim_reps).to_csv(out / "sim_precision.csv", index=False)

    from experiments import loaders
    ds = {}
    if a.hk:
        for col, src, extra in (("strain", "Perth", {}),
                                ("age.y", None, {"cut": 40, "labels": ["under40", "40plus"]})):
            df = loaders.hk_flu(a.hk, {"column": col, **extra})[0]
            ds[f"hk_{col}"] = (df, "marker_log2_hai_pre", src or "under40", "cluster_id")
    if a.png:
        png = getattr(loaders, "png_pv", None)
        for col, src, extra in (("village_group", "Ilaita", {}),
                                ("age_yr", None, {"cut": 2.5, "labels": ["young", "old"]})):
            sp = {"column": col, **extra}
            df = png(a.png, sp)[0] if png else load_png(a.png_loader, a.png, sp)
            ds[f"png_{col}"] = (df, [c for c in df.columns if c.startswith("marker_")], src or "young", None)
    if a.h505:
        for col, cut, labels in (("age", 30, ["age_18_29", "age_30_50"]),
                                 ("BMI", 25, ["bmi_lt25", "bmi_ge25"])):
            df = loaders.hvtn505(a.h505, {"column": col, "cut": cut, "labels": labels})[0]
            df = df[df["phase2"] == 1]
            ds[f"h505_{col}"] = (df, [c for c in df.columns if c.startswith("marker_")], labels[0], None)
    res = part_b(ds, rng, a.reps)
    res.to_csv(out / "real_splits.csv", index=False)

    primary = res[res["split"].isin(["hk_strain", "png_village_group", "h505_age"])]
    meta = dersimonian_laird(primary["attenuation"], primary["se_attenuation"])
    (out / "meta.json").write_text(json.dumps(meta, indent=2, default=float))
    print("C:", meta)


if __name__ == "__main__":
    main()
