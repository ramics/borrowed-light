"""Pool published per-study correlate coefficients and measure between-study
transfer on the attenuation scale (docs/power-and-framing.md, rank 3).

Input: a CSV with one row per study (or per study x stratum), columns
    family      group of comparable rows (same marker, same unit), e.g. covid_id50
    study       label
    log_hr      log hazard (or odds) ratio per marker unit (protective < 0)
    se          its standard error
    + any descriptive columns (platform, region, n_cases, source, ...)

For each family:
  1. Random-effects pooling: DerSimonian-Laird and REML, Q-profile CI for tau,
     Knapp-Hartung CI for the mean, and a grid posterior for (mu, tau) with a
     half-normal prior on tau (useful when k is 3-6).
  2. Leave-one-study-out transfer, the summary-level version of the project's
     operation: source = pool of the other studies, target = the held-out
     study. Attenuation = sign(beta_S) * (beta_S - beta_T), the same sign as
     PR 10's primary metric (positive = the correlate is weaker in the target),
     with SE sqrt(se_S^2 + se_T^2) where se_S includes tau.
  3. Null check, the summary-level stand-in for the random-split check: draw
     every study from N(mu_fixed, se_i^2), i.e. a world where studies differ
     only by their own sampling error, and recompute tau and the LOO
     attenuations. The observed values are compared with that null.

Usage: python spikes/summary_pool.py data/summary_estimates.csv --out DIR
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats


def dl(y, s):
    w = 1 / s**2
    fixed = np.sum(w * y) / w.sum()
    Q = float(np.sum(w * (y - fixed) ** 2))
    k = len(y)
    c = w.sum() - np.sum(w**2) / w.sum()
    return max(0.0, (Q - (k - 1)) / c), fixed, Q


def reml_tau2(y, s):
    def nll(t2):
        v = s**2 + t2
        w = 1 / v
        mu = np.sum(w * y) / w.sum()
        return 0.5 * (np.sum(np.log(v)) + np.log(w.sum()) + np.sum(w * (y - mu) ** 2))
    if len(y) < 3:
        return dl(y, s)[0]
    r = optimize.minimize_scalar(nll, bounds=(0, 25), method="bounded")
    return float(max(0.0, r.x))


def q_profile_ci(y, s, alpha=0.05):
    k = len(y)
    def q_at(t2):
        w = 1 / (s**2 + t2)
        m = np.sum(w * y) / w.sum()
        return np.sum(w * (y - m) ** 2)
    def solve(target):
        if q_at(0.0) < target:
            return 0.0
        lo, hi = 0.0, 100.0
        for _ in range(200):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if q_at(mid) > target else (lo, mid)
        return (lo + hi) / 2
    return (np.sqrt(solve(stats.chi2.ppf(1 - alpha / 2, k - 1))),
            np.sqrt(solve(stats.chi2.ppf(alpha / 2, k - 1))))


def pooled(y, s, t2):
    w = 1 / (s**2 + t2)
    mu = np.sum(w * y) / w.sum()
    return mu, np.sqrt(1 / w.sum())


def knapp_hartung(y, s, t2, alpha=0.05):
    k = len(y)
    w = 1 / (s**2 + t2)
    mu = np.sum(w * y) / w.sum()
    q = max(1.0, np.sum(w * (y - mu) ** 2) / (k - 1))
    se = np.sqrt(q / w.sum())
    t = stats.t.ppf(1 - alpha / 2, k - 1)
    return mu - t * se, mu + t * se


def grid_posterior(y, s, tau_prior_sd=0.5, n=600):
    """Flat prior on mu, half-normal(tau_prior_sd) on tau; mu integrated out
    analytically, tau on a grid. Returns tau median, 95% interval, and the
    posterior predictive SD for a new study (sqrt(tau^2 + var(mu)))."""
    taus = np.linspace(0, max(3.0, 6 * tau_prior_sd), n)
    logp = []
    for t in taus:
        v = s**2 + t**2
        w = 1 / v
        mu = np.sum(w * y) / w.sum()
        ll = -0.5 * (np.sum(np.log(v)) + np.log(w.sum()) + np.sum(w * (y - mu) ** 2))
        logp.append(ll + stats.halfnorm.logpdf(t, scale=tau_prior_sd))
    p = np.exp(np.array(logp) - max(logp))
    p /= p.sum()
    cdf = np.cumsum(p)
    q = lambda a: float(taus[np.searchsorted(cdf, a)])
    return {"tau_median": q(0.5), "tau_95": [q(0.025), q(0.975)],
            "tau_prior_sd": tau_prior_sd}


def loo(y, s, labels):
    rows = []
    for i in range(len(y)):
        m = np.ones(len(y), bool)
        m[i] = False
        if m.sum() >= 3:
            t2 = reml_tau2(y[m], s[m])
        else:
            t2 = 0.0
        mu, se_mu = pooled(y[m], s[m], t2)
        # Predicting a new study: source uncertainty includes tau.
        se_src = np.sqrt(se_mu**2 + t2)
        att = np.sign(mu) * (mu - y[i])
        se_att = np.sqrt(se_src**2 + s[i] ** 2)
        rows.append({"held_out": labels[i], "beta_source": mu, "se_source_pred": se_src,
                     "beta_target": y[i], "se_target": s[i],
                     "attenuation": att, "se_attenuation": se_att,
                     "z_fixed_world": np.sign(mu) * (mu - y[i]) / np.sqrt(se_mu**2 + s[i] ** 2),
                     "relative_attenuation": att / abs(mu) if mu != 0 else np.nan})
    return pd.DataFrame(rows)


def null_check(y, s, labels, reps, rng):
    """World with no between-study heterogeneity: y_i ~ N(mu_fixed, se_i^2)."""
    _, fixed, Q_obs = dl(y, s)
    tau_obs = np.sqrt(reml_tau2(y, s))
    loo_obs = loo(y, s, labels)
    max_abs_obs = float(np.max(np.abs(loo_obs["z_fixed_world"])))
    taus, Qs, maxz = [], [], []
    for _ in range(reps):
        yy = rng.normal(fixed, s)
        taus.append(np.sqrt(reml_tau2(yy, s)))
        Qs.append(dl(yy, s)[2])
        lz = loo(yy, s, labels)["z_fixed_world"]
        maxz.append(np.max(np.abs(lz)))
    taus, Qs, maxz = map(np.array, (taus, Qs, maxz))
    return {"tau_obs": tau_obs, "tau_null_95": float(np.quantile(taus, 0.95)),
            "p_tau": float(np.mean(taus >= tau_obs - 1e-12)),
            "Q_obs": Q_obs, "p_Q_sim": float(np.mean(Qs >= Q_obs)),
            "p_Q_chi2": float(stats.chi2.sf(Q_obs, len(y) - 1)),
            "max_abs_loo_z_obs": max_abs_obs,
            "p_max_loo_z": float(np.mean(maxz >= max_abs_obs))}


def analyse_family(df, reps, rng):
    y = df["log_hr"].to_numpy(float)
    s = df["se"].to_numpy(float)
    labels = df["study"].tolist()
    k = len(y)
    t2_dl, fixed, Q = dl(y, s)
    t2_reml = reml_tau2(y, s)
    mu, se_mu = pooled(y, s, t2_reml)
    out = {"k": k, "fixed": fixed, "Q": Q, "I2": max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0,
           "tau_dl": np.sqrt(t2_dl), "tau_reml": np.sqrt(t2_reml),
           "tau_qprofile_95": list(q_profile_ci(y, s)) if k >= 2 else None,
           "mu_reml": mu, "se_mu": se_mu,
           "mu_kh_95": list(knapp_hartung(y, s, t2_reml)) if k >= 2 else None}
    t = stats.t.ppf(0.975, max(k - 2, 1))
    out["prediction_95"] = [mu - t * np.sqrt(t2_reml + se_mu**2), mu + t * np.sqrt(t2_reml + se_mu**2)]
    out["bayes"] = {f"prior_{p}": grid_posterior(y, s, p) for p in (0.25, 0.5, 1.0)}
    out["null_check"] = null_check(y, s, labels, reps, rng) if k >= 3 else None
    return out, loo(y, s, labels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261122)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(a.csv, comment="#")
    if "log_hr" not in d:
        d["log_hr"] = np.log(d["hr"])
        d["se"] = (np.log(d["hi"]) - np.log(d["lo"])) / (2 * stats.norm.ppf(0.975))
    d = d[d["include"].fillna(1).astype(int) == 1] if "include" in d else d
    res, loos = {}, []
    for fam, g in d.groupby("family", sort=False):
        r, l = analyse_family(g.reset_index(drop=True), a.reps, rng)
        res[fam] = r
        loos.append(l.assign(family=fam))
        print(fam, json.dumps({k: v for k, v in r.items() if k != "bayes"}, default=float)[:900])
    (out / "pooled.json").write_text(json.dumps(res, indent=2, default=float))
    pd.concat(loos).to_csv(out / "loo_transfer.csv", index=False)


if __name__ == "__main__":
    main()
