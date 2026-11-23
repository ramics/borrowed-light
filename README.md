# borrowed-light

Measures how much a correlate-of-protection model loses when it is borrowed across a gap:
another study, platform, population, or a much smaller cohort.

```
fit        m = model(source_train)
baseline   b = evaluate(m, source_heldout)   # what the model is worth at home
transfer   t = evaluate(m, target)           # what it is worth across the gap
penalty    = b - t                           # per metric, with a bootstrap interval
```

## Install

```
pip install -e ".[test]"
pytest
```

## Run

```
python -m experiments.e2_within_study experiments/configs/e2_hvtn505_age.yaml \
    --data /path/outside/repo
```

Each run writes `penalties.csv`, `bootstrap_draws.csv`, `report.md` and a manifest.
Participant data stays outside the checkout, and a run reports aggregate numbers only.

## More

- [`SPEC.md`](SPEC.md): the design
- [`docs/literature-review.md`](docs/literature-review.md): the background
- [`docs/data-access.md`](docs/data-access.md): which datasets are public and how to request the rest
- `penalty/` holds the estimator, `sim/` the simulator, `experiments/` the runs and their configs
