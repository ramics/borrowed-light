# borrowed-light

Estimates how much a correlate-of-protection model loses when it is borrowed
across a gap: another study, platform, population or a much smaller cohort.
See [SPEC.md](SPEC.md) for the design and
[docs/literature-review.md](docs/literature-review.md) for the background.

```
pip install -e ".[test]"
pytest
```

| Module | Role |
| --- | --- |
| `penalty/schema.py` | tidy-frame validation, sampling weights, two-phase rows, marker mapping |
| `penalty/metrics.py` | weighted AUC, Brier, calibration-in-the-large, calibration slope |
| `penalty/models.py` | weighted ridge logistic (CV within source), best single marker |
| `penalty/transfer.py` | fit / baseline / transfer / penalty |
| `sim/generate.py` | dial-driven simulator and case-cohort sampler |

No real participant data is ever committed here (AC7). Real data is analysed
where it lives; `.gitignore` blocks the common data formats as a backstop.
