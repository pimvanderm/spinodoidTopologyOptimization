# spinodoid-surrogate

A PyTorch port of the MATLAB neural-network surrogate model for **inverse-designed
spinodoid metamaterials**, based on Kumar et al. (2020), *npj Computational Materials*
([doi:10.1038/s41524-020-0341-6](https://doi.org/10.1038/s41524-020-0341-6)).

It replicates the two-network design from the original `trainNeuralNetworks.m`:

| Model | Direction | Architecture | Loss |
|-------|-----------|--------------|------|
| **Forward** (`ForwardNet`) | structure → stiffness (4 → 9) | FC `128,128,64,64,32,32` + ReLU, z-score input layer | MSE |
| **Inverse** (`InverseNet`) | stiffness → structure (9 → 4) | FC `6×100` + ReLU | `recon + λ·param` through a frozen forward net |

* **Inputs** `[density, θ₁, θ₂, θ₃]` — min-max scaled to `[0,1]` with fixed bounds `[0,0,0,0]…[1,90,90,90]`.
* **Outputs** — the 9 stiffness components `C11,C12,C13,C22,C23,C33,C44,C55,C66`, min-max scaled using training-set bounds.
* The inverse model is trained against a **frozen** forward model so that
  reconstructed stiffness matches the queried stiffness; the parameter-matching
  term `λ` anneals from `0.5` to `0` after epoch 40.

## Install

```bash
cd spinodoid_surrogate
pip install -e ".[test]"
```

## Train

```bash
spinodoid-train \
    --train /path/to/spinodoid-train.csv \
    --test  /path/to/spinodoid-test.csv \
    --save-dir ./weights
```

Or from Python:

```python
from spinodoid_surrogate import (
    Normalizer, SpinodoidDataset, load_csv, train_forward, train_inverse,
)

theta_tr, stiff_tr = load_csv("spinodoid-train.csv")
norm = Normalizer.fit(theta_tr, stiff_tr)
train_ds = SpinodoidDataset(theta_tr, stiff_tr, norm)

forward_model, _ = train_forward(train_ds, norm)
inverse_model, _ = train_inverse(train_ds, forward_model)
```

## Test

```bash
pytest                       # synthetic-data tests only
SPINODOID_DATA_DIR=/path/to/csvs pytest   # also exercises the real CSVs
```

## Package layout

```
src/spinodoid_surrogate/
  config.py   # ForwardConfig / InverseConfig (hyperparameters)
  data.py     # load_csv, Normalizer, SpinodoidDataset
  models.py   # ForwardNet, InverseNet
  train.py    # train_forward, train_inverse
  metrics.py  # r2_score (per-component, matches MATLAB)
  cli.py      # `spinodoid-train` entry point
tests/
  test_surrogate.py
```

## License

[MIT](LICENSE).
