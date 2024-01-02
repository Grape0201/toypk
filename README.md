# toypk
__this is my personal project for learning Rust.__

toypk is a toy project for point-kernel photon shielding analyses

Python bindings are ready, powered by [PyO3](https://github.com/PyO3/pyo3)

```python
import toypk

toypk.run({
    "primitives": [
        {"shape": "Sphere", "data": [0.0, 0.0, 0.0, 1.0]},
        {"shape": "Sphere", "data": [0.0, 0.0, 0.0, 2.0]},
    ],
    "cells": [
        {"name": "inside", "material_index": 0, "csg_operations": [-1]},
        {"name": "middle", "material_index": 1, "csg_operations": [1, -2]},
        {"name": "outside", "material_index": 2, "csg_operations": [2]},
    ],
    "materials": [
        {"name": "fuel", "linear_attenuation_coefficient": [0.1]},
        {"name": "iron", "linear_attenuation_coefficient": [0.3]},
        {"name": "air", "linear_attenuation_coefficient": [0.01]},
    ],
    "source": {
        "factor_by_group": [1.0],
        "source_type": "test",
        "fdata": [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5],
        "udata": [10, 10, 10],
    },
    "tally_points": [[2.0, 0.0, 0.0]],
    "buildup": {
        "conversion_factor": [1.0, 1.0],
        "form": "test",
        "data": [[2.0], [2.0]],
    },
})
```

# TODOs
- [x] read build up factors from python dict
- [ ] new source types
- [x] standard linear coefficient libraries
- [ ] standard buildup factor library
