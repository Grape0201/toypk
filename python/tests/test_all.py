import toypk


def test_simple():
    toypk.run(
        {
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
                "form": "test",
                "data": [[2.0], [2.0]],
                "conversion_factor": [1.0, 1.0],
            },
        }
    )
