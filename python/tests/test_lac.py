from toypk.lac import get_linear_attenuation_coefficient


def test_H():
    lac = get_linear_attenuation_coefficient(0.001, {"H": 1})
    assert lac == 7.217
