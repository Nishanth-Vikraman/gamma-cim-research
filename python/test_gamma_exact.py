import random
from gamma_exact_model import exact_value, normalize, encode_int, multiply_ints, dot_product_ints


def test_local_identity_exact():
    assert exact_value({0:3}) == exact_value({-1:1, 1:1})


def test_normalize_preserves_value():
    for case in [{0:3},{0:-3},{0:26},{-2:7,0:-11,3:19},{-5:12,-1:-8,4:15}]:
        before = exact_value(case)
        after, st = normalize(case)
        assert exact_value(after) == before
        assert all(abs(v) <= 2 for v in after.values())


def test_full_int8_encoding():
    for n in range(-128,128):
        d = encode_int(n)
        assert exact_value(d) == (n,0)
        assert all(abs(v) <= 2 for v in d.values())


def test_small_exhaustive_multiply():
    for a in range(-16,17):
        for b in range(-16,17):
            result, st, raw, da, db = multiply_ints(a,b)
            assert exact_value(result) == (a*b,0)


def test_random_int8_multiply():
    rng = random.Random(0)
    for _ in range(1000):
        a = rng.randint(-128,127)
        b = rng.randint(-128,127)
        result, st, raw, da, db = multiply_ints(a,b)
        assert exact_value(result) == (a*b,0)


def test_random_dot_products():
    rng = random.Random(1)
    for interval in [1,4,16,64]:
        for _ in range(10):
            xs = [rng.randint(-128,127) for _ in range(64)]
            ws = [rng.randint(-128,127) for _ in range(64)]
            result, st = dot_product_ints(xs,ws,normalize_every=interval)
            assert exact_value(result) == (sum(x*w for x,w in zip(xs,ws)),0)
