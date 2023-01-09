from math import sqrt
import random
import string
import numpy as np
from fractions import Fraction


def _random_weight(semiring, **kwargs):
    from rayuela.base.semiring import (
        Count,
        Integer,
        Real,
        Rational,
        Tropical,
        Boolean,
        MaxPlus,
        String,
    )

    if semiring is String:
        str_len = int(random.random() * 8 + 1)
        return semiring(
            "".join(random.choice(string.ascii_lowercase) for _ in range(str_len))
        )

    elif semiring is Boolean:
        return semiring(True)

    elif semiring is Real:
        tol = 1e-3
        s = kwargs.get("divide_by", 6)
        random_weight = round(random.random() / s, 3)
        while random_weight < sqrt(tol):
            random_weight = round(random.random() / s, 3)
        return semiring(random_weight)

    elif semiring is Rational:
        return semiring(Fraction(f"{random.randint(1, 1)}/{random.randint(10, 15)}"))

    elif semiring is Tropical:
        return semiring(random.randint(0, 50))

    elif semiring is Integer:
        return semiring(random.randint(1, 10))

    elif semiring is MaxPlus:
        return semiring(random.randint(-10, -1))

    elif semiring is Count:
        return semiring(1.0)


def random_weight_negative(semiring):
    from rayuela.base.semiring import Real, Rational, Tropical, Boolean, MaxPlus, String

    if semiring is Tropical:
        return semiring(random.randint(-50, 50))
    else:
        raise AssertionError("Unsupported Semiring")

