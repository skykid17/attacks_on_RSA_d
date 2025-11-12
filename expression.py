from poly import Poly
from poly_mat import PolyMat


# Simple variables for ease of use.
x = Poly(1, (1, 0, 0))
y = Poly(1, (0, 1, 0))
z = Poly(1, (0, 0, 1))


def remove_x_factor(eqn: Poly) -> Poly:
    assert eqn.highest_y_power() == 0
    assert eqn.highest_z_power() == 0
    result = Poly()
    lowest_x: int | None = None
    for x in range(eqn.highest_x_power()+1):
        coeff = eqn[x, 0, 0]
        if coeff == 0:
            continue
        if lowest_x is None:
            lowest_x = x
        result[x-lowest_x, 0, 0] = coeff
    return result


def resultant(a: Poly, b: Poly) -> Poly:
    self_y = a.highest_y_power()
    other_y = b.highest_y_power()
    if self_y + other_y == 0:
        return Poly(0)
    sylvester_matrix = PolyMat(self_y + other_y)
    for i in range(other_y):
        for y_pow in range(self_y+1):
            sylvester_matrix[i, i+y_pow] = a.univariate_y(self_y - y_pow)
    for i in range(self_y):
        for y_pow in range(other_y+1):
            sylvester_matrix[i+other_y, i +
                             y_pow] = b.univariate_y(other_y - y_pow)
    return sylvester_matrix.det()
