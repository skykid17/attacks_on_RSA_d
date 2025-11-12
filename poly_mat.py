from poly import Poly
from typing import overload


class PolyMat:
    def __init__(self, dim: int) -> None:
        assert dim > 0
        self._dim = dim
        self._val = [Poly()] * (dim*dim)

    def det(self) -> Poly:
        if self._dim == 1:
            return self._val[0]
        if self._dim == 2:
            return self[0, 0] * self[1, 1] - self[1, 0]*self[0, 1]
        # Bareiss algorithm.
        mat = PolyMat(self._dim)  # Make a copy.
        for i in range(self._dim):
            mat[i] = self[i]
        prev_pivot = Poly(1)
        sign = 1
        for k in range(0, self._dim-1):
            pivot_row = k
            while mat[pivot_row, k].is_zero():
                pivot_row += 1
                if pivot_row == self._dim:
                    # No non-zero pivot, the determinant is zero.
                    return Poly(0)
            if pivot_row != k:
                # Need to do the swap.
                # Switching rows force us to switch signs.
                mat[pivot_row], mat[k] = mat[k], mat[pivot_row]
                sign = -sign
            pivot = mat[k, k]
            for i in range(k+1, self._dim):
                for j in range(k+1, self._dim):
                    num = mat[k, k]*mat[i, j] - mat[i, k]*mat[k, j]
                    mat[i, j] = num / prev_pivot
                mat[i, k] = Poly(0)
            prev_pivot = pivot
        return Poly(sign) * mat[self._dim-1, self._dim-1]

    @overload
    def __getitem__(self, key: int) -> list[Poly]:
        pass

    @overload
    def __getitem__(self, key: tuple[int, int]) -> Poly:
        pass

    def __getitem__(self, key: int | tuple[int, int]):
        if isinstance(key, tuple):
            assert key[0] < self._dim
            assert key[1] < self._dim
            return self._val[key[0]*self._dim + key[1]]
        assert key < self._dim
        return self._val[key*self._dim:(key+1)*self._dim]

    @overload
    def __setitem__(self, key: int, val: list[Poly]) -> None:
        pass

    @overload
    def __setitem__(self, key: tuple[int, int], val: Poly) -> None:
        pass

    def __setitem__(self,
                    key: int | tuple[int, int],
                    val: list[Poly] | Poly) -> None:
        if isinstance(key, tuple) and isinstance(val, Poly):
            assert key[0] >= 0 and key[0] < self._dim
            assert key[1] >= 0 and key[1] < self._dim
            idx = key[0] * self._dim + key[1]
            self._val[idx] = val
            return
        if isinstance(key, int) and isinstance(val, list):
            assert key >= 0 and key < self._dim
            assert len(val) == self._dim
            for i, v in enumerate(val):
                self[key, i] = v
            return
        raise AssertionError("invalid types")

    def __str__(self) -> str:
        result = "["
        for i in range(self._dim):
            result += f"\n\t{self[i]}"
        result += "\n]"
        return result
