from __future__ import annotations
from decimal import Context, Decimal, getcontext, setcontext
from math import atan, pi, sin, cos
from typing import Optional


class ComplexContext:
    _ctx: Context
    _in_ctx: bool
    _prev_ctx: Optional[Context]

    def __init__(self, prec: int, Emax: int):
        self._ctx = Context(prec=prec, Emax=Emax)
        self._in_ctx = False
        self._prev_ctx = None

    def __enter__(self):
        self._prev_ctx = getcontext()
        setcontext(self._ctx)
        self._in_ctx = True

    def __exit__(self, _typ, _val, _tb):
        setcontext(self._prev_ctx)
        self._prev_ctx = None
        self._in_ctx = False

    def create(self,
               real: int | float | Decimal,
               imag: int | float | Decimal = 0) -> Complex:
        assert self._in_ctx
        return Complex(self.create_decimal(real), self.create_decimal(imag))

    def create_decimal(self, abs: int | float | Decimal) -> Decimal:
        assert self._in_ctx
        return self._ctx.create_decimal(abs)

    def create_unit(self, arg: float) -> Complex:
        assert self._in_ctx
        return Complex(
            self.create_decimal(cos(arg)),
            self.create_decimal(sin(arg)),
        )


EPSILON = 1e-9
TWO_PI = Decimal(2 * pi)


class Complex:
    real: Decimal
    imag: Decimal

    def __init__(self, real: Decimal, imag: Decimal):
        self.real = real
        self.imag = imag

    def __add__(self, other: Complex) -> Complex:
        return Complex(self.real + other.real, self.imag + other.imag)

    def __sub__(self, other: Complex) -> Complex:
        return Complex(self.real - other.real, self.imag - other.imag)

    def __mul__(self, other: Complex) -> Complex:
        return Complex(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + other.real * self.imag,
        )

    def __truediv__(self, other: Complex) -> Complex:
        denom = other.real**2 + other.imag**2
        return Complex(
            (self.real*other.real + self.imag*other.imag) / denom,
            (self.imag*other.real - self.real*other.imag) / denom,
        )

    def __pow__(self, pow: float) -> Complex:
        new_abs = abs(self) ** Decimal(pow)
        new_arg = (self.arg * pow) % (2 * pi)
        return Complex(
            new_abs * Decimal(cos(new_arg)),
            new_abs * Decimal(sin(new_arg)),
        )

    def __abs__(self) -> Decimal:
        return (self.real**2 + self.imag**2).sqrt()

    def __str__(self) -> str:
        abs_val = abs(self)
        arg = self.arg
        if abs(arg) < EPSILON:
            return str(abs_val)
        return f"{abs_val}∠{arg:.2f}"

    def __repr__(self):
        return f"Complex({repr(self.real)}, {self.imag})"

    @property
    def arg(self) -> float:
        return self._calc_angle(self.real, self.imag)

    @staticmethod
    def _calc_angle(x: Decimal, y: Decimal) -> float:
        if abs(x) < EPSILON:
            # Handle (0, y).
            if y > 0:
                return pi/2
            return -pi/2
        # It's (x, y) where x != 0.
        # angle returns -pi/2 to pi/2.
        angle = atan(y/x)
        if x < 0:
            # The range should actually be pi/2 to 3pi/2.
            angle += pi
        return angle
