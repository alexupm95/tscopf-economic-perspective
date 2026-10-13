"""
Printed values: rounding and LaTeX formatting shared by the table scripts.

Money is printed with two decimals, every money cell is the rounded product of the printed factors,
and sums and prices quoted in the text are built from printed values.
"""

from decimal import ROUND_HALF_UP, Decimal


def q(x):
    """Round a Decimal to cents, half up."""
    return x.quantize(Decimal("0.01"), ROUND_HALF_UP)


def r1(v):
    """Round a float to the 0.1 table resolution, half up, as a Decimal."""
    return Decimal(repr(float(v))).quantize(Decimal("0.1"), ROUND_HALF_UP)


def tex_rows(path):
    """Data rows (after the first \\midrule) of a booktabs table, as lists of cell strings."""
    body = path.read_text(encoding="utf-8").split(r"\midrule", 1)[1]
    rows = []
    for line in body.split(r"\\"):
        cells = [c.strip() for c in line.replace(r"\midrule", "").replace(r"\bottomrule", "").split("&")]
        if len(cells) > 1:
            rows.append(cells)
    return rows


def tex_num(s):
    """'$-$1{,}160.10' -> Decimal('-1160.10')."""
    return Decimal(s.replace("$-$", "-").replace("$+$", "+").replace("{,}", "").strip())


def tex_money(x, signed=False):
    """Decimal -> LaTeX money string with two decimals, {,} thousands and $-$ / $+$ signs."""
    s = f"{abs(x):,.2f}".replace(",", "{,}")
    if x < 0:
        return f"$-${s}"
    return f"$+${s}" if signed and x > 0 else s


def tex_mw(x):
    """Decimal -> MW string without trailing zeros (187.5, 135, 122.9)."""
    return format(x.normalize(), "f").replace(",", "{,}")
