#!/usr/bin/env sage
"""Verify the explicit affine local models used in Paper 4.

The calculations are exact over QQ. They reconstruct the Hilbert bases and
toric ideals of the dP6 and dP7 cones, verify Altmann's dP6 versal equations,
identify all three reduced smoothing-component families with Cayley-cone
toric varieties, and check representative smooth and discriminant fibres by
the Jacobian criterion through Singular.

Run:  sage local_models.sage
"""

from sage.interfaces.singular import singular


DP6_VERTICES = (
    (0, 0),
    (1, 0),
    (2, 1),
    (2, 2),
    (1, 2),
    (0, 1),
)

DP7_VERTICES = (
    (0, 0),
    (1, 1),
    (0, 2),
    (-1, 2),
    (-1, 1),
)


# Hilbert-basis order used in the research note. The final element in each
# tuple is the Gorenstein-height character T.
DP6_HILBERT_BASIS = (
    (0, 1, 0),       # z1
    (1, 0, 0),       # z2
    (1, -1, 1),      # z3
    (0, -1, 2),      # z4
    (-1, 0, 2),      # z5
    (-1, 1, 1),      # z6
    (0, 0, 1),       # T
)

DP7_HILBERT_BASIS = (
    (1, 0, 1),       # y1
    (1, 1, 0),       # y2
    (0, -1, 2),      # y3
    (-1, -1, 2),     # y4
    (-1, 1, 0),      # y5
    (0, 1, 0),       # y6
    (-1, 0, 1),      # y7
    (0, 0, 1),       # T
)


# Columns for the three Cayley cones, ordered by the named total-space
# coordinates used below. Pure height characters come last.
DP6_LINE_CAYLEY = (
    (0, 1, 0, 0),       # z1
    (1, 0, 0, 0),       # z2
    (1, -1, 0, 1),      # z3
    (0, -1, 1, 1),      # z4
    (-1, 0, 1, 1),      # z5
    (-1, 1, 1, 0),      # z6
    (0, 0, 1, 0),       # A: odd-edge triangle
    (0, 0, 0, 1),       # B: even-edge triangle
)

DP6_PLANE_CAYLEY = (
    (0, 1, 0, 0, 0),       # z1
    (1, 0, 0, 0, 0),       # z2
    (1, -1, 0, 1, 0),      # z3
    (0, -1, 0, 1, 1),      # z4
    (-1, 0, 1, 0, 1),      # z5
    (-1, 1, 1, 0, 0),      # z6
    (0, 0, 1, 0, 0),       # A0: first Cayley height
    (0, 0, 0, 1, 0),       # A1: second Cayley height
    (0, 0, 0, 0, 1),       # A2: third Cayley height
)

DP7_LINE_CAYLEY = (
    (1, 0, 0, 1),       # y1
    (1, 1, 0, 0),       # y2
    (0, -1, 1, 1),      # y3
    (-1, -1, 2, 0),     # y4
    (-1, 1, 0, 0),      # y5
    (0, 1, 0, 0),       # y6
    (-1, 0, 1, 0),      # y7
    (0, 0, 1, 0),       # A: triangle
    (0, 0, 0, 1),       # B: segment
)


def column_matrix(columns):
    """Return the integral matrix whose columns are the displayed tuples."""
    return matrix(
        ZZ,
        len(columns[0]),
        len(columns),
        lambda row, column: columns[column][row],
    )


def mapped_toric_ideal(columns, target_ring, images):
    """Map the toric ideal of a column configuration into target_ring."""
    source_ideal = ToricIdeal(column_matrix(columns))
    source_ring = source_ideal.ring()
    morphism = source_ring.hom(images, target_ring)
    return target_ring.ideal(
        [morphism(generator) for generator in source_ideal.gens()]
    )


def cone_hilbert_basis(vertices):
    """Compute the Hilbert basis of the dual cone over a lattice polygon."""
    cone = Cone(rays=[tuple(vertex) + (1,) for vertex in vertices])
    return frozenset(tuple(element) for element in cone.dual().Hilbert_basis())


singular.lib("sing.lib")


def is_smooth(ideal):
    """Test affine smoothness using Singular's Jacobian singular-locus ideal."""
    singular_locus = singular.slocus(singular(ideal)).sage()
    return singular_locus.is_one()


def has_singular_tangent_at(ideal, point):
    """Verify that a rational point has tangent dimension above dimension."""
    ring = ideal.ring()
    assert len(point) == ring.ngens()
    assert all(equation(*point) == 0 for equation in ideal.gens())
    jacobian = matrix(
        QQ,
        [
            [equation.derivative(variable)(*point) for variable in ring.gens()]
            for equation in ideal.gens()
        ],
    )
    tangent_dimension = ring.ngens() - jacobian.rank()
    return tangent_dimension > ideal.dimension()


# ---------------------------------------------------------------------------
# dP6: central cone and full three-parameter Altmann family
# ---------------------------------------------------------------------------

assert cone_hilbert_basis(DP6_VERTICES) == frozenset(DP6_HILBERT_BASIS)

R6 = PolynomialRing(QQ, names=("z1", "z2", "z3", "z4", "z5", "z6", "T"), order="dp")
z1, z2, z3, z4, z5, z6, T = R6.gens()

DP6_CENTRAL_EQUATIONS = (
    z1 * T - z6 * z2,
    z2 * T - z1 * z3,
    z3 * T - z2 * z4,
    z4 * T - z3 * z5,
    z5 * T - z4 * z6,
    z6 * T - z5 * z1,
    T**2 - z1 * z4,
    T**2 - z2 * z5,
    T**2 - z3 * z6,
)
DP6_CENTRAL_IDEAL = R6.ideal(DP6_CENTRAL_EQUATIONS)
assert DP6_CENTRAL_IDEAL == mapped_toric_ideal(
    DP6_HILBERT_BASIS,
    R6,
    R6.gens(),
)
assert DP6_CENTRAL_IDEAL.dimension() == 3
assert has_singular_tangent_at(DP6_CENTRAL_IDEAL, (0,) * 7)


def dp6_equations(coordinates, s1_value, s2_value, s3_value):
    """Construct Altmann's nine dP6 equations in one polynomial ring."""
    ez1, ez2, ez3, ez4, ez5, ez6, eT = coordinates
    t1 = eT
    t2 = eT - s2_value - s3_value
    t3 = eT - s1_value
    t4 = eT - s3_value
    t5 = eT - s2_value
    t6 = eT - s1_value - s3_value
    equations = (
        ez1 * t1 - ez6 * ez2,
        ez2 * t2 - ez1 * ez3,
        ez3 * t3 - ez2 * ez4,
        ez4 * t4 - ez3 * ez5,
        ez5 * t5 - ez4 * ez6,
        ez6 * t6 - ez5 * ez1,
        t5 * t6 - ez1 * ez4,
        t3 * t4 - ez2 * ez5,
        t1 * t2 - ez3 * ez6,
    )
    return equations


def dp6_fibre(s1_value, s2_value, s3_value):
    """Return Altmann's dP6 fibre at the displayed base point."""
    return R6.ideal(
        dp6_equations(R6.gens(), s1_value, s2_value, s3_value)
    )


assert dp6_fibre(0, 0, 0) == DP6_CENTRAL_IDEAL

# The line branch is s1=s2=0. Its Cayley heights are A=T and B=T-s3.
R6L = PolynomialRing(
    QQ,
    names=("z1", "z2", "z3", "z4", "z5", "z6", "T", "s3"),
    order="dp",
)
lz1, lz2, lz3, lz4, lz5, lz6, lT, ls3 = R6L.gens()


def dp6_line_equations():
    return dp6_equations(R6L.gens()[:7], 0, 0, ls3)


DP6_LINE_IDEAL = R6L.ideal(dp6_line_equations())
DP6_LINE_TORIC_IDEAL = mapped_toric_ideal(
    DP6_LINE_CAYLEY,
    R6L,
    (lz1, lz2, lz3, lz4, lz5, lz6, lT, lT - ls3),
)
assert DP6_LINE_IDEAL == DP6_LINE_TORIC_IDEAL
assert DP6_LINE_IDEAL.dimension() == 4

# The plane branch is s3=0. Its Cayley heights are
# A0=T, A1=T-s2, and A2=T-s1 in the fixed Hilbert-coordinate ordering.
R6P = PolynomialRing(
    QQ,
    names=("z1", "z2", "z3", "z4", "z5", "z6", "T", "s1", "s2"),
    order="dp",
)
pz1, pz2, pz3, pz4, pz5, pz6, pT, ps1, ps2 = R6P.gens()


def dp6_plane_equations():
    return dp6_equations(R6P.gens()[:7], ps1, ps2, 0)


DP6_PLANE_IDEAL = R6P.ideal(dp6_plane_equations())
DP6_PLANE_TORIC_IDEAL = mapped_toric_ideal(
    DP6_PLANE_CAYLEY,
    R6P,
    (pz1, pz2, pz3, pz4, pz5, pz6, pT, pT - ps2, pT - ps1),
)
assert DP6_PLANE_IDEAL == DP6_PLANE_TORIC_IDEAL
assert DP6_PLANE_IDEAL.dimension() == 5

# Exact sample fibres locate the reduced discriminant on the plane branch.
print("checking smooth dP6 plane sample...", flush=True)
assert is_smooth(dp6_fibre(1, 2, 0))
assert has_singular_tangent_at(dp6_fibre(0, 1, 0), (0, 0, 0, 0, 0, 0, 0))
assert has_singular_tangent_at(dp6_fibre(1, 0, 0), (0, 0, 0, 0, 0, 0, 0))
assert has_singular_tangent_at(dp6_fibre(1, 1, 0), (0, 0, 0, 0, 0, 0, 1))
print("checking smooth dP6 line sample...", flush=True)
assert is_smooth(dp6_fibre(0, 0, 1))


# ---------------------------------------------------------------------------
# dP7: central cone and its unique reduced smoothing component
# ---------------------------------------------------------------------------

assert cone_hilbert_basis(DP7_VERTICES) == frozenset(DP7_HILBERT_BASIS)

R7 = PolynomialRing(
    QQ,
    names=("y1", "y2", "y3", "y4", "y5", "y6", "y7", "T"),
    order="dp",
)
y1, y2, y3, y4, y5, y6, y7, uT = R7.gens()
DP7_CENTRAL_IDEAL = mapped_toric_ideal(DP7_HILBERT_BASIS, R7, R7.gens())
assert DP7_CENTRAL_IDEAL.dimension() == 3
assert has_singular_tangent_at(DP7_CENTRAL_IDEAL, (0,) * 8)


def dp7_reduced_equations(coordinates, beta_value):
    """Construct the dP7 reduced-component equations in one ring."""
    ey1, ey2, ey3, ey4, ey5, ey6, ey7, eT = coordinates
    A = eT
    B = eT - beta_value
    equations = (
        ey2 * ey3 - ey1 * A,
        ey2 * B - ey1 * ey6,
        ey2 * ey5 - ey6**2,
        ey1 * ey5 - B * ey6,
        ey3 * ey5 - B * ey7,
        ey2 * ey4 - A**2,
        ey1 * ey4 - ey3 * A,
        B * ey4 - ey3 * ey7,
        ey5 * ey4 - ey7**2,
        ey2 * ey7 - ey6 * A,
        ey1 * ey7 - B * A,
        ey3 * ey6 - B * A,
        ey4 * ey6 - ey7 * A,
        ey7 * ey6 - ey5 * A,
    )
    return equations


def dp7_reduced_fibre(beta_value):
    """Return the dP7 Cayley fibre on the reduced branch alpha=0."""
    return R7.ideal(dp7_reduced_equations(R7.gens(), beta_value))


assert dp7_reduced_fibre(0) == DP7_CENTRAL_IDEAL
print("checking smooth dP7 line sample...", flush=True)
assert is_smooth(dp7_reduced_fibre(1))

R7L = PolynomialRing(
    QQ,
    names=("y1", "y2", "y3", "y4", "y5", "y6", "y7", "T", "beta"),
    order="dp",
)
ly1, ly2, ly3, ly4, ly5, ly6, ly7, luT, lbeta = R7L.gens()
DP7_LINE_IDEAL = R7L.ideal(
    dp7_reduced_equations(R7L.gens()[:8], lbeta)
)
DP7_LINE_TORIC_IDEAL = mapped_toric_ideal(
    DP7_LINE_CAYLEY,
    R7L,
    (ly1, ly2, ly3, ly4, ly5, ly6, ly7, luT, luT - lbeta),
)
assert DP7_LINE_IDEAL == DP7_LINE_TORIC_IDEAL
assert DP7_LINE_IDEAL.dimension() == 4


print("local_models: all exact assertions passed")
print("dP6 central cone: 7 Hilbert generators, 9 quadratic equations")
print("dP6 line branch: two-triangle Cayley family; sample s3=1 smooth")
print(
    "dP6 plane branch: three-segment Cayley family; sample (s1,s2)=(1,2) "
    "smooth"
)
print("dP6 plane discriminant samples: s1=0, s2=0, and s1=s2 are singular")
print("dP7 central cone: 8 Hilbert generators")
print("dP7 reduced branch: triangle-plus-segment Cayley family; beta=1 smooth")
