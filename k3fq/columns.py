"""
Columns of the k3_fq_isog table (must agree with schemas/k3_fq_isog.md;
tests/test_schema.py checks this), and which stage produces each.
"""

ISOG_COLUMNS = [
    ("label", "text"),
    ("q", "integer"),
    ("p", "integer"),
    ("r", "smallint"),
    ("sign", "smallint"),
    ("poly", "bigint[]"),
    ("poly_str", "text"),
    ("real_poly", "bigint[]"),
    ("angles", "double precision[]"),
    ("angle_rank", "smallint"),
    ("angle_corank", "smallint"),
    ("slopes", "text[]"),
    ("height", "smallint"),
    ("has_k3_newton_polygon", "boolean"),
    ("newton_elevation", "smallint"),
    ("newton_coelevation", "smallint"),
    ("is_ordinary", "boolean"),
    ("is_supersingular", "boolean"),
    ("cyclotomic_factors", "integer[]"),
    ("picard_rank", "smallint"),
    ("geometric_picard_rank", "smallint"),
    ("geometric_extension_degree", "integer"),
    ("extension_degrees", "integer[]"),
    ("picard_ranks", "smallint[]"),
    ("artin_tate_values", "numeric[]"),
    ("artin_tate_value", "numeric"),
    ("ns_disc_class", "bigint"),
    ("ej_condition", "boolean"),
    ("transcendental_poly", "bigint[]"),
    ("transcendental_factors", "jsonb"),
    ("transcendental_irreducible", "boolean"),
    ("number_fields", "text[]"),
    ("number_field_coeffs", "jsonb"),
    ("galois_groups", "text[]"),
    ("k3_counts", "numeric[]"),
    ("k3_count", "bigint"),
    ("cubic_counts", "numeric[]"),
    ("cubic_count", "bigint"),
    ("gm_counts", "numeric[]"),
    ("gm_count", "bigint"),
    ("verra_counts", "numeric[]"),
    ("verra_count", "bigint"),
    ("has_nonnegative_counts", "boolean"),
    ("has_count_growth", "boolean"),
    ("is_pointless", "boolean"),
    ("is_primitive", "boolean"),
    ("primitive_models", "text[]"),
    ("geometric_class", "text"),
    ("twist_count", "integer"),
    ("twists", "jsonb"),
    ("max_twist_degree", "integer"),
    ("kummer_av_classes", "text[]"),
    ("realizable_k3", "smallint"),
    ("realizable_nck3", "smallint"),
    ("realizable_cubic", "smallint"),
    ("realizable_nccubic", "smallint"),
    ("realizable_gm", "smallint"),
    ("realizable_ncgm", "smallint"),
    ("realizable_verra", "smallint"),
    ("realizable_ncverra", "smallint"),
    ("realizable_kummer", "smallint"),
    ("polarization_degrees", "integer[]"),
    ("polarization_degrees_complete", "boolean"),
    ("num_lattices", "integer"),
    ("min_ample_degree", "integer"),
    ("min_quasi_ample_degree", "integer"),
]

ISOG_TYPES = dict(ISOG_COLUMNS)

# Columns produced by the later stages; everything else comes from the
# basic stage (K3Poly.compute).
FIELDS_COLUMNS = ["number_fields", "number_field_coeffs", "galois_groups"]
TWIST_COLUMNS = ["geometric_class", "twist_count", "twists", "max_twist_degree"]
BASECHANGE_COLUMNS = ["is_primitive", "primitive_models"]
KUMMER_COLUMNS = ["kummer_av_classes", "realizable_kummer"]
# The census stage overrides realizable_k3 and sets these (null for other classes)
POLARIZATION_COLUMNS = ["polarization_degrees", "polarization_degrees_complete"]
CENSUS_COLUMNS = ["realizable_k3"] + POLARIZATION_COLUMNS
# Not computed by this repository yet (exported as null)
FUTURE_COLUMNS = ["num_lattices", "min_ample_degree", "min_quasi_ample_degree"]

BASIC_COLUMNS = [c for c, _ in ISOG_COLUMNS
                 if c not in TWIST_COLUMNS + BASECHANGE_COLUMNS + KUMMER_COLUMNS + POLARIZATION_COLUMNS + FUTURE_COLUMNS]

# Internal columns carried between stages but not exported
GEOMETRIC_KEY = ("geometric_key", "text")
