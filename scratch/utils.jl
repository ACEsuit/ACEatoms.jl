# --------------------- Radial Basis Construction

Rn_basis(; r0 = 2.5,
            trans = PolyTransform(2, r0),
            maxn = 10,
            rcut = 5.0,
            rin = 0.5 * r0,
            pcut = 2,
            pin = 0 ) =
   transformed_jacobi(maxn, trans, rcut, rin; pcut=pcut, pin=pin)


# RnYlm_basis(; species = :X, D = SparsePSHDegree(), maxdeg = 8, kwargs...)  =
#       RnYlm1pBasis(
#             radial_basis(maxn = get_maxn(D, maxdeg, species), kwargs...);
#             species = species, D = D )


# ----------------
