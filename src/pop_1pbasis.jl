import ACE: B1pComponent
import ACE: evaluate, evaluate_d, evaluate_ed, GetVal

export Pop1PBasis


@doc raw"""
`Pop1PBasis():` Generates a one-particle basis component similar to Scal1pBasis, but it has always length 1 such that P(x) = x. 

TODO: This is really a multiplier rather than a basis and should be implemented as such. For now we keep this for the sake of compatibility. 
"""
function Pop1PBasis(xsym = :pop, isym::Symbol = :P, label = "Pop")
  spec = [ NamedTuple{(isym,)}((1,)), ]
  degrees = [ 0, ]
  Bin = PopX()
  B = B1pComponent(Bin, GetVal{xsym}(), spec, degrees, "Pop")
end

struct PopX end

evaluate(::PopX, x) = [x,]

evaluate_d(::PopX, x::T) where {T} = [ one(T), ]

evaluate_ed(::PopX, x::T) where {T} = [x,], [ one(T), ]

write_dict(V::PopX) = Dict( "__id__" => "ACE_PopX")

read_dict(::Val{:ACE_PopX}, D::Dict) = Pop1PBasis()
