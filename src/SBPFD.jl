module SBPFD

using StaticArrays

export FirstDerivativeDiagonalSBP

"""
    FirstDerivativeDiagonalSBP{PI, PB}

Storage for coefficients for diagonal first derivative SBP operators of interior
order `PI` and boundary order `PB`.
"""
struct FirstDerivativeDiagonalSBP{PI, PB, T, B, C, BC, D}
  "boundary quadrature matrix"
  HB::SArray{Tuple{B}, T, 1, B}
  "boundary derivative matrix"
  DB::SArray{Tuple{B, C}, T, 2, BC}
  "interior derivative matrix"
  DI::SArray{Tuple{D}, T, 1, D}

  function FirstDerivativeDiagonalSBP{PI, PB}(
    HB::SVector{B, T},
    DB::SMatrix{B, C, T},
    DI::SVector{D, T},
  ) where {T, B, C, D, PI, PB}
    new{PI, PB, T, B, C, B * C, D}(HB, DB, DI)
  end
end

#=
@ARTICLE{Strand94,
  author = {B. Strand},
  title = {Summation by parts for finite difference approximations for $d/dx$},
  journal = {Journal of Computational Physics},
  year = {1994},
  volume = {110},
  pages = {47--67},
  number = {1},
  doi = {10.1006/jcph.1994.1005}
}
=#
FirstDerivativeDiagonalSBP{T, 2}() where {T} =
  FirstDerivativeDiagonalSBP{T, 2, 1}()
function FirstDerivativeDiagonalSBP{
  T,
  2,
  1,
}() where {T <: Union{AbstractFloat, Rational}}
  HB = SVector{1, T}([1 // 2])
  DB = SMatrix{1, 2, T}([-1 1])
  DI = SVector{3, T}([-1 // 2, 0, 1 // 2])
  FirstDerivativeDiagonalSBP{2, 1}(HB, DB, DI)
end

#=
@ARTICLE{Strand94,
  author = {B. Strand},
  title = {Summation by parts for finite difference approximations for $d/dx$},
  journal = {Journal of Computational Physics},
  year = {1994},
  volume = {110},
  pages = {47--67},
  number = {1},
  doi = {10.1006/jcph.1994.1005}
}
=#
FirstDerivativeDiagonalSBP{T, 4}() where {T} =
  FirstDerivativeDiagonalSBP{T, 4, 2}()
function FirstDerivativeDiagonalSBP{
  T,
  4,
  2,
}() where {T <: Union{AbstractFloat, Rational}}
  #! format: off
  HB = SVector{4, T}([17 // 48, 59 // 48, 43 // 48, 49 // 48])
  DB = SMatrix{4, 6, T}([ -24 // 17  59 // 34  -4 // 17 -3 // 34  0        0
                           -1 // 2    0         1 // 2   0        0        0
                            4 // 43 -59 // 86   0       59 // 86 -4 // 43  0
                            3 // 98   0       -59 // 98  0       32 // 49 -4 // 49
                           ])
  DI = SVector{5, T}([1 // 12 -2 // 3 0 2 // 3 -1 // 12])
  #! format: on
  FirstDerivativeDiagonalSBP{4, 2}(HB, DB, DI)
end

end # module
