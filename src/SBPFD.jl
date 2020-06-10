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

end # module
