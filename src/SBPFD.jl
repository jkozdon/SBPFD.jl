module SBPFD

using StaticArrays

export FirstDerivativeDiagonalSBP, SecondDerivativeDiagonalSBP

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

"""
    SecondDerivativeDiagonalSBP{PI, PB}

Storage for coefficients for diagonal second derivative SBP operators of
interior order `PI` and boundary order `PB`.
"""
struct SecondDerivativeDiagonalSBP{PI, PB, T, B, C, BWB, WWB, SWB, WWI, SWI}
  "boundary quadrature matrix"
  HB::SArray{Tuple{B}, T, 1, B}
  "boundary first derivative matrix"
  BS::SArray{Tuple{C}, T, 1, C}
  "Weights for boundary coefficients"
  WB::SArray{Tuple{BWB, WWB, WWB}, T, 3, SWB}
  "Weights for interior coefficients"
  WI::SArray{Tuple{WWI, WWI}, T, 2, SWI}

  function SecondDerivativeDiagonalSBP{PI, PB}(
    HB::SArray{Tuple{B}, T, 1, B},
    BS::SArray{Tuple{C}, T, 1, C},
    WB::SArray{Tuple{BWB, WWB, WWB}, T, 3, SWB},
    WI::SArray{Tuple{WWI, WWI}, T, 2, SWI},
  ) where {PI, PB, T, B, C, BWB, WWB, SWB, WWI, SWI}
    new{PI, PB, T, B, C, BWB, WWB, SWB, WWI, SWI}(HB, BS, WB, WI)
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

SecondDerivativeDiagonalSBP{T, 2}() where {T} =
  SecondDerivativeDiagonalSBP{T, 2, 1}()
function SecondDerivativeDiagonalSBP{
  T,
  2,
  1,
}() where {T <: Union{AbstractFloat, Rational}}
  HB = SVector{1, T}([1 // 2])
  BS = SVector{3, T}([3 // 2 -2 1 // 2])
  WB = SArray{Tuple{2, 1, 1}, T}(1 // 2, 1 // 2)
  WI = SArray{Tuple{3, 3}, T}(
    -1 // 2,
    -1 // 2,
    0,
    1 // 2,
    1,
    1 // 2,
    0,
    -1 // 2,
    -1 // 2,
  )
  SecondDerivativeDiagonalSBP{2, 1}(HB, BS, WB, WI)
end

SecondDerivativeDiagonalSBP{T, 4}() where {T} =
  SecondDerivativeDiagonalSBP{T, 4, 2}()
function SecondDerivativeDiagonalSBP{
  T,
  4,
  2,
}() where {T <: Union{AbstractFloat, Rational}}
  HB = SVector{4, T}([17 // 48, 59 // 48, 43 // 48, 49 // 48])
  BS = SVector{4, T}([11 // 6 -3 3 // 2 -1 // 3])
  WB = zeros(T, 8, 6, 6)
  WB[1, 1, 1] = 12 // 17
  WB[2, 1, 1] = 59 // 192
  WB[3, 1, 1] = 27010400129 // 345067064608
  WB[4, 1, 1] = 69462376031 // 2070402387648

  WB[1, 1, 2] = WB[1, 2, 1] = -59 // 68
  WB[3, 1, 2] = WB[3, 2, 1] = -6025413881 // 21126554976
  WB[4, 1, 2] = WB[4, 2, 1] = -537416663 // 7042184992

  WB[1, 1, 3] = WB[1, 3, 1] = 2 // 17
  WB[2, 1, 3] = WB[2, 3, 1] = -59 // 192
  WB[3, 1, 3] = WB[3, 3, 1] = 2083938599 // 8024815456
  WB[4, 1, 3] = WB[4, 3, 1] = 213318005 // 16049630912

  WB[1, 1, 4] = WB[1, 4, 1] = 3 // 68
  WB[3, 1, 4] = WB[3, 4, 1] = -1244724001 // 21126554976
  WB[4, 1, 4] = WB[4, 4, 1] = 752806667 // 21126554976

  WB[3, 1, 5] = WB[3, 5, 1] = 49579087 // 10149031312
  WB[4, 1, 5] = WB[4, 5, 1] = -49579087 // 10149031312

  WB[3, 1, 6] = WB[3, 6, 1] = 1 // 784
  WB[4, 1, 6] = WB[4, 6, 1] = -1 // 784

  WB[1, 2, 2] = 3481 // 3264
  WB[3, 2, 2] = 9258282831623875 // 7669235228057664
  WB[4, 2, 2] = 236024329996203 // 1278205871342944

  WB[1, 2, 3] = WB[1, 3, 2] = -59 // 408
  WB[3, 2, 3] = WB[3, 3, 2] = -29294615794607 // 29725717938208
  WB[4, 2, 3] = WB[4, 3, 2] = -2944673881023 // 29725717938208

  WB[1, 2, 4] = WB[1, 4, 2] = -59 // 1088
  WB[3, 2, 4] = WB[3, 4, 2] = 260297319232891 // 2556411742685888
  WB[4, 2, 4] = WB[4, 4, 2] = -60834186813841 // 1278205871342944

  WB[3, 2, 5] = WB[3, 5, 2] = -1328188692663 // 37594290333616
  WB[4, 2, 5] = WB[4, 5, 2] = 1328188692663 // 37594290333616

  WB[3, 2, 6] = WB[3, 6, 2] = -8673 // 2904112
  WB[4, 2, 6] = WB[4, 6, 2] = 8673 // 2904112

  WB[1, 3, 3] = 1 // 51
  WB[2, 3, 3] = 59 // 192
  WB[3, 3, 3] = 378288882302546512209 // 270764341349677687456
  WB[4, 3, 3] = 13777050223300597 // 26218083221499456
  WB[5, 3, 3] = 564461 // 13384296

  WB[1, 3, 4] = WB[1, 4, 3] = 1 // 136
  WB[3, 3, 4] = WB[3, 4, 3] = -4836340090442187227 // 5525802884687299744
  WB[4, 3, 4] = WB[4, 4, 3] = -17220493277981 // 89177153814624
  WB[5, 3, 4] = WB[5, 4, 3] = -125059 // 743572

  WB[3, 3, 5] = WB[3, 5, 3] = 1613976761032884305 // 7963657098519931984
  WB[4, 3, 5] = WB[4, 5, 3] = -10532412077335 // 42840005263888
  WB[5, 3, 5] = WB[5, 5, 3] = 564461 // 4461432

  WB[3, 3, 6] = WB[3, 6, 3] = 33235054191 // 26452850508784
  WB[4, 3, 6] = WB[4, 6, 3] = -960119 // 1280713392
  WB[5, 3, 6] = WB[5, 6, 3] = -3391 // 6692148

  WB[1, 4, 4] = 3 // 1088
  WB[3, 4, 4] = 507284006600757858213 // 475219048083107777984
  WB[4, 4, 4] = 1950062198436997 // 3834617614028832
  WB[5, 4, 4] = 1869103 // 2230716
  WB[6, 4, 4] = 1 // 24

  WB[3, 4, 5] = WB[3, 5, 4] = -4959271814984644613 // 20965546238960637264
  WB[4, 4, 5] = WB[4, 5, 4] = -15998714909649 // 37594290333616
  WB[5, 4, 5] = WB[5, 5, 4] = -375177 // 743572
  WB[6, 4, 5] = WB[6, 5, 4] = -1 // 6

  WB[3, 4, 6] = WB[3, 6, 4] = 752806667 // 539854092016
  WB[4, 4, 6] = WB[4, 6, 4] = 1063649 // 8712336
  WB[5, 4, 6] = WB[5, 6, 4] = -368395 // 2230716
  WB[6, 4, 6] = WB[6, 6, 4] = 1 // 8

  WB[3, 5, 5] = 8386761355510099813 // 128413970713633903242
  WB[4, 5, 5] = 2224717261773437 // 2763180339520776
  WB[5, 5, 5] = 280535 // 371786
  WB[6, 5, 5] = 5 // 6
  WB[7, 5, 5] = 1 // 24

  WB[3, 5, 6] = WB[3, 6, 5] = -13091810925 // 13226425254392
  WB[4, 5, 6] = WB[4, 6, 5] = -35039615 // 213452232
  WB[5, 5, 6] = WB[5, 6, 5] = -1118749 // 2230716
  WB[6, 5, 6] = WB[6, 6, 5] = -1 // 2
  WB[7, 5, 6] = WB[7, 6, 5] = -1 // 6

  WB[3, 6, 6] = 660204843 // 13226425254392
  WB[4, 6, 6] = 3290636 // 80044587
  WB[5, 6, 6] = 5580181 // 6692148
  WB[6, 6, 6] = 3 // 4
  WB[7, 6, 6] = 5 // 6
  WB[8, 6, 6] = 1 // 24

  WB = SArray{Tuple{8, 6, 6}, T}(WB)

  WI = zeros(T, 5, 5)
  WI[1, 1] = WI[5, 5] = 1 // 8
  WI[1, 2] = WI[5, 4] = -1 // 6
  WI[1, 3] = WI[5, 3] = 1 // 8

  WI[2, 1] = WI[4, 4] = -1 // 6
  WI[2, 2] = WI[4, 3] = -1 // 3
  WI[2, 3] = WI[4, 2] = -1 // 2
  WI[2, 4] = WI[4, 1] = -1 // 6

  WI[3, 1] = 1 // 24
  WI[3, 2] = 5 // 6
  WI[3, 3] = 3 // 4
  WI[3, 4] = 5 // 6
  WI[3, 5] = 1 // 24

  WI = SArray{Tuple{5, 5}, T}(WI)
  SecondDerivativeDiagonalSBP{4, 2}(HB, BS, WB, WI)
end

end # module
