include("physconsts.jl")

using Lints
using LinearAlgebra
using TensorOperations
using Test

using .PhysConsts

leri4 = false
lDFdebug = false
molfile="mol.xyz"

struct molecule
  nat::Int
  nam::String
  coords::Matrix{Float64}
  atchrg::Vector{Int}
end

#Read in molecule and set up basis sets

println("Reading in the molecule")
if isfile(molfile)
  lines = readlines(molfile)
else
end

nat = parse(Int,lines[1])

pos   = zeros(Float64,nat,3)
charg = Array{Int}(undef,nat)

coord = [0.0, 0.0, 0.0]
for i = 1:nat
  splitline = split(lines[i+2]," ")
  charg[i] = get(PhysConsts.atlist,splitline[1],0)
  popat!(splitline,1)
  coord .= parse.(Float64,splitline)
  pos[i,:] = collect(coord)
end
nocc = sum(charg)
mol = molecule(nat,lines[2],pos,charg)


@lints begin

  t_coord = Array{Array{Float64,1}}
  t_coord = [mol.coords[i,:] for i in 1:size(mol.coords,1)]

  lintmol = Lints.Molecule(mol.atchrg,t_coord)
  println("Done reading molecule")

  println(" ")
  println("Constructing basis")
  bas = Lints.BasisSet("CC-PVDZ",lintmol)
  bas_df = Lints.BasisSet("CC-PVDZ-RIFIT",lintmol)
  println("Successfully constructed basis")

  #generate AO integtals and AO-DF integrals (P|Q), (P|mn), (mn|kl)
  PQ  = Lints.make_ERI2(bas_df)
  Pmn = Lints.make_ERI3(bas,bas_df)
  if leri4 == true
    println("Generating eri4\n")
    mnkl = Lints.make_ERI4(bas)
  end #if leri4

  println("\n Dimensions of (P|Q):",size(PQ))
  println("Dimensions of (P|mn):",size(Pmn))
  if leri4 == true
    println("Size of (mn|kl):",size(mnkl))
  end #if leri4
  println("")

  #build inverse  (P|Q)^{-1/2}
  println("Building (P|Q)^{-1/2}...")
  PQh = PQ^(-1/2)
  println("done.\n")

  println("Building 3idx matrices (P|mn)...")
  @tensor Bmn[Q,m,n] := Pmn[P,m,n] * PQh[P,Q]
  println("done.\n")
  if lDFdebug == true
    println("Building 4idx matrices (mn|kl) from (P|mn)...")
    @tensor eri2[m,n,k,l] := Bmn[Q,m,n]*Bmn[Q,k,l]
    println("done.\n")
  end


  #CORE GUESS
  #build S, T, V, and I in AO Basis (for core guess)
  println("Building Overlap matrix S...")
  S = Lints.make_S(bas)
  println("built S")
  println("Building kinetic matrix T...")
  T = Lints.make_T(bas)
  println("built T")
  println("Building Vne matrix V...")
  V = Lints.make_V(bas)
  println("built V\n")
  if leri4 == true
    println("Building 4 index integral matrix...")
    I = Lints.make_ERI4(bas)
    println("done.\n")
  end #if leri4

  nmo = size(S,1)
  nvir = nmo - nocc
  println("\n nMO:",nmo)
  println("nocc:",nocc)
  println("nvir:",nvir)
  println("")

  #build inverse of overlap
  println("Building S^{-1/2}...")
  Sh = S^(-1/2)
  println("done.")
  H  = T + V
  F  = Array{Float64,2}(undef, nmo,nmo)
  F .= H
  Ft = Sh*F*transpose(Sh) 

  #diagonalize Ft
  eps,Ct = eigen(Hermitian(Ft))

  #transform Ct with Sh to get MO-coefficients
  C = Sh*Ct


end #lints
