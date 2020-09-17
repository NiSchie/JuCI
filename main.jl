using Lints
using LinearAlgebra
using TensorOperations
using Test

leri4 = true

@lints begin

  #Read in molecule and set up basis sets
  println("Reading in the molecule")
  mol = Lints.Molecule("h2o.xyz")
  println("Done reading molecule")
  println(" ")
  println("Constructing basis")
  bas = Lints.BasisSet("CC-PVDZ",mol)
  bas_df = Lints.BasisSet("CC-PVDZ-RIFIT",mol)
  println("Successfully constructed basis")

  #generate AO integtals and AO-DF integrals (P|Q), (P|mn), (mn|kl)
  PQ  = Lints.make_ERI2(bas_df)
  Pmn = Lints.make_ERI3(bas,bas_df)
  if leri4 == true
    println("Generating eri4\n")
    mnkl = Lints.make_ERI4(bas)
  end #if leri4

  println("Size of eri2:",size(PQ))
  println("Size of eri3:",size(Pmn))
  println("Size of eri4:",size(mnkl))

  #build inverse  (P|Q)^{-1/2}
  PQh = PQ^(-1/2)

  @tensor Bmn[Q,m,n] := Pmn[P,m,n] * PQh[P,Q]
  @tensor eri2[m,n,k,l] := B[Q,m,n]*B[Q,k,l]


end #lints
