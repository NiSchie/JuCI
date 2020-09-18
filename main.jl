include("physconsts.jl")

using Lints
using LinearAlgebra
using TensorOperations
using LoopVectorization
using Test

using .PhysConsts

leri4 = false
lDFdebug = false
molfile="mol.xyz"

Ethr = 1.0E-6
Dthr = 1.0E-6

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
println("Coordinates")
display(pos)

nocc = Int(sum(charg)/2)
mol = molecule(nat,lines[2],pos,charg)

Vnuc = 0.0
for i in 1:nat, j in 1:(i-1)
  global Vnuc
  x = mol.coords[i,1] - mol.coords[j,1]
  y = mol.coords[i,2] - mol.coords[j,2]
  z = mol.coords[i,3] - mol.coords[j,3]
  Vnuc += (mol.atchrg[i]*mol.atchrg[j])/(sqrt(x*x + y*y + z*z)*PhysConsts.angstrom_to_bohr)
end
println("Vnuc: ",Vnuc)

@lints begin

  t_coord = Array{Array{Float64,1}}
  t_coord = [mol.coords[i,:] for i in 1:size(mol.coords,1)]

  lintmol = Lints.Molecule(mol.atchrg,t_coord)
  println("Done reading molecule")

  println(" ")
  println("Constructing basis")
  bas = Lints.BasisSet("cc-pVDZ",lintmol)
  bas_df = Lints.BasisSet("cc-pVDZ-RIFIT",lintmol)
  println("Successfully constructed basis")

  #generate AO integtals and AO-DF integrals (P|Q), (P|mn), (mn|kl)
  PQ  = Lints.make_ERI2(bas_df)
  Pmn = Lints.make_ERI3(bas,bas_df)
  if leri4 == true
    println("Generating eri4\n")
    mnkl = Lints.make_ERI4(bas)
    println("done.\n")
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

end #lints

  nao = size(S,1)
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
  F  .= H 
  Ft = Sh*F*transpose(Sh) 
  Ftinit = deepcopy(Ft)

  #diagonalize Ft
  eps,Ctinit = eigen(Hermitian(Ft))

  #transform Ct with Sh to get MO-coefficients
  Cinit = Sh*Ctinit

  Cocc = Cinit[:,1:nocc]
  @tensor Dao[m,n] := Cocc[m,p] * Cocc[n,p] 
  Dold = deepcopy(Dao)

  Daoold   = deepcopy(Dao)
  epsold = deepcopy(eps)
  Cold   = deepcopy(Cinit)

  ite = 1
  E = 0.0
  dE = 1.0
  Eold = 1.0E20
  Drms = 1.0
  #diis = false
  #damp = 0.0
  converged = false

  maxit = 300
  #RHF LOOP
  while ite < maxit
    println("\nHF Iteration : $ite")
    global Ft,Dold,Eold,dE,Drms,Vnuc

    #Get new orb energies and coeff
    eps,Ct = eigen(Hermitian(real.(Ft)))
    C = Sh * Ct

    #new density matrix
    Cocc = C[:,1:nocc]
    @tensor Dao[m,n] := Cocc[m,p] * Cocc[n,p]
    if ite > 1
      dD = Dao - Dold
      Drms = sqrt(sum(dD.^2))
      println("Drms: ",Drms)
    end
    Dold = deepcopy(Dao)

    #Build the Fock matrix
    #println("Building the new Fock matrix...")
    @tensor F[m,n] = H[m,n] + Dao[k,l]*( (Bmn[P,m,n]*Bmn[P,k,l]) - 0.5*(Bmn[Q,m,l]*Bmn[Q,k,n]) )
    #@tensor begin
    #  #F[m,n] = H[m,n]
    #  G[m,n] := Dao[k,l] * (Bmn[Q,m,n]*Bmn[Q,k,l])
    #  G[m,n] = G[m,n] - 0.5 * (Dao[k,l] * (Bmn[Q,m,l]*Bmn[Q,k,n]))
    #  F[m,n] = H[m,n] + G[m,n]
    #end
    
    #println("done.\n")
    Ft = Sh*F*transpose(Sh)


    #E = 0.5*sum(Dao .* (T+V .+ F))
    @tensor E = Dao[m,n] * ( T[n,m]+V[n,m] + F[n,m] )
    if ite > 1
      dE = E - Eold
    end
    Eold = E
    println("HF Energy: ",E+Vnuc," Difference: ",dE, "Electronic Energy: ",E)

    if (abs(dE) < Ethr) & (Drms < Dthr) & (ite > 5)
      converged = true
      break
    end

    global ite +=1
  end #ite < maxit

  println("\nTM reference energy: -76.01678545283")

