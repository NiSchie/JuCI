include("physconsts.jl")
include("misc.jl")
include("scf.jl")
include("settings.jl")

using Lints
using LinearAlgebra
using TensorOperations
using LoopVectorization
using Test
using Printf
using Plots

using .PhysConsts
using .Misc
using .Settings


  #setting standard values if no config file is present
  sett = deepcopy(Settings.sett_init)
  display(sett)
  @printf("\n")

  lDFdebug = false
  Ethr = 1.0 * 10^(-parse(Float64,sett["econv"]))
  Dthr = 1.0 * 10^(-parse(Float64,sett["denconv"]))

  @printf("\n")
  @printf("SCF Econv:   %.1e\n",Ethr)
  @printf("SCF Denconv: %.1e\n",Ethr)
  @printf("\n")

  
  struct molecule
    nat::Int
    nam::String
    coords::Matrix{Float64}
    atchrg::Vector{Int}
  end
  
  AOint = Dict()
  
  #Read in molecule and set up basis sets
  
  @printf("Reading in the molecule\n")
  if isfile(sett["molfile"])
    lines = readlines(sett["molfile"])
  else
    throw("Molecule file not found")
  end
  
  nat = parse(Int,lines[1])
  
  pos   = zeros(Float64,3,nat)
  charg = Array{Int}(undef,nat)

  coord = [0.0, 0.0, 0.0]
  for i = 1:nat
    splitline = split(lines[i+2]," ")
    charg[i] = get(PhysConsts.atlist,splitline[1],0)
    popat!(splitline,1)
    filter!(Misc.isempty,splitline)
    coord .= parse.(Float64,splitline)
    pos[:,i] = collect(coord) * PhysConsts.angstrom_to_bohr
  end
  @printf("Coordinates in bohr:\n")
  display(pos)
  
  nocc = Int(sum(charg)/2)
  mol = molecule(nat,lines[2],pos,charg)
  
  Vnuc = 0.0
  for i in 1:nat, j in 1:(i-1)
    global Vnuc
    x = mol.coords[1,i] - mol.coords[1,j]
    y = mol.coords[2,i] - mol.coords[2,j]
    z = mol.coords[3,i] - mol.coords[3,j]
    Vnuc += (mol.atchrg[i]*mol.atchrg[j])/(sqrt(x*x + y*y + z*z))
  end
  @printf("\n\nVnuc: %.5f\n\n",Vnuc)


@lints begin

  using Printf

  t_coord = Array{Array{Float64,1}}
  t_coord = [mol.coords[:,i]*PhysConsts.bohr_to_angstrom for i in 1:size(mol.coords,2)]

  lintmol = Lints.Molecule(mol.atchrg,t_coord)

  timingstring=@elapsed bas    = Lints.BasisSet("cc-pVDZ",lintmol)
  @printf("Time needed to construct basis:         %.4f s\n",timingstring)
  
  #build S, T and V in AO Basis (for core guess)
  S = Lints.make_S(bas)
  T = Lints.make_T(bas)
  V = Lints.make_V(bas)

  nao = size(S,1)
  nmo = size(S,1)
  nvir = nmo - nocc

  timingstring=@elapsed bas_df = Lints.BasisSet("cc-pVDZ-JKFIT",lintmol)
  @printf("Time needed to construct DF basis:      %.4f s\n",timingstring)

  #generate AO integtals and AO-DF integrals (P|Q), (P|mn), (mn|kl)
  timingstring=@elapsed PQ  = Lints.make_ERI2(bas_df)
  @printf("Time needed to construct ERI2:          %.4f s\n",timingstring)
  timingstring=@elapsed Pmn = Lints.make_ERI3(bas,bas_df)
  @printf("Time needed to construct ERI3:          %.4f s\n",timingstring)

  if get(sett,"DF","false") == "false"
    timingstring=@elapsed mnkl = Lints.make_ERI4(bas)
    @printf("Time needed to construct ERI4:          %.4f s\n",timingstring)
  end #if leri4

  naux = size(Pmn,1)

  #build inverse  (P|Q)^{-1/2}
  timingstring=@elapsed PQh = PQ^(-1/2)
  @printf("Time needed to construct (P|Q)^{-1/2}:  %.4f s\n",timingstring)
  timingstring=@elapsed @tensor Bmn[Q,m,n] := Pmn[P,m,n] * PQh[P,Q]
  @printf("Time needed to construct (P|mn):        %.4f s\n",timingstring)
  if lDFdebug == true
    timingstring=@elapsed @tensor eri4[m,n,k,l] := Bmn[Q,m,n]*Bmn[Q,k,l]
    @printf("Time needed to construct (mn|kl) from (P|mn):  %.4f s\n",timingstring)
  end
  #
  #Fill AOintegrals Dict with Pointers to the ao integrals
  push!(AOint,"S" => S)
  push!(AOint,"T" => T)
  push!(AOint,"V" => V)
  push!(AOint,"B" => PQ)
  if (sett["DF"] == "false")
    push!(AOint,"ERI4" => mnkl)
  end

end #lints


  @printf("\nnMO: %d\n",nmo)
  @printf("nocc: %d\n",nocc)
  @printf("nvir: %d\n",nvir)
  @printf("naux: %d\n",naux)

  #build inverse of overlap
  Sh = AOint["S"]^(-1/2)
  H  = AOint["T"] .+ AOint["V"]
  push!(AOint,"H" => H)
  F  = deepcopy(H) 
  Ftinit = Sh*F*transpose(Sh) 

  #diagonalize Ft
  eps,Ctinit = eigen(Hermitian(Ftinit))

  #transform Ct with Sh to get MO-coefficients
  Cinit = Sh*Ctinit

  Cocc = Cinit[:,1:nocc]
  @tensor Dold[m,n] := 2.0 * Cocc[m,p] * Cocc[n,p] 

  Eold = 0.5*sum(Dold .* ( H .+ F ))

  scf(Cinit,Bmn,Sh)

  dE = 0.0
  Drms = 0.0
  #diis = false
  #damp = 0.0
  converged = false

  F   = zeros(nao,nao)
  Jti = zeros(naux)
  J   = zeros(nao,nao)
  K1  = zeros(naux,nao,nocc)
  K   = zeros(nao,nao)
  C   = deepcopy(Cinit)

  maxit = 100
  #RHF LOOP
  hfenfile=open("hf_energy","a+")
  Ftit = zeros(nao,nao)
  SCFEN = Array{Float64,1}()
  @printf("\n")
  @printf("HF It.  |  HF-Energy  |abs(dE)   |abs(dD)\n")
  @printf("-----------------------------------------\n")
  @printf(" GUESS  |  %.5f | %.2e | %.2e\n",Eold+Vnuc,dE,Drms)
  for ite = 1:maxit
    global Dold,Eold,dE,Drms,Vnuc,Ftit,C
    global eps,Cocc,F,converged,F,Eold,SCFEN


    #new density matrix
    Cocc = C[:,1:nocc]
    @tensor Dao[m,n] := 2.0 * Cocc[m,p] * Cocc[n,p]
    dD = Dao .- Dold
    Drms = sqrt(sum(dD.^2))
    Dold = deepcopy(Dao)

    #Build the Fock matrix
    #@printf("Building the new Fock matrix...")
    tstring=@elapsed @tensor begin
      Jti[P]    = Dao[k,l] * Bmn[P,k,l]
      J[m,n]    = Bmn[P,m,n] * Jti[P]
      K1[Q,m,c] = Cocc[r,c] * Bmn[Q,m,r]
      K[k,l]    = K1[Q,k,p] * K1[Q,l,p] 
      F[m,n]    = H[m,n] + J[m,n] - K[m,n]
    end
    Ftit = Sh*F*transpose(Sh)
    
    #Get new orb energies and coeff
    eps,Ct = eigen(Hermitian(real.(Ftit)))
    C = Sh * Ct

    E = 0.5*sum(Dao .* ( H .+ F ))
    dE = E - Eold
    Eold = E
    push!(SCFEN,E+Vnuc)

    @printf("  %3d   |  %.5f | %.2e | %.2e\n",ite,E+Vnuc,abs(dE),Drms)
    s = @sprintf("  %3d    %.10f  %.15f  %.15f\n",ite,E+Vnuc,abs(dE),Drms)
    write(hfenfile,s)

    if (abs(dE) < Ethr) & (Drms < Dthr) & (ite > 5)
      converged = true
      break
    end

  end #ite < maxit

  if(false)
    scfenplot = plot(SCFEN)
    png(scfenplot,"scfenergy.png")
  end 

  close(hfenfile)
