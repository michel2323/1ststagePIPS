using PyPlot


# Read dumped matrix
# Format is m, n integers setting the sizeof
# and m*n entries
function openMat(filename)
  fp=open(filename)
  m=read(fp, Int32)
  n=read(fp, Int32)
  A=read(fp, Float64, m, n)
  close(fp)
  println("Reading matrix from ", filename, " ", m, "x", n)
  return A
end

# Read first stage matrix and compute solution
function firststage()
  println("Reading matrix M")
  A=openMat("1ststageM.dmp")
  # show structure
  spy(A)
  show()
  # Reading RHS  
  println("Reading RHS")
  rhs=openMat("1ststageRHS.dmp")
  # Reading solution computed in PIPS  
  println("Reading solution")
  sol=openMat("1ststageSol.dmp")
  # Solve system in Julia  
  println("Solving system")
  x=\(A,rhs)
  
  # Compare both solutions  
  println("Solution from file:")
  println(norm(sol))
  println("Computed solution:")
  println(norm(x))
  
  println("Condition number: ", cond(A))
end

# Read 1st stage block of the global matrix
# This is the lower right matrix in the notes (all '0' entries)
# Iteration defines which iteration should be read.
function buildM(iter)
  child=0
  filename="globalKKT" * "$child" * "_" * "$iter" * ".dmp"
  W0=openMat(filename)
  return W0
end

# This reads Q0, B0, D0 and the diagonals, then assembles 
# the 1st stage matrix. This should do exactly the same as buildM(),
# however it doesn't. Use buildM().
# See the notes for how exactly the assembling is done
function buildM2(iter)
  # reading in block
 child=0
  filename="globalQ0" * "$child" * "_" * "$iter" * ".dmp"
  println(filename)
  Q=openMat(filename)
  filename="globalB0" * "$child" * "_" * "$iter" * ".dmp"
  B=openMat(filename)
  filename="globalDss0" * "$child" * "_" * "$iter" * ".dmp"
  Ds=openMat(filename)
  filename="globalDsx0" * "$child" * "_" * "$iter" * ".dmp"
  Dx=openMat(filename)
  filename="globalDsy0" * "$child" * "_" * "$iter" * ".dmp"
  Dy=openMat(filename)
  filename="globalDsz0" * "$child" * "_" * "$iter" * ".dmp"
  Dz=openMat(filename)
  filename="globalD0" * "$child" * "_" * "$iter" * ".dmp"
  D=openMat(filename)
  
  # computing block size
  @assert(size(Q,1)==size(B,2))
  @assert(size(Q,1)==size(D,2))
  @assert(size(B,1)==size(Dy,1))
  @assert(size(D,1)==size(Dz,1))
  n=size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)
  println(n)
  W=zeros(n,n)
  
  # move blocks to the right places
  # Read Dx
  # Dsx
  T=zeros(size(Dx,1),size(Dx,1))
  for i=1:size(Dx,1)
    T[i,i]=Dx[i,1]
  end
  
  # Compute \bar{Q}=Q*Dx
  Q=Q+T
  # start with Qs
  W[1:size(Q,1),1:size(Q,2)]=Q
  
  # Bs
  W[1+size(Q,1)+size(Ds,1):size(Q,1)+size(Ds,1)+size(B,1),1:size(B,2)]=B
  # transpose symmetric
  W[1:size(B,2),1+size(Q,1)+size(Ds,1):size(Q,1)+size(Ds,1)+size(B,1)]=transpose(B)
  
  # Ds
  W[1+size(Q,1)+size(Ds,1)+size(B,1):size(Q,1)+size(Ds,1)+size(B,1)+size(D,1),1:size(D,2)]=D
  # transpose symmetric
  W[1:size(D,2),1+size(Q,1)+size(Ds,1)+size(B,1):size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)]=transpose(D)
  
  # All the diagonals D
  
  # Dss
  fromy=size(Q,1)+1
  toy=size(Q,1)+size(Ds,1)
  fromx=size(Q,2)+1
  tox=size(Q,2)+size(Ds,1)
  for i=1:size(Ds,1)-1
    W[fromy+i,fromx+i]=Ds[i,1]
  end
  
  # Dsy
  fromy=size(Q,1)+size(Ds,1)+1
  toy=size(Q,1)+size(Dy,1)+size(Ds,1)
  fromx=size(Q,2)+size(Ds,1)+1
  tox=size(Q,2)+size(Dy,1)+size(Ds,1)
  for i=1:size(Dy,1)-1
    W[fromy+i,fromx+i]=Dy[i,1]
  end
  
  # Dsz
  fromy=size(Q,1)+size(Dy,1)+size(Ds,1)+1
  toy=size(Q,1)+size(Dz,1)+size(Dy,1)+size(Ds,1)
  fromx=size(Q,2)+size(Dy,1)+size(Ds,1)+1
  tox=size(Q,2)+size(Dz,1)+size(Dy,1)+size(Ds,1)
  for i=1:size(Dz,1)-1
    W[fromy+i,fromx+i]=Dz[i,1]
  end
  
  # Identity matrix entries
  id1=eye(size(D,1),size(Ds,1))
  id2=eye(size(Ds,1),size(D,1))
  
  fromy=1+size(Q,1)+size(Ds,1)+size(B,1)
  toy=size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)
  fromx=1+size(Q,1)
  tox=size(Q,1)+size(Ds,1)
  W[fromy:toy,fromx:tox]=-id1
  
  fromy=1+size(Q,1)+size(Ds,1)+size(B,1)
  toy=size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)
  fromx=1+size(Q,1)
  tox=size(Q,1)+size(Ds,1)
  W[fromx:tox,fromy:toy]=-id2
  
  return W
end

# This Reads the block diagonal W and off diagonal blocks O and returns them. 
# Each scenario has one of these.
# Iteration defines which iteration should be read.
function buildW(child,iter)
  filename="globalWs" * "$child" * "_" * "$iter" * ".dmp"
  W=openMat(filename)
  
# Now offdiagonal matrices
  
  filename="globalRs" * "$child" * "_" * "$iter" * ".dmp"
  R=openMat(filename)
  filename="globalAs" * "$child" * "_" * "$iter" * ".dmp"
  A=openMat(filename)
  filename="globalCs" * "$child" * "_" * "$iter" * ".dmp"
  C=openMat(filename)
  
  @assert(size(R,2)==size(A,2))
  @assert(size(R,2)==size(C,2))
  # @assert(size(R,1)==size(Q,2))
  # @assert(size(A,1)==size(B,1))
  # @assert(size(C,1)==size(D,1))
  # O=zeros(size(R,2),size(R,1)+size(Ds,1)+size(A,1)+size(C,1))
  O=zeros(size(R,2),size(W,2))
  sDs=size(W,2)-(size(R,1)+size(A,1)+size(C,1))
  
  O[1:size(R,2),1:size(R,1)]=transpose(R)
  O[1:size(R,2),1+size(R,1)+sDs:size(R,1)+sDs+size(A,1)]=transpose(A)
  O[1:size(R,2),1+size(R,1)+sDs+size(A,1):size(R,1)+sDs+size(A,1)+size(C,1)]=transpose(C)

  @assert(size(W,2)==size(O,2))
  
  return W,O
end

# This should do exactly the same as buildW, except that the block diagonal
# matrices are assmebled by reading Q, B, D and the diagonals.
function buildW2(child,iter)
  # reading in block
  filename="globalQs" * "$child" * "_" * "$iter" * ".dmp"
  Q=openMat(filename)
  filename="globalBs" * "$child" * "_" * "$iter" * ".dmp"
  B=openMat(filename)
  filename="globalDsss" * "$child" * "_" * "$iter" * ".dmp"
  Ds=openMat(filename)
  filename="globalDsxs" * "$child" * "_" * "$iter" * ".dmp"
  Dx=openMat(filename)
  filename="globalDsys" * "$child" * "_" * "$iter" * ".dmp"
  Dy=openMat(filename)
  filename="globalDszs" * "$child" * "_" * "$iter" * ".dmp"
  Dz=openMat(filename)
  filename="globalDs" * "$child" * "_" * "$iter" * ".dmp"
  D=openMat(filename)
  
  # computing block size
  @assert(size(Q,1)==size(B,2))
  @assert(size(Q,1)==size(D,2))
  @assert(size(B,1)==size(Dy,1))
  @assert(size(D,1)==size(Dz,1))
  n=size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)
  W=zeros(n,n)
  
  # move blocks to the right places
  # Read Dx
  # Dsx
  T=zeros(size(Dx,1),size(Dx,1))
  for i=1:size(Dx,1)
    T[i,i]=Dx[i,1]
  end
  
  # Compute \bar{Q}=Q*Dx
  Q=Q+T
  # start with Qs
  W[1:size(Q,1),1:size(Q,2)]=Q
  
  # Bs
  W[1+size(Q,1)+size(Ds,1):size(Q,1)+size(Ds,1)+size(B,1),1:size(B,2)]=B
  # transpose symmetric
  W[1:size(B,2),1+size(Q,1)+size(Ds,1):size(Q,1)+size(Ds,1)+size(B,1)]=transpose(B)
  
  # Ds
  W[1+size(Q,1)+size(Ds,1)+size(B,1):size(Q,1)+size(Ds,1)+size(B,1)+size(D,1),1:size(D,2)]=D
  # transpose symmetric
  W[1:size(D,2),1+size(Q,1)+size(Ds,1)+size(B,1):size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)]=transpose(D)
  
  # All the diagonals D
  
  # Dss
  fromy=size(Q,1)+1
  toy=size(Q,1)+size(Ds,1)
  fromx=size(Q,2)+1
  tox=size(Q,2)+size(Ds,1)
  for i=1:size(Ds,1)-1
    W[fromy+i,fromx+i]=Ds[i,1]
  end
  
  # Dsy
  fromy=size(Q,1)+size(Ds,1)+1
  toy=size(Q,1)+size(Dy,1)+size(Ds,1)
  fromx=size(Q,2)+size(Ds,1)+1
  tox=size(Q,2)+size(Dy,1)+size(Ds,1)
  for i=1:size(Dy,1)-1
    W[fromy+i,fromx+i]=Dy[i,1]
  end
  
  # Dsz
  fromy=size(Q,1)+size(Dy,1)+size(Ds,1)+1
  toy=size(Q,1)+size(Dz,1)+size(Dy,1)+size(Ds,1)
  fromx=size(Q,2)+size(Dy,1)+size(Ds,1)+1
  tox=size(Q,2)+size(Dz,1)+size(Dy,1)+size(Ds,1)
  for i=1:size(Dz,1)-1
    W[fromy+i,fromx+i]=Dz[i,1]
  end
  
  # Identity matrix entries
  id1=eye(size(D,1),size(Ds,1))
  id2=eye(size(Ds,1),size(D,1))
  
  fromy=1+size(Q,1)+size(Ds,1)+size(B,1)
  toy=size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)
  fromx=1+size(Q,1)
  tox=size(Q,1)+size(Ds,1)
  W[fromy:toy,fromx:tox]=-id1
  
  fromy=1+size(Q,1)+size(Ds,1)+size(B,1)
  toy=size(Q,1)+size(Ds,1)+size(B,1)+size(D,1)
  fromx=1+size(Q,1)
  tox=size(Q,1)+size(Ds,1)
  W[fromx:tox,fromy:toy]=-id2
  
  
# Now offdiagonal matrices
  
  filename="globalRs" * "$child" * "_" * "$iter" * ".dmp"
  R=openMat(filename)
  filename="globalAs" * "$child" * "_" * "$iter" * ".dmp"
  A=openMat(filename)
  filename="globalCs" * "$child" * "_" * "$iter" * ".dmp"
  C=openMat(filename)
  
  @assert(size(R,2)==size(A,2))
  @assert(size(R,2)==size(C,2))
  @assert(size(R,1)==size(Q,2))
  @assert(size(A,1)==size(B,1))
  @assert(size(C,1)==size(D,1))
  O=zeros(size(R,2),size(R,1)+size(Ds,1)+size(A,1)+size(C,1))
  
  O[1:size(R,2),1:size(R,1)]=transpose(R)
  O[1:size(R,2),1+size(R,1)+size(Ds,1):size(R,1)+size(Ds,1)+size(A,1)]=transpose(A)
  O[1:size(R,2),1+size(R,1)+size(Ds,1)+size(A,1):size(R,1)+size(Ds,1)+size(A,1)+size(C,1)]=transpose(C)

  @assert(size(W,2)==size(O,2))
  
  return W,O
  
end

function readRHS(iter)
  filename="globalRHS0_" * "$iter" * ".dmp"
  return RHS=openMat(filename)
end

function readSOL(iter)
  filename="globalSOL0_" * "$iter" * ".dmp"
  return SOL=openMat(filename)
end

ioff()
iter=1
if (size(ARGS,1) < 4)
  println("Not enough arguments.")
  println("Usage: solve.jl [#scenarios] [iteration] [assemble master] [assemble scenarios]")
  exit()
end
scenarios=parse(Int32,ARGS[1])
# Read which iteration
iter=parse(Int32,ARGS[2])
# Which method of assmembling. Either buildW, buildW2 or buildM, buildM2.
massemble=parse(Int32,ARGS[3])
sassemble=parse(Int32,ARGS[4])
if (massemble==0)
  M=buildM(iter)
else
  M=buildM2(iter)
end
if (!issym(M))
  println("M is not symmetric")
  exit()
else
  println("M is symmetric")
end

if (sassemble==0)
  W,O=buildW(0,iter)
else
  W,O=buildW2(0,iter)
end
if (!issym(W))
  println("W is not symmetric")
  exit()
else
  println("W is symmetric")
end

# Put the matrices at the right spot. See notes for that.
A=zeros(size(M,1)+scenarios*size(W,1),size(M,2)+scenarios*size(W,2))

A[1+scenarios*size(W,1):size(M,1)+scenarios*size(W,1),1+scenarios*size(W,2):size(M,2)+scenarios*size(W,2)]=M
A[1:size(W,1),1:size(W,2)]=W
A[1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1),1:size(O,2)]=O
A[1:size(O,2),1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)]=transpose(O)

if (!issym(A))
  println("A is not symmetric 1")
  exit()
end

for i=1:scenarios-1
  W,O=buildW(1,iter)
  A[1+i*size(W,1):(i+1)*size(W,1),1+i*size(W,2):(i+1)*size(W,2)]=W
  A[1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1),1+i*size(O,2):(i+1)*size(O,2)]=O
  A[1+i*size(O,2):(i+1)*size(O,2),1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)]=transpose(O)
end
# A[10:10]=10
if (!issym(A))
  println("A is not symmetric 2")
  exit()
end
# exit()

# Show structure
spy(A)
show()

# Read right-hand side and solution from PIPS. Do the same computation PIPS and compare
RHS=readRHS(iter)
SOL=readSOL(iter)
println("A: ", size(A,1), "x", size(A,2))
A_cond=cond(A)
println("Condition number: ", A_cond)
println("Solving system")
x=\(A,RHS)
r_j=A*x-RHS
r_p=A*SOL-RHS
println("Residual Julia: ", norm(r_j,Inf))
println("Residual PIPS: ", norm(r_p,Inf))

ro_j=norm(r_j,Inf)/(norm(A,Inf)*norm(x,Inf))
ro_p=norm(r_p,Inf)/(norm(A,Inf)*norm(SOL,Inf))
println("Weighted residual Julia:", ro_j)
println("Weighed residual PIPS", ro_p)

println("Error bound Julia: ", norm(x-SOL,Inf)/norm(x,Inf), " ", A_cond*ro_j)
println("Error bound PIPS: ", norm(SOL-x,Inf)/norm(SOL,Inf), " ", A_cond*ro_p)

println("Solution from file:")
println(norm(SOL,Inf))
println("Computed solution:")
println(norm(x,Inf))
println("RHS from file:")
println(norm(RHS,Inf))
println("Computed RHS:")
println(norm(A*x,Inf))




