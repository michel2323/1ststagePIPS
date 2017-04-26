using PyPlot

function openMat(filename)
  fp=open(filename)
  m=read(fp, Int32)
  n=read(fp, Int32)
  A=read(fp, Float64, m, n)
  close(fp)
  println("Reading matrix from ", filename, " ", m, "x", n)
  return A
end

function firststage()
  println("Reading matrix M")
  A=openMat("1ststageM.dmp")
  spy(A)
  show()
  
  println("Reading RHS")
  rhs=openMat("1ststageRHS.dmp")
  
  println("Reading solution")
  sol=openMat("1ststageSol.dmp")
  
  println("Solving system")
  x=\(A,rhs)
  
  println("Solution from file:")
  println(norm(sol))
  println("Computed solution:")
  println(norm(x))
  
  println("Condition number: ", cond(A))
end

function buildMaster()
  Q=openMat("globalQ00_1.dmp")
  spy(Q)
  show()
  B=openMat("globalB00_1.dmp")
  spy(B)
  show()
  C=openMat("globalC00_1.dmp")
  spy(B)
  show()
end
function buildM(iter)
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
  Q=Q*T
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
function buildW(child,iter)
  # reading in block
  filename="globalQs" * "$child" * "_" * "$iter" * ".dmp"
  println(filename)
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
  Q=Q*T
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

function readOffDiag(child,iter)
  
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
iter=2
scenarios=2
M=buildM(iter)
println(size(M,1))
W,O=buildW(0,iter)
A=zeros(size(M,1)+scenarios*size(W,1),size(M,2)+scenarios*size(W,2))
println("A: ", size(A,1), " ", size(A,2))

A[1+scenarios*size(W,1):size(M,1)+scenarios*size(W,1),1+scenarios*size(W,2):size(M,2)+scenarios*size(W,2)]=M
A[1:size(W,1),1:size(W,2)]=W
A[1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1),1:size(O,2)]=O
A[1:size(O,2),1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)]=transpose(O)

for i=1:scenarios-1
  W,O=buildW(1,iter)
  A[1+i*size(W,1):(i+1)*size(W,1),1+i*size(W,2):(i+1)*size(W,2)]=W
  A[1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1),1+i*size(O,2):(i+1)*size(O,2)]=O
  A[1+i*size(O,2):(i+1)*size(O,2),1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)]=transpose(O)
end
spy(A)
show()

RHS=readRHS(iter)
SOL=readSOL(iter)
println("Solving system")
x=\(A,RHS)

println("Solution from file:")
println(norm(SOL))
println("Computed solution:")
println(norm(x))

println("Condition number: ", cond(A))

