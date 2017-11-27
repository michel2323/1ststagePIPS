using PyPlot

# declare triplet type for output to file
type triplet
    i::Int64
    j::Int64
    v::Float64
end

function symmetry(A)
    issym(A)
end

# Read dumped matrix
# Format is m, n integers setting the sizeof
# and m*n entries
function openMat(filename)
  if (isfile(filename))
    fp=open(filename)
    m=read(fp, Int32)
    n=read(fp, Int32)
    A=read(fp, Float64, m, n)
    close(fp)
    println("Reading matrix from ", filename, " ", m, "x", n)
  else
    println("Reading empty matrix from ", filename)
    A=Array{Float64}(0,2)
  end
  return A
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
function buildSelfAssembleM(iter)
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
function buildSelfAssembleW(child,iter)
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

function main_()
ioff()
iter=1
if (size(ARGS,1) < 5)
  println("Not enough arguments.")
  println("Usage: solve.jl [#scenarios] [iteration] [assemble master] [assemble scenarios] [dump triplet format file]")
  exit()
end
scenarios=parse(Int32,ARGS[1])
# Read which iteration
iter=parse(Int32,ARGS[2])
# Which method of assembling. Either buildW, buildW2 or buildM, buildM2.
massemble=parse(Int32,ARGS[3])
sassemble=parse(Int32,ARGS[4])
dump_triplet=parse(Int32,ARGS[5])
if (massemble==0)
  M=buildM(iter)
else
  M=buildSelfAssembleM(iter)
end
if (!symmetry(M))
  println("M is not symmetric")
  exit()
else
  println("M is symmetric")
end

if (sassemble==0)
  W,O=buildW(0,iter)
else
  W,O=buildSelfAssembleW(0,iter)
end
if (!symmetry(W))
  println("W is not symmetric")
  exit()
else
  println("W is symmetric")
end
#declare list to store triplets
list=triplet[]
    
# build triplet list for output to file
j_=0
i_=0
if dump_triplet==1
    # master problem
    i_=1
    for i in 1+scenarios*size(W,1):size(M,1)+scenarios*size(W,1)
        j_=1
        for j in 1+scenarios*size(W,2):size(M,2)+scenarios*size(W,2)
            tmp=triplet(i,j,M[i_,j_])
            if tmp.v!=0.0
                push!(list,tmp)
            end
            j_=j_+1
        end
        i_=i_+1
    end
    i_=1
    for i in 1:size(W,1)
        j_=1
        for j in 1:size(W,2)
            tmp=triplet(i,j,W[i_,j_])
            if tmp.v!=0.0
                push!(list,tmp)
            end
            j_=j_+1
        end
        i_=i_+1
    end
    i_=1
    for i in 1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)
        j_=1
        for j in 1:size(O,2)
            tmp=triplet(i,j,O[i_,j_])
            if tmp.v!=0.0
                push!(list,tmp)
            end
            j_=j_+1
        end
        i_=i_+1
    end
    Ot=transpose(O)
    i_=1
    for i in 1:size(O,2)
        j_=1
        for j in 1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)
            tmp=triplet(i,j,Ot[i_,j_])
            if tmp.v!=0.0
                push!(list,tmp)
            end
            j_=j_+1
        end
        i_=i_+1
    end
else
    # Put the matrices at the right spot. See notes for that.
    A=zeros(size(M,1)+scenarios*size(W,1),size(M,2)+scenarios*size(W,2))
    
    A[1+scenarios*size(W,1):size(M,1)+scenarios*size(W,1),1+scenarios*size(W,2):size(M,2)+scenarios*size(W,2)]=M
    A[1:size(W,1),1:size(W,2)]=W
    A[1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1),1:size(O,2)]=O
    A[1:size(O,2),1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)]=transpose(O)
end

if dump_triplet!=1
    if (!symmetry(A))
        println("A is not symmetric 1")
        exit()
    end
end

for s=1:scenarios-1
    if (sassemble==0)
      W,O=buildW(s,iter)
    else
      W,O=buildSelfAssembleW(s,iter)
    end
  # W,O=buildW(2,iter)
  println(s)
  
  if dump_triplet==1
      i_=1
      for i2 in 1+s*size(W,1):(s+1)*size(W,1)
          j_=1
          for j in 1+s*size(W,2):(s+1)*size(W,2)
              tmp=triplet(i2,j,W[i_,j_])
              if tmp.v!=0.0
                  push!(list,tmp)
              end
              j_=j_+1
          end
          i_=i_+1
      end
      i_=1
      for i2 in 1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)
          j_=1
          for j in 1+s*size(O,2):(s+1)*size(O,2)
              tmp=triplet(i2,j,O[i_,j_])
              if tmp.v!=0.0
                  push!(list,tmp)
              end
              j_=j_+1
          end
          i_=i_+1
      end
      Ot=transpose(O)
      i_=1
      for i2 in 1+s*size(O,2):(s+1)*size(O,2)
          j_=1
          for j in 1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)
              tmp=triplet(i2,j,Ot[i_,j_])
              if tmp.v!=0.0
                  push!(list,tmp)
              end
              j_=j_+1
          end
          i_=i_+1
      end
  else
      A[1+s*size(W,1):(s+1)*size(W,1),1+s*size(W,2):(s+1)*size(W,2)]=W
      A[1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1),1+s*size(O,2):(s+1)*size(O,2)]=O
      A[1+s*size(O,2):(s+1)*size(O,2),1+scenarios*size(W,1):scenarios*size(W,1)+size(O,1)]=transpose(O)
  end
end
println("Size: ", size(list))
if dump_triplet!=0 
    f = open("globalmat.mtx", "w");
    for ln in list
        n1=ln.i
        n2=ln.j
        n3=ln.v
        write(f,"$n1 $n2 $n3\n")
    end
    close(f)
    exit(1)
end
# A[10:10]=10
if (!symmetry(A))
  println("A is not symmetric 2")
  exit()
end
# exit()

# Show structure
 # spy(A)
 # show()

# Read right-hand side and solution from PIPS. Do the same computation PIPS and compare
RHS=readRHS(iter)
SOL=readSOL(iter)
println("A: ", size(A,1), "x", size(A,2))
A_cond=cond(A,Inf)
println("Condition number: ", A_cond)
println("Solving system")
x=\(A,RHS)
r_j=A*x-RHS
r_p=A*SOL-RHS
println("Residual Julia: ", norm(r_j,Inf))
println("Residual PIPS: ", norm(r_p,Inf))

ro_j=norm(r_j,Inf)/(norm(A,Inf)*norm(x,Inf))
ro_p=norm(r_p,Inf)/(norm(A,Inf)*norm(SOL,Inf))
println("Weighted residual Julia: ", ro_j)
println("Weighed residual PIPS: ", ro_p)

n0=norm(SOL,Inf)
n1=norm(x-SOL,Inf)
n2=norm(x,Inf)
println("Cosmin term 3: ", n0, " ", n1, " ", n2)
norminf=n1/(1+n2)
println("Cosmin term 1: ", norminf)
# println("Cosmin term 1: ", norm(x-SOL,Inf))
n0=norm(SOL,2)
n1=norm(x-SOL,2)
n2=norm(x,2)
norm2=n1/(1+n2)
println("Cosmin term 4: ", n0, " ", n1, " ", n2)
println("Cosmin term 2: ", norm2)
# println("Cosmin term 2: ", norm(x-SOL,2))

println("Error bound Julia: ", norm(x-SOL,2)/norm(x,2), " ", A_cond*ro_j)
println("Error bound PIPS: ", norm(SOL-x,2)/norm(SOL,2), " ", A_cond*ro_p)

# println("Solution from file: ", norm(SOL,Inf))
println("Solution from file: ", norm(SOL,2))
# println("Computed solution: ", norm(x,Inf))
println("Computed solution: ", norm(x,2))
println("RHS from file: ", norm(RHS,Inf))
println("Computed RHS: ", norm(A*x,Inf))
println("Diff RHS: ", norm(RHS-A*x,Inf))
end
main_()




