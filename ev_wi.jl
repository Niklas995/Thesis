

import Pkg
Pkg.add("OffsetArrays")
Pkg.add("Plots")
Pkg.add("ArgParse")


#!/usr/local/bin/julia

# arrays with arbitrary starting index
using OffsetArrays
# statistics package
using Statistics
# plotting
using Plots
# delimited fields file output
using DelimitedFiles
# command line argument parsing
using ArgParse
# LAPACK
using LinearAlgebra


β= 3.2
NX=4
NT=4
NTHERM=10000
NMEAS=10
OFS=100 #ndisc
#outdir=out
#topupd=true  

# helper functions for more readable modulo operations
T(t)=mod(t,NT)
X(x)=mod(x,NX)

# define the staple (i.e. A^†) at point (t,x) in direction μ (exercise 1a)
stap(u,t,x,μ)=(exp(im*(-u[t,x,1-μ]-u[T(t+μ),X(x+1-μ),μ]+u[T(t+1-μ),X(x+μ),1-μ]))
              +exp(im*(u[T(t-μ),X(x+μ-1),1-μ]-u[T(t-μ),X(x+μ-1),μ]-u[T(t+1-2*μ),X(x+2*μ-1),1-μ])))

# define the plaquette
plaq(u)=exp.(im*(u[:,:,0]+circshift(u[:,:,1],(-1,0))-circshift(u[:,:,0],(0,-1))-u[:,:,1]))

# define the action
act(u,β)=-β*sum(cos.((u[:,:,0]+circshift(u[:,:,1],(-1,0))-circshift(u[:,:,0],(0,-1))-u[:,:,1])).-1)

# applying a gauge transformation
gt(u,g)=OffsetArray(cat(g+u[:,:,0]-circshift(g,(-1,0)),g+u[:,:,1]-circshift(g,(0,-1)),dims=3),0:NT-1,0:NX-1,0:1)

# the topological charge
topch(u)=imag(sum(log.(plaq(u))))/(2*π)

# constructing a gauge transformation for temporal gauge
function tempgauge(u)
# return a config that is gauge equivalent
# to u and has temporal gauge

# define the gauge transformation
      g=OffsetArray(zeros(NT,NX),0:NT-1,0:NX-1)
# line by line except the last one
      for t=0:NT-2
# gauge transformation at that line from phase of previous line
            g[t+1,:]=g[t,:]+u[t,:,0]
      end

# apply the constructed gauge transformation
      return gt(u,g)
end

# construct an instanton of charge ch
function topconf(ch)
# start with a trivial gauge config
      u=OffsetArray(zeros(NT,NX,0:1),0:NT-1,0:NX-1,0:1)
# the field per Plaquette
      F=2*π*ch/(NX*NT)
# set field on every timeslice but the last
      for t=1:NT-1
            u[t,:,1].=t*F
      end
# field on last time slice
      for x=1:NX-1
            u[NT-1,x,0]=-x*F*NT
      end
# return the final construct
      return u
end

# Metropolis update
function metro(link,staple,β)
# metropolis update link given the stable and β
# proposed new link
      prop=link+(rand()-rand())
# compute the action change
      ΔS=β*real(staple*(exp(im*link)-exp(im*prop)))
# metropolis step
      if exp(-ΔS)>rand()
            #print(" accept! (",ΔS,") ")
            return prop
      else
            #print(" reject! (",ΔS,") ")
            return link
      end
end

# checkerboard loop over updates
function update(u,β)
# outermost loop: loop over directions
      for μ=0:1
# then even-odd
            for eo=0:1
# then loop over space
                  for x=0:NX-1
# compute even-odd offset
                        ofs=(x&1)⊻eo
# loop over time
                        for t=ofs:2:NT-1
# perform metropolis update
                              u[t,x,μ]=metro(u[t,x,μ],stap(u,t,x,μ),β)
                        end
                  end
            end
      end
# local update is done
      if topupd
# produce new MC suggestion by multiplying with a topological config
            um=u.+topconf(rand(0:1)*2-1)
# compute the action difference
            ΔS=act(um,β)-act(u,β)
# and perform the MC step
            if exp(-ΔS)>rand()
                  u[:,:,:]=um[:,:,:]
            end
      end
end

# exercise 11, fermion operators

# the Euclidean γ matrices

# the Euclidean γ0 is equal to the Minkovski γ0
γ0=[1.0 0.0 ; 0.0 -1.0]
# the Euclidean γ1 is equal to -im*γ1 the Minkovski γ1
γ1=[0.0 -im ;  im  0.0]
# th 2*2 identity matrix
Id=Matrix{Float64}(I,2,2)



# the Naive Dirac operator
function Dnaive(U,m)
# Gauge field U, mass m
# first zero all elements
      D=OffsetArray(complex(zeros(NT,NX,2,NT,NX,2)),0:NT-1,0:NX-1,1:2,0:NT-1,0:NX-1,1:2)
# loop over all points
      for x in 0:NX-1, t in 0:NT-1
# add the diagonal mass term
            D[t,x,:,t,x,:]=Id*m
# add the forward in time neighbor
            D[t,x,:,T(t+1),x,:]=exp(im*U[t,x,0])*γ0/2.0
# backward neighbor in time
            D[t,x,:,T(t-1),x,:]=exp(-im*U[T(t-1),x,0])*(-γ0)/2.0
# forward neighbor in space
            D[t,x,:,t,X(x+1),:]=exp(im*U[t,x,1])*γ1/2.0
# backward neighbor in space
            D[t,x,:,t,X(x-1),:]=exp(-im*U[t,X(x-1),1])*(-γ1)/2.0
      end

      return D
end

# the Wilson Dirac operator
function Dwilson(U,m)
# Gauge field U, mass m
# first zero all elements
      D=OffsetArray(complex(zeros(NT,NX,2,NT,NX,2)),0:NT-1,0:NX-1,1:2,0:NT-1,0:NX-1,1:2)
# loop over all points
      for x in 0:NX-1, t in 0:NT-1
# add the diagonal mass term
            D[t,x,:,t,x,:]=Id*(m+2)
# add the forward in time neighbor
            D[t,x,:,T(t+1),x,:]=exp(im*U[t,x,0])*(-Id+γ0)/2.0
# backward neighbor in time
            D[t,x,:,T(t-1),x,:]=exp(-im*U[T(t-1),x,0])*(-Id-γ0)/2.0
# forward neighbor in space
            D[t,x,:,t,X(x+1),:]=exp(im*U[t,x,1])*(-Id+γ1)/2.0
# backward neighbor in space
            D[t,x,:,t,X(x-1),:]=exp(-im*U[t,X(x-1),1])*(-Id-γ1)/2.0
      end

      return D
end

# the Overlap Dirac operator
function Doverlap(U,rho,m)
# Gauge field U, projection point rho, mass m
# first get the wilson operator
    Dw=Dwilson(U,-rho)
# make an SVD
    S=svd(reshape(Dw,(2*NT*NX,2*NT*NX)))
# and produce the overlap operator from it
    D=(S.U*S.Vt+I)*rho
    
    return D
end
    
# Eigenvalues of the Naive operator at mass m on gauge field U
Evn(U,m)=eigvals(reshape(Dnaive(U,0.0),(2*NT*NX,2*NT*NX)))
# Eigenvalues of the Wilson operator at mass m on gauge field U
Evw(U,m)=eigvals(reshape(Dwilson(U,0.0),(2*NT*NX,2*NT*NX)))
# Eigenvalues of the Overlap operator with rho=1 at mass m on gauge field U
Evo(U,m)=eigvals(Doverlap(U,1.0,0.0))


##### the main program #####

U0=OffsetArray(zeros(NT,NX,2),0:NT-1,0:NX-1,0:1)

ev0=Evw(U0,0.0)
print(real(ev0))
print(imag(ev0))

#savefig(plot(real(ev0),imag(ev0), seriestype = :scatter, title = "Eigenvalues", fmt = :png))


U=OffsetArray(2*π*rand(NT,NX,2),0:NT-1,0:NX-1,0:1)

ev=Evw(U,0.0)


#plot(real(ev), imag(ev), seriestype = :scatter, title = "Eigenvalues", fmt = :png))
print(real(ev))
print(imag(ev))
