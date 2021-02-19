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

function decode()
    s = ArgParseSettings()

    @add_arg_table s begin
          "--beta","-b"
          help="gauge coupling β, default 1.8"
          arg_type=Float64
          default=1.8
	  "--rho", "-r"
	  help="rho parameter ρ, default 1"
	  arg_type=Float64
	  default=1.
	  "--Nf", "-f"
	  help = "Flavour number Nf, default 1 or 2"
	  arg_type=Float64
	  default=1.
          "--nx","-x"
          help="Spatial extent NX, default 12"
          arg_type=Int
          default=10
	  "--nt","-t"
	  help="Temporal extent NT, default 12"
          arg_type=Int
          default = 10
          "--bssteps","-s"
	  help = "number of bootstrap steps, default 100"
	  arg_type=Int
	  default= 20
	  "--ntherm","-u"
          help="Thermalization steps, default 10000"
          arg_type=Int
          default=500
          "--nmeas","-m"
	  help="Measurments, default 100"
          arg_type=Int
          default=20
          "--ndisc","-d"
          help="Updates between measurmants, default 500"
          arg_type=Int
          default=100
          "--outdir","-o"
          help="Output directory name, default out"
          arg_type=String
          default="out"
          "--topupdate","-c"
          help="Use topological update, default false"
          arg_type=Bool
          default=true
	  

    end
    arg=parse_args(s)
    global β=arg["beta"]
    global ρ=arg["rho"]
    global Nf=arg["Nf"]
    global NX=arg["nx"]
    global NT=arg["nt"]
    global NTHERM=arg["ntherm"]
    global NMEAS=arg["nmeas"]
    global OFS=arg["ndisc"]
    global outdir=arg["outdir"]
    global topupd=arg["topupdate"]
    global BS=arg["bssteps"]
end


# helper functions for more readable modulo operations
T(t)=mod(t,NT)
X(x)=mod(x,NX)

# define the staple (i.e. A^†) at point (t,x) in direction μ (exercise 1a)
stap(u,t,x,μ)=(exp(im*(-u[t,x,1-μ]-u[T(t+μ),X(x+1-μ),μ]+u[T(t+1-μ),X(x+μ),1-μ]))
              +exp(im*(u[T(t-μ),X(x+μ-1),1-μ]-u[T(t-μ),X(x+μ-1),μ]-u[T(t+1-2*μ),X(x+2*μ-1),1-μ])))

# define the plaquette
plaq(u)=exp.(im*(u[:,:,0]+circshift(u[:,:,1],(-1,0))-circshift(u[:,:,0],(0,-1))-u[:,:,1]))

# define the action
act(u,β)=-β*sum(cos.(u[:,:,0]+circshift(u[:,:,1],(-1,0))-circshift(u[:,:,0],(0,-1))-u[:,:,1]).-1)

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
	#println(prop)
# compute the action change
      ΔS=β*real(staple*(exp(im*link)-exp(im*prop)))   
      	#println("delta ", ΔS )
# metropolis step
      if exp(-ΔS)> rand()
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
			      u[t, x, μ] = metro(u[t, x, μ] , stap(u,t,x,μ), β)
                        end
                  end
            end
      end
      #println(U)
      #println("plaquettes", plaq(u))
      
# local update is done
      if topupd
# produce new MC suggestion by multiplying with a topological config
            um=u.+topconf(rand(0:1)*2-1)
# compute the action difference
            ΔS=act(um,β)-act(u,β)
            #println("delta S", exp(-ΔS))
# and perform the MC step
            if exp(-ΔS)>rand()
		  #println("great success topo multiplication in" )
                  u[:,:,:]=um[:,:,:]
	    
            end
      end
      #println("charge: ",topch(u))
      return u
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


# write all observables of a config to disc
function Confout(U,ncfg)
# config U, config index ncfg

# the base configuration name
      basname=string(outdir,"/",string(β),"_",string(NX),"_",string(NT),"_",string(ncfg))
# name for gauge observables file
      gfname=string(basname,".ptch")
# name for the naive spectrum file
      nsname=string(basname,".spectN")
# name for the Wilson spectrum file
      wsname=string(basname,".spectW")

# write gauge observables
      open(gfname,"w") do io
# write the plaquette averages
            println(io,mean(plaq(U)))
# write the topological charge
            println(io,topch(U))
      end

# compute the naive eigenmodes
      ev=Evn(U,0.0)
# print them to naive spectrum file
      open(nsname,"w") do io
# write the eigenvalues in a useful manner
            for i in 1:2*NT*NX
                  println(io,real(ev[i]),"   ",imag(ev[i]))
            end
      end

# compute the Wilson eigenmodes
      ev=Evw(U,0.0)
# print them to Wilson spectrum file
      open(wsname,"w") do io
# write the eigenvalues in a useful manner
            for i in 1:2*NT*NX
                  println(io,real(ev[i]),"   ",imag(ev[i]))
            end
      end
end



##### Bootstrap ###########

function Bootstrap( U_saved)
println("start Bootstrap algorithm")
m =  [-0.1:0.001:0.1;]
U_new_config = OffsetArray(zeros(NMEAS,NT,NX,2),0:NMEAS-1,0:NT-1,0:NX-1,0:1)
ov_condensates =  OffsetArray(zeros(BS, length(m)),0:BS-1, 1:length(m))
st_condensates =  OffsetArray(zeros(BS, length(m)),0:BS-1, 1:length(m))

	for b=0:BS-1
		println("BOOOOOOOTstrap", b)
		for k=0:NMEAS-1

			r = rand(0:NMEAS-1)
			U_new_config[k,:,:,:] = U_saved[r,:,:,:]
			
		end

		ov_condensates[b,:], st_condensates[b,:] = condensate(U_new_config, m)
	end
	return ov_condensates, st_condensates, m
end


function condensate(U_config, m)
ov_numerator = OffsetArray(zeros(length(m)),0:length(m)-1)
ov_denominator = OffsetArray(zeros(length(m)),0:length(m)-1)
ov_condensate = OffsetArray(zeros(length(m)),0:length(m)-1)
st_numerator = OffsetArray(zeros(length(m)),0:length(m)-1)
st_denominator = OffsetArray(zeros(length(m)),0:length(m)-1)
st_condensate = OffsetArray(zeros(length(m)),0:length(m)-1)

	for i = 0:NMEAS-1
		
		evo = Evo(U_config[i,:,:,:], 0.0)
		evw = Evw(U_config[i,:,:,:], 0.0)
		for j=0:length(m)-1
			numerator, denominator =  ov_scal(evo, m[j+1])
			ov_numerator[j] += numerator
			ov_denominator[j] += denominator
			
			numerator, denominator = st_scal(evw, m[j+1])
			st_numerator[j] += numerator
			st_denominator[j] += denominator
		end
	end
	for j = 0:length(m)-1
		ov_condensate[j] =(ov_numerator[j]/ov_denominator[j]) *sqrt(β)/(NT*NX)^2
		st_condensate[j] =(st_numerator[j]/st_denominator[j]) *sqrt(β)/(2*(NT*NX)^2)
	end
	return ov_condensate, st_condensate
end

function ov_scal(evo, m)	#overlap scalar condensate
	primed_sum = 0.0
	det = 1
	numerator, denominator = 0.0,0.0
	if abs(m) > 10^(-2)
		for i = 0: NT*NX*2-1
			#primed sum
			if abs(evo[i+1]-2*ρ) > 10^(-2)		#for complex conjugated pairs
				lambda_tilde =  2*ρ*evo[i+1]/(2*ρ-evo[i+1])
				if abs(imag(evo[i+1])) > 10^(-4)
					primed_sum += (real(lambda_tilde) + m) / (abs(lambda_tilde)^2 + m^2 + 2*m*real(lambda_tilde))
					
				end

				if abs(imag(evo[i+1])) < 10^(-4)	#for topological eigenvalues
					primed_sum += 1/(lambda_tilde +m)
					
				end
				
			end
			#determinant
			det *= ((1-m / (2*ρ)) * evo[i+1] + m)
		end
		numerator = real(primed_sum) * real(det^Nf)/NMEAS
		denominator = real(det^Nf)/NMEAS
		
	end	
	return numerator, denominator
end


function st_scal(evw, m)	#staggered scalar condensate
	sum = 0.0
	det = 1
	numerator, denominator = 0.0 ,0.0
	
	if abs(m) > 10^(-2)
		for i = 0:NT*NX*2-1
			#sum
			if abs(imag(evw[i+1])) > 10^(-4)   	#for complex conjugated pairs
				sum += 1(real(evw[i+1]) + m) / (abs(evw[i+1])^2 + m^2 + 2*m*real(evw[i+1]))
			end
		
			if abs(evw[i+1]) < 10^(-4) 		#for topological eigenvalues
				sum += 1/(evw[i+1]+ m)
			end
			#determinant
			det *= (evw[i+1]+ m)
		end
		numerator = real(sum) * real(det^(Nf/2))/NMEAS
		denominator = real(det^(Nf/2))/NMEAS
	end
	return numerator, denominator
end
		


##### the main program #####


# parse command line arguments
decode()
#println("Lattice setup: β=",β," NX=",NX," NT=",NT)
#println("Discarding first ",NTHERM," updates")
#println("Performing ",NMEAS," measurments separated by ",OFS," updates")

# create the results directory
#mkpath(outdir)
# the trivial gauge config
#U=OffsetArray(zeros(NT,NX,2),0:NT-1,0:NX-1,0:1)



U_saved = OffsetArray(zeros(NMEAS, NT,NX,2),0:NMEAS-1,0:NT-1,0:NX-1,0:1)

# hot start: a random gauge field
U= OffsetArray(2*π*rand(NT,NX,2),0:NT-1,0:NX-1,0:1)

# thermalization loop
println("start metropolis algorithm")
for i=1:NTHERM
	
      update(U,β)
end

# measurment loop
for m=0:NMEAS-1
	U_saved[m,:,:,:] = U[:,:,:]
	# skip measurments to decorrelate
	#println("step", m, "U", U)
     	for i=1:OFS
		update(U,β)
    	end
end

ov_condensate, st_condensate, m =Bootstrap(U_saved)
y_ov = OffsetArray(zeros(length(m)), 1:length(m))
y_ov_err = OffsetArray(zeros(length(m)), 1:length(m))
y_st = OffsetArray(zeros(length(m)), 1:length(m))
y_st_err = OffsetArray(zeros(length(m)), 1:length(m))
for i=1:length(m)
	y_ov[i] = mean(ov_condensate[:,i])
	y_ov_err[i] = std(ov_condensate[:,i])
	y_st[i] = mean(st_condensate[:,i])
	y_st_err[i] = std(st_condensate[:,i])
	
end
p1 = plot(m,ov_condensate[1,:], grid=false,ribbon=y_ov_err,fillalpha=0.5) 
p2 = plot(m,y_st,grid=false,ribbon=y_st_err,fillalpha=0.5) 
plot(p1,p2,layout=(1,2), legend=false)

