#= This file executes numerically solves the dynamics of a freely-deforming
Joukowski foil shedding vortices from its trailing edge =#
cd("/Users/blake/Dropbox/CMU/julia/swim")

include("generateStructs.jl")
include("joukowskiDynamics.jl")
# include("animateJoukowski.jl")

using .generateStructs, .joukowskiDynamics
using DifferentialEquations
using Debugger
using Plots

# Initialize vortex parameters
maxNumOfVortices = 1000;
vortexFlag = 0;
mvortex = zeros(maxNumOfVortices, 3);
mvortex0 = zeros(maxNumOfVortices, 6);
timeOfVortexSheddingEvent=0;

vp = vortexParams(maxNumOfVortices, timeOfVortexSheddingEvent, vortexFlag, mvortex, mvortex0);

# Initialize foil parameters
x = 0.0; y = 0.0; theta = 0.0; U = 0.0; V = 0.0; Omega = 0.0; gammac = 0.0; rc=1.0; ua = -0.3; uf = 2.0;
fp = foilParams(x,y,theta,U,V,Omega,gammac,rc,ua,uf);

# Initialize generic system parameters
n = 3;
N = n + 2*maxNumOfVortices;
Pinital = [0, 0, 0];
npv = 0; # number of initially-present vortices in the flow
sp = sysParams(n,N,Pinital,npv)

# npp = 0;
# U0 = 0; V0 = 0; Omega0 = 0; gammac = 0; npv = 0; Pic = [0 0 0];
# U=0; V=0; Omega=0; Pvortex0=zeros(3,1);

# Initialize solver vars
T = 4.0;
deltaT = 0.001;
Tspan = (0.0,T);

x0 = zeros(N,1);

#= Compute initial conditions corresponding to momentum conservation if at time
t = 0 there are already vortices or other flow features =#
p = vp, fp, sp;
prob = ODEProblem(joukowskiDyn,x0,Tspan,p);
# break_on(:error)
alg = ExplicitRK(tableau=constructDormandPrince())
# alg = Tsit5();
sol = solve(prob,alg,saveat=deltaT);

# traj = plot(sol.t,[sol[1,:] sol[2,:] sol[3,:]],linewidth=2,label=["x" "y" "θ"],xaxis = ("t", (0.0,T)))
# savefig(traj,"traj.pdf")
