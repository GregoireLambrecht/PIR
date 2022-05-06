include("tools.jl")
using Plots

const ħ=6.63e-34
const m = 9.1e-32  #Un dixième de la masse de l'éléctron libre
const eV = 1.6e-19 #Un électron Volt
const nm = 1e-9 #Le nanométre

#La hauteur des potentiels dans les composants sont généralement
#Entre 0.2 et 0.5 eV
#Les largeurs barrières et puits sont environ de 2 à 5 nm

mutable struct Barr
	L
	V
	ϵ
	k
	K
	r
	t
	T
	R
	N
end

function init_Barr(ϵ,L,V)
	k = sqrt(2m*ϵ)/ħ
	K = sqrt(2m*(V-ϵ))/ħ
	κ1 = k*k - K*K
	κ2 = k*k + K*K
	r = -(sinh(K*L)*κ2)/(2im*k*K*cosh(K*L) + sinh(K*L)*κ1)
	t = (2im*k*K/exp(1im*K*L)) * (1/(2im*K*k*cosh(K*L) + sinh(K*L)*κ1))
	T = real(t't)
	R = real(r'r)
	N = zeros(ComplexF64,2,2)
	N[1,1] = 1/conj(t)
	N[1,2] = r/t
	N[2,1] = conj(r)/conj(t)
	N[2,2] = 1/t
	Barr(L,V,ϵ,k,K,r,t,T,R,N)
end

function TraceT(n)
	E = [1/n*eV*i for i=1:(n-1)]
	T = zeros(n-1)
	for k=1:(n-1)
		P = init_Barr(nm,eV,E[k])
		T[k] = P.T
	end
	plot(E./eV,T)
	(T,E./eV)
end

function ParamDB(ϵ,V1,V2,L1,L2,Lp)
	Bar1 = init_Barr(ϵ,L1,V1)
	Bar2 = init_Barr(ϵ,L2,V2)
	N = Bar1.N * Bar2.N
	k = Bar1.k
	t = 1/conj(N[1,1])
	r = N[1,2]*t
	A1 = 1 + r
	B1 =1im*k*(1-r)/Bar1.K
	α = -(1im*exp(-1im*k*L1)/(2*k))*((1im*k*A1 + B1*Bar1.K)*cosh(L1*Bar1.K) + (1im*k*B1 + A1*Bar1.K)*sinh(L1*Bar1.K))
	β = (1im*exp(1im*k*L1)/(2*k))*((-1im*k*A1 + B1*Bar1.K)*cosh(L1*Bar1.K) + (-1im*k*B1 + A1*Bar1.K)*sinh(L1*Bar1.K))
	L = L1 + L2 + Lp
	A2 = (exp(1im*k*L)*t/Bar2.K)*(cosh(L*Bar2.K)*Bar2.K - 1im*k*sinh(L*Bar2.K))
	B2 = (exp(1im*k*L)*t/Bar2.K)*(cosh(L*Bar2.K)*1im*k - Bar2.K*sinh(L*Bar2.K))
	(t,r,A1,B1,α,β,A2,B2)
end

function doubleBarrT(n)
	V = 0.4
	l = 3
	E = [1/n*V*eV*i for i=1:(n-1)]
	T = zeros(n-1)
	#Premiere situation
	for k=1:(n-1)
		(t,r,A1,B1,α,β,A2,B2) = ParamDB(E[k],V*eV,V*eV,l*nm,l*nm,l*nm)
		T[k] = real.(t't)
	end
	ϵres = argmax(T) * V/n #L'energie de resonnance en eV
	ϵres = (round(ϵres * 10000))/10000 #On arrondie
	Γ = Largeur(T,0.001)/n #La largeur à mi-hauteur
	pl1 = plot(E./eV,T,xlabel = "Energie en eV, ϵres = $ϵres eV, Γ = $Γ eV" ,ylabel = "T(E)", title = "T(E) avec L1 = L2 = Lp = $l nm et V1 = V2 = $V eV")
	savefig(pl1,"situation1.pdf")
	#Situation 2 V1 = 2 V2
	for k=1:(n-1)
		(t,r,A1,B1,α,β,A2,B2) = ParamDB(E[k],2*V*eV,V*eV,l*nm,l*nm,l*nm)
		T[k] = real.(t't)
	end
	pl1 = plot(E./eV,T,xlabel = "Energie en eV",ylabel = "T(E)", title = "T(E) avec L1 = L2 = Lp = $l nm et V2 = V1/2 = $V eV")
	savefig(pl1,"situation2.pdf")
	#Situation 3, V2 varie, ϵ = 0.4 eV
	V2 = [(eV*V/2) *(1 + 2*i/n) for i = 1:(n-1)]
	for k=1:(n-1)
		(t,r,A1,B1,α,β,A2,B2) = ParamDB(ϵres*eV,V*eV,V2[k],l*nm,l*nm,l*nm)
		T[k] = real.(t't)
	end
	pl1 = plot(V2./eV,T,xlabel = "V2 en eV",ylabel = "T(V2)", title = "T(V2), L1=L2=Lp=$l nm, V1=$V eV, E=$ϵres eV")
	savefig(pl1,"situation3.pdf")
end



doubleBarrT(10000)
