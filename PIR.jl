include("tools.jl")
using Plots
using QuadGK

const ħ=6.63e-34
const m = 9.1e-32  #Un dixième de la masse de l'éléctron libre
const eV = 1.6e-19 #Un électron Volt
const nm = 1e-9 #Le nanométre

#La hauteur des potentiels dans les composants sont généralement
#Entre 0.2 et 0.5 eV
#Les largeurs barrières et puits sont environ de 2 à 5 nm

mutable struct Barr
	L
	Vg
	V
	Vd
	ϵ
	kg
	K
	kd
	r
	t
	T
	R
	N
end


#Une barrière simple avec Vd = Vg = 0
function init_BarrSimple(ϵ,L,V)
	k = sqrt(2m*ϵ)/ħ
	K = sqrt(2m*(V-ϵ))/ħ
	κ1 = k*k - K*K
	κ2 = k*k + K*K
	r = -(sinh(K*L)*κ2)/(2im*k*K*cosh(K*L) + sinh(K*L)*κ1)
	#Attention
	t = (2im*k*K/exp(1im*K*L)) * (1/(2im*K*k*cosh(K*L) + sinh(K*L)*κ1))
	T = real(t't)
	R = real(r'r)
	N = zeros(ComplexF64,2,2)
	N[1,1] = 1/conj(t)
	N[1,2] = r/t
	N[2,1] = conj(r)/conj(t)
	N[2,2] = 1/t
	Barr(L,0,V,0,ϵ,k,K,k,r,t,T,R,N)
end

#Vd >= 0
#Vg >= 0
function init_Barr(ϵ,L,Vg,V,Vd)
	kd = sqrt(2m*(ϵ + Vd))/ħ
	kg = sqrt(2m*(ϵ + Vg))/ħ
	K = sqrt(2m*(V-ϵ))/ħ
	Ξ = (1im*kd * cosh(K*L) - K *sinh(K*L))/(K*cosh(K*L) - 1im * kd *sinh(K*L))
	r = ((-1im * kg/K) + Ξ)/((1im * kg/K) + Ξ)
	#t = (1 - r)*K*exp(-1im*kd*L)/(K*cosh(K*L) - 1im*kd*sinh(K*L))
	#Attention
	t = (2*K*exp(-1im*K*L))/(K*cosh(K*L) - 1im*kd*sinh(K*L) + (K/(1im*kg))*(1im*kd*cosh(K*L)-K*sinh(K*L)))
	T = real(t't)
	R = real(r'r)
	N = zeros(ComplexF64,2,2)
	N[1,1] = 1/conj(t)
	N[1,2] = -conj(r)/t
	N[2,1] = -r/conj(t)
	N[2,2] = 1/t
	Barr(L,Vg,V,Vd,ϵ,kg,K,kd,r,t,T,R,N)
end


#Vp >= 0
function init_Puits(ϵ,Lp,Vp)
	k = sqrt(2m*ϵ)/ħ
	K = sqrt(2m*(Vp+ϵ))/ħ
	r = 0
	t = exp(-1im*K*Lp)
	T = 1
	R = 0
	N = zeros(ComplexF64, 2,2)
	N[1,1] = 1/conj(t)
	N[2,2] = 1/t
	Barr(Lp,0,Vp,0,ϵ,k,K,k,r,t,T,R,N)
end

function ParamDBsimple(ϵ,V1,V2,L1,L2,Lp)
	Bar1 = init_BarrSimple(ϵ,L1,V1)
	Bar2 = init_BarrSimple(ϵ,L2,V2)
	Puits = init_Puits(ϵ,Lp,0)
	N = Bar1.N*Puits.N*Bar2.N
	k = Bar1.kd
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

function ParamDB(ϵ,V1,Vp,V2,Vd,L1,Lp,L2)
	Bar1 = init_Barr(ϵ,L1,0,V1,Vp)
	Puits = init_Puits(ϵ,Lp,Vp)
	Bar2 = init_Barr(ϵ,L2,Vp,V2,Vd)
	N = Bar1.N*Puits.N*Bar2.N
	t = 1/conj(N[1,1])
	r = -conj(N[1,2]*t)
	(t,r)
end

function doubleBarrTsimpleBis(n)
	V = 0.4
	l = 3
	E = [(i/n)*V*eV for i=1:(n-1)]
	T = zeros(n-1)
	Lp = l/3
	#Premiere situation
	for k=1:(n-1)
		(t,r) = ParamDB(E[k],V*eV,0,V*eV,0,l*nm,Lp*nm,l*nm)
		T[k] = real.(t't)
	end
	ϵres = argmax(T) * V/n #L'energie de resonnance en eV
	ϵres = (round(ϵres * 10000))/10000 #On arrondie
	Γ = Largeur(T,0.01)*(V/n) #La largeur à mi-hauteur
	Γ = (round(Γ * 100000))/100000 #On arrondie
	pl1 = plot(E./eV,T,xlabel = "Energie en eV, ϵres = $ϵres eV, Γ = $Γ eV" ,ylabel = "T(E)", title = "T(E), L1=L2=$l nm Lp=$Lp nm, V1=V2=$V eV")
	savefig(pl1,"situation1Bis.pdf")
	#Situation 2 V1 = 2 V2
	for k=1:(n-1)
		(t,r) = ParamDB(E[k],2*V*eV,0,V*eV,0,l*nm,Lp*nm,l*nm)
		T[k] = real.(t't)
	end
	pl1 = plot(E./eV,T,xlabel = "Energie en eV",ylabel = "T(E)", title = "T(E), L1=L2=$l nm Lp=$Lp nm, V2=V1/2=$V eV")
	savefig(pl1,"situation2Bis.pdf")
	#Situation 3, V2 varie, ϵ = 0.4 eV
	V2 = [(eV*V/2) *(1 + 2*i/n) for i = 1:(n-1)]
	for k=1:(n-1)
		(t,r) = ParamDB(ϵres*eV,V*eV,0,V2[k],0,l*nm,Lp*nm,l*nm)
		T[k] = real.(t't)
	end
	pl1 = plot(V2./eV,T,xlabel = "V2 en eV",ylabel = "T(V2)", title = "T(V2), L1=L2=Lp=$l nm, V1=$V eV, E=$ϵres eV")
	savefig(pl1,"situation3Bis.pdf")
end



function doubleBarrTsimple(n)
	V = 0.4
	l = 3
	E = [(i/n)*V*eV for i=1:(n-1)]
	T = zeros(n-1)
	Lp = l/3
	#Premiere situation
	for k=1:(n-1)
		(t,r,A1,B1,α,β,A2,B2) = ParamDBsimple(E[k],V*eV,V*eV,l*nm,l*nm,Lp*nm)
		T[k] = real.(t't)
	end
	ϵres = argmax(T) * V/n #L'energie de resonnance en eV
	ϵres = (round(ϵres * 10000))/10000 #On arrondie
	Γ = Largeur(T,0.01)*(V/n) #La largeur à mi-hauteur
	Γ = (round(Γ * 100000))/100000 #On arrondie
	pl1 = plot(E./eV,T,xlabel = "Energie en eV, ϵres = $ϵres eV, Γ = $Γ eV" ,ylabel = "T(E)", title = "T(E), L1=L2=$l nm Lp=$Lp nm, V1=V2=$V eV")
	savefig(pl1,"situation1.pdf")
	#Situation 2 V1 = 2 V2
	for k=1:(n-1)
		(t,r,A1,B1,α,β,A2,B2) = ParamDBsimple(E[k],2*V*eV,V*eV,l*nm,l*nm,Lp*nm)
		T[k] = real.(t't)
	end
	pl1 = plot(E./eV,T,xlabel = "Energie en eV",ylabel = "T(E)", title = "T(E), L1=L2=$l nm Lp=$Lp nm, V2=V1/2=$V eV")
	savefig(pl1,"situation2.pdf")
	#Situation 3, V2 varie, ϵ = 0.4 eV
	V2 = [(eV*V/2) *(1 + 2*i/n) for i = 1:(n-1)]
	for k=1:(n-1)
		(t,r,A1,B1,α,β,A2,B2) = ParamDBsimple(ϵres*eV,V*eV,V2[k],l*nm,l*nm,Lp*nm)
		T[k] = real.(t't)
	end
	pl1 = plot(V2./eV,T,xlabel = "V2 en eV",ylabel = "T(V2)", title = "T(V2), L1=L2=Lp=$l nm, V1=$V eV, E=$ϵres eV")
	savefig(pl1,"situation3.pdf")
end


#doubleBarrTsimpleBis(10000)


function elemMatrice(B1,P,B2,A)
	(a,b,c,d,α1,β1,e,f) = ParamDBsimple(B1.ϵ,B1.V,B2.V,B1.L,B2.L,P.L)
	(a,b,c,d,α2,β2,e,f) = ParamDBsimple(A.ϵ,A.V,B2.V,A.L,B2.L,P.L)
    quadgk(x-> conj(α1*exp(1im*P.kd*x) + β1*exp(-1im*P.kd*x))*x*(α2*exp(1im*P.kd*x)+β2*exp(-1im*P.kd*x)),B1.L,B1.L + B1.L)
end
