using Plots

const ħ=6.63e-34
const m = 9.1e-31
const eV = 1.6e-19
const nm = 1e-9

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

function doubleBarrT(n,ϵ,V1,V2,L1,L2,Lp)
	Bar1 = init_Barr(ϵ,L1,V1)
	Bar2 = init_Barr(ϵ,L2,V2)
	N = Bar1.N * Bar2.N
	k = Bar1.k
	t = 1/conj(N[1,1])
	r = N[1,2]*t
	A1 = 1 + r
	B1 =1im*k*(1-r)/Bar1.K
	α = -(1im*exp(-1im*k*L1)/2*k)*((1im*k*A1 + B1*Bar1.K)*cosh(L1*Bar1.K) + (1im*k*B1 + A1*Bar1.K)*sinh(L1*Bar1.K))
	β = (1im*exp(1im*k*L1)/2*k)*((-1im*k*A1 + B1*Bar1.K)*cosh(L1*Bar1.K) + (-1im*k*B1 + A1*Bar1.K)*sinh(L1*Bar1.K))
	L = L1 + L2 + Lp
	A2 = (exp(1im*k*L)*t/Bar2.K)*(cosh(L*Bar2.K)*Bar2.K - 1im*k*sinh(L*Bar2.K))
	B2 = (exp(1im*k*L)*t/Bar2.K)*(cosh(L*Bar2.K)*1im*k - Bar2.K*sinh(L*Bar2.K))
	z = [i/n for i=1:n]
	dl = 2*L/n
	ψ = zeros(ComplexF64,n)
	for i = 1:Int64(n/4)
		ψ[i] = exp(1im*k*2L*z[i]) + r*exp(-1im*k*2L*z[i])
	end
	for i = Int64(n/4):Int64(round(((L/2)+L1)/dl))
		ψ[i] = A1*cosh(Bar1.K*2L*z[i]) + B1*sinh(Bar1.K*2L*z[i])
	end
	for i = Int64(round(((L/2)+L1)/dl)):Int64(round(((L/2)+L1 + Lp)/dl))
		ψ[i] = α*exp(1im*k*2L*z[i]) + β*exp(-1im*k*2L*z[i])
	end
	for i = Int64(round(((L/2)+L1 + Lp)/dl)):Int64(round(((L/2)+L)/dl))
		ψ[i] = A2*cosh(Bar1.K*2L*z[i]) + B2*sinh(Bar1.K*2L*z[i])
	end
	for i = Int64(round(((L/2)+L)/dl)):n
		ψ[i] = t*exp(1im*k*2L*z[i])
	end

	pl = plot(z,abs2.(ψ))
	savefig(pl,"solution.pdf")
end

doubleBarrT(100,0.3*eV,eV,eV,nm,nm,nm)
