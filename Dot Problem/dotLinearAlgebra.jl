using LinearAlgebra: dot

foo(n) = randn(n) + im*randn(n)

dot(foo(10^7),foo(10^7))
