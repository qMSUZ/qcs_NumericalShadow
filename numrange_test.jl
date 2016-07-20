
# x=randn(10) + im*randn(10)
# x=convert(Array{Complex64}, x)

A=[0.5 0.25 ; 0.25 0.5]
A=convert(Array{Complex64}, A)

wRe=[0 0 0 0];
wIm=[0 0 0 0];
wRe=convert(Array{Float32}, wRe)
wIm=convert(Array{Float32}, wIm)

pts=zeros(Complex64, (629,1))

ccall( (:init_enviroment, "./libNumRange.so"), Int64, ())

ccall( (:calc_bounding_box, "./libNumRange.so"), Int64, (Ptr{Complex64}, Int64, Ref{Float32}, Ref{Float32}), A, 2, wRe, wIm )

ccall( (:calc_numerical_range, "./libNumRange.so"), Int64, (Ptr{Complex64}, Int64, Float32, Float32, Int64, Ref{Complex64}), A, 2, 0, 0.01, 629, pts )

ccall( (:finalize_environment, "./libNumRange.so"), Int64, ())

pts

