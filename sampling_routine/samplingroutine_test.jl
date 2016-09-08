# x=randn(10) + im*randn(10)
# x=convert(Array{Complex64}, x)

A=[  0.640405-1.57554im   0.234606+2.42184im   -0.450655-1.04569im    0.931771+1.44655im   	;
	-0.602501-1.16399im   2.65645-2.60768im    -1.52335-0.69291im    -1.84434+2.64027im 	;
	 0.715993-1.06037im   0.323488-1.03543im   -1.24544+0.0185366im  -0.422736+0.596736im	;
	 1.39541-0.618205im  -0.735714-0.345376im   0.17286-2.1607im      0.643438+0.85441im ]

Aorg =  A;

pts = zeros(Complex64, (4,1))

A=convert(Array{Complex64}, A)
A=reshape(A,(16,1))

ccall( (:init_enviroment, "./libSamplingRoutine.so"), Int64, ())

ccall( (:sampling_routine_simpleComplex_float, "./libSamplingRoutine.so"), Int64, (Ptr{Complex64}, Ref{Complex64}, Int64, Int64, Int64), A, pts, 4, 4, 0 )

ccall( (:finalize_environment, "./libSamplingRoutine.so"), Int64, ())

pts

