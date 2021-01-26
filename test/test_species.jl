@testset "species" begin
	amon = Monomer(; id = 1, vol = 2.0)
	bmon = Monomer(; id = 2, vol = 1.0)

	point = Point(amon)
	@test point.Nref ≈ 2.0

	chain = Homopolymer(bmon, 100, 1.0, 200)
	@test chain.Nref ≈ 100.0

	chain = Diblock(amon, bmon, 100, 0.25, 1.0, 1.0, 200)
	@test chain.f_block[1] ≈ 0.4
	@test chain.Nref_block[2] ≈ 75.0

	chain = Multiblock([amon, amon, bmon, amon, bmon, bmon], [10, 10, 10, 10, 10, 10], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 60)
	@test monomer_fraction(chain, 1) ≈ 2/3
end