@testset "system" begin
	dims = (32, 1, 1)
	cell = UnitCell(1, :lamellar, 10.0)

	amon = Monomer(; id = 1, vol = 1.0, charge = 0.0)
	bmon = Monomer(; id = 2, vol = 1.0, charge = 0.0)
	cmon = Monomer(; id = 3, vol = 1.0, charge = 0.0)
	point = Point(cmon)
	chain = Diblock(amon, bmon, 100, 0.5, 1.0, 1.0, 100)

	itx = FloryInteraction()
	set_interaction!(itx, amon.id, bmon.id, 0.2)

	sys = FieldSystem(dims, cell; monomers = [amon, bmon, cmon], ensemble = Canonical)
	add_species!(sys, point, 0.75)
	add_species!(sys, chain, 0.25)
	add_interaction!(sys, itx)

	# Checks many initialization fields
	@test_throws ErrorException validate(sys)

	# Check actual system usage
	fieldinit!(sys; scale = 0.1, seed = 12345)

	# Species density operators
	density!(point)
	density!(chain)

	@test point.Q ≈ 1.004177455788938
	@test chain.Q ≈ 8.24300700515971

	@test free_energy(sys) ≈ -0.1681580628744596
	@test free_energy_bulk(sys) ≈ -0.9686022902416352
	@test scfstress(sys)[1] ≈ -0.0004272750367130805
end