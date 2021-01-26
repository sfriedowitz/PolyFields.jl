@testset "updaters" begin
	dims = (32, 1, 1)
	cell = UnitCell(1, :lamellar, 10.0)

	amon = Monomer(; id = 1, vol = 1.0, charge = 0.0)
	bmon = Monomer(; id = 2, vol = 1.0, charge = 0.0)
	chain = Diblock(amon, bmon, 100, 0.5, 1.0, 1.0, 100)

	itx = FloryInteraction()
	set_interaction!(itx, amon.id, bmon.id, 0.2)

	sys = FieldSystem(dims, cell; monomers = [amon, bmon], ensemble = Canonical)
	add_species!(sys, chain, 1.0)
	add_interaction!(sys, itx)

	# Create some updaters and see what happens
	fieldinit!(sys; scale = 0.1, seed = 12345)

	u1 = updater_with_system(EulerUpdater, sys; lam = 0.005, method = :PECE)
	u2 = updater_with_system(MomentumUpdater, sys; method = :NAG, lam = 0.005)
	u3 = updater_with_system(AndersonUpdater, sys; nhist = 5)

	for i = 1:10; step!(u1); end
	@test free_energy(sys) ≈ 0.0626842750325202

	for i = 1:10; step!(u2); end
	# @test free_energy(sys) ≈ 0.06172669754976512

	for i = 1:10; step!(u3); end
	@test free_energy(sys) ≈ 0.12798422498635653
end