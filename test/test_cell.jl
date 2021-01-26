@testset "cell" begin
	# 1D crystals
	cell = UnitCell(1, :lamellar, 10.0)
	@test cell.volume ≈ 10.0

	# 2D crystals
	cell = UnitCell(2, :square, 10.0)
	@test cell.volume ≈ 100.0

	cell = UnitCell(2, :rectangular, 10.0, 5.0)
	@test cell.volume ≈ 50.0

	cell = UnitCell(2, :hexagonal, 10.0)
	@test cell.volume ≈ 86.60254037844385

	cell = UnitCell(2, :oblique, 10.0, 5.0, pi/4)
	@test cell.volume ≈ 35.35533905932737

	# 3D crystals
	cell = UnitCell(3, :cubic, 10.0)
	@test cell.volume ≈ 1000.0

	cell = cell = UnitCell(3, :tetragonal, 10.0, 5.0)
	@test cell.volume ≈ 500.0

	cell = UnitCell(3, :orthorhombic, 10.0, 5.0, 2.0)
	@test cell.volume ≈ 100.0

	cell = UnitCell(3, :hexagonal, 10.0, 5.0)
	@test cell.volume ≈ 433.01270189221924

	cell = UnitCell(3, :monoclinic, 10.0, 5.0, 2.0, pi/4)
	@test cell.volume ≈ 70.71067811865474

	cell = UnitCell(3, :triclinic, 10.0, 5.0, 2.0, pi/2, pi/2, pi/2)
	@test cell.volume ≈ 100.0	
end