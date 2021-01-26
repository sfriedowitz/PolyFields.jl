@testset "interactions" begin
	fh = FloryInteraction()
	set_interaction!(fh, 1, 2, 0.5)
	set_interaction!(fh, 2, 1, 0.5)
	@test length(fh.chis) == 1
	@test_throws AssertionError set_interaction!(fh, 1, 1, 0.5) 
end