#==============================================================================#
# SCFT Simulation
#==============================================================================#

"""
	struct SCFTOptions

	SCFTOptions(steps = 1000, nsout = 0, nsave = 0, ftol = 1e-5, stol = 1e-5, domain = false, savepath = "")

Input options for an SCFT calculation.
"""
@with_kw struct SCFTOptions
	nsteps    :: Int = 100
	nsout     :: Int = 0
	nsave	  :: Int = 0
	ftol      :: Float64 = 1e-5
	stol      :: Float64 = 1e-5
	domain    :: Bool = false
	savepath  :: String = ""
end

"""
	struct SCFTResults

Container for output summary from an SCFT calculation.
"""
@with_kw struct SCFTResults
	time      :: Float64
	nsteps    :: Int
	converged :: Bool
	ferr      :: Float64
	serr      :: Float64
	fhelm     :: Float64
end

#==============================================================================#

"""
	scft!(sys, updater; opts = SCFTOptions())

Run an SCFT calculation to relax the system configuration with the provided field updater.
The system is checked with the `validate` method prior to the calculation
to ensure all fields are properly initialized.

Return a summary of the numerical scheme as a result.
"""
function scft!(sys::FieldSystem, updater::AbstractFieldUpdater; opts::SCFTOptions = SCFTOptions())
	sim_time = @elapsed begin
		# Validate and initial setup
		validate(sys)
		residuals!(sys)
		setup!(updater, sys)
		#if !isnothing(domain); setup!(domain, sys); end

		# Iterate SCFT
		step = 0
		ferr = field_error(sys)
		fhelm = free_energy(sys)

		converged = false
		while step < opts.nsteps && !converged
			if opts.nsout > 0 && step % opts.nsout == 0
				@printf("Step (%-5d / %5d): ferr = %8e, fhelm = %.5e\n", step, opts.nsteps, ferr, fhelm)
			end

			# Check for convergence of fields and cell
			fconverged = false
			sconverged = false
			if ferr < opts.ftol
				fconverged = true
				converged = fconverged

				# if !isnothing(domain)
				# 	if stress_error(domain) < opts.stol
				# 		cell_converged = true
				# 	end
				# else
				# 	cell_converged = true
				# end
				# converged = fields_converged && cell_converged
				if converged
					@printf("Converged @ step = %d: ferr = %8e, fhelm = %.5e\n", step, ferr, fhelm)
					break
				end
			end

			# If not, take a step with updater
			step!(updater)
			ferr = field_error(sys)
			fhelm = free_energy(sys)

			# # Update the cell if necessary
			# if !isnothing(domain)
			# 	if step_ready(domain, step, err) && stress_error(sys) > opts.stol
			# 		relax_cell!(domain)
			# 		@printf("Updating cell dimensions: serr = %.4e", stress_error(sys))
			# 		println(", parameters = [" * parameter_string(sys.cell) * "]")
			# 	end
			# end

			if isnan(ferr)
				error("Diverging field updater -- breaking.") # Better way to check for this?
			end

			# Check if output should be saved
			if opts.nsave > 0 && step % opts.nsave == 0
				if !isempty(opts.savepath)
					savefields(sys, opts.savepath)
				end
			end

			step += 1
		end

		# Make sure system state is current
		residuals!(sys)
		fhelm = free_energy(sys)
		if !isempty(opts.savepath)
			save_fields(sys, opts.field_out)
		end
	end
	
	return SCFTResults(sim_time, step, converged, ferr, 0.0, fhelm)
end