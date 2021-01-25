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

function isconverged(ferr, ftol, serr, stol, vcell = false)
	fconverged = false
	sconverged = true
	if ferr < ftol
		fconverged = true
		if vcell && serr < stol
			sconverged = true
		else
			sconverged = false
		end
	end
	return fconverged && sconverged
end

"""
	scft!(sys, updater; domain = nothing, opts = SCFTOptions())

Run an SCFT calculation to relax the system configuration with the provided field updater.
The system is checked with the `validate` method prior to the calculation
to ensure all fields are properly initialized.

Returns a summary in the form of an `SCFTResults` object.
"""
function scft!(sys::FieldSystem, updater::AbstractFieldUpdater; 
	domain::Option{DomainUpdater} = nothing, opts::SCFTOptions = SCFTOptions())
	sim_time = @elapsed begin
		# Validate and initial setup
		validate(sys)
		residuals!(sys)
		setup!(updater, sys)

		# Iterate SCFT
		step = 0
		ferr = field_error(sys)
		fhelm = free_energy(sys)

		serr = 0.0
		vcell = !isnothing(domain)
		if vcell
			setup!(domain, sys)
			serr = stress_error(domain)
		end

		converged = isconverged(ferr, opts.ftol, serr, opts.stol, vcell)
		while step < opts.nsteps && !converged
			if opts.nsout > 0 && step % opts.nsout == 0
				@printf("Step (%-5d / %5d): ferr = %8e, fhelm = %.5e\n", step, opts.nsteps, ferr, fhelm)
			end

			# Update the fields
			step!(updater)

			# Update the cell if necessary
			if vcell && step % domain.nskip == 0
				scfstress!(domain)
				if stress_error(domain) > opts.stol
					step!(domain; nsteps = domain.nper, tol = opts.stol)
					serr = stress_error(domain)
					@printf("Updating cell dimensions: serr = %.5e\n", serr)
				end
			end

			# Track error
			ferr = field_error(sys)
			fhelm = free_energy(sys)

			if isnan(ferr)
				error("Diverging field updater -- breaking.") # Better way to check for this?
			end

			# Check for convergence of fields and cell
			converged = isconverged(ferr, opts.ftol, serr, opts.stol, vcell)
			if converged
				@printf("Converged @ step = %d: ferr = %8e, fhelm = %.5e\n", step, ferr, fhelm)
				if vcell
					@printf("Converged cell parameters: %s\n", sys.cell.params)
				end
				break
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
		ferr = field_error(sys)
		fhelm = free_energy(sys)
		
		if vcell
			serr = stress_error(domain)
		end

		if !isempty(opts.savepath)
			savefields(sys, opts.field_out)
		end
	end
	
	return SCFTResults(sim_time, step, converged, ferr, serr, fhelm)
end