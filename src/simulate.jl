#==============================================================================#
# Field initialization
#==============================================================================#

"""
	init_fields!(sys; seed = -1, scale = 1.0)
	init_fields!(sys, fields; interpolate = true)
	init_fields!(sys, fpath; interpolate = true)

Initialize fields for all monomer types in the system.
Fields can be initialized randomly with no provided data,
loaded from a field dictionary, or loaded from file.
"""
function init_fields!(sys::FieldSystem; seed::Integer = -1, scale::Real = 1.0)
	if seed >= 0
		Random.seed!(seed)
	end
	for (mid, omega) in sys.fields
		omega .= scale * randn(sys.dims)
	end
	make_residuals!(sys)
	return nothing
end

function init_fields!(sys::FieldSystem, fields::Dict; interpolate = true)
	for (mid, omega) in fields
		if haskey(sys.monomers, mid)
			if size(omega) == sys.dims
				sys.fields[mid] = omega
			elseif interpolate
				omega_itp = interpolate_grid(omega, sys.dims)
				sys.fields[mid] = omega_itp
			else
				error("Input field dimensions $(size(omega)) does not match system grid $(sys.dims).")
			end
		else
			@warn "Monomer id = $(mid) not present in system -- skipping."
		end
	end	
	make_residuals!(sys)
	return nothing
end

function init_fields!(sys::FieldSystem, fpath::AbstractString; interpolate::Bool = true)
	fields = load_fields(fpath)
	init_fields!(sys, fields; interpolate = interpolate)
	return nothing
end

#==============================================================================#
# Core SCFT simulation routine
#==============================================================================#

"""
	struct SCFTOptions

Input options for an SCFT calculation.

### Constructors
```
SCFTOptions(; nsteps = 1000, nsout = 0, nsave = -1, domain = true, rtol = 1e-5, stol = 1e-5, field_out = "")
```
"""
@with_kw struct SCFTOptions
	nsteps    :: Int = 1000
	nsout     :: Int = 0
	nsave	  :: Int = -1
	field_out :: String = ""
	rtol      :: Float64 = 1e-5
	stol      :: Float64 = 1e-5    
end

"""
	struct SCFTResults

Container for a summary of output from an SCFT calculation.
"""
@with_kw struct SCFTResults
	time      :: Float64
	nsteps    :: Int
	converged :: Bool
	rerr      :: Float64
	serr      :: Float64
	fhelm     :: Float64
end

#==============================================================================#

"""
	scft!(sys, updater; opts = SCFTOptions())

Run an SCFT calculation to relax the configuration of `sys` with the field updater `updater`.
The system is checked with the `validate` method prior to the calculation
to ensure all fields are properly initialized.

Return a summary of the numerical optimization scheme as an `SCFTResults` struct
after the computation has finished.
"""
function scft!(
	sys::FieldSystem, 
	updater::AbstractFieldUpdater,
	domain::Option{DomainUpdater} = nothing
	; 
	opts::SCFTOptions = SCFTOptions()
)
	sim_time = @elapsed begin
		# Validate and initial setup
		validate(sys)
		make_residuals!(sys)
		setup!(updater, sys)
		if !isnothing(domain); setup!(domain, sys); end

		# Iterate SCFT
		step = 0
		rerr = field_error(sys)
		fhelm = free_energy(sys)

		converged = false
		while step < opts.nsteps && !converged
			if opts.nsout > 0 && step % opts.nsout == 0
				@printf("Step (%-5d / %5d): rerr = %8e, fhelm = %.5e\n", step, opts.nsteps, rerr, fhelm)
			end

			# Check for convergence of fields and cell
			fields_converged = false
			cell_converged = false
			if rerr < opts.rtol
				fields_converged = true
				if !isnothing(domain)
					if stress_error(domain) < opts.stol
						cell_converged = true
					end
				else
					cell_converged = true
				end
				converged = fields_converged && cell_converged
				if converged
					@printf("Converged @ step = %d: rerr = %8e, fhelm = %.5e\n", step, rerr, fhelm)
					break
				end
			end

			# If not, take a step with updater
			update_state!(sys, updater)
			rerr = field_error(sys)
			fhelm = free_energy(sys)

			# Mean-subtract the system fields if canonical ensemble
			if sys.ensemble == Canonical
				for omega in values(sys.fields)
					omega .-= mean(omega)
				end
			end

			# Update the cell if necessary
			if !isnothing(domain)
				if step_ready(domain, step, err) && stress_error(sys) > opts.stol
					relax_cell!(domain)
					@printf("Updating cell dimensions: serr = %.4e", stress_error(sys))
					println(", parameters = [" * parameter_string(sys.cell) * "]")
				end
			end

			if isnan(rerr)
				error("Diverging field updater -- breaking.") # Better way to check for this?
			end

			# Check if output should be saved
			if opts.nsave > 0 && step % opts.nsave == 0
				if !isempty(opts.field_out); save_fields(sys, opts.field_out); end
			end

			step += 1
		end

		# Make sure system state is current
		make_residuals!(sys)
		fhelm = free_energy(sys)
		if !isempty(opts.field_out); save_fields(sys, opts.field_out); end
	end
	
	return SCFTResults(sim_time, step, converged, rerr, stress_error(sys), fhelm)
end