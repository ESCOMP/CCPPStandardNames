# CCPP Standard Name Library
#### Table of Contents
* [dimensions](#dimensions)
* [constants](#constants)
* [coordinates](#coordinates)
* [state_variables](#state_variables)
* [diagnostics](#diagnostics)
* [constituents](#constituents)
* [standard_variables](#standard_variables)
* [non-category](#non-category)

## dimensions
Dimension standard names may come in sets of six related standard names for each dimension:
```
[dim_name]_dimension -- The full dimension size
[dim_name]_loop_extent -- Size of dim for current call
[dim_name]_begin - Start index for dimension
[dim_name]_end - End index for dimension
[dim_name]_index - Single index for dimension
[dim_name]_selection - Array of selected indices for dimension
```
Note that the cap generator may substitute among standard names in this category in order to properly call suite parts and individual schemes. In the substitutions below, the name on the left is the standard_name in the dimensions field of the caller while the name(s) on the right is (are) the standard name(s) of the callee (in the form used in the subroutine call).
```
[dim_name]_dimension ==> 1:[dim_name]_loop_extent
[dim_name]_loop_extent ==> 1:[dim_name]_loop_extent
[dim_name]_begin:[dim_name]_end ==> 1:[dim_name]_loop_extent
[dim_name]_begin:[dim_name]_end ==> 1:[dim_name]_dimension
```
Also note that horizontal_dimension should be used in xxx_[timestep_]init and xxx_[timestep_]final routines but not in xxx_run routines.
Currently, the only dimension which supports all six dimension types is horizontal_dimension. This and other supported dimension standard names are listed below.
* `horizontal_dimension`: Size horizontal dimension
    * `integer`: units = count
* `vertical_dimension`: number of vertical levels
    * `integer`: units = count
* `vertical_dimension_plus_one`: number of vertical levels + 1
    * `integer`: units = count
* `vertical_layer_index`: index of a particular vertical level
    * `integer`: units = count
* `vertical_interface_index`: index of a particular vertical interface
    * `integer`: units = count
* `vertical_index_at_surface_adjacent_layer`: Vertical index at surface adjacent layer
    * `integer`: units = index
* `vertical_index_at_top_adjacent_layer`: Vertical index at top adjacent layer
    * `integer`: units = index
* `vertical_index_at_surface_interface`: Vertical index at surface interface
    * `integer`: units = index
* `vertical_index_at_top_interface`: Vertical index at top interface
    * `integer`: units = index
* `number_of_openmp_threads`: Total number of thread blocks which the host model may use to call CCPP physics run groups during the CCPP run phase.
    * `integer`: units = none
* `ccpp_thread_number`: Number of current thread block. This variable may only be used during CCPP run phase
    * `integer`: units = none
## constants
* `avogadro_number`: Avogadro number
    * `real(kind=kind_phys)`: units = molecules mole-1
* `reference_air_pressure_for_atmosphere_vertical_coordinate`: Reference air pressure for atmosphere vertical coordinate
    * `real(kind=kind_phys)`: units = Pa
* `boltzmann_constant`: Boltzmann constant
    * `real(kind=kind_phys)`: units = J K-1
* `gas_constant_of_dry_air`: Gas constant of dry air
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `seconds_in_calendar_day`: Seconds in calendar day
    * `integer(kind=kind_phys)`: units = s
* `specific_heat_of_dry_air_at_constant_pressure`: Specific heat of dry air at constant pressure
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `specific_heat_of_liquid_water_at_20c`: specific heat of liquid water at 20c
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `latent_heat_of_vaporization_of_water_at_0c`: latent heat of vaporization of water at 0c
    * `real(kind=kind_phys)`: units = J kg-1
* `dry_air_density_at_stp`: density of dry air at STP
    * `real(kind=kind_phys)`: units = kg m-3
* `fresh_liquid_water_density_at_0c`: density of liquid water at 0c
    * `real(kind=kind_phys)`: units = kg m-3
* `ratio_of_water_vapor_to_dry_air_gas_constants_minus_one`: (Rwv / Rdair) - 1.0
    * `real(kind=kind_phys)`: units = 1
## coordinates
* `latitude`: Latitude
    * `real(kind=kind_phys)`: units = degree_north
* `longitude`: Longitude
    * `real(kind=kind_phys)`: units = degree_east
* `gravitational_acceleration`: Gravitational acceleration
    * `real(kind=kind_phys)`: units = m s-2
* `cell_area`: Cell area
    * `real(kind=kind_phys)`: units = steradian
* `cell_weight`: Cell weight
    * `real(kind=kind_phys)`: units = none
## state_variables
Note that appending '_on_previous_timestep' to standard_names in this section yields another valid standard_name
* `physics_state_due_to_dynamics`: Physics state due to dynamics
    * `physics_state(kind=kind_phys)`: units = none
* `timestep_for_physics`: Timestep for physics
    * `integer(kind=kind_phys)`: units = s
* `total_tendency_of_physics`: Total tendency of physics
    * `physics_tend(kind=kind_phys)`: units = none
* `air_pressure_at_top_of_atmosphere_model`: Air pressure at top of atmosphere model
    * `real(kind=kind_phys)`: units = Pa
* `air_pressure_at_sea_level`: Air pressure at sea level
    * `real(kind=kind_phys)`: units = Pa
* `surface_air_pressure`: Surface air pressure
    * `real(kind=kind_phys)`: units = Pa
* `surface_dry_air_pressure`: Surface dry air pressure
    * `real(kind=kind_phys)`: units = Pa
* `surface_geopotential`: Surface geopotential
    * `real(kind=kind_phys)`: units = m2 s-2
* `air_temperature`: Air temperature
    * `real(kind=kind_phys)`: units = K
* `air_temperature_on_previous_timestep`: Air temperature on previous timestep
    * `real(kind=kind_phys)`: units = K
* `x_wind`: Horizontal wind in a direction perdendicular to y_wind
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind`: Horizontal wind in a direction perdendicular to x_wind
    * `real(kind=kind_phys)`: units = m s-1
* `dry_static_energy`: Dry static energy Content of Atmosphere Layer
    * `real(kind=kind_phys)`: units = J kg-1
* `flag_for_lagrangian_vertical_coordinate`: flag indicating if vertical coordinate is lagrangian
    * `logical(kind=)`: units = flag
* `lagrangian_tendency_of_air_pressure`: Vertical pressure velocity
    * `real(kind=kind_phys)`: units = Pa s-1
* `air_pressure`: Midpoint air pressure
    * `real(kind=kind_phys)`: units = Pa
* `dry_air_pressure`: Dry midpoint pressure
    * `real(kind=kind_phys)`: units = Pa
* `air_pressure_thickness`: Air pressure thickness
    * `real(kind=kind_phys)`: units = Pa
* `dry_air_pressure_thickness`: Dry air pressure thickness
    * `real(kind=kind_phys)`: units = Pa
* `reciprocal_of_air_pressure_thickness`: Reciprocal of air pressure thickness
    * `real(kind=kind_phys)`: units = Pa-1
* `reciprocal_of_dry_air_pressure_thickness`: Reciprocal of dry air pressure thickness
    * `real(kind=kind_phys)`: units = Pa-1
* `ln_air_pressure`: Ln air pressure
    * `real(kind=kind_phys)`: units = 1
* `ln_dry_air_pressure`: Ln dry air pressure
    * `real(kind=kind_phys)`: units = 1
* `reciprocal_of_dimensionless_exner_function_wrt_surface_air_pressure`: inverse exner function w.r.t. surface pressure, (ps/p)^(R/cp)
    * `real(kind=kind_phys)`: units = 1
* `geopotential_height`: Geopotential height
    * `real(kind=kind_phys)`: units = m
* `constituent_mixing_ratio`: Constituent mixing ratio
    * `real(kind=kind_phys)`: units = kg/kg moist or dry air depending on type
* `air_pressure_at_interface`: Air pressure at interface
    * `real(kind=kind_phys)`: units = Pa
* `dry_air_pressure_at_interface`: Dry air pressure at interface
    * `real(kind=kind_phys)`: units = Pa
* `ln_air_pressure_at_interface`: Ln air pressure at interface
    * `real(kind=kind_phys)`: units = 1
* `ln_dry_air_pressure_at_interface`: Ln dry air pressure at interface
    * `real(kind=kind_phys)`: units = 1
* `largest_model_top_pressure_that_allows_molecular_diffusion`: Largest model top pressure that allows molecular diffusion
    * `real(kind=kind_phys)`: units = Pa
* `flag_for_molecular_diffusion`: Flag for molecular diffusion
    * `logical(kind=kind_phys)`: units = flag
* `flag_for_physics_grid_initialization`: Flag to indicate if physics grid is initialized
    * `logical(kind=kind_phys)`: units = flag
* `geopotential_height_at_interface`: Geopotential height at interface
    * `real(kind=kind_phys)`: units = m
* `vertically_integrated_total_energy_of_initial_state`: Vertically integrated total energy of initial state
    * `real(kind=kind_phys)`: units = J m-2
* `vertically_integrated_total_energy_of_current_state`: Vertically integrated total energy of current state
    * `real(kind=kind_phys)`: units = J m-2
* `vertically_integrated_total_water_of_initial_state`: Vertically integrated total water of initial state
    * `real(kind=kind_phys)`: units = kg m-2
* `vertically_integrated_total_water_of_current_state`: Vertically integrated total water of current state
    * `real(kind=kind_phys)`: units = kg m-2
* `tendency_of_air_temperature_due_to_X`: Change in temperature from a parameterization
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_air_temperature_due_to_model_physics`: Total change in temperature from a                               physics suite
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_air_potential_temperature_due_to_X`: Change in potential temperature from a parameterization
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_air_potential_temperature_due_to_model_physics`: Tendency of air potential temperature due to model physics
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_x_wind_due_to_X`: Change in x wind from a parameterization
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_x_wind_due_to_model_physics`: Tendency of x wind due to model physics
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_wind_due_to_X`: Change in y wind from a parameterization
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_wind_due_to_model_physics`: Tendency of y wind due to model physics
    * `real(kind=kind_phys)`: units = m s-2
* `surface_upward_heat_flux_in_air`: Surface upward heat flux in air
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_boundary_flux_of_total_energy`: Cumulative boundary flux of total energy
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_boundary_flux_of_total_water`: Cumulative boundary flux of total water
    * `real(kind=kind_phys)`: units = W m-2
* `reference_air_pressure`: reference pressure 
    * `real(kind=kind_phys)`: units = Pa
* `US_standard_air_pressure_at_sea_level`: US Standard Atmospheric pressure at sea level
    * `real(kind=kind_phys)`: units = Pa
* `surface_reference_air_pressure`: reference pressure at surface
    * `real(kind=kind_phys)`: units = Pa
* `reference_air_pressure_normalized_by_surface_air_pressure`: reference pressure normalized by surface pressure
    * `real(kind=kind_phys)`: units = 1
* `dimensionless_exner_function`: exner function
    * `real(kind=kind_phys)`: units = 1
* `air_potential_temperature`: air potential temperature
    * `real(kind=kind_phys)`: units = K
* `air_potential_temperature_on_previous_timestep`: air potential temperature on previous timestep
    * `real(kind=kind_phys)`: units = K
* `pressure_dependent_gas_constant_of_dry_air`: Pressure dependent gas constant of dry air
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `pressure_dependent_ratio_of_dry_air_to_water_vapor_gas_constants_minus_one`: (Rwv / Rdair) - 1.0
    * `real(kind=kind_phys)`: units = 1
## diagnostics
* `total_precipitation_rate_at_surface`: Total precipitation rate at surface
    * `real(kind=kind_phys)`: units = m s-1
## constituents
* `number_of_chemical_species`: Number of chemical species
    * `integer(kind=kind_phys)`: units = count
* `number_of_tracers`: Number of tracers
    * `integer(kind=kind_phys)`: units = count
* `specific_humidity`: Specific humidity
    * `real(kind=kind_phys)`: units = kg kg-1
* `mole_fraction_of_water_vapor`: Mole fraction of water vapor
    * `real(kind=kind_phys)`: units = mole mole-1
* `cloud_liquid_water_mixing_ratio_of_moist_air`: Cloud liquid water mixing ratio of moist air
    * `real(kind=kind_phys)`: units = kg kg-1
* `cloud_liquid_water_mixing_ratio`: Cloud liquid water mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `cloud_ice_mixing_ratio`: Ratio of the mass of ice to the mass of dry air
    * `real(kind=kind_phys)`: units = kg kg-1
* `rain_water_mixing_ratio`: Rain water mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_ch4`: CH4 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_co`: CO volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_co2`: CO2 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_ccl4`: CCL4 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_cfc11`: CFC11 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_cfc12`: CFC12 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_cfc113`: CFC113 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_cfc22`: CFC22 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_o2`: O2 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_of_n2o`: N2O volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
## standard_variables
Standard / required CCPP variables
* `ccpp_error_message`: Error message for error handling in CCPP
    * `character(kind=len=512)`: units = 1
* `ccpp_error_flag`: Error flag for error handling in CCPP
    * `integer(kind=)`: units = flag
## non-category
* `control_for_scale_aware_tke_moist_edmf_pbl_scheme`: Control for scale aware tke moist edmf pbl scheme
    * `integer`: units = none
* `identifier_for_2018_scale_aware_tke_moist_edmf_pbl`: Identifier for 2018 scale aware tke moist edmf pbl
    * `integer`: units = none
* `horizontal_loop_extent`: Horizontal loop extent
    * `integer`: units = count
* `number_of_vertical_diffusion_tracers`: Number of vertical diffusion tracers
    * `integer`: units = count
* `index_of_cloud_liquid_water_mixing_ratio_in_tracer_concentration_array`: Index of cloud liquid water mixing ratio in tracer concentration array
    * `integer`: units = index
* `index_for_ice_cloud_condensate_vertical_diffusion_tracer`: Index for ice cloud condensate vertical diffusion tracer
    * `integer`: units = index
* `index_for_turbulent_kinetic_energy_vertical_diffusion_tracer`: Index for turbulent kinetic energy vertical diffusion tracer
    * `integer`: units = index
* `vertical_index_at_top_of_atmosphere_boundary_layer`: Vertical index at top of atmosphere boundary layer
    * `integer`: units = index
* `flag_TKE_dissipation_heating`: Flag TKE dissipation heating
    * `logical`: units = flag
* `index_of_highest_temperature_inversion`: Index of highest temperature inversion
    * `integer`: units = index
* `flag_for_generic_tendency_due_to_plantary_boundary_layer`: Flag for generic tendency due to plantary boundary layer
    * `logical`: units = flag
* `flag_for_xyz_dimensioned_diagnostics`: Flag for xyz dimensioned diagnostics
    * `logical`: units = flag
* `flag_for_tracer_xyz_dimensioned_diagnostics`: Flag for tracer xyz dimensioned diagnostics
    * `logical`: units = flag
* `number_of_hydrometeors`: Number of hydrometeors
    * `integer`: units = count
* `flag_for_arakawa_wu_adjustment`: Flag for arakawa wu adjustment
    * `logical`: units = flag
* `number_of_tracers_plus_one`: Number of tracers plus one
    * `integer`: units = count
* `number_of_tracers_for_convective_transport`: Number of tracers for convective transport
    * `integer`: units = count
* `number_of_tracers_for_CS`: Number of tracers for CS
    * `integer`: units = count
* `number_of_cloud_types_in_chikira_sugiyama_deep_convection`: Number of cloud types in chikira sugiyama deep convection
    * `integer`: units = count
* `latitude_index_in_debug_printouts`: Latitude index in debug printouts
    * `integer`: units = index
* `index_of_time_step`: Index of time step
    * `integer`: units = index
* `mpi_rank`: Mpi rank
    * `integer`: units = index
* `flag_for_arakawa_wu_downdrafts_for_deep_convection`: Flag for arakawa wu downdrafts for deep convection
    * `logical`: units = flag
* `flag_for_flux_form_in_chikira_sugiyama_deep_convection_scheme`: Flag for flux form in chikira sugiyama deep convection scheme
    * `logical`: units = flag
* `flag_print`: Flag print
    * `logical`: units = flag
* `horizontal_index_of_printed_column`: Horizontal index of printed column
    * `integer`: units = index
* `flag_for_occurrence_of_deep_convection_in_column`: Flag for occurrence of deep convection in column
    * `integer`: units = flag
* `control_for_microphysics_scheme`: Control for microphysics scheme
    * `integer`: units = flag
* `number_of_perturbed_land_surface_variables`: Number of perturbed land surface variables
    * `integer`: units = count
* `control_for_stochastic_land_surface_perturbation`: Control for stochastic land surface perturbation
    * `integer`: units = index
* `flag_for_calling_shortwave_radiation`: Flag for calling shortwave radiation
    * `logical`: units = flag
* `number_of_daytime_points`: Number of daytime points
    * `integer`: units = count
* `daytime_points`: Daytime points
    * `integer`: units = index
* `control_for_land_surface_scheme`: Control for land surface scheme
    * `integer`: units = flag
* `identifier_for_noah_wrfv4_land_surface_scheme`: Identifier for noah wrfv4 land surface scheme
    * `integer`: units = flag
* `vertical_dimension_of_soil`: Vertical dimension of soil
    * `integer`: units = count
* `flag_for_noah_lsm_ua_extension`: Flag for noah lsm ua extension
    * `logical`: units = flag
* `control_for_flux_adjusting_surface_data_assimilation_system`: Control for flux adjusting surface data assimilation system
    * `integer`: units = flag
* `flag_for_restart`: Flag for restart
    * `logical`: units = flag
* `index_of_ice_vegetation_category`: Index of ice vegetation category
    * `integer`: units = index
* `flag_for_calling_land_surface_model`: Flag for calling land surface model
    * `logical`: units = flag
* `flag_for_calling_land_surface_model_glacier`: Flag for calling land surface model glacier
    * `logical`: units = flag
* `index_of_urban_vegetation_category`: Index of urban vegetation category
    * `integer`: units = index
* `flag_for_reading_leaf_area_index_from_input`: Flag for reading leaf area index from input
    * `logical`: units = flag
* `flag_for_reading_surface_diffused_shortwave_albedo_from_input`: Flag for reading surface diffused shortwave albedo from input
    * `logical`: units = flag
* `vegetation_type_classification`: Vegetation type classification
    * `integer`: units = index
* `soil_type_classification`: Soil type classification
    * `integer`: units = index
* `surface_slope_classification`: Surface slope classification
    * `integer`: units = index
* `control_for_land_surface_scheme_thermal_conductivity_option`: Control for land surface scheme thermal conductivity option
    * `integer`: units = index
* `mpi_root`: Mpi root
    * `integer`: units = index
* `iounit_of_namelist`: Iounit of namelist
    * `integer`: units = none
* `iounit_of_log`: Iounit of log
    * `integer`: units = none
* `number_of_equatorial_longitude_points`: Number of equatorial longitude points
    * `integer`: units = count
* `number_of_latitude_points`: Number of latitude points
    * `integer`: units = count
* `flag_for_unified_gravity_wave_physics_gravity_wave_drag_scheme`: Flag for unified gravity wave physics gravity wave drag scheme
    * `logical`: units = flag
* `number_of_statistical_measures_of_subgrid_orography`: Number of statistical measures of subgrid orography
    * `integer`: units = count
* `flag_for_turbulent_orographic_form_drag_in_unified_gravity_wave_physics_gravitiy_wave_drag_scheme`: Flag for turbulent orographic form drag in unified gravity wave physics gravitiy wave drag scheme
    * `logical`: units = flag
* `flag_for_unified_gravity_wave_physics_diagnostics`: Flag for unified gravity wave physics diagnostics
    * `logical`: units = flag
* `index_of_turbulent_kinetic_energy_in_tracer_concentration_array`: Index of turbulent kinetic energy in tracer concentration array
    * `integer`: units = index
* `flag_for_diagnostics`: Flag for diagnostics
    * `logical`: units = flag
* `flag_for_generic_tendency_due_to_gravity_wave_drag`: Flag for generic tendency due to gravity wave drag
    * `logical`: units = flag
* `flag_for_mellor_yamada_janjic_surface_layer_scheme`: Flag for mellor yamada janjic surface layer scheme
    * `logical`: units = flag
* `index_of_cloud_ice_mixing_ratio_in_tracer_concentration_array`: Index of cloud ice mixing ratio in tracer concentration array
    * `integer`: units = index
* `index_of_rain_mixing_ratio_in_tracer_concentration_array`: Index of rain mixing ratio in tracer concentration array
    * `integer`: units = index
* `index_of_snow_mixing_ratio_in_tracer_concentration_array`: Index of snow mixing ratio in tracer concentration array
    * `integer`: units = index
* `index_of_graupel_mixing_ratio_in_tracer_concentration_array`: Index of graupel mixing ratio in tracer concentration array
    * `integer`: units = index
* `control_for_lake_surface_scheme`: Control for lake surface scheme
    * `integer`: units = flag
* `flag_for_fractional_landmask`: Flag for fractional landmask
    * `logical`: units = flag
* `flag_for_cice`: Flag for cice
    * `logical`: units = flag
* `flag_for_surface_flux_coupling`: Flag for surface flux coupling
    * `logical`: units = flag
* `flag_for_one_way_ocean_wave_coupling_to_atmosphere`: Flag for one way ocean wave coupling to atmosphere
    * `logical`: units = flag
* `flag_for_nonzero_land_surface_fraction`: Flag for nonzero land surface fraction
    * `logical`: units = flag
* `flag_for_nonzero_sea_ice_surface_fraction`: Flag for nonzero sea ice surface fraction
    * `logical`: units = flag
* `flag_for_nonzero_lake_surface_fraction`: Flag for nonzero lake surface fraction
    * `logical`: units = flag
* `flag_for_nonzero_ocean_surface_fraction`: Flag for nonzero ocean surface fraction
    * `logical`: units = flag
* `flag_for_nonzero_wet_surface_fraction`: Flag for nonzero wet surface fraction
    * `logical`: units = flag
* `sea_land_ice_mask`: Sea land ice mask
    * `integer`: units = flag
* `sea_land_ice_mask_cice`: Sea land ice mask cice
    * `integer`: units = flag
* `vertical_dimension_of_sea_ice`: Vertical dimension of sea ice
    * `integer`: units = count
* `flag_for_first_time_step`: Flag for first time step
    * `logical`: units = flag
* `flag_for_saturation_adjustment_for_microphysics_in_dynamics`: Flag for saturation adjustment for microphysics in dynamics
    * `logical`: units = none
* `top_layer_index_for_fast_physics`: Top layer index for fast physics
    * `integer`: units = index
* `number_of_water_species`: Number of water species
    * `integer`: units = count
* `number_of_gases_for_multi_gases_physics`: Number of gases for multi gases physics
    * `integer`: units = count
* `mpi_rank_for_fast_physics`: Mpi rank for fast physics
    * `integer`: units = index
* `mpi_root_for_fast_physics`: Mpi root for fast physics
    * `integer`: units = index
* `starting_x_direction_index`: Starting x direction index
    * `integer`: units = count
* `ending_x_direction_index`: Ending x direction index
    * `integer`: units = count
* `starting_x_direction_index_domain`: Starting x direction index domain
    * `integer`: units = count
* `ending_x_direction_index_domain`: Ending x direction index domain
    * `integer`: units = count
* `vertical_dimension_for_fast_physics`: Vertical dimension for fast physics
    * `integer`: units = count
* `vertical_dimension_for_thickness_at_Lagrangian_surface`: Vertical dimension for thickness at Lagrangian surface
    * `integer`: units = count
* `starting_y_direction_index`: Starting y direction index
    * `integer`: units = count
* `ending_y_direction_index`: Ending y direction index
    * `integer`: units = count
* `starting_y_direction_index_domain`: Starting y direction index domain
    * `integer`: units = count
* `ending_y_direction_index_domain`: Ending y direction index domain
    * `integer`: units = count
* `number_of_ghost_zones`: Number of ghost zones
    * `integer`: units = count
* `flag_for_hydrostatic_solver_for_fast_physics`: Flag for hydrostatic solver for fast physics
    * `logical`: units = flag
* `flag_for_fast_microphysics_energy_conservation`: Flag for fast microphysics energy conservation
    * `logical`: units = flag
* `flag_for_tendency_of_air_temperature_at_Lagrangian_surface`: Flag for tendency of air temperature at Lagrangian surface
    * `logical`: units = flag
* `flag_for_the_last_step_of_k_split_remapping`: Flag for the last step of k split remapping
    * `logical`: units = flag
* `flag_for_inline_cloud_fraction_calculation`: Flag for inline cloud fraction calculation
    * `logical`: units = flag
* `omp_threads_for_fast_physics`: Omp threads for fast physics
    * `integer`: units = count
* `flag_for_using_rrtmg_cloud_optics`: Flag for using rrtmg cloud optics
    * `logical`: units = flag
* `flag_for_using_rrtmgp_cloud_optics_with_pade_approximation`: Flag for using rrtmgp cloud optics with pade approximation
    * `logical`: units = flag
* `flag_for_using_rrtmgp_cloud_optics_look_up_table`: Flag for using rrtmgp cloud optics look up table
    * `logical`: units = flag
* `number_of_ice_roughness_categories`: Number of ice roughness categories
    * `integer`: units = count
* `mpi_communicator`: Mpi communicator
    * `integer`: units = index
* `flag_for_calling_longwave_radiation`: Flag for calling longwave radiation
    * `logical`: units = flag
* `flag_for_optical_property_for_liquid_clouds_for_longwave_radiation`: Flag for optical property for liquid clouds for longwave radiation
    * `integer`: units = flag
* `flag_for_optical_property_for_ice_clouds_for_longwave_radiation`: Flag for optical property for ice clouds for longwave radiation
    * `integer`: units = flag
* `number_of_longwave_bands`: Number of longwave bands
    * `integer`: units = count
* `longwave_optical_properties_for_cloudy_atmosphere_by_band`: Longwave optical properties for cloudy atmosphere by band
    * `ty_optical_props_2str`: units = DDT
* `longwave_optical_properties_for_precipitation_by_band`: Longwave optical properties for precipitation by band
    * `ty_optical_props_2str`: units = DDT
* `index_of_ozone_mixing_ratio_in_tracer_concentration_array`: Index of ozone mixing ratio in tracer concentration array
    * `integer`: units = index
* `flag_for_stratospheric_water_vapor_physics`: Flag for stratospheric water vapor physics
    * `logical`: units = flag
* `flag_for_prescribed_aerosols`: Flag for prescribed aerosols
    * `logical`: units = flag
* `control_for_ice_cloud_condensation_nuclei_forcing`: Control for ice cloud condensation nuclei forcing
    * `integer`: units = none
* `control_for_vertical_index_direction`: Control for vertical index direction
    * `integer`: units = flag
* `number_of_x_points_for_current_mpi_rank`: Number of x points for current mpi rank
    * `integer`: units = count
* `number_of_y_points_for_current_mpi_rank`: Number of y points for current mpi rank
    * `integer`: units = count
* `date_and_time_at_model_initialization_in_United_States_order`: Date and time at model initialization in United States order
    * `integer`: units = none
* `lower_latitude_index_of_ozone_forcing_for_interpolation`: Lower latitude index of ozone forcing for interpolation
    * `integer`: units = index
* `upper_latitude_index_of_ozone_forcing_for_interpolation`: Upper latitude index of ozone forcing for interpolation
    * `integer`: units = index
* `lower_latitude_index_of_stratospheric_water_vapor_forcing_for_interpolation`: Lower latitude index of stratospheric water vapor forcing for interpolation
    * `integer`: units = index
* `upper_latitude_index_of_stratospheric_water_vapor_forcing_for_interpolation`: Upper latitude index of stratospheric water vapor forcing for interpolation
    * `integer`: units = index
* `lower_latitude_index_of_aerosol_forcing_for_interpolation`: Lower latitude index of aerosol forcing for interpolation
    * `integer`: units = index
* `upper_latitude_index_of_aerosol_forcing_for_interpolation`: Upper latitude index of aerosol forcing for interpolation
    * `integer`: units = index
* `lower_longitude_index_of_aerosol_forcing_for_interpolation`: Lower longitude index of aerosol forcing for interpolation
    * `integer`: units = index
* `upper_longitude_index_of_aerosol_forcing_for_interpolation`: Upper longitude index of aerosol forcing for interpolation
    * `integer`: units = index
* `lower_latitude_index_of_cloud_nuclei_forcing_for_interpolation`: Lower latitude index of cloud nuclei forcing for interpolation
    * `integer`: units = index
* `upper_latitude_index_of_cloud_nuclei_forcing_for_interpolation`: Upper latitude index of cloud nuclei forcing for interpolation
    * `integer`: units = index
* `lower_longitude_index_of_cloud_nuclei_forcing_for_interpolation`: Lower longitude index of cloud nuclei forcing for interpolation
    * `integer`: units = index
* `upper_longitude_index_of_cloud_nuclei_forcing_for_interpolation`: Upper longitude index of cloud nuclei forcing for interpolation
    * `integer`: units = index
* `map_of_block_column_number_to_global_i_index`: Map of block column number to global i index
    * `integer`: units = none
* `map_of_block_column_number_to_global_j_index`: Map of block column number to global j index
    * `integer`: units = none
* `flag_for_ugwp_version_1`: Flag for ugwp version 1
    * `logical`: units = flag
* `lower_latitude_index_of_absolute_momentum_flux_due_to_nonorographic_gravity_wave_drag_for_interpolation`: Lower latitude index of absolute momentum flux due to nonorographic gravity wave drag for interpolation
    * `integer`: units = none
* `upper_latitude_index_of_absolute_momentum_flux_due_to_nonorographic_gravity_wave_drag_for_interpolation`: Upper latitude index of absolute momentum flux due to nonorographic gravity wave drag for interpolation
    * `integer`: units = none
* `control_for_soil_type_dataset`: Control for soil type dataset
    * `integer`: units = index
* `control_for_vegetation_dataset`: Control for vegetation dataset
    * `integer`: units = index
* `identifier_for_noahmp_land_surface_scheme`: Identifier for noahmp land surface scheme
    * `integer`: units = flag
* `identifier_for_ruc_land_surface_scheme`: Identifier for ruc land surface scheme
    * `integer`: units = flag
* `lower_bound_of_vertical_dimension_of_surface_snow`: Lower bound of vertical dimension of surface snow
    * `integer`: units = count
* `upper_bound_of_snow_vertical_dimension_for_land_surface_model`: Upper bound of snow vertical dimension for land surface model
    * `integer`: units = count
* `number_of_x_points_for_current_cubed_sphere_tile`: Number of x points for current cubed sphere tile
    * `integer`: units = count
* `number_of_y_points_for_current_cubed_sphere_tile`: Number of y points for current cubed sphere tile
    * `integer`: units = count
* `starting_x_index_for_current_mpi_rank`: Starting x index for current mpi rank
    * `integer`: units = count
* `starting_y_index_for_current_mpi_rank`: Starting y index for current mpi rank
    * `integer`: units = count
* `number_of_random_numbers`: Number of random numbers
    * `integer`: units = count
* `number_of_timesteps_between_shortwave_radiation_calls`: Number of timesteps between shortwave radiation calls
    * `integer`: units = 
* `control_for_deep_convection_scheme`: Control for deep convection scheme
    * `integer`: units = flag
* `flag_for_dominant_precipitation_type_partition`: Flag for dominant precipitation type partition
    * `logical`: units = flag
* `flag_for_random_clouds_in_relaxed_arakawa_schubert_deep_convection`: Flag for random clouds in relaxed arakawa schubert deep convection
    * `logical`: units = flag
* `number_of_timesteps_between_surface_cycling_calls`: Number of timesteps between surface cycling calls
    * `integer`: units = count
* `random_number_seed_for_deep_convection`: Random number seed for deep convection
    * `integer`: units = none
* `control_for_nsstm`: Control for nsstm
    * `integer`: units = flag
* `index_of_cubed_sphere_tile`: Index of cubed sphere tile
    * `integer`: units = none
* `vertical_dimension_of_soil_internal_to_land_surface_scheme`: Vertical dimension of soil internal to land surface scheme
    * `integer`: units = count
* `control_for_surface_albedo`: Control for surface albedo
    * `integer`: units = flag
* `flag_for_gcycle_surface_option`: Flag for gcycle surface option
    * `logical`: units = flag
* `flag_for_nsstm_analysis_in_gcycle`: Flag for nsstm analysis in gcycle
    * `logical`: units = flag
* `flag_for_output_of_tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep_assuming_clear_sky`: Flag for output of tendency of air temperature due to shortwave heating on radiation timestep assuming clear sky
    * `logical`: units = flag
* `model_layer_number_at_cloud_base`: Model layer number at cloud base
    * `integer`: units = index
* `model_layer_number_at_cloud_top`: Model layer number at cloud top
    * `integer`: units = index
* `surface_sw_fluxes_assuming_total_and_clear_sky_on_radiation_timestep`: Surface sw fluxes assuming total and clear sky on radiation timestep
    * `sfcfsw_type`: units = W m-2
* `toa_sw_fluxes_assuming_total_and_clear_sky_on_radiation_timestep`: Toa sw fluxes assuming total and clear sky on radiation timestep
    * `topfsw_type`: units = W m-2
* `RRTMGP_sw_fluxes`: RRTMGP sw fluxes
    * `profsw_type`: units = W m-2
* `components_of_surface_downward_shortwave_fluxes`: Components of surface downward shortwave fluxes
    * `cmpfsw_type`: units = W m-2
* `identifier_for_fer_hires_microphysics_scheme`: Identifier for fer hires microphysics scheme
    * `integer`: units = flag
* `identifier_for_gfdl_microphysics_scheme`: Identifier for gfdl microphysics scheme
    * `integer`: units = flag
* `identifier_for_thompson_microphysics_scheme`: Identifier for thompson microphysics scheme
    * `integer`: units = flag
* `identifier_for_wsm6_microphysics_scheme`: Identifier for wsm6 microphysics scheme
    * `integer`: units = flag
* `identifier_for_zhao_carr_microphysics_scheme`: Identifier for zhao carr microphysics scheme
    * `integer`: units = flag
* `identifier_for_zhao_carr_pdf_microphysics_scheme`: Identifier for zhao carr pdf microphysics scheme
    * `integer`: units = flag
* `identifier_for_morrison_gettelman_microphysics_scheme`: Identifier for morrison gettelman microphysics scheme
    * `integer `: units = flag
* `vertical_dimension_for_radiation`: Vertical dimension for radiation
    * `integer`: units = count
* `flag_for_initial_time_date_control`: Flag for initial time date control
    * `integer`: units = flag
* `control_for_solar_constant`: Control for solar constant
    * `integer`: units = flag
* `control_for_co2`: Control for co2
    * `integer`: units = flag
* `control_for_shortwave_radiation_aerosols`: Control for shortwave radiation aerosols
    * `integer`: units = flag
* `control_for_surface_emissivity`: Control for surface emissivity
    * `integer`: units = flag
* `number_of_microphysics_varaibles_in_xyz_dimensioned_restart_array`: Number of microphysics varaibles in xyz dimensioned restart array
    * `integer`: units = count
* `flag_for_cloud_overlap_method_for_radiation`: Flag for cloud overlap method for radiation
    * `integer`: units = flag
* `flag_for_sw_clouds_grid_approximation`: Flag for sw clouds grid approximation
    * `integer`: units = flag
* `flag_for_lw_clouds_sub_grid_approximation`: Flag for lw clouds sub grid approximation
    * `integer`: units = flag
* `control_for_shortwave_radiation_liquid_clouds`: Control for shortwave radiation liquid clouds
    * `integer`: units = flag
* `flag_for_crick_elimination`: Flag for crick elimination
    * `logical`: units = flag
* `flag_for_in_cloud_condensate`: Flag for in cloud condensate
    * `logical`: units = flag
* `flag_for_turning_off_precipitation_radiative_effect`: Flag for turning off precipitation radiative effect
    * `logical`: units = flag
* `date_and_time_at_model_initialization_in_ISO_order`: Date and time at model initialization in ISO order
    * `integer`: units = none
* `date_and_time_of_forecast_in_United_States_order`: Date and time of forecast in United States order
    * `integer`: units = none
* `flag_for_reset_maximum_hourly_fields`: Flag for reset maximum hourly fields
    * `logical`: units = flag
* `flag_for_radar_reflectivity`: Flag for radar reflectivity
    * `logical`: units = flag
* `kind_dyn`: Kind dyn
    * `integer`: units = none
* `kind_grid`: Kind grid
    * `integer`: units = none
* `kind_phys`: Kind phys
    * `integer`: units = none
* `kind_LOGICAL`: Kind LOGICAL
    * `integer`: units = none
* `kind_INTEGER`: Kind INTEGER
    * `integer`: units = none
* `number_of_aerosol_tracers_for_convection`: Number of aerosol tracers for convection
    * `integer`: units = count
* `number_of_chemical_tracers`: Number of chemical tracers
    * `integer`: units = count
* `index_for_turbulent_kinetic_energy_convective_transport_tracer`: Index for turbulent kinetic energy convective transport tracer
    * `integer`: units = index
* `number_of_tracers_for_samf`: Number of tracers for samf
    * `integer`: units = count
* `vertical_index_at_cloud_base`: Vertical index at cloud base
    * `integer`: units = index
* `vertical_index_at_cloud_top`: Vertical index at cloud top
    * `integer`: units = index
* `flag_for_hurricane_specific_code_in_scale_aware_mass_flux_shallow_convection`: Flag for hurricane specific code in scale aware mass flux shallow convection
    * `logical`: units = flag
* `flag_for_canopy_heat_storage_in_land_surface_scheme`: Flag for canopy heat storage in land surface scheme
    * `logical`: units = flag
* `flag_for_integrated_dynamics_through_earths_atmosphere`: Flag for integrated dynamics through earths atmosphere
    * `logical`: units = flag
* `number_of_plumes`: Number of plumes
    * `integer`: units = count
* `k_level_of_highest_plume`: K level of highest plume
    * `integer`: units = count
* `control_for_gravitational_settling_of_cloud_droplets`: Control for gravitational settling of cloud droplets
    * `integer`: units = flag
* `control_for_tke_budget_output`: Control for tke budget output
    * `integer`: units = flag
* `flag_for_tke_advection`: Flag for tke advection
    * `logical`: units = flag
* `control_for_cloud_pdf_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for cloud pdf in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_mixing_length_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for mixing length in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_edmf_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for edmf in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_edmf_momentum_transport_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for edmf momentum transport in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_edmf_tke_transport_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for edmf tke transport in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_edmf_partitioning_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for edmf partitioning in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_cloud_species_mixing_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for cloud species mixing in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_total_water_mixing_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for total water mixing in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_additional_diagnostics_in_mellor_yamada_nakanishi_niino_pbl_scheme`: Control for additional diagnostics in mellor yamada nakanishi niino pbl scheme
    * `integer`: units = flag
* `control_for_sgs_cloud_radiation_coupling_in_mellor_yamamda_nakanishi_niino_pbl_scheme`: Control for sgs cloud radiation coupling in mellor yamamda nakanishi niino pbl scheme
    * `integer`: units = flag
* `flag_for_mellor_yamada_nakanishi_niino_surface_layer_scheme`: Flag for mellor yamada nakanishi niino surface layer scheme
    * `logical`: units = flag
* `flag_for_aerosol_physics`: Flag for aerosol physics
    * `logical`: units = flag
* `number_of_timesteps_between_longwave_radiation_calls`: Number of timesteps between longwave radiation calls
    * `integer`: units = 
* `flag_for_debug_output`: Flag for debug output
    * `logical`: units = flag
* `number_of_days_in_current_year`: Number of days in current year
    * `integer`: units = days
* `index_of_horizontal_gridpoint_for_debug_output`: Index of horizontal gridpoint for debug output
    * `integer`: units = index 
* `number_of_radiatively_active_gases`: Number of radiatively active gases
    * `integer  `: units = count
* `Gas_concentrations_for_RRTMGP_suite`: Gas concentrations for RRTMGP suite
    * `ty_gas_concs`: units = DDT
* `flag_for_nrl_2006_ozone_scheme`: Flag for nrl 2006 ozone scheme
    * `logical`: units = flag
* `vertical_dimension_of_ozone_forcing_data`: Vertical dimension of ozone forcing data
    * `integer`: units = count
* `number_of_coefficients_in_ozone_forcing_data`: Number of coefficients in ozone forcing data
    * `integer`: units = index
* `adjusted_vertical_layer_dimension_for_radiation`: Adjusted vertical layer dimension for radiation
    * `integer`: units = count
* `adjusted_vertical_level_dimension_for_radiation`: Adjusted vertical level dimension for radiation
    * `integer`: units = count
* `identifier_for_grell_freitas_deep_convection`: Identifier for grell freitas deep convection
    * `integer`: units = flag
* `number_of_condensate_species`: Number of condensate species
    * `integer`: units = count
* `number_of_pdf_based_variables_in_xyz_dimensioned_restart_array`: Number of pdf based variables in xyz dimensioned restart array
    * `integer`: units = count
* `number_of_convective_cloud_variables_in_xyz_dimensioned_restart_array`: Number of convective cloud variables in xyz dimensioned restart array
    * `integer`: units = count
* `index_of_specific_humidity_in_tracer_concentration_array`: Index of specific humidity in tracer concentration array
    * `integer`: units = index
* `index_of_mass_number_concentration_of_cloud_droplets_in_tracer_concentration_array`: Index of mass number concentration of cloud droplets in tracer concentration array
    * `integer`: units = index
* `index_of_mass_number_concentration_of_cloud_ice_in_tracer_concentration_array`: Index of mass number concentration of cloud ice in tracer concentration array
    * `integer`: units = index
* `index_of_mass_number_concentration_of_hygroscopic_aerosols_in_tracer_concentration_array`: Index of mass number concentration of hygroscopic aerosols in tracer concentration array
    * `integer`: units = index
* `index_of_cloud_area_fraction_in_atmosphere_layer_in_tracer_concentration_array`: Index of cloud area fraction in atmosphere layer in tracer concentration array
    * `integer`: units = index
* `index_of_cloud_liquid_water_effective_radius_in_xyz_dimensioned_restart_array`: Index of cloud liquid water effective radius in xyz dimensioned restart array
    * `integer`: units = index
* `index_of_cloud_ice_effective_radius_in_xyz_dimensioned_restart_array`: Index of cloud ice effective radius in xyz dimensioned restart array
    * `integer`: units = index
* `index_of_snow_effective_radius_in_xyz_dimensioned_restart_array`: Index of snow effective radius in xyz dimensioned restart array
    * `integer`: units = index
* `flag_for_gfdl_microphysics_radiation_interaction`: Flag for gfdl microphysics radiation interaction
    * `logical`: units = flag
* `flag_for_shoc_cloud_area_fraction_for_radiation`: Flag for shoc cloud area fraction for radiation
    * `logical`: units = flag
* `flag_for_cloud_effective_radii`: Flag for cloud effective radii
    * `logical`: units = flag
* `flag_for_mellor_yamada_nakanishi_niino_pbl_scheme`: Flag for mellor yamada nakanishi niino pbl scheme
    * `logical`: units = flag
* `flag_for_cloud_area_fraction_option_for_radiation`: Flag for cloud area fraction option for radiation
    * `logical`: units = flag
* `flag_for_scale_aware_mass_flux_convection`: Flag for scale aware mass flux convection
    * `logical`: units = flag
* `flag_for_stochastic_cloud_fraction_perturbations`: Flag for stochastic cloud fraction perturbations
    * `logical`: units = flag
* `cloud_effect_to_optical_depth_and_cloud_fraction`: Cloud effect to optical depth and cloud fraction
    * `integer`: units = flag
* `vertical_index_difference_between_inout_and_local`: Vertical index difference between inout and local
    * `integer`: units = index
* `vertical_index_difference_between_layer_and_upper_bound`: Vertical index difference between layer and upper bound
    * `integer`: units = index
* `vertical_index_difference_between_layer_and_lower_bound`: Vertical index difference between layer and lower bound
    * `integer`: units = index
* `random_number_seed_for_mcica_longwave`: Random number seed for mcica longwave
    * `integer`: units = none
* `toa_lw_fluxes_assuming_total_and_clear_sky_on_radiation_timestep`: Toa lw fluxes assuming total and clear sky on radiation timestep
    * `topflw_type`: units = W m-2
* `surface_lw_fluxes_assuming_total_and_clear_sky_on_radiation_timestep`: Surface lw fluxes assuming total and clear sky on radiation timestep
    * `sfcflw_type`: units = W m-2
* `counter_for_GF`: Counter for GF
    * `integer`: units = none
* `control_for_shallow_convection_scheme`: Control for shallow convection scheme
    * `integer`: units = flag
* `flag_for_generic_tendency_due_to_shallow_convection`: Flag for generic tendency due to shallow convection
    * `logical`: units = flag
* `flag_for_generic_tendency_due_to_deep_convection`: Flag for generic tendency due to deep convection
    * `logical`: units = flag
* `GFS_control_type_instance`: GFS control type instance
    * `GFS_control_type`: units = DDT
* `GFS_data_type_instance_all_blocks`: GFS data type instance all blocks
    * `GFS_data_type`: units = DDT
* `GFS_interstitial_type_instance_all_threads`: GFS interstitial type instance all threads
    * `GFS_interstitial_type`: units = DDT
* `GFS_statein_type_instance`: GFS statein type instance
    * `GFS_statein_type`: units = DDT
* `GFS_stateout_type_instance`: GFS stateout type instance
    * `GFS_stateout_type`: units = DDT
* `GFS_sfcprop_type_instance`: GFS sfcprop type instance
    * `GFS_sfcprop_type`: units = DDT
* `GFS_coupling_type_instance`: GFS coupling type instance
    * `GFS_coupling_type`: units = DDT
* `GFS_grid_type_instance`: GFS grid type instance
    * `GFS_grid_type`: units = DDT
* `GFS_tbd_type_instance`: GFS tbd type instance
    * `GFS_tbd_type`: units = DDT
* `GFS_cldprop_type_instance`: GFS cldprop type instance
    * `GFS_cldprop_type`: units = DDT
* `GFS_radtend_type_instance`: GFS radtend type instance
    * `GFS_radtend_type`: units = DDT
* `GFS_diag_type_instance`: GFS diag type instance
    * `GFS_diag_type`: units = DDT
* `GFS_interstitial_type_instance`: GFS interstitial type instance
    * `GFS_interstitial_type`: units = DDT
* `ccpp_block_number`: Ccpp block number
    * `integer`: units = index
* `ccpp_loop_counter`: Ccpp loop counter
    * `integer`: units = index
* `flag_for_iteration`: Flag for iteration
    * `logical`: units = flag
* `flag_for_guess_run`: Flag for guess run
    * `logical`: units = flag
* `flag_for_nrl_2015_ozone_scheme`: Flag for nrl 2015 ozone scheme
    * `logical`: units = flag
* `extra_top_layer`: Extra top layer
    * `integer`: units = none
* `flag_for_output_of_tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep_assuming_clear_sky`: Flag for output of tendency of air temperature due to longwave heating on radiation timestep assuming clear sky
    * `logical`: units = flag
* `number_of_frozen_precipitation_species`: Number of frozen precipitation species
    * `integer`: units = count
* `flag_for_uniform_subcolumns`: Flag for uniform subcolumns
    * `logical`: units = flag
* `flag_for_cloud_ice_processes`: Flag for cloud ice processes
    * `logical`: units = flag
* `flag_for_heterogeneous_nucleation`: Flag for heterogeneous nucleation
    * `logical`: units = flag
* `flag_for_allowance_of_supersaturation_after_sedimentation`: Flag for allowance of supersaturation after sedimentation
    * `logical`: units = flag
* `flag_for_seifert_and_beheng_2001_autoconversion`: Flag for seifert and beheng 2001 autoconversion
    * `logical`: units = flag
* `flag_for_hail_instead_of_graupel`: Flag for hail instead of graupel
    * `logical`: units = flag
* `flag_for_graupel_instead_of_hail`: Flag for graupel instead of hail
    * `logical`: units = flag
* `flag_for_prescribed_cloud_droplet_number_concentration`: Flag for prescribed cloud droplet number concentration
    * `logical`: units = flag
* `flag_for_prescribed_cloud_ice_number_concentration`: Flag for prescribed cloud ice number concentration
    * `logical`: units = flag
* `flag_for_prescribed_graupel_number_concentration`: Flag for prescribed graupel number concentration
    * `logical`: units = flag
* `flag_for_gmao_autoconversion_to_snow`: Flag for gmao autoconversion to snow
    * `logical`: units = flag
* `flag_for_liu_autoconversion_to_rain`: Flag for liu autoconversion to rain
    * `logical`: units = flag
* `flag_flip`: Flag flip
    * `logical`: units = flag
* `flag_for_skip_cloud_macrophysics_in_MG`: Flag for skip cloud macrophysics in MG
    * `logical`: units = flag
* `control_for_pdf_shape_for_microphysics`: Control for pdf shape for microphysics
    * `integer`: units = flag
* `number_of_tracers_for_cloud_condensate`: Number of tracers for cloud condensate
    * `integer`: units = count
* `control_for_drag_suite_gravity_wave_drag`: Control for drag suite gravity wave drag
    * `integer`: units = flag
* `flag_for_gsl_drag_suite_large_scale_orographic_and_blocking_drag`: Flag for gsl drag suite large scale orographic and blocking drag
    * `logical`: units = flag
* `flag_for_gsl_drag_suite_small_scale_orographic_drag`: Flag for gsl drag suite small scale orographic drag
    * `logical`: units = flag
* `flag_for_gsl_drag_suite_turbulent_orographic_form_drag`: Flag for gsl drag suite turbulent orographic form drag
    * `logical`: units = flag
* `flag_for_optical_property_for_ice_clouds_for_shortwave_radiation`: Flag for optical property for ice clouds for shortwave radiation
    * `integer`: units = flag
* `number_of_shortwave_bands`: Number of shortwave bands
    * `integer`: units = count
* `shortwave_optical_properties_for_cloudy_atmosphere_by_band`: Shortwave optical properties for cloudy atmosphere by band
    * `ty_optical_props_2str`: units = DDT
* `shortwave_optical_properties_for_precipitation_by_band`: Shortwave optical properties for precipitation by band
    * `ty_optical_props_2str`: units = DDT
* `identifier_for_simplified_arakawa_schubert_deep_convection`: Identifier for simplified arakawa schubert deep convection
    * `integer`: units = flag
* `number_of_spectral_wave_truncation_for_simplified_arakawa_schubert_convection`: Number of spectral wave truncation for simplified arakawa schubert convection
    * `integer`: units = count
* `identifier_for_2019_scale_aware_tke_moist_edmf_pbl`: Identifier for 2019 scale aware tke moist edmf pbl
    * `integer`: units = none
* `flag_for_convective_gravity_wave_drag`: Flag for convective gravity wave drag
    * `logical`: units = flag
* `random_number_seed_for_mcica_shortwave`: Random number seed for mcica shortwave
    * `integer`: units = none
* `flag_for_chemistry_coupling`: Flag for chemistry coupling
    * `logical`: units = flag
* `flag_for_relaxed_arakawa_schubert_deep_convection`: Flag for relaxed arakawa schubert deep convection
    * `logical`: units = flag
* `flag_for_chikira_sugiyama_deep_convection_scheme`: Flag for chikira sugiyama deep convection scheme
    * `logical`: units = flag
* `flag_for_stochastic_physics_perturbations`: Flag for stochastic physics perturbations
    * `logical`: units = flag
* `flag_for_stochastic_microphysics_perturbations`: Flag for stochastic microphysics perturbations
    * `logical`: units = flag
* `flag_for_mountain_blocking_for_sppt`: Flag for mountain blocking for sppt
    * `logical`: units = flag
* `flag_for_stochastic_shum_option`: Flag for stochastic shum option
    * `logical`: units = flag
* `flag_for_stochastic_skeb_option`: Flag for stochastic skeb option
    * `logical`: units = flag
* `flag_for_cellular_automata`: Flag for cellular automata
    * `logical`: units = flag
* `flag_for_global_cellular_automata`: Flag for global cellular automata
    * `logical`: units = flag
* `forecast_month`: Forecast month
    * `integer`: units = none
* `flag_for_shoc`: Flag for shoc
    * `logical`: units = flag
* `flag_for_mg3_as_mg2`: Flag for mg3 as mg2
    * `logical`: units = flag
* `index_of_mass_number_concentration_of_rain_in_tracer_concentration_array`: Index of mass number concentration of rain in tracer concentration array
    * `integer`: units = index
* `index_of_mass_number_concentration_of_snow_in_tracer_concentration_array`: Index of mass number concentration of snow in tracer concentration array
    * `integer`: units = index
* `index_of_mass_number_concentration_of_graupel_in_tracer_concentration_array`: Index of mass number concentration of graupel in tracer concentration array
    * `integer`: units = index
* `index_of_mass_number_concentration_of_nonhygroscopic_ice_nucleating_aerosols_in_tracer_concentration_array`: Index of mass number concentration of nonhygroscopic ice nucleating aerosols in tracer concentration array
    * `integer`: units = index
* `index_of_mass_weighted_rime_factor_in_tracer_concentration_array`: Index of mass weighted rime factor in tracer concentration array
    * `integer`: units = index
* `flag_for_aerosol_convective_transport_and_PBL_diffusion`: Flag for aerosol convective transport and PBL diffusion
    * `logical`: units = flag
* `index_of_first_chemical_tracer_in_tracer_concentration_array`: Index of first chemical tracer in tracer concentration array
    * `integer`: units = index
* `flag_for_hybrid_edmf_pbl_scheme`: Flag for hybrid edmf pbl scheme
    * `logical`: units = flag
* `flag_for_scale_aware_TKE_moist_EDMF_PBL`: Flag for scale aware TKE moist EDMF PBL
    * `logical`: units = flag
* `flag_for_scale_aware_shin_hong_pbl_scheme`: Flag for scale aware shin hong pbl scheme
    * `logical`: units = flag
* `flag_for_ysu_pbl_scheme`: Flag for ysu pbl scheme
    * `logical`: units = flag
* `flag_for_hydrostatic_solver`: Flag for hydrostatic solver
    * `logical`: units = flag
* `flag_for_hydrostatic_heating_from_physics`: Flag for hydrostatic heating from physics
    * `logical`: units = flag
* `flag_for_hurricane_specific_code_in_scale_aware_mass_flux_deep_convection`: Flag for hurricane specific code in scale aware mass flux deep convection
    * `logical`: units = flag
* `flag_for_global_cellular_automata_closure`: Flag for global cellular automata closure
    * `logical`: units = flag
* `flag_for_global_cellular_automata_deep_convective_entrainment`: Flag for global cellular automata deep convective entrainment
    * `logical`: units = flag
* `flag_for_global_cellular_automata_trigger`: Flag for global cellular automata trigger
    * `logical`: units = flag
* `vertical_dimension_of_water_vapor_forcing_data`: Vertical dimension of water vapor forcing data
    * `integer`: units = count
* `number_of_coefficients_in_water_vapor_forcing_data`: Number of coefficients in water vapor forcing data
    * `integer`: units = index
* `flag_for_ugwp_version_0`: Flag for ugwp version 0
    * `logical`: units = flag
* `flag_for_ugwp_version_0_orographic_gwd`: Flag for ugwp version 0 orographic gwd
    * `logical`: units = flag
* `flag_for_ugwp_version_1_orographic_gwd`: Flag for ugwp version 1 orographic gwd
    * `logical`: units = flag
* `flag_for_ugwp_version_1_nonorographic_gwd`: Flag for ugwp version 1 nonorographic gwd
    * `logical`: units = flag
* `flag_for_rrtmgp_longwave_jacobian`: Flag for rrtmgp longwave jacobian
    * `logical `: units = flag
* `number_of_gaussian_quadrature_angles_for_radiation`: Number of gaussian quadrature angles for radiation
    * `integer`: units = count
* `longwave_optical_properties_for_clear_sky`: Longwave optical properties for clear sky
    * `ty_optical_props_1scl`: units = DDT
* `longwave_optical_properties_for_cloudy_atmosphere`: Longwave optical properties for cloudy atmosphere
    * `ty_optical_props_2str`: units = DDT
* `longwave_optical_properties_for_aerosols`: Longwave optical properties for aerosols
    * `ty_optical_props_1scl`: units = DDT
* `longwave_source_function`: Longwave source function
    * `ty_source_func_lw`: units = DDT
* `flag_for_simplified_arakawa_schubert_shallow_convection`: Flag for simplified arakawa schubert shallow convection
    * `logical`: units = flag
* `flag_for_old_PBL_scheme`: Flag for old PBL scheme
    * `logical`: units = flag
* `flag_for_moorthi_stratus`: Flag for moorthi stratus
    * `logical`: units = flag
* `flag_for_convective_transport_of_tracers`: Flag for convective transport of tracers
    * `logical`: units = flag
* `number_of_total_tracers`: Number of total tracers
    * `integer`: units = count
* `flag_for_converting_hydrometeors_from_moist_to_dry_air`: Flag for converting hydrometeors from moist to dry air
    * `logical`: units = flag
* `number_of_aerosol_tracers_MG`: Number of aerosol tracers MG
    * `integer`: units = count
* `shortwave_optical_properties_for_aerosols`: Shortwave optical properties for aerosols
    * `ty_optical_props_2str`: units = DDT
* `flag_for_limited_surface_roughness_length_over_ocean`: Flag for limited surface roughness length over ocean
    * `logical`: units = flag
* `flag_for_surface_roughness_option_over_water`: Flag for surface roughness option over water
    * `integer`: units = flag
* `flag_for_decorrelation_length_method`: Flag for decorrelation length method
    * `integer`: units = flag
* `flag_for_constant_decorrelation_length_method`: Flag for constant decorrelation length method
    * `integer`: units = flag
* `flag_for_hogan_decorrelation_length_method`: Flag for hogan decorrelation length method
    * `integer`: units = flag
* `flag_for_oreopoulos_decorrelation_length_method`: Flag for oreopoulos decorrelation length method
    * `integer`: units = flag
* `flag_shallow_convective_cloud`: Flag shallow convective cloud
    * `logical`: units = 
* `identifier_for_simplified_arakawa_schubert_shallow_convection`: Identifier for simplified arakawa schubert shallow convection
    * `integer`: units = flag
* `identifier_for_scale_aware_mass_flux_shallow_convection`: Identifier for scale aware mass flux shallow convection
    * `integer`: units = flag
* `lower_bound_for_depth_of_sea_temperature_for_nsstm`: Lower bound for depth of sea temperature for nsstm
    * `integer`: units = mm
* `upper_bound_for_depth_of_sea_temperature_for_nsstm`: Upper bound for depth of sea temperature for nsstm
    * `integer`: units = mm
* `RRTMGP_lw_fluxes`: RRTMGP lw fluxes
    * `proflw_type`: units = W m-2
* `topflw_type`: Topflw type
    * `topflw_type`: units = DDT
* `sfcflw_type`: Sfcflw type
    * `sfcflw_type`: units = DDT
* `proflw_type`: Proflw type
    * `proflw_type`: units = DDT
* `flag_for_hurricane_specific_code_in_hybrid_edmf_pbl_scheme`: Flag for hurricane specific code in hybrid edmf pbl scheme
    * `logical`: units = flag
* `initial_permutation_seed_lw`: Initial permutation seed lw
    * `integer`: units = none
* `flag_for_maximum_random_cloud_overlap_method`: Flag for maximum random cloud overlap method
    * `integer`: units = flag
* `flag_for_random_cloud_overlap_method`: Flag for random cloud overlap method
    * `integer`: units = flag
* `flag_for_maximum_cloud_overlap_method`: Flag for maximum cloud overlap method
    * `integer`: units = flag
* `longwave_optical_properties_for_precipitation`: Longwave optical properties for precipitation
    * `ty_optical_props_2str`: units = DDT
* `topfsw_type`: Topfsw type
    * `topfsw_type`: units = DDT
* `sfcfsw_type`: Sfcfsw type
    * `sfcfsw_type`: units = DDT
* `cmpfsw_type`: Cmpfsw type
    * `cmpfsw_type`: units = DDT
* `profsw_type`: Profsw type
    * `profsw_type`: units = DDT
* `flag_for_resetting_radar_reflectivity_calculation`: Flag for resetting radar reflectivity calculation
    * `logical`: units = flag
* `flag_for_stochastic_radiative_heating_perturbations`: Flag for stochastic radiative heating perturbations
    * `logical`: units = flag
* `flag_for_ocean_wave_coupling`: Flag for ocean wave coupling
    * `logical`: units = flag
* `flag_for_surface_layer_scheme_ocean_currents`: Flag for surface layer scheme ocean currents
    * `logical`: units = flag
* `flag_for_surface_layer_scheme_surface_drag_coefficient_for_momentum_in_air_perturbations`: Flag for surface layer scheme surface drag coefficient for momentum in air perturbations
    * `logical`: units = flag
* `control_for_surface_layer_scheme_skin_temperature_update`: Control for surface layer scheme skin temperature update
    * `integer`: units = flag
* `identifier_for_noah_land_surface_scheme`: Identifier for noah land surface scheme
    * `integer`: units = flag
* `initial_permutation_seed_sw`: Initial permutation seed sw
    * `integer`: units = none
* `shortwave_optical_properties_for_cloudy_atmosphere`: Shortwave optical properties for cloudy atmosphere
    * `ty_optical_props_2str`: units = DDT
* `shortwave_optical_properties_for_precipitation`: Shortwave optical properties for precipitation
    * `ty_optical_props_2str    `: units = DDT
* `control_for_land_surface_scheme_dynamic_vegetation`: Control for land surface scheme dynamic vegetation
    * `integer`: units = index
* `control_for_land_surface_scheme_canopy_stomatal_resistance`: Control for land surface scheme canopy stomatal resistance
    * `integer`: units = index
* `control_for_land_surface_scheme_soil_moisture_factor_stomatal_resistance`: Control for land surface scheme soil moisture factor stomatal resistance
    * `integer`: units = index
* `control_for_land_surface_scheme_runoff_and_groundwater`: Control for land surface scheme runoff and groundwater
    * `integer`: units = index
* `control_for_land_surface_scheme_surface_layer_drag_coefficient`: Control for land surface scheme surface layer drag coefficient
    * `integer`: units = index
* `control_for_land_surface_scheme_supercooled_liquid_water`: Control for land surface scheme supercooled liquid water
    * `integer`: units = index
* `control_for_land_surface_scheme_frozen_soil_permeability`: Control for land surface scheme frozen soil permeability
    * `integer`: units = index
* `control_for_land_surface_scheme_radiative_transfer`: Control for land surface scheme radiative transfer
    * `integer`: units = index
* `control_for_land_surface_scheme_surface_snow_albedo`: Control for land surface scheme surface snow albedo
    * `integer`: units = index
* `control_for_land_surface_scheme_precipitation_type_partition`: Control for land surface scheme precipitation type partition
    * `integer`: units = index
* `control_for_land_surface_scheme_lower_boundary_soil_temperature`: Control for land surface scheme lower boundary soil temperature
    * `integer`: units = index
* `control_for_land_surface_scheme_soil_and_snow_temperature_time_integration`: Control for land surface scheme soil and snow temperature time integration
    * `integer`: units = index
* `shortwave_optical_properties_for_clear_sky`: Shortwave optical properties for clear sky
    * `ty_optical_props_2str`: units = DDT
* `flag_for_ugwp_version_0_nonorographic_gwd`: Flag for ugwp version 0 nonorographic gwd
    * `logical`: units = flag
* `index_of_water_vegetation_category`: Index of water vegetation category
    * `integer`: units = index
* `number_of_species_for_aerosol_optical_depth`: Number of species for aerosol optical depth
    * `integer`: units = count
* `number_of_diagnostics_variables_for_radiation`: Number of diagnostics variables for radiation
    * `integer`: units = count
* `number_of_shortwave_spectral_points`: Number of shortwave spectral points
    * `integer`: units = count
* `flag_for_separate_advection_of_condensate_species`: Flag for separate advection of condensate species
    * `logical`: units = flag
* `area_type`: Area type
    * `real`: units = flag
* `number_of_timesteps_for_concurrent_radiation_and_remainder_physics_calls_after_model_initialization`: Number of timesteps for concurrent radiation and remainder physics calls after model initialization
    * `integer`: units = count
* `gravitational_acceleration`: Gravitational acceleration
    * `real(kind=kind_phys)`: units = m s-2
* `gas_constant_dry_air`: Gas constant dry air
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `gas_constant_water_vapor`: Gas constant water vapor
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `latent_heat_of_fusion_of_water_at_0C`: Latent heat of fusion of water at 0C
    * `real(kind=kind_phys)`: units = J kg-1
* `ratio_of_vapor_to_dry_air_gas_constants_minus_one`: Ratio of vapor to dry air gas constants minus one
    * `real(kind=kind_phys)`: units = none
* `ratio_of_dry_air_to_water_vapor_gas_constants`: Ratio of dry air to water vapor gas constants
    * `real(kind=kind_phys)`: units = none
* `ratio_of_dry_air_to_water_vapor_gas_constants_minus_one`: Ratio of dry air to water vapor gas constants minus one
    * `real(kind=kind_phys)`: units = none
* `tendency_of_y_wind_due_to_model_physics`: Tendency of y wind due to model physics
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_x_wind_due_to_model_physics`: Tendency of x wind due to model physics
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_air_temperature_due_to_model_physics`: Tendency of air temperature due to model physics
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_vertically_diffused_tracer_concentration`: Tendency of vertically diffused tracer concentration
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `vertically_diffused_tracer_concentration`: Vertically diffused tracer concentration
    * `real(kind=kind_phys)`: units = kg kg-1
* `tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep`: Tendency of air temperature due to shortwave heating on radiation timestep
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep`: Tendency of air temperature due to longwave heating on radiation timestep
    * `real(kind=kind_phys)`: units = K s-1
* `zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes`: Zenith angle temporal adjustment factor for shortwave fluxes
    * `real(kind=kind_phys)`: units = none
* `surface_dimensionless_exner_function`: Surface dimensionless exner function
    * `real(kind=kind_phys)`: units = none
* `bulk_richardson_number_at_lowest_model_level`: Bulk richardson number at lowest model level
    * `real(kind=kind_phys)`: units = none
* `surface_roughness_length`: Surface roughness length
    * `real(kind=kind_phys)`: units = cm
* `x_wind_at_10m`: X wind at 10m
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind_at_10m`: Y wind at 10m
    * `real(kind=kind_phys)`: units = m s-1
* `monin_obukhov_similarity_function_for_momentum`: Monin obukhov similarity function for momentum
    * `real(kind=kind_phys)`: units = none
* `monin_obukhov_similarity_function_for_heat`: Monin obukhov similarity function for heat
    * `real(kind=kind_phys)`: units = none
* `surface_skin_temperature`: Surface skin temperature
    * `real(kind=kind_phys)`: units = K
* `kinematic_surface_upward_sensible_heat_flux_reduced_by_surface_roughness`: Kinematic surface upward sensible heat flux reduced by surface roughness
    * `real(kind=kind_phys)`: units = K m s-1
* `kinematic_surface_upward_latent_heat_flux_reduced_by_surface_roughness`: Kinematic surface upward latent heat flux reduced by surface roughness
    * `real(kind=kind_phys)`: units = kg kg-1 m s-1
* `surface_wind_stress`: Surface wind stress
    * `real(kind=kind_phys)`: units = m2 s-2
* `wind_speed_at_lowest_model_layer`: Wind speed at lowest model layer
    * `real(kind=kind_phys)`: units = m s-1
* `air_pressure_difference_between_midlayers`: Air pressure difference between midlayers
    * `real(kind=kind_phys)`: units = Pa
* `geopotential_at_interface`: Geopotential at interface
    * `real(kind=kind_phys)`: units = m2 s-2
* `geopotential`: Geopotential
    * `real(kind=kind_phys)`: units = m2 s-2
* `instantaneous_surface_x_momentum_flux`: Instantaneous surface x momentum flux
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_surface_y_momentum_flux`: Instantaneous surface y momentum flux
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_surface_upward_sensible_heat_flux`: Instantaneous surface upward sensible heat flux
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_upward_latent_heat_flux`: Instantaneous surface upward latent heat flux
    * `real(kind=kind_phys)`: units = W m-2
* `atmosphere_boundary_layer_thickness`: Atmosphere boundary layer thickness
    * `real(kind=kind_phys)`: units = m
* `atmosphere_momentum_diffusivity_due_to_background`: Atmosphere momentum diffusivity due to background
    * `real(kind=kind_phys)`: units = m2 s-1
* `atmosphere_heat_diffusivity_due_to_background`: Atmosphere heat diffusivity due to background
    * `real(kind=kind_phys)`: units = m2 s-1
* `sigma_pressure_level_at_upper_extent_of_background_diffusion`: Sigma pressure level at upper extent of background diffusion
    * `real(kind=kind_phys)`: units = none
* `cumulative_change_in_temperature_due_to_PBL`: Cumulative change in temperature due to PBL
    * `real(kind=kind_phys)`: units = K
* `cumulative_change_in_x_wind_due_to_PBL`: Cumulative change in x wind due to PBL
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_y_wind_due_to_PBL`: Cumulative change in y wind due to PBL
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_water_vapor_specific_humidity_due_to_PBL`: Cumulative change in water vapor specific humidity due to PBL
    * `real(kind=kind_phys)`: units = kg kg-1
* `cumulative_change_in_ozone_mixing_ratio_due_to_PBL`: Cumulative change in ozone mixing ratio due to PBL
    * `real(kind=kind_phys)`: units = kg kg-1
* `specific_humidity_of_new_state`: Specific humidity of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `ice_water_mixing_ratio_convective_transport_tracer`: Ice water mixing ratio convective transport tracer
    * `real(kind=kind_phys)`: units = kg kg-1
* `cloud_condensed_water_mixing_ratio_convective_transport_tracer`: Cloud condensed water mixing ratio convective transport tracer
    * `real(kind=kind_phys)`: units = kg kg-1
* `grid_size_related_coefficient_used_in_scale_sensitive_schemes`: Grid size related coefficient used in scale sensitive schemes
    * `real(kind=kind_phys)`: units = none
* `grid_size_related_coefficient_used_in_scale_sensitive_schemes_complement`: Grid size related coefficient used in scale sensitive schemes complement
    * `real(kind=kind_phys)`: units = none
* `tunable_parameter_1_for_maximum_cloud_base_updraft_velocity_in_chikira_sugiyama_deep_convection`: Tunable parameter 1 for maximum cloud base updraft velocity in chikira sugiyama deep convection
    * `real(kind=kind_phys)`: units = m s-1
* `tunable_parameter_2_for_maximum_cloud_base_updraft_velocity_in_chikira_sugiyama_deep_convection`: Tunable parameter 2 for maximum cloud base updraft velocity in chikira sugiyama deep convection
    * `real(kind=kind_phys)`: units = m s-1
* `maximum_updraft_velocity_at_cloud_base`: Maximum updraft velocity at cloud base
    * `real(kind=kind_phys)`: units = m s-1
* `fraction_of_cloud_top_water_scavenged`: Fraction of cloud top water scavenged
    * `real(kind=kind_phys)`: units = km-1
* `fraction_of_tracer_scavenged`: Fraction of tracer scavenged
    * `real(kind=kind_phys)`: units = km-1
* `water_vapor_specific_humidity_save`: Water vapor specific humidity save
    * `real(kind=kind_phys)`: units = kg kg-1
* `cloud_condensed_water_mixing_ratio_save`: Cloud condensed water mixing ratio save
    * `real(kind=kind_phys)`: units = kg kg-1
* `ice_water_mixing_ratio_save`: Ice water mixing ratio save
    * `real(kind=kind_phys)`: units = kg kg-1
* `convective_updraft_area_fraction_at_model_interfaces`: Convective updraft area fraction at model interfaces
    * `real(kind=kind_phys)`: units = frac
* `convective_updraft_area_fraction`: Convective updraft area fraction
    * `real(kind=kind_phys)`: units = frac
* `air_temperature_of_new_state`: Air temperature of new state
    * `real(kind=kind_phys)`: units = K
* `lwe_thickness_of_deep_convective_precipitation_amount`: Lwe thickness of deep convective precipitation amount
    * `real(kind=kind_phys)`: units = m
* `convective_transportable_tracers`: Convective transportable tracers
    * `real(kind=kind_phys)`: units = kg kg-1
* `timestep_for_dynamics`: Timestep for dynamics
    * `real(kind=kind_phys)`: units = s
* `instantaneous_atmosphere_updraft_convective_mass_flux`: Instantaneous atmosphere updraft convective mass flux
    * `real(kind=kind_phys)`: units = kg m-2
* `atmosphere_downdraft_convective_mass_flux_integrated_over_physics_timestep`: Atmosphere downdraft convective mass flux integrated over physics timestep
    * `real(kind=kind_phys)`: units = kg m-2
* `atmosphere_detrainment_convective_mass_flux`: Atmosphere detrainment convective mass flux
    * `real(kind=kind_phys)`: units = kg m-2
* `x_wind_of_new_state`: X wind of new state
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind_of_new_state`: Y wind of new state
    * `real(kind=kind_phys)`: units = m s-1
* `atmosphere_updraft_convective_mass_flux_at_cloud_base_by_cloud_type`: Atmosphere updraft convective mass flux at cloud base by cloud type
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `tunable_parameter_1_for_detrainment_and_precipitation_partitioning_in_chikira_sugiyama_deep_convection`: Tunable parameter 1 for detrainment and precipitation partitioning in chikira sugiyama deep convection
    * `real(kind=kind_phys)`: units = m
* `tunable_parameter_2_for_detrainment_and_precipitation_partitioning_in_chikira_sugiyama_deep_convection`: Tunable parameter 2 for detrainment and precipitation partitioning in chikira sugiyama deep convection
    * `real(kind=kind_phys)`: units = m
* `tunable_parameter_for_entrainment_efficiency_in_chikira_sugiyama_deep_convection`: Tunable parameter for entrainment efficiency in chikira sugiyama deep convection
    * `real(kind=kind_phys)`: units = none
* `mass_fraction_of_convective_cloud_liquid_water`: Mass fraction of convective cloud liquid water
    * `real(kind=kind_phys)`: units = kg kg-1
* `mass_fraction_of_convective_cloud_ice`: Mass fraction of convective cloud ice
    * `real(kind=kind_phys)`: units = kg kg-1
* `vertical_velocity_for_updraft`: Vertical velocity for updraft
    * `real(kind=kind_phys)`: units = m s-1
* `convective_cloud_fraction_for_microphysics`: Convective cloud fraction for microphysics
    * `real(kind=kind_phys)`: units = frac
* `atmosphere_detrained_convective_mass_flux`: Atmosphere detrained convective mass flux
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `tendency_of_cloud_liquid_water_due_to_convective_microphysics`: Tendency of cloud liquid water due to convective microphysics
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `convective_cloud_volume_fraction`: Convective cloud volume fraction
    * `real(kind=kind_phys)`: units = frac
* `ice_mass_fraction_in_convective_tower`: Ice mass fraction in convective tower
    * `real(kind=kind_phys)`: units = frac
* `number_concentration_of_cloud_liquid_water_particles_for_detrainment`: Number concentration of cloud liquid water particles for detrainment
    * `real(kind=kind_phys)`: units = m-3
* `number_concentration_of_ice_crystals_for_detrainment`: Number concentration of ice crystals for detrainment
    * `real(kind=kind_phys)`: units = m-3
* `land_surface_perturbation_magnitudes`: Land surface perturbation magnitudes
    * `real(kind=kind_phys)`: units = variable
* `forecast_utc_hour`: Forecast utc hour
    * `real(kind=kind_phys)`: units = h
* `longitude_in_radians`: Longitude in radians
    * `real(kind=kind_phys)`: units = radian
* `cosine_of_latitude`: Cosine of latitude
    * `real(kind=kind_phys)`: units = none
* `sine_of_latitude`: Sine of latitude
    * `real(kind=kind_phys)`: units = none
* `lwe_surface_snow`: Lwe surface snow
    * `real(kind=kind_phys)`: units = mm
* `surface_snow_area_fraction_over_land`: Surface snow area fraction over land
    * `real(kind=kind_phys)`: units = frac
* `upper_bound_of_max_albedo_assuming_deep_snow`: Upper bound of max albedo assuming deep snow
    * `real(kind=kind_phys)`: units = frac
* `surface_ground_temperature_for_radiation`: Surface ground temperature for radiation
    * `real(kind=kind_phys)`: units = K
* `surface_air_temperature_for_radiation`: Surface air temperature for radiation
    * `real(kind=kind_phys)`: units = K
* `standard_deviation_of_subgrid_orography`: Standard deviation of subgrid orography
    * `real(kind=kind_phys)`: units = m
* `vis_albedo_strong_cosz`: Vis albedo strong cosz
    * `real(kind=kind_phys)`: units = frac
* `nir_albedo_strong_cosz`: Nir albedo strong cosz
    * `real(kind=kind_phys)`: units = frac
* `vis_albedo_weak_cosz`: Vis albedo weak cosz
    * `real(kind=kind_phys)`: units = frac
* `nir_albedo_weak_cosz`: Nir albedo weak cosz
    * `real(kind=kind_phys)`: units = frac
* `strong_cosz_area_fraction`: Strong cosz area fraction
    * `real(kind=kind_phys)`: units = frac
* `weak_cosz_area_fraction`: Weak cosz area fraction
    * `real(kind=kind_phys)`: units = frac
* `sea_ice_area_fraction_in_sea_water`: Sea ice area fraction in sea water
    * `real(kind=kind_phys)`: units = frac
* `sea_ice_temperature`: Sea ice temperature
    * `real(kind=kind_phys)`: units = K
* `surface_albedo_direct_visible`: Surface albedo direct visible
    * `real(kind=kind_phys)`: units = frac
* `surface_albedo_direct_NIR`: Surface albedo direct NIR
    * `real(kind=kind_phys)`: units = frac
* `surface_albedo_diffuse_visible`: Surface albedo diffuse visible
    * `real(kind=kind_phys)`: units = frac
* `surface_albedo_diffuse_NIR`: Surface albedo diffuse NIR
    * `real(kind=kind_phys)`: units = frac
* `area_type`: Area type
    * `real(kind=kind_phys)`: units = flag
* `surface_stochastic_weights_from_coupled_process`: Surface stochastic weights from coupled process
    * `real(kind=kind_phys)`: units = none
* `air_pressure_at_layer_for_RRTMGP_in_hPa`: Air pressure at layer for RRTMGP in hPa
    * `real(kind=kind_phys)`: units = hPa
* `virtual_temperature`: Virtual temperature
    * `real(kind=kind_phys)`: units = K
* `relative_humidity`: Relative humidity
    * `real(kind=kind_phys)`: units = frac
* `air_pressure_at_interface_for_RRTMGP_in_hPa`: Air pressure at interface for RRTMGP in hPa
    * `real(kind=kind_phys)`: units = hPa
* `cosine_of_solar_zenith_angle_for_daytime_points_on_radiation_timestep`: Cosine of solar zenith angle for daytime points on radiation timestep
    * `real(kind=kind_phys)`: units = none
* `cosine_of_solar_zenith_angle_on_radiation_timestep`: Cosine of solar zenith angle on radiation timestep
    * `real(kind=kind_phys)`: units = none
* `surface_albedo_nearIR_direct`: Surface albedo nearIR direct
    * `real(kind=kind_phys)`: units = none
* `surface_albedo_nearIR_diffuse`: Surface albedo nearIR diffuse
    * `real(kind=kind_phys)`: units = none
* ` surface_albedo_uvvis_dir`:  surface albedo uvvis dir
    * `real(kind=kind_phys)`: units = none
* ` surface_albedo_uvvis_dif`:  surface albedo uvvis dif
    * `real(kind=kind_phys)`: units = none
* `surface_albedo_for_diffused_shortwave_on_radiation_timestep`: Surface albedo for diffused shortwave on radiation timestep
    * `real(kind=kind_phys)`: units = frac
* `precipitation_type`: Precipitation type
    * `real(kind=kind_phys)`: units = flag
* `multiplicative_tuning_parameter_for_potential_evaporation`: Multiplicative tuning parameter for potential evaporation
    * `real(kind=kind_phys)`: units = none
* `height_above_ground_at_lowest_model_layer`: Height above ground at lowest model layer
    * `real(kind=kind_phys)`: units = m
* `thickness_of_soil_levels_for_land_surface_model`: Thickness of soil levels for land surface model
    * `real(kind=kind_phys)`: units = m
* `surface_downwelling_longwave_flux_absorbed_by_surface_over_land`: Surface downwelling longwave flux absorbed by surface over land
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_shortwave_flux`: Surface downwelling shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_net_downwelling_shortwave_flux`: Surface net downwelling shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `air_pressure_at_surface_adjacent_layer`: Air pressure at surface adjacent layer
    * `real(kind=kind_phys)`: units = Pa
* `total_precipitation_rate_on_dynamics_timestep_over_land`: Total precipitation rate on dynamics timestep over land
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `air_temperature_at_surface_adjacent_layer`: Air temperature at surface adjacent layer
    * `real(kind=kind_phys)`: units = K
* `bounded_specific_humidity_at_lowest_model_layer_over_land`: Bounded specific humidity at lowest model layer over land
    * `real(kind=kind_phys)`: units = kg kg-1
* `potential_temperature_at_lowest_model_layer`: Potential temperature at lowest model layer
    * `real(kind=kind_phys)`: units = K
* `saturation_specific_humidity_at_lowest_model_layer`: Saturation specific humidity at lowest model layer
    * `real(kind=kind_phys)`: units = kg kg-1
* `saturation_specific_humidity_slope`: Saturation specific humidity slope
    * `real(kind=kind_phys)`: units = K-1
* `bounded_vegetation_area_fraction`: Bounded vegetation area fraction
    * `real(kind=kind_phys)`: units = frac
* `min_vegetation_area_fraction`: Min vegetation area fraction
    * `real(kind=kind_phys)`: units = frac
* `max_vegetation_area_fraction`: Max vegetation area fraction
    * `real(kind=kind_phys)`: units = frac
* `deep_soil_temperature`: Deep soil temperature
    * `real(kind=kind_phys)`: units = K
* `baseline_surface_roughness_length`: Baseline surface roughness length
    * `real(kind=kind_phys)`: units = m
* `surface_roughness_length_over_land`: Surface roughness length over land
    * `real(kind=kind_phys)`: units = m
* `surface_longwave_emissivity_over_land_interstitial`: Surface longwave emissivity over land interstitial
    * `real(kind=kind_phys)`: units = frac
* `baseline_surface_longwave_emissivity`: Baseline surface longwave emissivity
    * `real(kind=kind_phys)`: units = frac
* `canopy_water_amount_in_m`: Canopy water amount in m
    * `real(kind=kind_phys)`: units = m
* `surface_skin_temperature_after_iteration_over_land`: Surface skin temperature after iteration over land
    * `real(kind=kind_phys)`: units = K
* `soil_temperature`: Soil temperature
    * `real(kind=kind_phys)`: units = K
* `volume_fraction_of_condensed_water_in_soil`: Volume fraction of condensed water in soil
    * `real(kind=kind_phys)`: units = frac
* `volume_fraction_of_unfrozen_water_in_soil`: Volume fraction of unfrozen water in soil
    * `real(kind=kind_phys)`: units = frac
* `actual_snow_depth`: Actual snow depth
    * `real(kind=kind_phys)`: units = m
* `water_equivalent_accumulated_snow_depth_over_land`: Water equivalent accumulated snow depth over land
    * `real(kind=kind_phys)`: units = m
* `surface_conductance_for_heat_and_moisture_in_air_over_land`: Surface conductance for heat and moisture in air over land
    * `real(kind=kind_phys)`: units = m s-1
* `stefan_boltzmann_constant`: Stefan boltzmann constant
    * `real(kind=kind_phys)`: units = W m-2 K-4
* `specific_heat_of_liquid_water_at_constant_pressure`: Specific heat of liquid water at constant pressure
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `specific_heat_of_ice_at_constant_pressure`: Specific heat of ice at constant pressure
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `upward_latent_heat_flux_in_canopy`: Upward latent heat flux in canopy
    * `real(kind=kind_phys)`: units = W m-2
* `upward_latent_heat_flux_in_soil`: Upward latent heat flux in soil
    * `real(kind=kind_phys)`: units = W m-2
* `transpiration_flux`: Transpiration flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upward_latent_heat_flux_due_to_snow_deposition_sublimation`: Surface upward latent heat flux due to snow deposition sublimation
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upward_potential_latent_heat_flux_over_land`: Surface upward potential latent heat flux over land
    * `real(kind=kind_phys)`: units = W m-2
* `upward_heat_flux_in_soil_over_land`: Upward heat flux in soil over land
    * `real(kind=kind_phys)`: units = W m-2
* `latent_heat_flux_from_precipitating_snow`: Latent heat flux from precipitating snow
    * `real(kind=kind_phys)`: units = W m-2
* `latent_heat_flux_from_freezing_rain`: Latent heat flux from freezing rain
    * `real(kind=kind_phys)`: units = W m-2
* `latent_heat_flux_due_to_snowmelt`: Latent heat flux due to snowmelt
    * `real(kind=kind_phys)`: units = W m-2
* `surface_runoff_flux_in_m_sm1`: Surface runoff flux in m sm1
    * `real(kind=kind_phys)`: units = m s-1
* `subsurface_runoff_flux_in_m_sm1`: Subsurface runoff flux in m sm1
    * `real(kind=kind_phys)`: units = m s-1
* `soil_moisture_content_in_m`: Soil moisture content in m
    * `real(kind=kind_phys)`: units = m
* `surface_specific_humidity_over_land`: Surface specific humidity over land
    * `real(kind=kind_phys)`: units = kg kg-1
* `bulk_richardson_number_at_lowest_model_level_over_land`: Bulk richardson number at lowest model level over land
    * `real(kind=kind_phys)`: units = none
* `volume_fraction_of_condensed_water_in_soil_at_wilting_point`: Volume fraction of condensed water in soil at wilting point
    * `real(kind=kind_phys)`: units = frac
* `threshold_volume_fraction_of_condensed_water_in_soil`: Threshold volume fraction of condensed water in soil
    * `real(kind=kind_phys)`: units = frac
* `soil_porosity`: Soil porosity
    * `real(kind=kind_phys)`: units = frac
* `time_since_last_snowfall`: Time since last snowfall
    * `real(kind=kind_phys)`: units = s
* `sigma_pressure_hybrid_coordinate_a_coefficient`: Sigma pressure hybrid coordinate a coefficient
    * `real(kind=kind_phys)`: units = Pa
* `sigma_pressure_hybrid_coordinate_b_coefficient`: Sigma pressure hybrid coordinate b coefficient
    * `real(kind=kind_phys)`: units = none
* `multiplicative_tunable_parameters_for_mountain_blocking_and_orographic_gravity_wave_drag`: Multiplicative tunable parameters for mountain blocking and orographic gravity wave drag
    * `real(kind=kind_phys)`: units = none
* `tunable_parameters_for_convective_gravity_wave_drag`: Tunable parameters for convective gravity wave drag
    * `real(kind=kind_phys)`: units = none
* `air_pressure_at_bottom_extent_of_rayleigh_damping`: Air pressure at bottom extent of rayleigh damping
    * `real(kind=kind_phys)`: units = Pa
* `timescale_for_rayleigh_damping`: Timescale for rayleigh damping
    * `real(kind=kind_phys)`: units = d
* `standard_atmospheric_pressure`: Standard atmospheric pressure
    * `real(kind=kind_phys)`: units = Pa
* `height_above_mean_sea_level`: Height above mean sea level
    * `real(kind=kind_phys)`: units = m
* `unfiltered_height_above_mean_sea_level`: Unfiltered height above mean sea level
    * `real(kind=kind_phys)`: units = m
* `convexity_of_subgrid_orography`: Convexity of subgrid orography
    * `real(kind=kind_phys)`: units = none
* `angle_from_east_of_maximum_subgrid_orographic_variations`: Angle from east of maximum subgrid orographic variations
    * `real(kind=kind_phys)`: units = degree
* `slope_of_subgrid_orography`: Slope of subgrid orography
    * `real(kind=kind_phys)`: units = none
* `anisotropy_of_subgrid_orography`: Anisotropy of subgrid orography
    * `real(kind=kind_phys)`: units = none
* `maximum_subgrid_orography`: Maximum subgrid orography
    * `real(kind=kind_phys)`: units = m
* `fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height`: Fraction of grid box with subgrid orography higher than critical height
    * `real(kind=kind_phys)`: units = frac
* `asymmetry_of_subgrid_orography`: Asymmetry of subgrid orography
    * `real(kind=kind_phys)`: units = none
* `latitude_in_radians`: Latitude in radians
    * `real(kind=kind_phys)`: units = radian
* `tracer_concentration`: Tracer concentration
    * `real(kind=kind_phys)`: units = kg kg-1
* `x_stress_due_to_gravity_wave_drag`: X stress due to gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `y_stress_due_to_gravity_wave_drag`: Y stress due to gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `tendency_of_x_wind_due_to_gravity_wave_drag`: Tendency of x wind due to gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_wind_due_to_gravity_wave_drag`: Tendency of y wind due to gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_air_temperature_due_to_gravity_wave_drag`: Tendency of air temperature due to gravity wave drag
    * `real(kind=kind_phys)`: units = K s-1
* `atmosphere_momentum_diffusivity_due_to_gravity_wave_drag`: Atmosphere momentum diffusivity due to gravity wave drag
    * `real(kind=kind_phys)`: units = m2 s-1
* `instantaneous_momentum_flux_due_to_turbulent_orographic_form_drag`: Instantaneous momentum flux due to turbulent orographic form drag
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_momentum_flux_due_to_mountain_blocking_drag`: Instantaneous momentum flux due to mountain blocking drag
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_momentum_flux_due_to_orographic_gravity_wave_drag`: Instantaneous momentum flux due to orographic gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_momentum_flux_due_to_nonstationary_gravity_wave`: Instantaneous momentum flux due to nonstationary gravity wave
    * `real(kind=kind_phys)`: units = Pa
* `height_of_mountain_blocking`: Height of mountain blocking
    * `real(kind=kind_phys)`: units = m
* `height_of_low_level_wave_breaking`: Height of low level wave breaking
    * `real(kind=kind_phys)`: units = m
* `height_of_launch_level_of_orographic_gravity_wave`: Height of launch level of orographic gravity wave
    * `real(kind=kind_phys)`: units = m
* `instantaneous_change_in_x_wind_due_to_mountain_blocking_drag`: Instantaneous change in x wind due to mountain blocking drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_x_wind_due_to_mesoscale_orographic_gravity_wave_drag`: Tendency of x wind due to mesoscale orographic gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_x_wind_due_to_turbulent_orographic_form_drag`: Tendency of x wind due to turbulent orographic form drag
    * `real(kind=kind_phys)`: units = m s-2
* `time_integral_of_change_in_x_wind_due_to_mountain_blocking_drag`: Time integral of change in x wind due to mountain blocking drag
    * `real(kind=kind_phys)`: units = m s-2
* `time_integral_of_change_in_x_wind_due_to_orographic_gravity_wave_drag`: Time integral of change in x wind due to orographic gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `time_integral_of_change_in_x_wind_due_to_turbulent_orographic_form_drag`: Time integral of change in x wind due to turbulent orographic form drag
    * `real(kind=kind_phys)`: units = m s-2
* `level_of_dividing_streamline`: Level of dividing streamline
    * `real(kind=kind_phys)`: units = none
* `pi`: Pi
    * `real(kind=kind_phys)`: units = none
* `angular_velocity_of_earth`: Angular velocity of earth
    * `real(kind=kind_phys)`: units = s-1
* `lwe_thickness_of_precipitation_amount_on_dynamics_timestep`: Lwe thickness of precipitation amount on dynamics timestep
    * `real(kind=kind_phys)`: units = m
* `turbulent_kinetic_energy`: Turbulent kinetic energy
    * `real(kind=kind_phys)`: units = J
* `cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag`: Cumulative change in x wind due to orographic gravity wave drag
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag`: Cumulative change in y wind due to orographic gravity wave drag
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_temperature_due_to_orographic_gravity_wave_drag`: Cumulative change in temperature due to orographic gravity wave drag
    * `real(kind=kind_phys)`: units = K
* `cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag`: Cumulative change in x wind due to convective gravity wave drag
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag`: Cumulative change in y wind due to convective gravity wave drag
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_temperature_due_to_convective_gravity_wave_drag`: Cumulative change in temperature due to convective gravity wave drag
    * `real(kind=kind_phys)`: units = K
* `dimensionless_exner_function_at_surface_adjacent_layer`: Dimensionless exner function at surface adjacent layer
    * `real(kind=kind_phys)`: units = none
* `ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer`: Ratio of exner function between midlayer and interface at lowest model layer
    * `real(kind=kind_phys)`: units = ratio
* `surface_specific_humidity`: Surface specific humidity
    * `real(kind=kind_phys)`: units = kg kg-1
* `surface_specific_humidity_for_MYJ_schemes`: Surface specific humidity for MYJ schemes
    * `real(kind=kind_phys)`: units = kg kg-1
* `air_potential_temperature_at_top_of_viscous_sublayer`: Air potential temperature at top of viscous sublayer
    * `real(kind=kind_phys)`: units = K
* `specific_humidity_at_top_of_viscous_sublayer`: Specific humidity at top of viscous sublayer
    * `real(kind=kind_phys)`: units = kg kg-1
* `x_wind_at_top_of_viscous_sublayer`: X wind at top of viscous sublayer
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind_at_top_of_viscous_sublayer`: Y wind at top of viscous sublayer
    * `real(kind=kind_phys)`: units = m s-1
* `heat_exchange_coefficient_for_MYJ_schemes`: Heat exchange coefficient for MYJ schemes
    * `real(kind=kind_phys)`: units = m s-1
* `momentum_exchange_coefficient_for_MYJ_schemes`: Momentum exchange coefficient for MYJ schemes
    * `real(kind=kind_phys)`: units = m s-1
* `control_for_surface_layer_evaporation`: Control for surface layer evaporation
    * `real(kind=kind_phys)`: units = none
* `kinematic_surface_latent_heat_flux`: Kinematic surface latent heat flux
    * `real(kind=kind_phys)`: units = m s-1 kg kg-1
* `weight_for_momentum_at_top_of_viscous_sublayer`: Weight for momentum at top of viscous sublayer
    * `real(kind=kind_phys)`: units = none
* `weight_for_potental_temperature_at_top_of_viscous_sublayer`: Weight for potental temperature at top of viscous sublayer
    * `real(kind=kind_phys)`: units = none
* `weight_for_specific_humidity_at_top_of_viscous_sublayer`: Weight for specific humidity at top of viscous sublayer
    * `real(kind=kind_phys)`: units = none
* `surface_friction_velocity`: Surface friction velocity
    * `real(kind=kind_phys)`: units = m s-1
* `surface_drag_coefficient_for_momentum_in_air`: Surface drag coefficient for momentum in air
    * `real(kind=kind_phys)`: units = none
* `surface_drag_coefficient_for_heat_and_moisture_in_air`: Surface drag coefficient for heat and moisture in air
    * `real(kind=kind_phys)`: units = none
* `atmosphere_heat_diffusivity`: Atmosphere heat diffusivity
    * `real(kind=kind_phys)`: units = m2 s-1
* `countergradient_mixing_term_for_temperature`: Countergradient mixing term for temperature
    * `real(kind=kind_phys)`: units = K
* `countergradient_mixing_term_for_water_vapor`: Countergradient mixing term for water vapor
    * `real(kind=kind_phys)`: units = kg kg-1
* `land_area_fraction`: Land area fraction
    * `real(kind=kind_phys)`: units = frac
* `lake_area_fraction`: Lake area fraction
    * `real(kind=kind_phys)`: units = frac
* `lake_depth`: Lake depth
    * `real(kind=kind_phys)`: units = m
* `sea_area_fraction`: Sea area fraction
    * `real(kind=kind_phys)`: units = frac
* `land_area_fraction_for_microphysics`: Land area fraction for microphysics
    * `real(kind=kind_phys)`: units = frac
* `sea_ice_thickness`: Sea ice thickness
    * `real(kind=kind_phys)`: units = m
* `surface_snow_thickness_water_equivalent_over_water`: Surface snow thickness water equivalent over water
    * `real(kind=kind_phys)`: units = mm
* `surface_snow_thickness_water_equivalent_over_land`: Surface snow thickness water equivalent over land
    * `real(kind=kind_phys)`: units = mm
* `surface_snow_thickness_water_equivalent_over_ice`: Surface snow thickness water equivalent over ice
    * `real(kind=kind_phys)`: units = mm
* `nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep`: Nonnegative lwe thickness of precipitation amount on dynamics timestep
    * `real(kind=kind_phys)`: units = m
* `nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_water`: Nonnegative lwe thickness of precipitation amount on dynamics timestep over water
    * `real(kind=kind_phys)`: units = m
* `nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_land`: Nonnegative lwe thickness of precipitation amount on dynamics timestep over land
    * `real(kind=kind_phys)`: units = m
* `nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_ice`: Nonnegative lwe thickness of precipitation amount on dynamics timestep over ice
    * `real(kind=kind_phys)`: units = m
* `surface_friction_velocity_over_water`: Surface friction velocity over water
    * `real(kind=kind_phys)`: units = m s-1
* `surface_friction_velocity_over_land`: Surface friction velocity over land
    * `real(kind=kind_phys)`: units = m s-1
* `surface_friction_velocity_over_ice`: Surface friction velocity over ice
    * `real(kind=kind_phys)`: units = m s-1
* `lwe_thickness_of_surface_snow_amount`: Lwe thickness of surface snow amount
    * `real(kind=kind_phys)`: units = mm
* `water_equivalent_accumulated_snow_depth_over_water`: Water equivalent accumulated snow depth over water
    * `real(kind=kind_phys)`: units = mm
* `water_equivalent_accumulated_snow_depth_over_ice`: Water equivalent accumulated snow depth over ice
    * `real(kind=kind_phys)`: units = mm
* `surface_upward_potential_latent_heat_flux_over_ice`: Surface upward potential latent heat flux over ice
    * `real(kind=kind_phys)`: units = W m-2
* `sea_surface_temperature`: Sea surface temperature
    * `real(kind=kind_phys)`: units = K
* `surface_skin_temperature_over_land`: Surface skin temperature over land
    * `real(kind=kind_phys)`: units = K
* `surface_skin_temperature_over_water_interstitial`: Surface skin temperature over water interstitial
    * `real(kind=kind_phys)`: units = K
* `surface_skin_temperature_over_land_interstitial`: Surface skin temperature over land interstitial
    * `real(kind=kind_phys)`: units = K
* `surface_skin_temperature_over_ice_interstitial`: Surface skin temperature over ice interstitial
    * `real(kind=kind_phys)`: units = K
* `sea_ice_temperature_interstitial`: Sea ice temperature interstitial
    * `real(kind=kind_phys)`: units = K
* `surface_skin_temperature_after_iteration`: Surface skin temperature after iteration
    * `real(kind=kind_phys)`: units = K
* `surface_skin_temperature_after_iteration_over_water`: Surface skin temperature after iteration over water
    * `real(kind=kind_phys)`: units = K
* `surface_skin_temperature_after_iteration_over_ice`: Surface skin temperature after iteration over ice
    * `real(kind=kind_phys)`: units = K
* `upward_heat_flux_in_soil_over_ice`: Upward heat flux in soil over ice
    * `real(kind=kind_phys)`: units = W m-2
* `freezing_point_temperature_of_seawater`: Freezing point temperature of seawater
    * `real(kind=kind_phys)`: units = K
* `surface_longwave_emissivity`: Surface longwave emissivity
    * `real(kind=kind_phys)`: units = frac
* `surface_longwave_emissivity_over_water_interstitial`: Surface longwave emissivity over water interstitial
    * `real(kind=kind_phys)`: units = frac
* `surface_longwave_emissivity_over_ice_interstitial`: Surface longwave emissivity over ice interstitial
    * `real(kind=kind_phys)`: units = frac
* `surface_specific_humidity_over_water`: Surface specific humidity over water
    * `real(kind=kind_phys)`: units = kg kg-1
* `surface_specific_humidity_over_ice`: Surface specific humidity over ice
    * `real(kind=kind_phys)`: units = kg kg-1
* `kinematic_surface_upward_sensible_heat_flux`: Kinematic surface upward sensible heat flux
    * `real(kind=kind_phys)`: units = K m s-1
* `kinematic_surface_upward_sensible_heat_flux_over_water`: Kinematic surface upward sensible heat flux over water
    * `real(kind=kind_phys)`: units = K m s-1
* `kinematic_surface_upward_sensible_heat_flux_over_land`: Kinematic surface upward sensible heat flux over land
    * `real(kind=kind_phys)`: units = K m s-1
* `kinematic_surface_upward_sensible_heat_flux_over_ice`: Kinematic surface upward sensible heat flux over ice
    * `real(kind=kind_phys)`: units = K m s-1
* `min_lake_ice_area_fraction`: Min lake ice area fraction
    * `real(kind=kind_phys)`: units = frac
* `min_sea_ice_area_fraction`: Min sea ice area fraction
    * `real(kind=kind_phys)`: units = frac
* `surface_roughness_length_over_water`: Surface roughness length over water
    * `real(kind=kind_phys)`: units = cm
* `surface_roughness_length_over_ice`: Surface roughness length over ice
    * `real(kind=kind_phys)`: units = cm
* `surface_downwelling_longwave_flux`: Surface downwelling longwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_longwave_flux_absorbed_by_surface_over_ice`: Surface downwelling longwave flux absorbed by surface over ice
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_longwave_flux_absorbed_by_ground_over_water`: Surface downwelling longwave flux absorbed by ground over water
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_shortwave_flux`: Surface upwelling shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_drag_coefficient_for_momentum_in_air_over_water`: Surface drag coefficient for momentum in air over water
    * `real(kind=kind_phys)`: units = none
* `surface_drag_coefficient_for_momentum_in_air_over_land`: Surface drag coefficient for momentum in air over land
    * `real(kind=kind_phys)`: units = none
* `surface_drag_coefficient_for_momentum_in_air_over_ice`: Surface drag coefficient for momentum in air over ice
    * `real(kind=kind_phys)`: units = none
* `surface_drag_coefficient_for_heat_and_moisture_in_air_over_water`: Surface drag coefficient for heat and moisture in air over water
    * `real(kind=kind_phys)`: units = none
* `surface_drag_coefficient_for_heat_and_moisture_in_air_over_land`: Surface drag coefficient for heat and moisture in air over land
    * `real(kind=kind_phys)`: units = none
* `surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice`: Surface drag coefficient for heat and moisture in air over ice
    * `real(kind=kind_phys)`: units = none
* `bulk_richardson_number_at_lowest_model_level_over_water`: Bulk richardson number at lowest model level over water
    * `real(kind=kind_phys)`: units = none
* `bulk_richardson_number_at_lowest_model_level_over_ice`: Bulk richardson number at lowest model level over ice
    * `real(kind=kind_phys)`: units = none
* `surface_wind_stress_over_water`: Surface wind stress over water
    * `real(kind=kind_phys)`: units = m2 s-2
* `surface_wind_stress_over_land`: Surface wind stress over land
    * `real(kind=kind_phys)`: units = m2 s-2
* `surface_wind_stress_over_ice`: Surface wind stress over ice
    * `real(kind=kind_phys)`: units = m2 s-2
* `Monin_Obukhov_similarity_function_for_momentum_over_water`: Monin Obukhov similarity function for momentum over water
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_momentum_over_land`: Monin Obukhov similarity function for momentum over land
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_momentum_over_ice`: Monin Obukhov similarity function for momentum over ice
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_heat_over_water`: Monin Obukhov similarity function for heat over water
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_heat_over_land`: Monin Obukhov similarity function for heat over land
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_heat_over_ice`: Monin Obukhov similarity function for heat over ice
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_momentum_at_10m`: Monin Obukhov similarity function for momentum at 10m
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_momentum_at_10m_over_water`: Monin Obukhov similarity function for momentum at 10m over water
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_momentum_at_10m_over_land`: Monin Obukhov similarity function for momentum at 10m over land
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_momentum_at_10m_over_ice`: Monin Obukhov similarity function for momentum at 10m over ice
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_heat_at_2m`: Monin Obukhov similarity function for heat at 2m
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_heat_at_2m_over_water`: Monin Obukhov similarity function for heat at 2m over water
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_heat_at_2m_over_land`: Monin Obukhov similarity function for heat at 2m over land
    * `real(kind=kind_phys)`: units = none
* `Monin_Obukhov_similarity_function_for_heat_at_2m_over_ice`: Monin Obukhov similarity function for heat at 2m over ice
    * `real(kind=kind_phys)`: units = none
* `surface_drag_wind_speed_for_momentum_in_air`: Surface drag wind speed for momentum in air
    * `real(kind=kind_phys)`: units = m s-1
* `surface_drag_wind_speed_for_momentum_in_air_over_water`: Surface drag wind speed for momentum in air over water
    * `real(kind=kind_phys)`: units = m s-1
* `surface_drag_wind_speed_for_momentum_in_air_over_land`: Surface drag wind speed for momentum in air over land
    * `real(kind=kind_phys)`: units = m s-1
* `surface_drag_wind_speed_for_momentum_in_air_over_ice`: Surface drag wind speed for momentum in air over ice
    * `real(kind=kind_phys)`: units = m s-1
* `surface_drag_mass_flux_for_heat_and_moisture_in_air`: Surface drag mass flux for heat and moisture in air
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `surface_drag_mass_flux_for_heat_and_moisture_in_air_over_water`: Surface drag mass flux for heat and moisture in air over water
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `surface_drag_mass_flux_for_heat_and_moisture_in_air_over_land`: Surface drag mass flux for heat and moisture in air over land
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ice`: Surface drag mass flux for heat and moisture in air over ice
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `upward_heat_flux_in_soil`: Upward heat flux in soil
    * `real(kind=kind_phys)`: units = W m-2
* `upward_heat_flux_in_soil_over_water`: Upward heat flux in soil over water
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upward_potential_latent_heat_flux`: Surface upward potential latent heat flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upward_potential_latent_heat_flux_over_water`: Surface upward potential latent heat flux over water
    * `real(kind=kind_phys)`: units = W m-2
* `kinematic_surface_upward_latent_heat_flux`: Kinematic surface upward latent heat flux
    * `real(kind=kind_phys)`: units = kg kg-1 m s-1
* `kinematic_surface_upward_latent_heat_flux_over_water`: Kinematic surface upward latent heat flux over water
    * `real(kind=kind_phys)`: units = kg kg-1 m s-1
* `kinematic_surface_upward_latent_heat_flux_over_land`: Kinematic surface upward latent heat flux over land
    * `real(kind=kind_phys)`: units = kg kg-1 m s-1
* `kinematic_surface_upward_latent_heat_flux_over_ice`: Kinematic surface upward latent heat flux over ice
    * `real(kind=kind_phys)`: units = kg kg-1 m s-1
* `temperature_in_ice_layer`: Temperature in ice layer
    * `real(kind=kind_phys)`: units = K
* `cloud_condensed_water_mixing_ratio`: Cloud condensed water mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `ice_water_mixing_ratio`: Ice water mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `cloud_area_fracton`: Cloud area fracton
    * `real(kind=kind_phys)`: units = frac
* `cloud_liquid_water_path`: Cloud liquid water path
    * `real(kind=kind_phys)`: units = g m-2
* `effective_radius_of_cloud_liquid_water_particles_in_um`: Effective radius of cloud liquid water particles in um
    * `real(kind=kind_phys)`: units = um
* `cloud_ice_water_path`: Cloud ice water path
    * `real(kind=kind_phys)`: units = g m-2
* `effective_radius_of_ice_liquid_water_particles_in_um`: Effective radius of ice liquid water particles in um
    * `real(kind=kind_phys)`: units = um
* `cloud_snow_water_path`: Cloud snow water path
    * `real(kind=kind_phys)`: units = g m-2
* `effective_radius_of_cloud_snow_particles_in_um`: Effective radius of cloud snow particles in um
    * `real(kind=kind_phys)`: units = um
* `cloud_rain_water_path`: Cloud rain water path
    * `real(kind=kind_phys)`: units = g m-2
* `effective_radius_of_cloud_rain_particles_in_um`: Effective radius of cloud rain particles in um
    * `real(kind=kind_phys)`: units = um
* `precipitation_fraction_by_layer`: Precipitation fraction by layer
    * `real(kind=kind_phys)`: units = frac
* `RRTMGP_cloud_optical_depth_layers_at_10mu_band`: RRTMGP cloud optical depth layers at 10mu band
    * `real(kind=kind_phys)`: units = none
* `latitude_interpolation_weight_for_ozone_forcing`: Latitude interpolation weight for ozone forcing
    * `real(kind=kind_phys)`: units = none
* `ozone_forcing`: Ozone forcing
    * `real(kind=kind_phys)`: units = various
* `latitude_interpolation_weight_for_stratospheric_water_vapor_forcing`: Latitude interpolation weight for stratospheric water vapor forcing
    * `real(kind=kind_phys)`: units = none
* `stratospheric_water_vapor_forcing`: Stratospheric water vapor forcing
    * `real(kind=kind_phys)`: units = various
* `latitude_interpolation_weight_for_aerosol_forcing`: Latitude interpolation weight for aerosol forcing
    * `real(kind=kind_phys)`: units = none
* `longitude_interpolation_weight_for_aerosol_forcing`: Longitude interpolation weight for aerosol forcing
    * `real(kind=kind_phys)`: units = none
* `mass_number_concentration_of_aerosol_from_gocart_climatology`: Mass number concentration of aerosol from gocart climatology
    * `real(kind=kind_phys)`: units = kg-1
* `latitude_interpolation_weight_for_cloud_nuclei_forcing`: Latitude interpolation weight for cloud nuclei forcing
    * `real(kind=kind_phys)`: units = none
* `longitude_interpolation_weight_for_cloud_nuclei_forcing`: Longitude interpolation weight for cloud nuclei forcing
    * `real(kind=kind_phys)`: units = none
* `latitude_interpolation_weight_complement_for_absolute_momentum_flux_due_to_nonorographic_gravity_wave_drag`: Latitude interpolation weight complement for absolute momentum flux due to nonorographic gravity wave drag
    * `real(kind=kind_phys)`: units = none
* `latitude_interpolation_weight_for_absolute_momentum_flux_due_to_nonorographic_gravity_wave_drag`: Latitude interpolation weight for absolute momentum flux due to nonorographic gravity wave drag
    * `real(kind=kind_phys)`: units = none
* `surface_snow_area_fraction_over_ice`: Surface snow area fraction over ice
    * `real(kind=kind_phys)`: units = frac
* `vegetation_type_classification_real`: Vegetation type classification real
    * `real(kind=kind_phys)`: units = index
* `depth_of_soil_levels_for_land_surface_model`: Depth of soil levels for land surface model
    * `real(kind=kind_phys)`: units = m
* `canopy_temperature`: Canopy temperature
    * `real(kind=kind_phys)`: units = K
* `ground_temperature_for_noahmp`: Ground temperature for noahmp
    * `real(kind=kind_phys)`: units = K
* `air_temperature_in_canopy`: Air temperature in canopy
    * `real(kind=kind_phys)`: units = K
* `canopy_intercepted_ice_mass`: Canopy intercepted ice mass
    * `real(kind=kind_phys)`: units = mm
* `canopy_intercepted_liquid_water`: Canopy intercepted liquid water
    * `real(kind=kind_phys)`: units = mm
* `air_vapor_pressure_in_canopy`: Air vapor pressure in canopy
    * `real(kind=kind_phys)`: units = Pa
* `surface_drag_coefficient_for_momentum_for_noahmp`: Surface drag coefficient for momentum for noahmp
    * `real(kind=kind_phys)`: units = none
* `surface_drag_coefficient_for_heat_and_moisture_for_noahmp`: Surface drag coefficient for heat and moisture for noahmp
    * `real(kind=kind_phys)`: units = none
* `wet_canopy_area_fraction`: Wet canopy area fraction
    * `real(kind=kind_phys)`: units = none
* `snow_mass_at_previous_time_step`: Snow mass at previous time step
    * `real(kind=kind_phys)`: units = mm
* `surface_albedo_assuming_deep_snow_on_previous_timestep`: Surface albedo assuming deep snow on previous timestep
    * `real(kind=kind_phys)`: units = frac
* `lwe_snowfall_rate`: Lwe snowfall rate
    * `real(kind=kind_phys)`: units = mm s-1
* `water_storage_in_lake`: Water storage in lake
    * `real(kind=kind_phys)`: units = mm
* `dimensionless_age_of_surface_snow`: Dimensionless age of surface snow
    * `real(kind=kind_phys)`: units = none
* `water_storage_in_aquifer`: Water storage in aquifer
    * `real(kind=kind_phys)`: units = mm
* `water_storage_in_aquifer_and_saturated_soil`: Water storage in aquifer and saturated soil
    * `real(kind=kind_phys)`: units = mm
* `water_table_depth`: Water table depth
    * `real(kind=kind_phys)`: units = m
* `leaf_area_index`: Leaf area index
    * `real(kind=kind_phys)`: units = none
* `stem_area_index`: Stem area index
    * `real(kind=kind_phys)`: units = none
* `leaf_mass_content`: Leaf mass content
    * `real(kind=kind_phys)`: units = g m-2
* `stem_mass_content`: Stem mass content
    * `real(kind=kind_phys)`: units = g m-2
* `fine_root_mass_content`: Fine root mass content
    * `real(kind=kind_phys)`: units = g m-2
* `wood_mass_content`: Wood mass content
    * `real(kind=kind_phys)`: units = g m-2
* `slow_soil_pool_mass_content_of_carbon`: Slow soil pool mass content of carbon
    * `real(kind=kind_phys)`: units = g m-2
* `fast_soil_pool_mass_content_of_carbon`: Fast soil pool mass content of carbon
    * `real(kind=kind_phys)`: units = g m-2
* `volumetric_soil_moisture_between_soil_bottom_and_water_table`: Volumetric soil moisture between soil bottom and water table
    * `real(kind=kind_phys)`: units = m3 m-3
* `water_table_recharge_assuming_deep`: Water table recharge assuming deep
    * `real(kind=kind_phys)`: units = m
* `water_table_recharge_assuming_shallow`: Water table recharge assuming shallow
    * `real(kind=kind_phys)`: units = m
* `surface_emissivity_lsm`: Surface emissivity lsm
    * `real(kind=kind_phys)`: units = frac
* `number_of_snow_layers`: Number of snow layers
    * `real(kind=kind_phys)`: units = count
* `snow_layer_ice`: Snow layer ice
    * `real(kind=kind_phys)`: units = mm
* `snow_layer_liquid_water`: Snow layer liquid water
    * `real(kind=kind_phys)`: units = mm
* `temperature_in_surface_snow`: Temperature in surface snow
    * `real(kind=kind_phys)`: units = K
* `volumetric_equilibrium_soil_moisture`: Volumetric equilibrium soil moisture
    * `real(kind=kind_phys)`: units = m3 m-3
* `depth_from_snow_surface_at_bottom_interface`: Depth from snow surface at bottom interface
    * `real(kind=kind_phys)`: units = m
* `canopy_water_amount`: Canopy water amount
    * `real(kind=kind_phys)`: units = kg m-2
* `soil_type_classification_real`: Soil type classification real
    * `real(kind=kind_phys)`: units = index
* `temperature_at_zero_celsius`: Temperature at zero celsius
    * `real(kind=kind_phys)`: units = K
* `period_of_shortwave_radiation_calls`: Period of shortwave radiation calls
    * `real(kind=kind_phys)`: units = s
* `forecast_time`: Forecast time
    * `real(kind=kind_phys)`: units = h
* `control_for_convective_cloud_diagnostics`: Control for convective cloud diagnostics
    * `real(kind=kind_phys)`: units = none
* `ice_nucleation_number_from_climatology`: Ice nucleation number from climatology
    * `real(kind=kind_phys)`: units = kg-1
* `tendency_of_activated_cloud_condensation_nuclei_from_climatology`: Tendency of activated cloud condensation nuclei from climatology
    * `real(kind=kind_phys)`: units = kg-1 s-1
* `random_number`: Random number
    * `real(kind=kind_phys)`: units = none
* `frequency_for_surface_cycling_calls`: Frequency for surface cycling calls
    * `real(kind=kind_phys)`: units = h
* `forecast_time_on_previous_timestep`: Forecast time on previous timestep
    * `real(kind=kind_phys)`: units = h
* `volume_fraction_of_soil_moisture_for_land_surface_model`: Volume fraction of soil moisture for land surface model
    * `real(kind=kind_phys)`: units = frac
* `volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model`: Volume fraction of unfrozen soil moisture for land surface model
    * `real(kind=kind_phys)`: units = frac
* `soil_temperature_for_land_surface_model`: Soil temperature for land surface model
    * `real(kind=kind_phys)`: units = K
* `reference_sea_surface_temperature`: Reference sea surface temperature
    * `real(kind=kind_phys)`: units = K
* `surface_slope_classification_real`: Surface slope classification real
    * `real(kind=kind_phys)`: units = index
* `vegetation_area_fraction`: Vegetation area fraction
    * `real(kind=kind_phys)`: units = frac
* `convective_cloud_area_fraction_between_sw_radiation_calls_from_cnvc90`: Convective cloud area fraction between sw radiation calls from cnvc90
    * `real(kind=kind_phys)`: units = frac
* `pressure_at_convective_cloud_base_between_sw_radiation_calls_from_cnvc90`: Pressure at convective cloud base between sw radiation calls from cnvc90
    * `real(kind=kind_phys)`: units = Pa
* `pressure_at_convective_cloud_top_between_sw_radiation_calls_from_cnvc90`: Pressure at convective cloud top between sw radiation calls from cnvc90
    * `real(kind=kind_phys)`: units = Pa
* `absolute_momentum_flux_due_to_nonorographic_gravity_wave_drag`: Absolute momentum flux due to nonorographic gravity wave drag
    * `real(kind=kind_phys)`: units = various
* `air_temperature_at_layer_for_RRTMGP`: Air temperature at layer for RRTMGP
    * `real(kind=kind_phys)`: units = K
* `RRTMGP_sw_flux_profile_upward_allsky`: RRTMGP sw flux profile upward allsky
    * `real(kind=kind_phys)`: units = W m-2
* `RRTMGP_sw_flux_profile_downward_allsky`: RRTMGP sw flux profile downward allsky
    * `real(kind=kind_phys)`: units = W m-2
* `RRTMGP_sw_flux_profile_upward_clrsky`: RRTMGP sw flux profile upward clrsky
    * `real(kind=kind_phys)`: units = W m-2
* `RRTMGP_sw_flux_profile_downward_clrsky`: RRTMGP sw flux profile downward clrsky
    * `real(kind=kind_phys)`: units = W m-2
* `timestep_for_radiation`: Timestep for radiation
    * `real(kind=kind_phys)`: units = s
* `atmosphere_optical_thickness_due_to_ambient_aerosol_particles`: Atmosphere optical thickness due to ambient aerosol particles
    * `real(kind=kind_phys)`: units = none
* `cloud_area_fraction_in_atmosphere_layer`: Cloud area fraction in atmosphere layer
    * `real(kind=kind_phys)`: units = frac
* `RRTMGP_cloud_optical_depth_layers_at_0_55mu_band`: RRTMGP cloud optical depth layers at 0 55mu band
    * `real(kind=kind_phys)`: units = none
* `diagnostics_for_shortwave_and_longwave_radiation`: Diagnostics for shortwave and longwave radiation
    * `real(kind=kind_phys)`: units = various
* `surface_downwelling_direct_nir_shortwave_flux_on_radiation_timestep`: Surface downwelling direct nir shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_diffuse_nir_shortwave_flux_on_radiation_timestep`: Surface downwelling diffuse nir shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_direct_uv_and_vis_shortwave_flux_on_radiation_timestep`: Surface downwelling direct uv and vis shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_diffuse_uv_and_vis_shortwave_flux_on_radiation_timestep`: Surface downwelling diffuse uv and vis shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_direct_nir_shortwave_flux_on_radiation_timestep`: Surface upwelling direct nir shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_diffuse_nir_shortwave_flux_on_radiation_timestep`: Surface upwelling diffuse nir shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_direct_uv_and_vis_shortwave_flux_on_radiation_timestep`: Surface upwelling direct uv and vis shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_diffuse_uv_and_vis_shortwave_flux_on_radiation_timestep`: Surface upwelling diffuse uv and vis shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_net_downwelling_shortwave_flux_on_radiation_timestep`: Surface net downwelling shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_shortwave_flux_on_radiation_timestep`: Surface downwelling shortwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_timestep`: Tendency of air temperature due to shortwave heating assuming clear sky on radiation timestep
    * `real(kind=kind_phys)`: units = K s-1
* `sigma_pressure_hybrid_vertical_coordinate`: Sigma pressure hybrid vertical coordinate
    * `real(kind=kind_phys)`: units = none
* `equation_of_time`: Equation of time
    * `real(kind=kind_phys)`: units = radian
* `sine_of_solar_declination_angle`: Sine of solar declination angle
    * `real(kind=kind_phys)`: units = none
* `cosine_of_solar_declination_angle`: Cosine of solar declination angle
    * `real(kind=kind_phys)`: units = none
* `solar_constant`: Solar constant
    * `real(kind=kind_phys)`: units = W m-2
* `radar_reflectivity_10cm`: Radar reflectivity 10cm
    * `real(kind=kind_phys)`: units = dBZ
* `maximum_reflectivity_at_1km_agl_over_maximum_hourly_time_interval`: Maximum reflectivity at 1km agl over maximum hourly time interval
    * `real(kind=kind_phys)`: units = dBZ
* `maximum_reflectivity_at_minus10c_over_maximum_hourly_time_interval`: Maximum reflectivity at minus10c over maximum hourly time interval
    * `real(kind=kind_phys)`: units = dBZ
* `maximum_u_wind_at_10m_over_maximum_hourly_time_interval`: Maximum u wind at 10m over maximum hourly time interval
    * `real(kind=kind_phys)`: units = m s-1
* `maximum_v_wind_at_10m_over_maximum_hourly_time_interval`: Maximum v wind at 10m over maximum hourly time interval
    * `real(kind=kind_phys)`: units = m s-1
* `maximum_wind_at_10m_over_maximum_hourly_time_interval`: Maximum wind at 10m over maximum hourly time interval
    * `real(kind=kind_phys)`: units = m s-1
* `air_temperature_at_2m`: Air temperature at 2m
    * `real(kind=kind_phys)`: units = K
* `specific_humidity_at_2m`: Specific humidity at 2m
    * `real(kind=kind_phys)`: units = kg kg-1
* `maximum_temperature_at_2m_over_maximum_hourly_time_interval`: Maximum temperature at 2m over maximum hourly time interval
    * `real(kind=kind_phys)`: units = K
* `minimum_temperature_at_2m_over_maximum_hourly_time_interval`: Minimum temperature at 2m over maximum hourly time interval
    * `real(kind=kind_phys)`: units = K
* `maximum_relative_humidity_at_2m_over_maximum_hourly_time_interval`: Maximum relative humidity at 2m over maximum hourly time interval
    * `real(kind=kind_phys)`: units = %
* `minimum_relative_humidity_at_2m_over_maximum_hourly_time_interval`: Minimum relative humidity at 2m over maximum hourly time interval
    * `real(kind=kind_phys)`: units = %
* `specific_heat_of_water_vapor_at_constant_pressure`: Specific heat of water vapor at constant pressure
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `chemical_tracer_scavenging_fractions`: Chemical tracer scavenging fractions
    * `real(kind=kind_phys)`: units = none
* `lwe_thickness_of_shallow_convective_precipitation_amount`: Lwe thickness of shallow convective precipitation amount
    * `real(kind=kind_phys)`: units = m
* `convective_cloud_water_mixing_ratio`: Convective cloud water mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `convective_cloud_area_fraction`: Convective cloud area fraction
    * `real(kind=kind_phys)`: units = frac
* `entrainment_rate_coefficient_for_shallow_convection`: Entrainment rate coefficient for shallow convection
    * `real(kind=kind_phys)`: units = none
* `rain_conversion_parameter_for_shallow_convection`: Rain conversion parameter for shallow convection
    * `real(kind=kind_phys)`: units = m-1
* `detrainment_conversion_parameter_for_shallow_convection`: Detrainment conversion parameter for shallow convection
    * `real(kind=kind_phys)`: units = m-1
* `momentum_transport_reduction_factor_due_to_pressure_gradient_force_for_shallow_convection`: Momentum transport reduction factor due to pressure gradient force for shallow convection
    * `real(kind=kind_phys)`: units = frac
* `aerosol_aware_multiplicative_rain_conversion_parameter_for_shallow_convection`: Aerosol aware multiplicative rain conversion parameter for shallow convection
    * `real(kind=kind_phys)`: units = none
* `characteristic_grid_lengthscale`: Characteristic grid lengthscale
    * `real(kind=kind_phys)`: units = m
* `mass_number_concentration_of_cloud_droplets`: Mass number concentration of cloud droplets
    * `real(kind=kind_phys)`: units = kg-1
* `mass_number_concentration_of_cloud_ice`: Mass number concentration of cloud ice
    * `real(kind=kind_phys)`: units = kg-1
* `ozone_mixing_ratio`: Ozone mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `mass_number_concentration_of_hygroscopic_aerosols`: Mass number concentration of hygroscopic aerosols
    * `real(kind=kind_phys)`: units = kg-1
* `mass_number_concentration_of_nonhygroscopic_ice_nucleating_aerosols`: Mass number concentration of nonhygroscopic ice nucleating aerosols
    * `real(kind=kind_phys)`: units = kg-1
* `instantaneous_surface_x_momentum_flux_for_diag`: Instantaneous surface x momentum flux for diag
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_surface_y_momentum_flux_for_diag`: Instantaneous surface y momentum flux for diag
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_surface_upward_sensible_heat_flux_for_diag`: Instantaneous surface upward sensible heat flux for diag
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_upward_latent_heat_flux_for_diag`: Instantaneous surface upward latent heat flux for diag
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep`: Cumulative surface x momentum flux for diag multiplied by timestep
    * `real(kind=kind_phys)`: units = Pa s
* `cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep`: Cumulative surface y momentum flux for diag multiplied by timestep
    * `real(kind=kind_phys)`: units = Pa s
* `cumulative_surface_upward_sensible_heat_flux_for_diag_multiplied_by_timestep`: Cumulative surface upward sensible heat flux for diag multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_upward_latent_heat_flux_for_diag_multiplied_by_timestep`: Cumulative surface upward latent heat flux for diag multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `surface_x_momentum_flux_from_coupled_process`: Surface x momentum flux from coupled process
    * `real(kind=kind_phys)`: units = Pa
* `surface_y_momentum_flux_from_coupled_process`: Surface y momentum flux from coupled process
    * `real(kind=kind_phys)`: units = Pa
* `surface_upward_sensible_heat_flux_from_coupled_process`: Surface upward sensible heat flux from coupled process
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upward_latent_heat_flux_from_coupled_process`: Surface upward latent heat flux from coupled process
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_x_momentum_flux_for_coupling`: Instantaneous surface x momentum flux for coupling
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_surface_y_momentum_flux_for_coupling`: Instantaneous surface y momentum flux for coupling
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_surface_upward_sensible_heat_flux_for_coupling`: Instantaneous surface upward sensible heat flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_upward_latent_heat_flux_for_coupling`: Instantaneous surface upward latent heat flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_surface_x_momentum_flux_for_coupling_multiplied_by_timestep`: Cumulative surface x momentum flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = Pa s
* `cumulative_surface_y_momentum_flux_for_coupling_multiplied_by_timestep`: Cumulative surface y momentum flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = Pa s
* `cumulative_surface_upward_sensible_heat_flux_for_coupling_multiplied_by_timestep`: Cumulative surface upward sensible heat flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_upward_latent_heat_flux_for_coupling_multiplied_by_timestep`: Cumulative surface upward latent heat flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `reciprocal_of_obukhov_length`: Reciprocal of obukhov length
    * `real(kind=kind_phys)`: units = m-1
* `nonadvected_turbulent_kinetic_energy_multiplied_by_2`: Nonadvected turbulent kinetic energy multiplied by 2
    * `real(kind=kind_phys)`: units = m2 s-2
* `variance_of_air_temperature`: Variance of air temperature
    * `real(kind=kind_phys)`: units = K2
* `variance_of_specific_humidity`: Variance of specific humidity
    * `real(kind=kind_phys)`: units = kg2 kg-2
* `covariance_of_air_temperature_and_specific_humidity`: Covariance of air temperature and specific humidity
    * `real(kind=kind_phys)`: units = K kg kg-1
* `turbulent_mixing_length`: Turbulent mixing length
    * `real(kind=kind_phys)`: units = m
* `stability_function_for_heat`: Stability function for heat
    * `real(kind=kind_phys)`: units = none
* `atmosphere_heat_diffusivity_for_mynnpbl`: Atmosphere heat diffusivity for mynnpbl
    * `real(kind=kind_phys)`: units = m2 s-1
* `atmosphere_momentum_diffusivity_for_mynnpbl`: Atmosphere momentum diffusivity for mynnpbl
    * `real(kind=kind_phys)`: units = m2 s-1
* `subgrid_cloud_water_mixing_ratio_pbl`: Subgrid cloud water mixing ratio pbl
    * `real(kind=kind_phys)`: units = kg kg-1
* `subgrid_scale_cloud_ice_mixing_ratio`: Subgrid scale cloud ice mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `subgrid_scale_cloud_area_fraction_in_atmosphere_layer`: Subgrid scale cloud area fraction in atmosphere layer
    * `real(kind=kind_phys)`: units = frac
* `emdf_updraft_area`: Emdf updraft area
    * `real(kind=kind_phys)`: units = frac
* `emdf_updraft_vertical_velocity`: Emdf updraft vertical velocity
    * `real(kind=kind_phys)`: units = m s-1
* `emdf_updraft_total_water`: Emdf updraft total water
    * `real(kind=kind_phys)`: units = kg kg-1
* `emdf_updraft_theta_l`: Emdf updraft theta l
    * `real(kind=kind_phys)`: units = K
* `emdf_updraft_entrainment_rate`: Emdf updraft entrainment rate
    * `real(kind=kind_phys)`: units = s-1
* `emdf_updraft_cloud_water`: Emdf updraft cloud water
    * `real(kind=kind_phys)`: units = kg kg-1
* `theta_subsidence_tendency`: Theta subsidence tendency
    * `real(kind=kind_phys)`: units = K s-1
* `water_vapor_subsidence_tendency`: Water vapor subsidence tendency
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `theta_detrainment_tendency`: Theta detrainment tendency
    * `real(kind=kind_phys)`: units = K s-1
* `water_vapor_detrainment_tendency`: Water vapor detrainment tendency
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `maximum_mass_flux`: Maximum mass flux
    * `real(kind=kind_phys)`: units = m s-1
* `tendency_of_water_vapor_specific_humidity_due_to_model_physics`: Tendency of water vapor specific humidity due to model physics
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `tendency_of_cloud_liquid_water_mixing_ratio_due_to_model_physics`: Tendency of cloud liquid water mixing ratio due to model physics
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `tendency_of_cloud_ice_mixing_ratio_due_to_model_physics`: Tendency of cloud ice mixing ratio due to model physics
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `tendency_of_ozone_mixing_ratio_due_to_model_physics`: Tendency of ozone mixing ratio due to model physics
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `tendency_of_mass_number_concentration_of_cloud_liquid_water_particles_in_air_due_to_model_physics`: Tendency of mass number concentration of cloud liquid water particles in air due to model physics
    * `real(kind=kind_phys)`: units = kg-1 s-1
* `tendency_of_mass_number_concentration_of_ice_crystals_in_air_due_to_model_physics`: Tendency of mass number concentration of ice crystals in air due to model physics
    * `real(kind=kind_phys)`: units = kg-1 s-1
* `tendency_of_water_friendly_aerosol_number_concentration_due_to_model_physics`: Tendency of water friendly aerosol number concentration due to model physics
    * `real(kind=kind_phys)`: units = kg-1 s-1
* `tendency_of_ice_friendly_aerosol_mass_number_concentration_due_to_model_physics`: Tendency of ice friendly aerosol mass number concentration due to model physics
    * `real(kind=kind_phys)`: units = kg-1 s-1
* `forecast_time_in_seconds`: Forecast time in seconds
    * `real(kind=kind_phys)`: units = s
* `time_elapsed_since_diagnostics_reset`: Time elapsed since diagnostics reset
    * `real(kind=kind_phys)`: units = h
* `forecast_julian_day`: Forecast julian day
    * `real(kind=kind_phys)`: units = days
* `max_tendency_of_air_potential_temperature_due_to_large_scale_precipitation`: Max tendency of air potential temperature due to large scale precipitation
    * `real(kind=kind_phys)`: units = K s-1
* `air_temperature_save`: Air temperature save
    * `real(kind=kind_phys)`: units = K
* `period_of_longwave_radiation_calls`: Period of longwave radiation calls
    * `real(kind=kind_phys)`: units = s
* `minimum_value_of_saturation_mixing_ratio`: Minimum value of saturation mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `minimum_pressure_in_RRTMGP`: Minimum pressure in RRTMGP
    * `real(kind=kind_phys)`: units = Pa
* `minimum_temperature_in_RRTMGP`: Minimum temperature in RRTMGP
    * `real(kind=kind_phys)`: units = K
* `air_temperature_at_interface_for_RRTMGP`: Air temperature at interface for RRTMGP
    * `real(kind=kind_phys)`: units = K
* `saturation_vapor_pressure`: Saturation vapor pressure
    * `real(kind=kind_phys)`: units = Pa
* `water_vapor_mixing_ratio`: Water vapor mixing ratio
    * `real(kind=kind_phys)`: units = kg/kg
* `chemical_tracers`: Chemical tracers
    * `real(kind=kind_phys)`: units = g g-1
* `lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep`: Lwe thickness of convective precipitation amount on dynamics timestep
    * `real(kind=kind_phys)`: units = m
* `cumulative_lwe_thickness_of_convective_precipitation_amount_between_sw_radiation_calls`: Cumulative lwe thickness of convective precipitation amount between sw radiation calls
    * `real(kind=kind_phys)`: units = m
* `cumulative_min_vertical_index_at_cloud_base_between_sw_radiation_calls`: Cumulative min vertical index at cloud base between sw radiation calls
    * `real(kind=kind_phys)`: units = index
* `cumulative_max_vertical_index_at_cloud_base_between_sw_radiation_calls`: Cumulative max vertical index at cloud base between sw radiation calls
    * `real(kind=kind_phys)`: units = index
* `ozone_concentration_of_new_state`: Ozone concentration of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `natural_log_of_ozone_forcing_data_pressure_levels`: Natural log of ozone forcing data pressure levels
    * `real(kind=kind_phys)`: units = log(Pa)
* `cumulative_change_in_ozone_concentration_due_to_production_and_loss_rate`: Cumulative change in ozone concentration due to production and loss rate
    * `real(kind=kind_phys)`: units = kg kg-1
* `cumulative_change_in_ozone_concentration_due_to_ozone_mixing_ratio`: Cumulative change in ozone concentration due to ozone mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `cumulative_change_in_ozone_concentration_due_to_temperature`: Cumulative change in ozone concentration due to temperature
    * `real(kind=kind_phys)`: units = kg kg-1
* `cumulative_change_in_ozone_concentration_due_to_overhead_ozone_column`: Cumulative change in ozone concentration due to overhead ozone column
    * `real(kind=kind_phys)`: units = kg kg-1
* `x_wind_at_surface_adjacent_layer`: X wind at surface adjacent layer
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind_at_surface_adjacent_layer`: Y wind at surface adjacent layer
    * `real(kind=kind_phys)`: units = m s-1
* `specific_humidity_at_surface_adjacent_layer`: Specific humidity at surface adjacent layer
    * `real(kind=kind_phys)`: units = kg kg-1
* `specified_kinematic_surface_upward_sensible_heat_flux`: Specified kinematic surface upward sensible heat flux
    * `real(kind=kind_phys)`: units = K m s-1
* `specified_kinematic_surface_upward_latent_heat_flux`: Specified kinematic surface upward latent heat flux
    * `real(kind=kind_phys)`: units = kg kg-1 m s-1
* `vonKarman_constant`: VonKarman constant
    * `real(kind=kind_phys)`: units = none
* `tunable_parameter_for_ice_supersaturation`: Tunable parameter for ice supersaturation
    * `real(kind=kind_phys)`: units = none
* `ratio_of_gas_constant_dry_air_to_gravitational_acceleration`: Ratio of gas constant dry air to gravitational acceleration
    * `real(kind=kind_phys)`: units = J s2 K-1 kg-1 m-1
* `ratio_of_gas_constant_dry_air_to_specific_heat_of_dry_air_at_constant_pressure`: Ratio of gas constant dry air to specific heat of dry air at constant pressure
    * `real(kind=kind_phys)`: units = none
* `cloud_fraction_for_MG`: Cloud fraction for MG
    * `real(kind=kind_phys)`: units = frac
* `effective_radius_of_stratiform_cloud_rain_particle`: Effective radius of stratiform cloud rain particle
    * `real(kind=kind_phys)`: units = um
* `sppt_weights_from_coupled_process`: Sppt weights from coupled process
    * `real(kind=kind_phys)`: units = none
* `total_ampltiude_of_sppt_perturbation`: Total ampltiude of sppt perturbation
    * `real(kind=kind_phys)`: units = none
* `convective_cloud_water_mixing_ratio_in_xyz_dimensioned_restart_array`: Convective cloud water mixing ratio in xyz dimensioned restart array
    * `real(kind=kind_phys)`: units = kg kg-1
* `convective_cloud_area_fraction_in_xyz_dimensioned_restart_array`: Convective cloud area fraction in xyz dimensioned restart array
    * `real(kind=kind_phys)`: units = frac
* `effective_radius_of_stratiform_cloud_liquid_water_particle`: Effective radius of stratiform cloud liquid water particle
    * `real(kind=kind_phys)`: units = um
* `effective_radius_of_stratiform_cloud_ice_particle`: Effective radius of stratiform cloud ice particle
    * `real(kind=kind_phys)`: units = um
* `effective_radius_of_stratiform_cloud_snow_particle`: Effective radius of stratiform cloud snow particle
    * `real(kind=kind_phys)`: units = um
* `cloud_decorrelation_length`: Cloud decorrelation length
    * `real(kind=kind_phys)`: units = km
* `perturbation_of_surface_albedo`: Perturbation of surface albedo
    * `real(kind=kind_phys)`: units = frac
* `layer_pressure_thickness`: Layer pressure thickness
    * `real(kind=kind_phys)`: units = hPa
* `layer_thickness_for_radiation`: Layer thickness for radiation
    * `real(kind=kind_phys)`: units = km
* `air_pressure_at_interface_in_hPa`: Air pressure at interface in hPa
    * `real(kind=kind_phys)`: units = hPa
* `air_pressure_in_hPa`: Air pressure in hPa
    * `real(kind=kind_phys)`: units = hPa
* `air_temperature_at_interface_for_radiation`: Air temperature at interface for radiation
    * `real(kind=kind_phys)`: units = K
* `air_temperature_at_layer_for_radiation`: Air temperature at layer for radiation
    * `real(kind=kind_phys)`: units = K
* `water_vapor_specific_humidity_at_layer_for_radiation`: Water vapor specific humidity at layer for radiation
    * `real(kind=kind_phys)`: units = kg kg-1
* `ozone_concentration`: Ozone concentration
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_co2`: Volume mixing ratio co2
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_n2o`: Volume mixing ratio n2o
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_ch4`: Volume mixing ratio ch4
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_o2`: Volume mixing ratio o2
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_co`: Volume mixing ratio co
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc11`: Volume mixing ratio cfc11
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc12`: Volume mixing ratio cfc12
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc22`: Volume mixing ratio cfc22
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_ccl4`: Volume mixing ratio ccl4
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc113`: Volume mixing ratio cfc113
    * `real(kind=kind_phys)`: units = kg kg-1
* `instantaneous_3d_cloud_fraction`: Instantaneous 3d cloud fraction
    * `real(kind=kind_phys)`: units = frac
* `aerosol_optical_depth_for_shortwave_bands_01_16`: Aerosol optical depth for shortwave bands 01 16
    * `real(kind=kind_phys)`: units = none
* `aerosol_single_scattering_albedo_for_shortwave_bands_01_16`: Aerosol single scattering albedo for shortwave bands 01 16
    * `real(kind=kind_phys)`: units = frac
* `aerosol_asymmetry_parameter_for_shortwave_bands_01_16`: Aerosol asymmetry parameter for shortwave bands 01 16
    * `real(kind=kind_phys)`: units = none
* `aerosol_optical_depth_for_longwave_bands_01_16`: Aerosol optical depth for longwave bands 01 16
    * `real(kind=kind_phys)`: units = none
* `aerosol_single_scattering_albedo_for_longwave_bands_01_16`: Aerosol single scattering albedo for longwave bands 01 16
    * `real(kind=kind_phys)`: units = frac
* `aerosol_asymmetry_parameter_for_longwave_bands_01_16`: Aerosol asymmetry parameter for longwave bands 01 16
    * `real(kind=kind_phys)`: units = none
* `cloud_overlap_decorrelation_parameter`: Cloud overlap decorrelation parameter
    * `real(kind=kind_phys)`: units = frac
* `vertically_integrated_x_momentum_flux_due_to_form_drag`: Vertically integrated x momentum flux due to form drag
    * `real(kind=kind_phys)`: units = Pa
* `vertically_integrated_x_momentum_flux_due_to_blocking_drag`: Vertically integrated x momentum flux due to blocking drag
    * `real(kind=kind_phys)`: units = Pa
* `tendency_of_x_momentum_due_to_blocking_drag`: Tendency of x momentum due to blocking drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_x_momentum_due_to_form_drag`: Tendency of x momentum due to form drag
    * `real(kind=kind_phys)`: units = m s-2
* `time_integral_of_height_of_mountain_blocking`: Time integral of height of mountain blocking
    * `real(kind=kind_phys)`: units = m
* `time_integral_of_height_of_low_level_wave_breaking`: Time integral of height of low level wave breaking
    * `real(kind=kind_phys)`: units = m
* `time_integral_of_height_of_launch_level_of_orographic_gravity_wave`: Time integral of height of launch level of orographic gravity wave
    * `real(kind=kind_phys)`: units = m
* `time_integral_of_momentum_flux_due_to_turbulent_orographic_form_drag`: Time integral of momentum flux due to turbulent orographic form drag
    * `real(kind=kind_phys)`: units = Pa
* `time_integral_of_momentum_flux_due_to_mountain_blocking_drag`: Time integral of momentum flux due to mountain blocking drag
    * `real(kind=kind_phys)`: units = Pa
* `time_integral_of_momentum_flux_due_to_orographic_gravity_wave_drag`: Time integral of momentum flux due to orographic gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `time_integral_of_momentum_flux_due_to_nonstationary_gravity_wave`: Time integral of momentum flux due to nonstationary gravity wave
    * `real(kind=kind_phys)`: units = Pa
* `time_integral_of_change_in_x_wind_due_to_nonstationary_gravity_wave`: Time integral of change in x wind due to nonstationary gravity wave
    * `real(kind=kind_phys)`: units = m s-2
* `time_integral_of_change_in_y_wind_due_to_nonstationary_gravity_wave`: Time integral of change in y wind due to nonstationary gravity wave
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step_and_radiation_levels`: Tendency of air temperature due to longwave heating on radiation time step and radiation levels
    * `real(kind=kind_phys)`: units = K s-1
* `atmosphere_optical_thickness_due_to_cloud_at_10mu_band`: Atmosphere optical thickness due to cloud at 10mu band
    * `real(kind=kind_phys)`: units = none
* `tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step_and_radiation_levels`: Tendency of air temperature due to longwave heating assuming clear sky on radiation time step and radiation levels
    * `real(kind=kind_phys)`: units = K s-1
* `cloud_liquid_water_mixing_ratio_of_new_state`: Cloud liquid water mixing ratio of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `air_temperature_two_timesteps_back`: Air temperature two timesteps back
    * `real(kind=kind_phys)`: units = K
* `water_vapor_specific_humidity_two_timesteps_back`: Water vapor specific humidity two timesteps back
    * `real(kind=kind_phys)`: units = kg kg-1
* `surface_air_pressure_two_timesteps_back`: Surface air pressure two timesteps back
    * `real(kind=kind_phys)`: units = Pa
* `water_vapor_specific_humidity_at_previous_timestep`: Water vapor specific humidity at previous timestep
    * `real(kind=kind_phys)`: units = kg kg-1
* `surface_air_pressure_at_previous_timestep`: Surface air pressure at previous timestep
    * `real(kind=kind_phys)`: units = Pa
* `critical_relative_humidity`: Critical relative humidity
    * `real(kind=kind_phys)`: units = frac
* `tendendy_of_specific_humidity_due_to_nonphysics`: Tendendy of specific humidity due to nonphysics
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `tendency_of_air_temperature_due_to_nonphysics`: Tendency of air temperature due to nonphysics
    * `real(kind=kind_phys)`: units = K s-1
* `statistical_measures_of_subgrid_orography_collection_array`: Statistical measures of subgrid orography collection array
    * `real(kind=kind_phys)`: units = various
* `standard_deviation_of_subgrid_orography_small_scale`: Standard deviation of subgrid orography small scale
    * `real(kind=kind_phys)`: units = m
* `convexity_of_subgrid_orography_small_scale`: Convexity of subgrid orography small scale
    * `real(kind=kind_phys)`: units = none
* `asymmetry_of_subgrid_orography_small_scale`: Asymmetry of subgrid orography small scale
    * `real(kind=kind_phys)`: units = none
* `fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height_small_scale`: Fraction of grid box with subgrid orography higher than critical height small scale
    * `real(kind=kind_phys)`: units = frac
* `time_integral_of_x_stress_due_to_gravity_wave_drag`: Time integral of x stress due to gravity wave drag
    * `real(kind=kind_phys)`: units = Pa s
* `time_integral_of_y_stress_due_to_gravity_wave_drag`: Time integral of y stress due to gravity wave drag
    * `real(kind=kind_phys)`: units = Pa s
* `cloud_work_function`: Cloud work function
    * `real(kind=kind_phys)`: units = m2 s-2
* `cumulative_change_in_x_wind_due_to_shallow_convection`: Cumulative change in x wind due to shallow convection
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_y_wind_due_to_shallow_convection`: Cumulative change in y wind due to shallow convection
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_temperature_due_to_shallow_convection`: Cumulative change in temperature due to shallow convection
    * `real(kind=kind_phys)`: units = K
* `cumulative_change_in_water_vapor_specific_humidity_due_to_shallow_convection`: Cumulative change in water vapor specific humidity due to shallow convection
    * `real(kind=kind_phys)`: units = kg kg-1
* `cumulative_change_in_x_wind_due_to_deep_convection`: Cumulative change in x wind due to deep convection
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_y_wind_due_to_deep_convection`: Cumulative change in y wind due to deep convection
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_temperature_due_to_deep_convection`: Cumulative change in temperature due to deep convection
    * `real(kind=kind_phys)`: units = K
* `cumulative_change_in_water_vapor_specific_humidity_due_to_deep_convection`: Cumulative change in water vapor specific humidity due to deep convection
    * `real(kind=kind_phys)`: units = kg kg-1
* `convective_cloud_condesate_after_rainout`: Convective cloud condesate after rainout
    * `real(kind=kind_phys)`: units = kg kg-1
* `air_temperature_at_surface_adjacent_layer_on_radiation_timestep`: Air temperature at surface adjacent layer on radiation timestep
    * `real(kind=kind_phys)`: units = K
* `surface_downwelling_longwave_flux_on_radiation_timestep`: Surface downwelling longwave flux on radiation timestep
    * `real(kind=kind_phys)`: units = W m-2
* `tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_timestep`: Tendency of air temperature due to longwave heating assuming clear sky on radiation timestep
    * `real(kind=kind_phys)`: units = K s-1
* `triple_point_temperature_of_water`: Triple point temperature of water
    * `real(kind=kind_phys)`: units = K
* `autoconverion_to_snow_size_threshold`: Autoconverion to snow size threshold
    * `real(kind=kind_phys)`: units = um
* `relative_variance_of_subgrid_cloud_condensate_distribution`: Relative variance of subgrid cloud condensate distribution
    * `real(kind=kind_phys)`: units = 
* `timescale_for_autoconversion_to_snow`: Timescale for autoconversion to snow
    * `real(kind=kind_phys)`: units = s
* `relative_humidity_threshold_for_ice_nucleation`: Relative humidity threshold for ice nucleation
    * `real(kind=kind_phys)`: units = none
* `bergeron_findeisen_process_efficiency_factor`: Bergeron findeisen process efficiency factor
    * `real(kind=kind_phys)`: units = frac
* `prescribed_cloud_droplet_number_concentration`: Prescribed cloud droplet number concentration
    * `real(kind=kind_phys)`: units = m-3
* `prescribed_cloud_ice_number_concentration`: Prescribed cloud ice number concentration
    * `real(kind=kind_phys)`: units = m-3
* `prescribed_graupel_number_concentration`: Prescribed graupel number concentration
    * `real(kind=kind_phys)`: units = m-3
* `cloud_ice_mixing_ratio_of_new_state`: Cloud ice mixing ratio of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `lwe_thickness_of_explicit_precipitation_amount`: Lwe thickness of explicit precipitation amount
    * `real(kind=kind_phys)`: units = m
* `ratio_of_snowfall_to_rainfall`: Ratio of snowfall to rainfall
    * `real(kind=kind_phys)`: units = frac
* `mass_number_concentration_of_cloud_droplets_of_new_state`: Mass number concentration of cloud droplets of new state
    * `real(kind=kind_phys)`: units = kg-1
* `mass_number_concentration_of_cloud_ice_of_new_state`: Mass number concentration of cloud ice of new state
    * `real(kind=kind_phys)`: units = kg-1
* `local_rain_water_mixing_ratio`: Local rain water mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `local_snow_water_mixing_ratio`: Local snow water mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `local_graupel_mixing_ratio`: Local graupel mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `local_mass_number_concentration_of_rain_in_air`: Local mass number concentration of rain in air
    * `real(kind=kind_phys)`: units = kg-1
* `local_mass_number_concentration_of_snow_in_air`: Local mass number concentration of snow in air
    * `real(kind=kind_phys)`: units = kg-1
* `mass_number_concentration_of_ice_crystals_in_air`: Mass number concentration of ice crystals in air
    * `real(kind=kind_phys)`: units = kg-1
* `effective_radius_of_stratiform_cloud_graupel_particle`: Effective radius of stratiform cloud graupel particle
    * `real(kind=kind_phys)`: units = um
* `alpha_tuning_coefficient_for_morrison_gettelman_microphysics_scheme`: Alpha tuning coefficient for morrison gettelman microphysics scheme
    * `real(kind=kind_phys)`: units = none
* `minimum_cloud_condensate_mixing_ratio_thresholds`: Minimum cloud condensate mixing ratio thresholds
    * `real(kind=kind_phys)`: units = kg kg-1
* `subgrid_scale_cloud_fraction_from_shoc`: Subgrid scale cloud fraction from shoc
    * `real(kind=kind_phys)`: units = frac
* `minimum_value_of_specific_humidity`: Minimum value of specific humidity
    * `real(kind=kind_phys)`: units = kg kg-1
* `layer_thickness`: Layer thickness
    * `real(kind=kind_phys)`: units = m
* `cloud_overlap_param`: Cloud overlap param
    * `real(kind=kind_phys)`: units = km
* `atmosphere_heat_diffusivity_from_shoc`: Atmosphere heat diffusivity from shoc
    * `real(kind=kind_phys)`: units = m2 s-1
* `prandtl_number`: Prandtl number
    * `real(kind=kind_phys)`: units = none
* `max_atmosphere_heat_diffusivity_due_to_background`: Max atmosphere heat diffusivity due to background
    * `real(kind=kind_phys)`: units = m2 s-1
* `tendency_of_y_wind_due_to_mesoscale_orographic_gravity_wave_drag`: Tendency of y wind due to mesoscale orographic gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_momentum_due_to_blocking_drag`: Tendency of y momentum due to blocking drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_x_momentum_due_to_small_scale_gravity_wave_drag`: Tendency of x momentum due to small scale gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_momentum_due_to_small_scale_gravity_wave_drag`: Tendency of y momentum due to small scale gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_momentum_due_to_form_drag`: Tendency of y momentum due to form drag
    * `real(kind=kind_phys)`: units = m s-2
* `vertically_integrated_x_momentum_flux_due_to_mesoscale_orographic_gravity_wave_drag`: Vertically integrated x momentum flux due to mesoscale orographic gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `vertically_integrated_y_momentum_flux_due_to_mesoscale_orographic_gravity_wave_drag`: Vertically integrated y momentum flux due to mesoscale orographic gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `vertically_integrated_y_momentum_flux_due_to_blocking_drag`: Vertically integrated y momentum flux due to blocking drag
    * `real(kind=kind_phys)`: units = Pa
* `vertically_integrated_x_momentum_flux_due_to_small_scale_gravity_wave_drag`: Vertically integrated x momentum flux due to small scale gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `vertically_integrated_y_momentum_flux_due_to_small_scale_gravity_wave_drag`: Vertically integrated y momentum flux due to small scale gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `vertically_integrated_y_momentum_flux_due_to_form_drag`: Vertically integrated y momentum flux due to form drag
    * `real(kind=kind_phys)`: units = Pa
* `cloud_phase_transition_threshold_temperature`: Cloud phase transition threshold temperature
    * `real(kind=kind_phys)`: units = K
* `reciprocal_of_cloud_phase_transition_temperature_range`: Reciprocal of cloud phase transition temperature range
    * `real(kind=kind_phys)`: units = K-1
* `pressure_threshold_for_increased_tke_dissipation`: Pressure threshold for increased tke dissipation
    * `real(kind=kind_phys)`: units = Pa
* `multiplicative_tunable_parameter_for_tke_dissipation`: Multiplicative tunable parameter for tke dissipation
    * `real(kind=kind_phys)`: units = none
* `multiplicative_tunable_parameter_for_tke_dissipation_at_surface_adjacent_layer`: Multiplicative tunable parameter for tke dissipation at surface adjacent layer
    * `real(kind=kind_phys)`: units = none
* `uncentering_coefficient_for_implicit_tke_integration`: Uncentering coefficient for implicit tke integration
    * `real(kind=kind_phys)`: units = none
* `control_for_tke_dissipation_method`: Control for tke dissipation method
    * `real(kind=kind_phys)`: units = none
* `tracer_concentration_of_new_state`: Tracer concentration of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `turbulent_kinetic_energy_convective_transport_tracer`: Turbulent kinetic energy convective transport tracer
    * `real(kind=kind_phys)`: units = m2 s-2
* `kinematic_buoyancy_flux`: Kinematic buoyancy flux
    * `real(kind=kind_phys)`: units = K m s-1
* `entrainment_rate_coefficient_for_deep_convection`: Entrainment rate coefficient for deep convection
    * `real(kind=kind_phys)`: units = none
* `rain_conversion_parameter_for_deep_convection`: Rain conversion parameter for deep convection
    * `real(kind=kind_phys)`: units = m-1
* `detrainment_conversion_parameter_for_deep_convection`: Detrainment conversion parameter for deep convection
    * `real(kind=kind_phys)`: units = m-1
* `downdraft_fraction_reaching_surface_over_land_for_deep_convection`: Downdraft fraction reaching surface over land for deep convection
    * `real(kind=kind_phys)`: units = frac
* `downdraft_fraction_reaching_surface_over_water_deep_convection`: Downdraft fraction reaching surface over water deep convection
    * `real(kind=kind_phys)`: units = frac
* `rain_evaporation_coefficient_over_ocean_for_deep_convection`: Rain evaporation coefficient over ocean for deep convection
    * `real(kind=kind_phys)`: units = frac
* `rain_evaporation_coefficient_over_land_for_deep_convection`: Rain evaporation coefficient over land for deep convection
    * `real(kind=kind_phys)`: units = frac
* `momentum_transport_reduction_factor_due_to_pressure_gradient_force_for_deep_convection`: Momentum transport reduction factor due to pressure gradient force for deep convection
    * `real(kind=kind_phys)`: units = frac
* `surface_snow_melt`: Surface snow melt
    * `real(kind=kind_phys)`: units = m
* `atmosphere_momentum_diffusivity`: Atmosphere momentum diffusivity
    * `real(kind=kind_phys)`: units = m2 s-1
* `multiplicative_tuning_parameter_for_tke_dissipative_heating`: Multiplicative tuning parameter for tke dissipative heating
    * `real(kind=kind_phys)`: units = none
* `updraft_area_fraction_in_scale_aware_tke_moist_edmf_pbl_scheme`: Updraft area fraction in scale aware tke moist edmf pbl scheme
    * `real(kind=kind_phys)`: units = none
* `downdraft_area_fraction_in_scale_aware_tke_moist_edmf_pbl_scheme`: Downdraft area fraction in scale aware tke moist edmf pbl scheme
    * `real(kind=kind_phys)`: units = none
* `characteristic_grid_length_scale`: Characteristic grid length scale
    * `real(kind=kind_phys)`: units = m
* `cloud_area_fraction`: Cloud area fraction
    * `real(kind=kind_phys)`: units = frac
* `maximum_column_heating_rate`: Maximum column heating rate
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_x_wind_due_to_convective_gravity_wave_drag`: Tendency of x wind due to convective gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_wind_due_to_convective_gravity_wave_drag`: Tendency of y wind due to convective gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `aerosol_optical_properties_for_longwave_bands_01_16`: Aerosol optical properties for longwave bands 01 16
    * `real(kind=kind_phys)`: units = various
* `aerosol_optical_properties_for_shortwave_bands_01_16`: Aerosol optical properties for shortwave bands 01 16
    * `real(kind=kind_phys)`: units = various
* `air_pressure_at_lowest_model_interface`: Air pressure at lowest model interface
    * `real(kind=kind_phys)`: units = Pa
* `tracer_concentration_save`: Tracer concentration save
    * `real(kind=kind_phys)`: units = kg kg-1
* `dynamics_to_physics_timestep_ratio`: Dynamics to physics timestep ratio
    * `real(kind=kind_phys)`: units = none
* `lwe_thickness_of_ice_amount_on_dynamics_timestep`: Lwe thickness of ice amount on dynamics timestep
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_snow_amount_on_dynamics_timestep`: Lwe thickness of snow amount on dynamics timestep
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_graupel_amount_on_dynamics_timestep`: Lwe thickness of graupel amount on dynamics timestep
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_explicit_rain_amount`: Lwe thickness of explicit rain amount
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_ice_amount`: Lwe thickness of ice amount
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_snow_amount`: Lwe thickness of snow amount
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_graupel_amount`: Lwe thickness of graupel amount
    * `real(kind=kind_phys)`: units = m
* `dominant_rain_type`: Dominant rain type
    * `real(kind=kind_phys)`: units = none
* `dominant_freezing_rain_type`: Dominant freezing rain type
    * `real(kind=kind_phys)`: units = none
* `dominant_sleet_type`: Dominant sleet type
    * `real(kind=kind_phys)`: units = none
* `dominant_snow_type`: Dominant snow type
    * `real(kind=kind_phys)`: units = none
* `cumulative_lwe_thickness_of_convective_precipitation_amount`: Cumulative lwe thickness of convective precipitation amount
    * `real(kind=kind_phys)`: units = m
* `accumulated_lwe_thickness_of_precipitation_amount`: Accumulated lwe thickness of precipitation amount
    * `real(kind=kind_phys)`: units = m
* `accumulated_lwe_thickness_of_ice_amount`: Accumulated lwe thickness of ice amount
    * `real(kind=kind_phys)`: units = kg m-2
* `accumulated_lwe_thickness_of_snow_amount`: Accumulated lwe thickness of snow amount
    * `real(kind=kind_phys)`: units = kg m-2
* `accumulated_lwe_thickness_of_graupel_amount`: Accumulated lwe thickness of graupel amount
    * `real(kind=kind_phys)`: units = kg m-2
* `cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket`: Cumulative lwe thickness of convective precipitation amount in bucket
    * `real(kind=kind_phys)`: units = m
* `accumulated_lwe_thickness_of_precipitation_amount_in_bucket`: Accumulated lwe thickness of precipitation amount in bucket
    * `real(kind=kind_phys)`: units = m
* `accumulated_lwe_thickness_of_ice_amount_in_bucket`: Accumulated lwe thickness of ice amount in bucket
    * `real(kind=kind_phys)`: units = kg m-2
* `accumulated_lwe_thickness_of_snow_amount_in_bucket`: Accumulated lwe thickness of snow amount in bucket
    * `real(kind=kind_phys)`: units = kg m-2
* `accumulated_lwe_thickness_of_graupel_amount_in_bucket`: Accumulated lwe thickness of graupel amount in bucket
    * `real(kind=kind_phys)`: units = kg m-2
* `cumulative_change_in_temperature_due_to_microphysics`: Cumulative change in temperature due to microphysics
    * `real(kind=kind_phys)`: units = K
* `cumulative_change_in_water_vapor_specific_humidity_due_to_microphysics`: Cumulative change in water vapor specific humidity due to microphysics
    * `real(kind=kind_phys)`: units = kg kg-1
* `cumulative_lwe_thickness_of_precipitation_amount_for_coupling`: Cumulative lwe thickness of precipitation amount for coupling
    * `real(kind=kind_phys)`: units = m
* `cumulative_lwe_thickness_of_convective_precipitation_amount_for_coupling`: Cumulative lwe thickness of convective precipitation amount for coupling
    * `real(kind=kind_phys)`: units = m
* `cumulative_lwe_thickness_of_snow_amount_for_coupling`: Cumulative lwe thickness of snow amount for coupling
    * `real(kind=kind_phys)`: units = m
* `column_precipitable_water`: Column precipitable water
    * `real(kind=kind_phys)`: units = kg m-2
* `lwe_thickness_of_rain_amount_on_dynamics_timestep_for_coupling`: Lwe thickness of rain amount on dynamics timestep for coupling
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_snowfall_amount_on_dynamics_timestep_for_coupling`: Lwe thickness of snowfall amount on dynamics timestep for coupling
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_convective_precipitation_amount_on_previous_timestep`: Lwe thickness of convective precipitation amount on previous timestep
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_explicit_precipitation_amount_on_previous_timestep`: Lwe thickness of explicit precipitation amount on previous timestep
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_ice_precipitation_amount_on_previous_timestep`: Lwe thickness of ice precipitation amount on previous timestep
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_snowfall_amount_on_previous_timestep`: Lwe thickness of snowfall amount on previous timestep
    * `real(kind=kind_phys)`: units = m
* `lwe_thickness_of_graupel_amount_on_previous_timestep`: Lwe thickness of graupel amount on previous timestep
    * `real(kind=kind_phys)`: units = m
* `convective_precipitation_rate_on_previous_timestep`: Convective precipitation rate on previous timestep
    * `real(kind=kind_phys)`: units = mm s-1
* `explicit_precipitation_rate_on_previous_timestep`: Explicit precipitation rate on previous timestep
    * `real(kind=kind_phys)`: units = mm s-1
* `ice_precipitation_rate_on_previous_timestep`: Ice precipitation rate on previous timestep
    * `real(kind=kind_phys)`: units = mm s-1
* `snowfall_rate_on_previous_timestep`: Snowfall rate on previous timestep
    * `real(kind=kind_phys)`: units = mm s-1
* `graupel_precipitation_rate_on_previous_timestep`: Graupel precipitation rate on previous timestep
    * `real(kind=kind_phys)`: units = mm s-1
* `precip_overlap_param`: Precip overlap param
    * `real(kind=kind_phys)`: units = km
* `maximum_soil_moisture_content_for_land_surface_model`: Maximum soil moisture content for land surface model
    * `real(kind=kind_phys)`: units = m
* `minimum_soil_moisture_content_for_land_surface_model`: Minimum soil moisture content for land surface model
    * `real(kind=kind_phys)`: units = m
* `perturbation_of_soil_type_b_parameter`: Perturbation of soil type b parameter
    * `real(kind=kind_phys)`: units = frac
* `perturbation_of_leaf_area_index`: Perturbation of leaf area index
    * `real(kind=kind_phys)`: units = frac
* `perturbation_of_vegetation_fraction`: Perturbation of vegetation fraction
    * `real(kind=kind_phys)`: units = frac
* `magnitude_of_perturbation_of_vegetation_fraction`: Magnitude of perturbation of vegetation fraction
    * `real(kind=kind_phys)`: units = frac
* `subsurface_runoff_flux`: Subsurface runoff flux
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `surface_runoff_flux`: Surface runoff flux
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `surface_snow_area_fraction`: Surface snow area fraction
    * `real(kind=kind_phys)`: units = frac
* `soil_moisture_content`: Soil moisture content
    * `real(kind=kind_phys)`: units = kg m-2
* `snow_freezing_rain_upward_latent_heat_flux`: Snow freezing rain upward latent heat flux
    * `real(kind=kind_phys)`: units = W m-2
* `normalized_soil_wetness`: Normalized soil wetness
    * `real(kind=kind_phys)`: units = frac
* `x_wind_save`: X wind save
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind_save`: Y wind save
    * `real(kind=kind_phys)`: units = m s-1
* `tendency_of_specific_humidity_due_to_moist_convection_for_coupling`: Tendency of specific humidity due to moist convection for coupling
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `cumulative_cloud_work_function`: Cumulative cloud work function
    * `real(kind=kind_phys)`: units = m2 s-1
* `cumulative_atmosphere_updraft_convective_mass_flux`: Cumulative atmosphere updraft convective mass flux
    * `real(kind=kind_phys)`: units = kg m-1 s-2
* `cumulative_atmosphere_downdraft_convective_mass_flux`: Cumulative atmosphere downdraft convective mass flux
    * `real(kind=kind_phys)`: units = kg m-1 s-2 
* `cumulative_atmosphere_detrainment_convective_mass_flux`: Cumulative atmosphere detrainment convective mass flux
    * `real(kind=kind_phys)`: units = kg m-1 s-2
* `cellular_automata_global_patter_from_coupled_process`: Cellular automata global patter from coupled process
    * `real(kind=kind_phys)`: units = flag
* `cellular_automata_vertical_weight`: Cellular automata vertical weight
    * `real(kind=kind_phys)`: units = frac
* `skeb_x_wind_weights_from_coupled_process`: Skeb x wind weights from coupled process
    * `real(kind=kind_phys)`: units = none
* `skeb_y_wind_weights_from_coupled_process`: Skeb y wind weights from coupled process
    * `real(kind=kind_phys)`: units = none
* `shum_weights_from_coupled_process`: Shum weights from coupled process
    * `real(kind=kind_phys)`: units = none
* `weights_for_stochastic_sppt_perturbation_flipped`: Weights for stochastic sppt perturbation flipped
    * `real(kind=kind_phys)`: units = none
* `weights_for_stochastic_skeb_perturbation_of_x_wind_flipped`: Weights for stochastic skeb perturbation of x wind flipped
    * `real(kind=kind_phys)`: units = none
* `weights_for_stochastic_skeb_perturbation_of_y_wind_flipped`: Weights for stochastic skeb perturbation of y wind flipped
    * `real(kind=kind_phys)`: units = none
* `weights_for_stochastic_shum_perturbation_flipped`: Weights for stochastic shum perturbation flipped
    * `real(kind=kind_phys)`: units = none
* `dissipation_estimate_of_air_temperature_at_model_layers`: Dissipation estimate of air temperature at model layers
    * `real(kind=kind_phys)`: units = K
* `rain_mixing_ratio`: Rain mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `snow_mixing_ratio`: Snow mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `graupel_mixing_ratio`: Graupel mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `tendency_of_air_temperature_to_withold_from_sppt`: Tendency of air temperature to withold from sppt
    * `real(kind=kind_phys)`: units = K s-1
* `rain_mixing_ratio_of_new_state`: Rain mixing ratio of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `snow_mixing_ratio_of_new_state`: Snow mixing ratio of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `graupel_mixing_ratio_of_new_state`: Graupel mixing ratio of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_fraction_of_frozen_soil_moisture_for_land_surface_model`: Volume fraction of frozen soil moisture for land surface model
    * `real(kind=kind_phys)`: units = frac
* `normalized_soil_wetness_for_land_surface_model`: Normalized soil wetness for land surface model
    * `real(kind=kind_phys)`: units = frac
* `cloud_liquid_water_mixing_ratio_at_surface_adjacent_layer`: Cloud liquid water mixing ratio at surface adjacent layer
    * `real(kind=kind_phys)`: units = kg kg-1
* `flag_for_frozen_soil_physics`: Flag for frozen soil physics
    * `real(kind=kind_phys)`: units = flag
* `snow_temperature_bottom_first_layer_over_land`: Snow temperature bottom first layer over land
    * `real(kind=kind_phys)`: units = K
* `cloud_condensed_water_mixing_ratio_at_surface_over_land`: Cloud condensed water mixing ratio at surface over land
    * `real(kind=kind_phys)`: units = kg kg-1
* `water_vapor_mixing_ratio_at_surface_over_land`: Water vapor mixing ratio at surface over land
    * `real(kind=kind_phys)`: units = kg kg-1
* `total_runoff`: Total runoff
    * `real(kind=kind_phys)`: units = kg m-2
* `cumulative_surface_runoff_amount`: Cumulative surface runoff amount
    * `real(kind=kind_phys)`: units = kg m-2
* `total_accumulated_snowfall_over_land`: Total accumulated snowfall over land
    * `real(kind=kind_phys)`: units = kg m-2
* `cloud_condensed_water_mixing_ratio_at_surface_over_ice`: Cloud condensed water mixing ratio at surface over ice
    * `real(kind=kind_phys)`: units = kg kg-1
* `water_vapor_mixing_ratio_at_surface_over_ice`: Water vapor mixing ratio at surface over ice
    * `real(kind=kind_phys)`: units = kg kg-1
* `snow_temperature_bottom_first_layer_over_ice`: Snow temperature bottom first layer over ice
    * `real(kind=kind_phys)`: units = K
* `total_accumulated_snowfall_over_ice`: Total accumulated snowfall over ice
    * `real(kind=kind_phys)`: units = kg m-2
* `frozen_precipitation_density`: Frozen precipitation density
    * `real(kind=kind_phys)`: units = kg m-3
* `mass_number_concentration_of_rain_of_new_state`: Mass number concentration of rain of new state
    * `real(kind=kind_phys)`: units = kg-1
* `mass_number_concentration_of_snow_of_new_state`: Mass number concentration of snow of new state
    * `real(kind=kind_phys)`: units = kg-1
* `mass_number_concentration_of_graupel_of_new_state`: Mass number concentration of graupel of new state
    * `real(kind=kind_phys)`: units = kg-1
* `tendency_of_rain_water_mixing_ratio_due_to_microphysics`: Tendency of rain water mixing ratio due to microphysics
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `autoconversion_to_snow_coefficient`: Autoconversion to snow coefficient
    * `real(kind=kind_phys)`: units = none
* `autoconversion_to_rain_coefficient`: Autoconversion to rain coefficient
    * `real(kind=kind_phys)`: units = none
* `precipitation_evaporation_coefficient`: Precipitation evaporation coefficient
    * `real(kind=kind_phys)`: units = none
* `cloud_condensate_autoconversion_threshold_coefficient`: Cloud condensate autoconversion threshold coefficient
    * `real(kind=kind_phys)`: units = none
* `tendency_of_tracers_due_to_model_physics`: Tendency of tracers due to model physics
    * `real(kind=kind_phys)`: units = kg kg-1 s-1
* `air_temperature_at_lowest_model_layer_for_diag`: Air temperature at lowest model layer for diag
    * `real(kind=kind_phys)`: units = K
* `water_vapor_specific_humidity_at_lowest_model_layer_for_diag`: Water vapor specific humidity at lowest model layer for diag
    * `real(kind=kind_phys)`: units = kg kg-1
* `instantaneous_surface_upward_sensible_heat_flux_for_chemistry_coupling`: Instantaneous surface upward sensible heat flux for chemistry coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_atmosphere_heat_diffusivity`: Instantaneous atmosphere heat diffusivity
    * `real(kind=kind_phys)`: units = m2 s-1
* `surface_upward_sensible_heat_flux_reduction_factor`: Surface upward sensible heat flux reduction factor
    * `real(kind=kind_phys)`: units = none
* `surface_upward_latent_heat_flux_reduction_factor`: Surface upward latent heat flux reduction factor
    * `real(kind=kind_phys)`: units = none
* `cloud_area_fraction_in_atmosphere_layer_of_new_state`: Cloud area fraction in atmosphere layer of new state
    * `real(kind=kind_phys)`: units = frac
* `aerosol_aware_multiplicative_rain_conversion_parameter_for_deep_convection`: Aerosol aware multiplicative rain conversion parameter for deep convection
    * `real(kind=kind_phys)`: units = none
* `cellular_automata_vertical_velocity_perturbation_threshold_for_deep_convection`: Cellular automata vertical velocity perturbation threshold for deep convection
    * `real(kind=kind_phys)`: units = m s-1
* `fraction_of_cellular_automata_for_deep_convection`: Fraction of cellular automata for deep convection
    * `real(kind=kind_phys)`: units = frac
* `physics_field_for_coupling`: Physics field for coupling
    * `real(kind=kind_phys)`: units = m2 s-2
* `surface_emissivity_in_each_RRTMGP_LW_band`: Surface emissivity in each RRTMGP LW band
    * `real(kind=kind_phys)`: units = none
* `temperature_at_2m_from_noahmp`: Temperature at 2m from noahmp
    * `real(kind=kind_phys)`: units = K
* `specific_humidity_at_2m_from_noahmp`: Specific humidity at 2m from noahmp
    * `real(kind=kind_phys)`: units = kg kg-1
* `minimum_temperature_at_2m`: Minimum temperature at 2m
    * `real(kind=kind_phys)`: units = K
* `maximum_temperature_at_2m`: Maximum temperature at 2m
    * `real(kind=kind_phys)`: units = K
* `minimum_specific_humidity_at_2m`: Minimum specific humidity at 2m
    * `real(kind=kind_phys)`: units = kg kg-1
* `maximum_specific_humidity_at_2m`: Maximum specific humidity at 2m
    * `real(kind=kind_phys)`: units = kg kg-1
* `maximum_wind_at_10m`: Maximum wind at 10m
    * `real(kind=kind_phys)`: units = m s-1
* `maximum_x_wind_at_10m`: Maximum x wind at 10m
    * `real(kind=kind_phys)`: units = m s-1
* `maximum_y_wind_at_10m`: Maximum y wind at 10m
    * `real(kind=kind_phys)`: units = m s-1
* `dewpoint_temperature_at_2m`: Dewpoint temperature at 2m
    * `real(kind=kind_phys)`: units = K
* `natural_log_of_water_vapor_forcing_data_at_pressure_levels`: Natural log of water vapor forcing data at pressure levels
    * `real(kind=kind_phys)`: units = log(Pa)
* `radius_of_earth`: Radius of earth
    * `real(kind=kind_phys)`: units = m
* `hours_between_clearing_of_diagnostic_buckets`: Hours between clearing of diagnostic buckets
    * `real(kind=kind_phys)`: units = h
* `tendency_of_x_wind_due_to_nonorographic_gravity_wave_drag`: Tendency of x wind due to nonorographic gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_wind_due_to_nonorographic_gravity_wave_drag`: Tendency of y wind due to nonorographic gravity wave drag
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_air_temperature_due_to_nonorographic_gravity_wave_drag`: Tendency of air temperature due to nonorographic gravity wave drag
    * `real(kind=kind_phys)`: units = K s-1
* `atmosphere_momentum_diffusivity_due_to_nonorographic_gravity_wave_drag`: Atmosphere momentum diffusivity due to nonorographic gravity wave drag
    * `real(kind=kind_phys)`: units = m2 s-1
* `momentum_flux_due_to_subgrid_scale_orographic_gravity_wave_drag`: Momentum flux due to subgrid scale orographic gravity wave drag
    * `real(kind=kind_phys)`: units = Pa
* `height_of_launch_level_of_nonorographic_gravity_waves`: Height of launch level of nonorographic gravity waves
    * `real(kind=kind_phys)`: units = m
* `cumulative_change_in_x_wind_due_to_rayleigh_damping`: Cumulative change in x wind due to rayleigh damping
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_y_wind_due_to_rayleigh_damping`: Cumulative change in y wind due to rayleigh damping
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_temperature_due_to_rayleigh_damping`: Cumulative change in temperature due to rayleigh damping
    * `real(kind=kind_phys)`: units = K
* `RRTMGP_lw_flux_profile_upward_allsky`: RRTMGP lw flux profile upward allsky
    * `real(kind=kind_phys)`: units = W m-2
* `RRTMGP_longwave_flux_downward_allsky`: RRTMGP longwave flux downward allsky
    * `real(kind=kind_phys)`: units = W m-2
* `RRTMGP_lw_flux_profile_upward_clrsky`: RRTMGP lw flux profile upward clrsky
    * `real(kind=kind_phys)`: units = W m-2
* `RRTMGP_lw_flux_profile_downward_clrsky`: RRTMGP lw flux profile downward clrsky
    * `real(kind=kind_phys)`: units = W m-2
* `RRTMGP_jacobian_of_lw_flux_upward`: RRTMGP jacobian of lw flux upward
    * `real(kind=kind_phys)`: units = W m-2 K-1
* `min_grid_scale`: Min grid scale
    * `real(kind=kind_phys)`: units = m2 rad-2
* `reciprocal_of_grid_scale_range`: Reciprocal of grid scale range
    * `real(kind=kind_phys)`: units = rad2 m-2
* `surface_air_pressure_diag`: Surface air pressure diag
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_cosine_of_zenith_angle`: Instantaneous cosine of zenith angle
    * `real(kind=kind_phys)`: units = none
* `surface_upwelling_longwave_flux_from_coupled_process`: Surface upwelling longwave flux from coupled process
    * `real(kind=kind_phys)`: units = W m-2
* `tendency_of_air_temperature_due_to_integrated_dynamics_through_earths_atmosphere`: Tendency of air temperature due to integrated dynamics through earths atmosphere
    * `real(kind=kind_phys)`: units = K s-1
* `tunable_parameter_for_critical_cloud_top_entrainment_instability_criteria`: Tunable parameter for critical cloud top entrainment instability criteria
    * `real(kind=kind_phys)`: units = none
* `duration_of_sunshine`: Duration of sunshine
    * `real(kind=kind_phys)`: units = s
* `updated_tendency_of_air_temperature_due_to_longwave_heating_on_physics_time_step`: Updated tendency of air temperature due to longwave heating on physics time step
    * `real(kind=kind_phys)`: units = K s-1
* `surface_upwelling_longwave_flux`: Surface upwelling longwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_longwave_flux_over_land`: Surface upwelling longwave flux over land
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_longwave_flux_over_ice`: Surface upwelling longwave flux over ice
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_longwave_flux_over_water_interstitial`: Surface upwelling longwave flux over water interstitial
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_surface_downwelling_longwave_flux_multiplied_by_timestep`: Cumulative surface downwelling longwave flux multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_upwelling_longwave_flux_multiplied_by_timestep`: Cumulative surface upwelling longwave flux multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_pressure_multiplied_by_timestep`: Cumulative surface pressure multiplied by timestep
    * `real(kind=kind_phys)`: units = Pa s
* `cumulative_change_in_temperature_due_to_longwave_radiation`: Cumulative change in temperature due to longwave radiation
    * `real(kind=kind_phys)`: units = K
* `cumulative_change_in_temperature_due_to_shortwave_radiation`: Cumulative change in temperature due to shortwave radiation
    * `real(kind=kind_phys)`: units = K
* `grid_sensitive_critical_cloud_top_entrainment_instability_criteria`: Grid sensitive critical cloud top entrainment instability criteria
    * `real(kind=kind_phys)`: units = none
* `cloud_top_entrainment_instability_value`: Cloud top entrainment instability value
    * `real(kind=kind_phys)`: units = none
* `netcdf_float_fillvalue`: Netcdf float fillvalue
    * `real(kind=kind_phys)`: units = none
* `critical_relative_humidity_at_surface`: Critical relative humidity at surface
    * `real(kind=kind_phys)`: units = frac
* `critical_relative_humidity_at_top_of_atmosphere_boundary_layer`: Critical relative humidity at top of atmosphere boundary layer
    * `real(kind=kind_phys)`: units = frac
* `critical_relative_humidity_at_toa`: Critical relative humidity at toa
    * `real(kind=kind_phys)`: units = frac
* `max_critical_relative_humidity`: Max critical relative humidity
    * `real(kind=kind_phys)`: units = frac
* `air_temperature_save_from_convective_parameterization`: Air temperature save from convective parameterization
    * `real(kind=kind_phys)`: units = K
* `surface_albedo_due_to_near_IR_direct`: Surface albedo due to near IR direct
    * `real(kind=kind_phys)`: units = frac
* `surface_albedo_due_to_near_IR_diffused`: Surface albedo due to near IR diffused
    * `real(kind=kind_phys)`: units = frac
* `surface_albedo_due_to_UV_and_VIS_direct`: Surface albedo due to UV and VIS direct
    * `real(kind=kind_phys)`: units = frac
* `surface_albedo_due_to_UV_and_VIS_diffused`: Surface albedo due to UV and VIS diffused
    * `real(kind=kind_phys)`: units = frac
* `tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step_and_radiation_levels`: Tendency of air temperature due to shortwave heating on radiation time step and radiation levels
    * `real(kind=kind_phys)`: units = K s-1
* `atmosphere_optical_thickness_due_to_cloud_at_0p55mu_band`: Atmosphere optical thickness due to cloud at 0.55mu band
    * `real(kind=kind_phys)`: units = none
* `tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step_and_radiation_levels`: Tendency of air temperature due to shortwave heating assuming clear sky on radiation time step and radiation levels
    * `real(kind=kind_phys)`: units = K s-1
* `air_temperature_on_previous_timestep`: Air temperature on previous timestep
    * `real(kind=kind_phys)`: units = K
* `specific_humidity_on_previous_timestep`: Specific humidity on previous timestep
    * `real(kind=kind_phys)`: units = kg kg-1
* `lwe_surface_snow_from_coupled_process`: Lwe surface snow from coupled process
    * `real(kind=kind_phys)`: units = m
* `perturbation_of_momentum_roughness_length`: Perturbation of momentum roughness length
    * `real(kind=kind_phys)`: units = frac
* `perturbation_of_heat_to_momentum_roughness_length_ratio`: Perturbation of heat to momentum roughness length ratio
    * `real(kind=kind_phys)`: units = frac
* `surface_roughness_length_from_wave_model`: Surface roughness length from wave model
    * `real(kind=kind_phys)`: units = cm
* `decorrelation_length_used_by_overlap_method`: Decorrelation length used by overlap method
    * `real(kind=kind_phys)`: units = km
* `geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature`: Geopotential difference between midlayers divided by midlayer virtual temperature
    * `real(kind=kind_phys)`: units = m2 s-2 K-1
* `joules_per_calorie_constant`: Joules per calorie constant
    * `real(kind=kind_phys)`: units = J cal-1
* `sea_water_reference_density`: Sea water reference density
    * `real(kind=kind_phys)`: units = kg m-3
* `surface_skin_temperature_for_nsst`: Surface skin temperature for nsst
    * `real(kind=kind_phys)`: units = K
* `heat_content_in_diurnal_thermocline`: Heat content in diurnal thermocline
    * `real(kind=kind_phys)`: units = K m
* `sea_water_salinity_in_diurnal_thermocline`: Sea water salinity in diurnal thermocline
    * `real(kind=kind_phys)`: units = ppt m
* `x_current_in_diurnal_thermocline`: X current in diurnal thermocline
    * `real(kind=kind_phys)`: units = m2 s-1
* `y_current_in_diurnal_thermocline`: Y current in diurnal thermocline
    * `real(kind=kind_phys)`: units = m2 s-1
* `diurnal_thermocline_layer_thickness`: Diurnal thermocline layer thickness
    * `real(kind=kind_phys)`: units = m
* `ocean_mixed_layer_thickness`: Ocean mixed layer thickness
    * `real(kind=kind_phys)`: units = m
* `derivative_of_heat_content_in_diurnal_thermocline_wrt_surface_skin_temperature`: Derivative of heat content in diurnal thermocline wrt surface skin temperature
    * `real(kind=kind_phys)`: units = m
* `derivative_of_diurnal_thermocline_layer_thickness_wrt_surface_skin_temperature`: Derivative of diurnal thermocline layer thickness wrt surface skin temperature
    * `real(kind=kind_phys)`: units = m K-1
* `molecular_sublayer_temperature_correction_in_sea_water`: Molecular sublayer temperature correction in sea water
    * `real(kind=kind_phys)`: units = K
* `molecular_sublayer_thickness_in_sea_water`: Molecular sublayer thickness in sea water
    * `real(kind=kind_phys)`: units = m
* `coefficient_c_0`: Coefficient c 0
    * `real(kind=kind_phys)`: units = none
* `coefficient_c_d`: Coefficient c d
    * `real(kind=kind_phys)`: units = none
* `coefficient_w_0`: Coefficient w 0
    * `real(kind=kind_phys)`: units = none
* `coefficient_w_d`: Coefficient w d
    * `real(kind=kind_phys)`: units = none
* `free_convection_layer_thickness_in_sea_water`: Free convection layer thickness in sea water
    * `real(kind=kind_phys)`: units = m
* `control_for_diurnal_thermocline_calculation`: Control for diurnal thermocline calculation
    * `real(kind=kind_phys)`: units = index
* `surface_sensible_heat_due_to_rainfall`: Surface sensible heat due to rainfall
    * `real(kind=kind_phys)`: units = W
* `air_temperature_lapse_rate_constant`: Air temperature lapse rate constant
    * `real(kind=kind_phys)`: units = K m-1
* `mean_change_over_depth_in_sea_water_temperature`: Mean change over depth in sea water temperature
    * `real(kind=kind_phys)`: units = K
* `surface_upwelling_longwave_flux_on_radiation_time_step`: Surface upwelling longwave flux on radiation time step
    * `real(kind=kind_phys)`: units = W m-2
* `multiplicative_tuning_parameter_for_atmosphere_diffusivity`: Multiplicative tuning parameter for atmosphere diffusivity
    * `real(kind=kind_phys)`: units = none
* `control_for_variable_bulk_richardson_number`: Control for variable bulk richardson number
    * `real(kind=kind_phys)`: units = flag
* `coefficient_for_variable_bulk_richardson_number_over_land`: Coefficient for variable bulk richardson number over land
    * `real(kind=kind_phys)`: units = none
* `coefficient_for_variable_bulk_richardson_number_over_water`: Coefficient for variable bulk richardson number over water
    * `real(kind=kind_phys)`: units = none
* `mass_number_concentration_of_rain_water_in_air`: Mass number concentration of rain water in air
    * `real(kind=kind_phys)`: units = kg-1
* `tendency_of_hygroscopic_aerosols_at_surface_adjacent_layer`: Tendency of hygroscopic aerosols at surface adjacent layer
    * `real(kind=kind_phys)`: units = kg-1 s-1
* `tendency_of_nonhygroscopic_ice_nucleating_aerosols_at_surface_adjacent_layer`: Tendency of nonhygroscopic ice nucleating aerosols at surface adjacent layer
    * `real(kind=kind_phys)`: units = kg-1 s-1
* `mass_number_concentration_of_hygroscopic_aerosols_of_new_state`: Mass number concentration of hygroscopic aerosols of new state
    * `real(kind=kind_phys)`: units = kg-1
* `mass_number_concentration_of_nonhygroscopic_ice_nucleating_aerosols_of_new_state`: Mass number concentration of nonhygroscopic ice nucleating aerosols of new state
    * `real(kind=kind_phys)`: units = kg-1
* `cumulative_change_in_x_wind_due_to_physics`: Cumulative change in x wind due to physics
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_y_wind_due_to_physics`: Cumulative change in y wind due to physics
    * `real(kind=kind_phys)`: units = m s-1
* `cumulative_change_in_temperature_due_to_physics`: Cumulative change in temperature due to physics
    * `real(kind=kind_phys)`: units = K
* `cumulative_change_in_water_vapor_specific_humidity_due_to_physics`: Cumulative change in water vapor specific humidity due to physics
    * `real(kind=kind_phys)`: units = kg kg-1
* `cumulative_change_in_ozone_concentration_due_to_physics`: Cumulative change in ozone concentration due to physics
    * `real(kind=kind_phys)`: units = kg kg-1
* `surface_upwelling_direct_near_infrared_shortwave_flux`: Surface upwelling direct near infrared shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_diffuse_near_infrared_shortwave_flux`: Surface upwelling diffuse near infrared shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux`: Surface upwelling direct ultraviolet and visible shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux`: Surface upwelling diffuse ultraviolet and visible shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_direct_near_infrared_shortwave_flux`: Surface downwelling direct near infrared shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_diffuse_near_infrared_shortwave_flux`: Surface downwelling diffuse near infrared shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux`: Surface downwelling direct ultraviolet and visible shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux`: Surface downwelling diffuse ultraviolet and visible shortwave flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_layer_scheme_enthalpy_flux_factor`: Surface layer scheme enthalpy flux factor
    * `real(kind=kind_phys)`: units = none
* `consecutive_calls_for_GF`: Consecutive calls for GF
    * `real(kind=kind_phys)`: units = none
* `weights_for_stochastic_surface_physics_perturbation_flipped`: Weights for stochastic surface physics perturbation flipped
    * `real(kind=kind_phys)`: units = none
* `area_type_from_coupled_process`: Area type from coupled process
    * `real(kind=kind_phys)`: units = flag
* `enhancement_to_wind_speed_at_surface_adjacent_layer_due_to_convection`: Enhancement to wind speed at surface adjacent layer due to convection
    * `real(kind=kind_phys)`: units = m s-1
* `instantaneous_surface_potential_evaporation`: Instantaneous surface potential evaporation
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_ground_heat_flux`: Instantaneous surface ground heat flux
    * `real(kind=kind_phys)`: units = W m-2
* `x_wind_at_lowest_model_layer_for_diag`: X wind at lowest model layer for diag
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind_at_lowest_model_layer_for_diag`: Y wind at lowest model layer for diag
    * `real(kind=kind_phys)`: units = m s-1
* `instantaneous_surface_downwelling_longwave_flux_for_coupling`: Instantaneous surface downwelling longwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_downwelling_shortwave_flux_for_coupling`: Instantaneous surface downwelling shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_surface_downwelling_longwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface downwelling longwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_downwelling_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface downwelling shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `instantaneous_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling`: Instantaneous surface downwelling direct near infrared shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling`: Instantaneous surface downwelling diffuse near infrared shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling`: Instantaneous surface downwelling direct ultraviolet and visible shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling`: Instantaneous surface downwelling diffuse ultraviolet and visible shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_surface_downwelling_direct_nir_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface downwelling direct nir shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_downwelling_diffuse_nir_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface downwelling diffuse nir shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_downwelling_direct_uv_and_vis_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface downwelling direct uv and vis shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_downwelling_diffuse_uv_and_vis_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface downwelling diffuse uv and vis shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `instantaneous_surface_net_downward_longwave_flux_for_coupling`: Instantaneous surface net downward longwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_surface_net_downwelling_longwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface net downwelling longwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `instantaneous_temperature_at_2m_for_coupling`: Instantaneous temperature at 2m for coupling
    * `real(kind=kind_phys)`: units = K
* `instantaneous_specific_humidity_at_2m_for_coupling`: Instantaneous specific humidity at 2m for coupling
    * `real(kind=kind_phys)`: units = kg kg-1
* `instantaneous_x_wind_at_10m_for_coupling`: Instantaneous x wind at 10m for coupling
    * `real(kind=kind_phys)`: units = m s-1
* `instantaneous_y_wind_at_10m_for_coupling`: Instantaneous y wind at 10m for coupling
    * `real(kind=kind_phys)`: units = m s-1
* `instantaneous_surface_skin_temperature_for_coupling`: Instantaneous surface skin temperature for coupling
    * `real(kind=kind_phys)`: units = K
* `instantaneous_surface_air_pressure_for_coupling`: Instantaneous surface air pressure for coupling
    * `real(kind=kind_phys)`: units = Pa
* `instantaneous_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling`: Instantaneous surface net downward direct near infrared shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling`: Instantaneous surface net downward diffuse near infrared shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling`: Instantaneous surface net downward direct ultraviolet and visible shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling`: Instantaneous surface net downward diffuse ultraviolet and visible shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `instantaneous_surface_net_downward_shortwave_flux_for_coupling`: Instantaneous surface net downward shortwave flux for coupling
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_surface_net_downwelling_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface net downwelling shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_net_downwelling_direct_nir_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface net downwelling direct nir shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_net_downwellling_diffuse_nir_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface net downwellling diffuse nir shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_net_downwelling_direct_uv_and_vis_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface net downwelling direct uv and vis shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_net_downwelling_diffuse_uv_and_vis_shortwave_flux_for_coupling_multiplied_by_timestep`: Cumulative surface net downwelling diffuse uv and vis shortwave flux for coupling multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_ground_heat_flux_multiplied_by_timestep`: Cumulative surface ground heat flux multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_upward_latent_heat_flux_in_soil_multiplied_by_timestep`: Cumulative upward latent heat flux in soil multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_canopy_upward_latent_heat_flu_multiplied_by_timestep`: Cumulative canopy upward latent heat flu multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_transpiration_flux_multiplied_by_timestep`: Cumulative transpiration flux multiplied by timestep
    * `real(kind=kind_phys)`: units = kg m-2
* `cumulative_snow_deposition_sublimation_upward_latent_heat_flux_multiplied_by_timestep`: Cumulative snow deposition sublimation upward latent heat flux multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_snow_area_fraction_multiplied_by_timestep`: Cumulative surface snow area fraction multiplied by timestep
    * `real(kind=kind_phys)`: units = s
* `cumulative_snow_freezing_rain_upward_latent_heat_flux_multiplied_by_timestep`: Cumulative snow freezing rain upward latent heat flux multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `cumulative_surface_upward_potential_latent_heat_flux_multiplied_by_timestep`: Cumulative surface upward potential latent heat flux multiplied by timestep
    * `real(kind=kind_phys)`: units = W m-2 s
* `multiplicative_tuning_parameter_for_reduced_surface_heat_fluxes_due_to_canopy_heat_storage`: Multiplicative tuning parameter for reduced surface heat fluxes due to canopy heat storage
    * `real(kind=kind_phys)`: units = none
* `multiplicative_tuning_parameter_for_reduced_latent_heat_flux_due_to_canopy_heat_storage`: Multiplicative tuning parameter for reduced latent heat flux due to canopy heat storage
    * `real(kind=kind_phys)`: units = none
* `tunable_parameter_for_critical_cloud_workfunction_in_relaxed_arakawa_schubert_deep_convection`: Tunable parameter for critical cloud workfunction in relaxed arakawa schubert deep convection
    * `real(kind=kind_phys)`: units = none
* `autoconversion_to_snow_coefficient_for_deep_convection`: Autoconversion to snow coefficient for deep convection
    * `real(kind=kind_phys)`: units = none
* `autoconversion_to_rain_coefficient_for_deep_convection`: Autoconversion to rain coefficient for deep convection
    * `real(kind=kind_phys)`: units = none
* `cloud_condensate_autoconversion_threshold_coefficient_for_deep_convection`: Cloud condensate autoconversion threshold coefficient for deep convection
    * `real(kind=kind_phys)`: units = none
* `cloud_condensate_detrainment_coefficient`: Cloud condensate detrainment coefficient
    * `real(kind=kind_phys)`: units = none
* `dimensionless_exner_function_at_interface`: Dimensionless exner function at interface
    * `real(kind=kind_phys)`: units = none
* `liquid_water_density`: Liquid water density
    * `real(kind=kind_phys)`: units = kg m-3
* `x_wind_of_new_state_at_surface_adjacent_layer`: X wind of new state at surface adjacent layer
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind_of_new_state_at_surface_adjacent_layer`: Y wind of new state at surface adjacent layer
    * `real(kind=kind_phys)`: units = m s-1
* `air_temperature_of_new_state_at_surface_adjacent_layer`: Air temperature of new state at surface adjacent layer
    * `real(kind=kind_phys)`: units = K
* `specific_humidity_of_new_state_at_surface_adjacent_layer`: Specific humidity of new state at surface adjacent layer
    * `real(kind=kind_phys)`: units = kg kg-1
* `ratio_of_wind_at_surface_adjacent_layer_to_wind_at_10m`: Ratio of wind at surface adjacent layer to wind at 10m
    * `real(kind=kind_phys)`: units = ratio
* `toa_incident_sw_flux_by_spectral_point`: Toa incident sw flux by spectral point
    * `real(kind=kind_phys)`: units = W m-2
* `air_density_at_lowest_model_layer`: Air density at lowest model layer
    * `real(kind=kind_phys)`: units = kg m-3
* `water_equivalent_accumulated_snow_depth_over_land_save`: Water equivalent accumulated snow depth over land save
    * `real(kind=kind_phys)`: units = mm
* `surface_snow_thickness_water_equivalent_over_land_save`: Surface snow thickness water equivalent over land save
    * `real(kind=kind_phys)`: units = mm
* `surface_skin_temperature_over_land_interstitial_save`: Surface skin temperature over land interstitial save
    * `real(kind=kind_phys)`: units = K
* `canopy_water_amount_save`: Canopy water amount save
    * `real(kind=kind_phys)`: units = kg m-2
* `volume_fraction_of_soil_moisture_save`: Volume fraction of soil moisture save
    * `real(kind=kind_phys)`: units = frac
* `soil_temperature_save`: Soil temperature save
    * `real(kind=kind_phys)`: units = K
* `volume_fraction_of_unfrozen_soil_moisture_save`: Volume fraction of unfrozen soil moisture save
    * `real(kind=kind_phys)`: units = frac
* `fraction_of_ice_water_cloud`: Fraction of ice water cloud
    * `real(kind=kind_phys)`: units = frac
* `fraction_of_rain_water_cloud`: Fraction of rain water cloud
    * `real(kind=kind_phys)`: units = frac
* `rime_factor`: Rime factor
    * `real(kind=kind_phys)`: units = frac
* `total_cloud_condensate_mixing_ratio_updated_by_physics`: Total cloud condensate mixing ratio updated by physics
    * `real(kind=kind_phys)`: units = kg kg-1
* `accumulated_change_of_air_temperature_due_to_FA_scheme`: Accumulated change of air temperature due to FA scheme
    * `real(kind=kind_phys)`: units = K
* `mass_weighted_rime_factor_of_new_state`: Mass weighted rime factor of new state
    * `real(kind=kind_phys)`: units = kg kg-1
* `relative_humidity_threshold_for_condensation`: Relative humidity threshold for condensation
    * `real(kind=kind_phys)`: units = none
* `surface_friction_velocity_for_momentum`: Surface friction velocity for momentum
    * `real(kind=kind_phys)`: units = m s-1
* `ratio_of_height_to_monin_obukhov_length`: Ratio of height to monin obukhov length
    * `real(kind=kind_phys)`: units = none
* `surface_temperature_scale`: Surface temperature scale
    * `real(kind=kind_phys)`: units = K
* `surface_upward_latent_heat_flux`: Surface upward latent heat flux
    * `real(kind=kind_phys)`: units = W m-2
* `surface_exchange_coefficient_for_heat`: Surface exchange coefficient for heat
    * `real(kind=kind_phys)`: units = W m-2 K-1
* `surface_exchange_coefficient_for_moisture`: Surface exchange coefficient for moisture
    * `real(kind=kind_phys)`: units = kg m-2 s-1
* `air_potential_temperature_at_2m`: Air potential temperature at 2m
    * `real(kind=kind_phys)`: units = K
* `surface_exchange_coefficient_for_heat_at_2m`: Surface exchange coefficient for heat at 2m
    * `real(kind=kind_phys)`: units = m s-1
* `surface_exchange_coefficient_for_moisture_at_2m`: Surface exchange coefficient for moisture at 2m
    * `real(kind=kind_phys)`: units = m s-1
