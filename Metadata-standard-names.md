# CCPP Standard Name Library
#### Table of Contents
* [dimensions](#dimensions)
* [constants](#constants)
* [coordinates](#coordinates)
* [state_variables](#state_variables)
* [diagnostics](#diagnostics)
* [constituents](#constituents)
* [standard_variables](#standard_variables)

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
* `vertical_layer_dimension`: number of vertical levels
    * `integer`: units = count
* `vertical_interface_dimension`: number of vertical interfaces
    * `integer`: units = count
* `vertical_layer_index`: index of a particular vertical level
    * `integer`: units = count
* `vertical_interface_index`: index of a particular vertical interface
    * `integer`: units = count
* `index_of_bottom_vertical_layer`: Index of bottom vertical layer
    * `integer`: units = index
* `index_of_top_vertical_layer`: Index of top vertical layer
    * `integer`: units = index
* `index_of_bottom_vertical_interface`: Index of bottom vertical interface
    * `integer`: units = index
* `index_of_top_vertical_interface`: Index of top vertical interface
    * `integer`: units = index
* `thread_block_dimension`: Total number of thread blocks which the host model may use to call CCPP physics run groups during the CCPP run phase.
    * `integer`: units = none
* `thread_block_index`: Number of current thread block. This variable may only be used during CCPP run phase
    * `integer`: units = none
## constants
* `avogadro_number`: Avogadro number
    * `real(kind=kind_phys)`: units = molecules mole-1
* `base_state_surface_pressure_for_hybrid_vertical_coordinate`: Base state surface pressure for hybrid vertical coordinate
    * `real(kind=kind_phys)`: units = Pa
* `boltzmann_constant`: Boltzmann constant
    * `real(kind=kind_phys)`: units = J K-1
* `gas_constant_of_dry_air`: Gas constant of dry air
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `seconds_in_calendar_day`: Seconds in calendar day
    * `integer(kind=kind_phys)`: units = s
* `specific_heat_of_dry_air_at_constant_pressure`: Specific heat of dry air at constant pressure
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `specific_heat_of_liquid_water_at_20c`: Specific heat of liquid water at 20c
    * `real(kind=kind_phys)`: units = J kg-1 K-1
* `latent_heat_of_vaporization_of_water_at_0c`: latent heat of vaporization of water at 0C
    * `real(kind=kind_phys)`: units = J kg-1
* `density_of_dry_air_at_stp`: Density of dry air at stp
    * `real(kind=kind_phys)`: units = kg m-3
* `density_of_liquid_water_at_0c`: density of liquid water at 0C
    * `real(kind=kind_phys)`: units = kg m-3
* `ratio_of_water_vapor_to_dry_air_gas_constants_minus_one`: (Rwv / Rdair) - 1.0
    * `real(kind=kind_phys)`: units = 1
## coordinates
* `latitude`: Latitude
    * `real(kind=kind_phys)`: units = radians
* `longitude`: Longitude
    * `real(kind=kind_phys)`: units = radians
* `gravitational_acceleration`: Gravitational acceleration
    * `real(kind=kind_phys)`: units = m s-2
* `cell_area`: Cell area
    * `real(kind=kind_phys)`: units = steradian
* `cell_weight`: Cell weight
    * `real(kind=kind_phys)`: units = none
## state_variables
Note that appending '_from_previous_timestep' to standard_names in this section yields another valid standard_name
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
* `surface_pressure_of_dry_air`: Surface pressure of dry air
    * `real(kind=kind_phys)`: units = Pa
* `geopotential_at_surface`: Geopotential at surface
    * `real(kind=kind_phys)`: units = m2 s-2
* `air_temperature`: Air temperature
    * `real(kind=kind_phys)`: units = K
* `air_temperature_on_previous_timestep`: Air temperature on previous timestep
    * `real(kind=kind_phys)`: units = K
* `equivalent_potential_temperature`: Equivalent potential temperature
    * `real(kind=kind_phys)`: units = K
* `x_wind`: Zonal wind
    * `real(kind=kind_phys)`: units = m s-1
* `y_wind`: Meridional wind
    * `real(kind=kind_phys)`: units = m s-1
* `dry_static_energy`: Dry static energy Content of Atmosphere Layer
    * `real(kind=kind_phys)`: units = J kg-1
* `lagrangian_vertical_coordinate`: flag indicating if vertical coordinate is lagrangian
    * `logical(kind=)`: units = flag
* `lagrangian_tendency_of_air_pressure`: Vertical pressure velocity
    * `real(kind=kind_phys)`: units = Pa s-1
* `air_pressure`: Midpoint air pressure
    * `real(kind=kind_phys)`: units = Pa
* `air_pressure_of_dry_air`: Dry midpoint pressure
    * `real(kind=kind_phys)`: units = Pa
* `pressure_thickness`: Pressure thickness
    * `real(kind=kind_phys)`: units = Pa
* `pressure_thickness_of_dry_air`: Pressure thickness of dry air
    * `real(kind=kind_phys)`: units = Pa
* `reciprocal_of_pressure_thickness`: Reciprocal of pressure thickness
    * `real(kind=kind_phys)`: units = Pa-1
* `reciprocal_of_pressure_thickness_of_dry_air`: Reciprocal of pressure thickness of dry air
    * `real(kind=kind_phys)`: units = Pa-1
* `ln_of_air_pressure`: Ln of air pressure
    * `real(kind=kind_phys)`: units = 1
* `log_of_air_pressure_of_dry_air`: Log of air pressure of dry air
    * `real(kind=kind_phys)`: units = Pa
* `inverse_dimensionless_exner_function_wrt_surface_pressure`: inverse exner function w.r.t. surface pressure, (ps/p)^(R/cp)
    * `real(kind=kind_phys)`: units = 1
* `geopotential_height`: Geopotential height
    * `real(kind=kind_phys)`: units = m
* `constituent_mixing_ratio`: Constituent mixing ratio
    * `real(kind=kind_phys)`: units = kg/kg moist or dry air depending on type
* `air_pressure_at_interface`: Air pressure at interface
    * `real(kind=kind_phys)`: units = Pa
* `air_pressure_of_dry_air_at_interface`: Air pressure of dry air at interface
    * `real(kind=kind_phys)`: units = Pa
* `ln_of_air_pressure_at_interface`: Ln of air pressure at interface
    * `real(kind=kind_phys)`: units = 1
* `ln_of_air_pressure_of_dry_air_at_interface`: Ln of air pressure of dry air at interface
    * `real(kind=kind_phys)`: units = 1
* `largest_model_top_pressure_that_allows_molecular_diffusion`: Largest model top pressure that allows molecular diffusion
    * `real(kind=kind_phys)`: units = Pa
* `flag_for_molecular_diffusion`: Flag for molecular diffusion
    * `logical(kind=kind_phys)`: units = flag
* `flag_for_physics_grid_initialization`: Flag for physics grid initialization
    * `logical(kind=kind_phys)`: units = flag
* `geopotential_height_at_interface`: Geopotential height at interface
    * `real(kind=kind_phys)`: units = m
* `column_integrated_total_kinetic_and_static_energy_of_initial_state`: Column integrated total kinetic and static energy of initial state
    * `real(kind=kind_phys)`: units = J m-2
* `column_integrated_total_kinetic_and_static_energy_of_current_state`: Column integrated total kinetic and static energy of current state
    * `real(kind=kind_phys)`: units = J m-2
* `column_integrated_total_water_of_initial_state`: Column integrated total water of initial state
    * `real(kind=kind_phys)`: units = kg m-2
* `column_integrated_total_water_of_new_state`: Column integrated total water of new state
    * `real(kind=kind_phys)`: units = kg m-2
* `tendency_of_temperature`: Change in temperature from a parameterization
    * `real(kind=kind_phys)`: units = K s-1
* `total_tendency_of_air_temperature`: Total change in temperature from a                               physics suite
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_potential_temperature`: Change in potential temperature from a parameterization
    * `real(kind=kind_phys)`: units = K s-1
* `total_tendency_of_potential_temperature`: Total tendency of potential temperature
    * `real(kind=kind_phys)`: units = K s-1
* `tendency_of_x_wind`: Change in zonal wind from a parameterization
    * `real(kind=kind_phys)`: units = m s-2
* `total_tendency_of_x_wind`: Total tendency of x wind
    * `real(kind=kind_phys)`: units = m s-2
* `tendency_of_y_wind`: Change in meridional from a parameterization
    * `real(kind=kind_phys)`: units = m s-2
* `total_tendency_of_y_wind`: Total tendency of y wind
    * `real(kind=kind_phys)`: units = m s-2
* `surface_energy_flux`: Surface energy flux
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_boundary_flux_of_total_energy`: Cumulative boundary flux of total energy
    * `real(kind=kind_phys)`: units = W m-2
* `cumulative_boundary_flux_of_total_water`: Cumulative boundary flux of total water
    * `real(kind=kind_phys)`: units = W m-2
* `reference_pressure`: reference pressure 
    * `real(kind=kind_phys)`: units = Pa
* `us_standard_atmospheric_pressure_at_sea_level`: US Standard Atmospheric pressure at sea level
    * `real(kind=kind_phys)`: units = Pa
* `reference_pressure_at_surface`: reference pressure at surface
    * `real(kind=kind_phys)`: units = Pa
* `reference_pressure_normalized_by_surface_pressure`: reference pressure normalized by surface pressure
    * `real(kind=kind_phys)`: units = none
* `exner_function`: exner function
    * `real(kind=kind_phys)`: units = none
* `potential_temperature`: potential temperature
    * `real(kind=kind_phys)`: units = K
* `potential_temperature_on_previous_timestep`: potential temperature on previous timestep
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
* `water_vapor_specific_humidity`: Water vapor specific humidity
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
* `volume_mixing_ratio_ch4`: CH4 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_co`: CO volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_co2`: CO2 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_ccl4`: CCL4 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc11`: CFC11 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc12`: CFC12 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc113`: CFC113 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_cfc22`: CFC22 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_o2`: O2 volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
* `volume_mixing_ratio_n2o`: N2O volume mixing ratio
    * `real(kind=kind_phys)`: units = kg kg-1
## standard_variables
Standard / required CCPP variables
* `ccpp_error_message`: Error message for error handling in CCPP
    * `character(kind=len=512)`: units = 1
* `ccpp_error_flag`: Error flag for error handling in CCPP
    * `integer(kind=)`: units = flag
