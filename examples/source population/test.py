event_priors_ = dict(
    merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018",
    mass_source_frame="binary_masses_BBH_popI_II_powerlaw_gaussian",
    spin="constant_values_n_size_m_params",
    zs=None,
    geocent_time="geocent_time_uniform",
    sky_position="sky_position_uniform_bilby",
    phase="coalescence_phase_uniform_bilby",
    psi="polarization_angle_uniform_bilby",
    iota="inclination_uniform_bilby",
)
event_prior_params_ = dict(
    merger_rate_density=dict(R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30),
    mass_source_frame=dict(
        mminbh=4.98,
        mmaxbh=112.5,
        alpha=3.78,
        mu_g=32.27,
        sigma_g=3.88,
        lambda_peak=0.03,
        delta_m=4.8,
        beta=0.81,
    ),
    spin=dict(value=0.0),
    zs=None,
    geocent_time=dict(
        start_time=1238166018, end_time=1238166018 + 31536000
    ),
    sky_position=None,
    phase=None,
    psi=None,
    iota=None,
)