from ler.rates import LeR

ler = LeR(
    npool=6,
    source_priors=dict(
        merger_rate_density='merger_rate_density_bbh_popI_II_oguri2018',
        source_frame_masses='binary_masses_BBH_popI_II_powerlaw_gaussian',
    ),
    source_priors_params=dict(
        merger_rate_density={'R0': 25e-09, 'b2': 1.6, 'b3': 2.0, 'b4': 30},
        source_frame_masses={'mminbh': 4.98, 'mmaxbh': 112.5, 'alpha': 3.78, 'mu_g': 32.27, 'sigma_g': 3.88, 'lambda_peak': 0.03, 'delta_m': 4.8, 'beta': 0.81},
    ),
    lens_functions=dict(
        strong_lensing_condition='rjs_with_cross_section',
        optical_depth='optical_depth_SIE_hemanta',
    ),
    lens_priors=dict(
        velocity_dispersion='velocity_dispersion_ewoud',
        axis_ratio='axis_ratio_rayleigh',
    ),
    lens_priors_params=dict(
        velocity_dispersion=None,
        axis_ratio=dict(q_min=0.2, q_max=1.),
    ),
    create_new_interpolator=dict(
        redshift_distribution=dict(create_new=False, resolution=500),
        velocity_dispersion=dict(create_new=False, resolution=100),
        optical_depth=dict(create_new=False, resolution=100),
    ),
)

ler.batch_size = 25000
unlensed_param = ler.unlensed_cbc_statistics(size=1000000, json_file="./unlensed_bbh_1M_new.json",resume=True)

lensed_param = ler.lensed_cbc_statistics(size=1000000, json_file="./lensed_bbh_1M_new.json", resume=True)
