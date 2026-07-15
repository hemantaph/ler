import os
import numpy as np
import contextlib
from .utils import append_json, get_param_from_json

class TrainingDataGenerator():

    def __init__(self,
        npool=4,
        z_min=0.0,
        z_max=5.0,
        verbose=True,
        **kwargs,  # ler and gwsnr arguments
    ):

        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.verbose = verbose

        pdet_kwargs = dict(
            snr_th=10.0,
            snr_th_net=10.0,
            pdet_type="boolean",
            distribution_type="noncentral_chi2",
            include_optimal_snr=True,
            include_observed_snr=False,
        )

        self.ler_init_args = dict(
            event_type="BBH",
            cosmology=None,
            ler_directory="./ler_data",
            spin_zero=False,
            spin_precession=True,
            # gwsnr args
            mtot_min=2*4.98, # 4.98 Mo is the minimum component mass of BBH systems in GWTC-3
            mtot_max=2*112.5+10.0, # 112.5 Mo is the maximum component mass of BBH systems in GWTC-3. 10.0 Mo is added to avoid edge effects.
            ratio_min=0.1,
            ratio_max=1.0,
            spin_max=0.99,
            mtot_resolution=200,
            ratio_resolution=20,
            spin_resolution=10,
            sampling_frequency=2048.0,
            waveform_approximant="IMRPhenomXPHM",
            minimum_frequency=20.0,
            snr_method="inner_product",
            snr_type="optimal_snr",
            pdet_kwargs=pdet_kwargs,
            psds=None,
            ifos=None,
            interpolator_dir="./interpolator_json",
            create_new_interpolator=False,
            gwsnr_verbose=True,
            multiprocessing_verbose=True,
            mtot_cut=False,
        )
        self.ler_init_args.update(kwargs)
        if "spin_precessing" in self.ler_init_args:
            self.ler_init_args["spin_precession"] = self.ler_init_args["spin_precessing"]
        pdet_kwargs.update(self.ler_init_args.get("pdet_kwargs") or {})
        self.ler_init_args["pdet_kwargs"] = pdet_kwargs

    @staticmethod
    def _snr_key(gw_param):
        for key in ("optimal_snr_net", "snr_net"):
            if key in gw_param:
                return key
        raise KeyError("GW parameter dictionary must contain 'optimal_snr_net' or 'snr_net'.")

    def gw_parameters_generator(self, 
        size, 
        batch_size=100000, 
        snr_recalculation=True,
        trim_to_size=False, 
        verbose=True, 
        replace=False, 
        data_distribution_range=None,
        output_jsonfile="gw_parameters.json",
    ):

        if data_distribution_range is None:
            data_distribution_range = [0, 2, 4, 6, 8, 10, 12, 14, 16, 100]

        args = self.ler_init_args.copy()
        if snr_recalculation:
            snr_method = 'interpolation_aligned_spins'
        else:
            snr_method = args['snr_method']

        # 
        from ler.rates import GWRATES

        # ler initialization
        ler = GWRATES(
            npool=self.npool,
            z_min=self.z_min,
            z_max=self.z_max,  # becareful with this value
            verbose=self.verbose,
            # ler
            event_type=args['event_type'],
            cosmology=args['cosmology'],
            ler_directory=args['ler_directory'],
            spin_zero=args['spin_zero'],
            spin_precession=args['spin_precession'],
            # gwsnr args
            mtot_min=args['mtot_min'],
            mtot_max=args['mtot_max'],
            ratio_min=args['ratio_min'],
            ratio_max=args['ratio_max'],
            spin_max=args['spin_max'],
            mtot_resolution=args['mtot_resolution'],
            ratio_resolution=args['ratio_resolution'],
            spin_resolution=args['spin_resolution'],
            sampling_frequency=args['sampling_frequency'],
            waveform_approximant=args['waveform_approximant'],
            minimum_frequency=args['minimum_frequency'],
            snr_method=snr_method,
            snr_type=args['snr_type'],
            pdet_kwargs=args['pdet_kwargs'],
            psds=args['psds'],
            ifos=args['ifos'],
            interpolator_dir=args['interpolator_dir'],
            create_new_interpolator=args['create_new_interpolator'],
            gwsnr_verbose=args['gwsnr_verbose'],
            multiprocessing_verbose=args['multiprocessing_verbose'],
            mtot_cut=args['mtot_cut'],
        )
        ler.batch_size = batch_size
        snr_recalculator = None
        if snr_recalculation:
            from gwsnr import GWSNR
            snr_recalculator = GWSNR(
                npool=self.npool,
                snr_method='inner_product',
                snr_type=args['snr_type'],
                pdet_kwargs=args['pdet_kwargs'],
                mtot_min=args['mtot_min'],
                mtot_max=args['mtot_max'],
                ratio_min=args['ratio_min'],
                ratio_max=args['ratio_max'],
                spin_max=args['spin_max'],
                mtot_resolution=args['mtot_resolution'],
                ratio_resolution=args['ratio_resolution'],
                spin_resolution=args['spin_resolution'],
                sampling_frequency=args['sampling_frequency'],
                waveform_approximant=args['waveform_approximant'],
                minimum_frequency=args['minimum_frequency'],
                psds=args['psds'],
                ifos=args['ifos'],
                interpolator_dir=args['interpolator_dir'],
                create_new_interpolator=args['create_new_interpolator'],
                gwsnr_verbose=args['gwsnr_verbose'],
                multiprocessing_verbose=args['multiprocessing_verbose'],
                mtot_cut=args['mtot_cut'],
            )

        # path to save parameters
        json_path = f"{args['ler_directory']}/{output_jsonfile}"
        if replace:
            if os.path.exists(json_path):
                os.remove(json_path)
            len_final = 0
        else:
            if os.path.exists(json_path):
                gw_param = get_param_from_json(json_path)
                len_final = len(gw_param[self._snr_key(gw_param)])
                print(f'current size of the json file: {len_final}\n')
            else:
                len_final = 0

        print(f'total event to collect: {size}\n')
        while len_final<size:
            with contextlib.redirect_stdout(None):
                gw_param = ler.gw_cbc_statistics(size=ler.batch_size, resume=False)

            if data_distribution_range is not None:
                gw_param = self.helper_data_distribution(gw_param, data_distribution_range)

                if gw_param is None:
                    continue

            if snr_recalculation:
                snrs = snr_recalculator.optimal_snr(gw_param_dict=gw_param)
                gw_param.update(snrs)

                if data_distribution_range is not None:
                    gw_param = self.helper_data_distribution(gw_param, data_distribution_range)

            if gw_param is None:
                print("No data in one of the given range")
                continue
            # save the parameters
            append_json(json_path, gw_param, replace=False);

            # print(f"Collected number of events: {len_}")
            len_final += len(gw_param[self._snr_key(gw_param)])
            if verbose:
                print(f"Collected number of events: {len_final}")

        if trim_to_size:
            gw_param = get_param_from_json(json_path)
            for key, value in gw_param.items():
                gw_param[key] = value[:size]
            append_json(json_path, gw_param, replace=True);
            len_final = len(gw_param[self._snr_key(gw_param)])

        print(f"final size: {len_final}\n")
        print(f"json file saved at: {json_path}\n")

    def helper_data_distribution(self, gw_param, data_distribution_range):
    # optimal SNR 
        snr = np.array(gw_param[self._snr_key(gw_param)])

        idx_arr = []
        snr_range = np.array(data_distribution_range)
        len_ = len(snr_range) 
        len_arr = []  # size of len_arr is len_-1  
        for j in range(len_-1):
            idx_ = np.argwhere((snr>=snr_range[j]) & (snr<snr_range[j+1])).flatten()
            idx_arr.append(idx_)
            len_arr.append(len(idx_))

        idx_arr = np.array(idx_arr, dtype=object)
        len_ref = min(len_arr)
        
        if len_ref==0:
            print("No data in one of the given range")
            return None
        else:
            gw_param_final = {}
            for j, len_ in enumerate(len_arr):  # loop over snr range
                idx_buffer = np.asarray(
                    np.random.choice(idx_arr[j], len_ref, replace=False),
                    dtype=int,
                )

                for key, value in gw_param.items(): 
                    
                    try:
                        buffer_ = np.array(value)[idx_buffer]
                    except IndexError:
                        print(f"IndexError")
                        print(f"key: {key}, len(value): {len(value)}, len(idx_buffer): {len(idx_buffer)}")
                        print(f"rerun the code again with: replace=False")
                        return None
                    if j==0:
                        gw_param_final[key] = buffer_
                    else:
                        gw_param_final[key] = np.concatenate([gw_param_final[key], buffer_])

            return gw_param_final

    def combine_dicts(self,
        file_name_list=None,
        path_list=None, 
        detector='L1',
        parameter_list=None,
        output_jsonfile="combined_data.json",
    ):

        if parameter_list is None:
            parameter_list = ['mass_1', 'mass_2', 'luminosity_distance', 'theta_jn', 'psi', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2']

        detector_snr_key = f"optimal_snr_{detector}"
        parameter_list = list(parameter_list) + [
            detector_snr_key,
            detector,
            "optimal_snr_net",
            "snr_net",
        ]
        combined_dict = {}

        if file_name_list is not None:
            path_list = [f"{self.ler_init_args['ler_directory']}/{file_name}" for file_name in file_name_list]
        elif path_list is None:
            print("Please provide either file_name_list or path_list")
            return None

        for path in path_list:
            data = get_param_from_json(path)
            for key, value in data.items():
                if key in parameter_list:
                    if key in combined_dict:
                        combined_dict[key] = np.concatenate([combined_dict[key], value])
                    else:
                        combined_dict[key] = value

        if "optimal_snr_net" not in combined_dict:
            if detector_snr_key in combined_dict:
                combined_dict["optimal_snr_net"] = combined_dict[detector_snr_key]
            elif detector in combined_dict:
                combined_dict["optimal_snr_net"] = combined_dict[detector]
            elif "snr_net" in combined_dict:
                combined_dict["optimal_snr_net"] = combined_dict["snr_net"]
            else:
                raise KeyError(
                    f"Could not find '{detector_snr_key}', '{detector}', or 'snr_net' "
                    "to build 'optimal_snr_net'."
                )

        json_path = f"{self.ler_init_args['ler_directory']}/{output_jsonfile}"
        print(f"json file saved at: {json_path}\n")
        append_json(json_path, combined_dict, replace=True);

    def delete_json_file(self, path_list):
        for path in path_list:
            if os.path.exists(path):
                os.remove(path)
