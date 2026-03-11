import os
import re

files_to_fix = [
    './ler/gw_source_population/cbc_source_parameter_distribution.py',
    './ler/utils/function_interpolation.py',
    './ler/utils/cosmological_coversions.py',
    './ler/lens_galaxy_population/lens_functions.py',
    './ler/lens_galaxy_population/cross_section_interpolator.py',
]

for filepath in files_to_fix:
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            content = f.read()
            
        new_content = re.sub(r',\s*inline="always"', '', content)
        if new_content != content:
            with open(filepath, 'w') as f:
                f.write(new_content)
            print(f"Updated {filepath}")
