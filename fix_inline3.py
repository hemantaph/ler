import os
import re

files_to_fix = [
    './ler/lens_galaxy_population/optical_depth.py',
    './ler/lens_galaxy_population/sampler_functions.py'
]

for filepath in files_to_fix:
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            content = f.read()
            
        # specifically target the lens sampler closures
        new_content = re.sub(r'@njit\(cache=True,\s*inline="always"\)\n\s*def (z[sl]_[a-z0-9_]+)', r'@njit(cache=True)\n        def \1', content)
        
        # for optical_depth closures
        new_content = re.sub(r'@njit\(cache=True,\s*inline="always"\)\n\s*def (theta_E_[a-z0-9_]+)', r'@njit(cache=True)\n        def \1', new_content)

        if new_content != content:
            with open(filepath, 'w') as f:
                f.write(new_content)
            print(f"Updated {filepath}")
