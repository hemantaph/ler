import os
import re


def process_file(filepath):
    with open(filepath, "r") as f:
        content = f.read()

    # Find parallel=True, cache=True, inline="always" and replace
    new_content = re.sub(
        r'parallel=True,\s*cache=True,\s*inline="always"',
        r"parallel=True, cache=True",
        content,
    )
    new_content = re.sub(
        r'cache=True,\s*parallel=True,\s*inline="always"',
        r"cache=True, parallel=True",
        new_content,
    )

    if new_content != content:
        with open(filepath, "w") as f:
            f.write(new_content)
        print(f"Updated {filepath}")


for root, _, files in os.walk("./ler"):
    for file in files:
        if file.endswith(".py"):
            process_file(os.path.join(root, file))
