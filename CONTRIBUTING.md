# Community Guidelines

Thank you for your interest in `ler`. Contributions, bug reports, questions,
and examples are welcome. These guidelines explain how to contribute, report
problems, and seek support.

## Ways to contribute

Useful contributions include:

- bug fixes
- documentation improvements
- examples and tutorials
- tests
- new features or scientific models
- performance improvements

For larger changes, please open an issue first so we can discuss the scope and
avoid duplicated work.

## Development setup

1. Fork the repository and clone your fork.
2. Create a Python environment with Python 3.10 or newer.
3. Install `ler` in editable mode with development tools:

```bash
pip install -e ".[dev]"
```

4. Create a branch for your work:

```bash
git checkout -b my-change
```

## Making changes

Please keep pull requests focused and easy to review. A good pull request:

- describes the problem and the proposed solution
- adds or updates tests when behavior changes
- updates documentation or examples when user-facing behavior changes
- keeps unrelated formatting or refactoring out of the change

Run the test suite before opening a pull request:

```bash
python -m pytest tests/ --tb=short
```

For changes that may affect slow numerical paths, also run:

```bash
python -m pytest tests/ --run-slow --tb=short
```


## Reporting issues or problems

Please report bugs through the GitHub issue tracker:

https://github.com/hemantaph/ler/issues

When possible, include:

- the `ler` version
- your Python version and operating system
- the command or code that produced the problem
- the full error message or traceback
- a small reproducible example
- any relevant input files, configuration, or random seeds

If the issue is about a scientific result, please also describe the model,
assumptions, and expected behavior.

## Seeking support

For questions about installation, usage, examples, or unexpected behavior,
please use the GitHub issue tracker:

https://github.com/hemantaph/ler/issues

Before opening a new issue, please check:

- the documentation: https://ler.hemantaph.com/
- existing issues: https://github.com/hemantaph/ler/issues
- the examples in the `docs/examples/` directory

## Code of conduct

Please be respectful and constructive in all project discussions. We welcome
contributions from people with different backgrounds, experience levels, and
research interests. Harassment, personal attacks, and other exclusionary
behavior are not acceptable.

## Citation

If `ler` supports your research, please cite the project as described in the
documentation and repository metadata.
