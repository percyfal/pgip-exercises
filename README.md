# Population Genomics in Practice exercises

Collection of exercises for course [Population Genomics in
Practice](https://uppsala.instructure.com/courses/52168)

**Under construction**

# Requirements

Create a conda environment called `pgip` using the environment file

	conda create --name pgip python=3.9
	conda activate pgip
	mamba env update -f=environment.yaml

and install additional python requirements with pip

	pip install -r requirements.txt

Install pgip and bash jupyter kernels

	python -m bash_kernel.install
	python -m ipykernel install --user --name pgip --display-name "Population Genomics in Practice (Python)"

Edit MyST markdown files and build the docs with

	make

# Development

## Test / development data

NB: once acceptable test data has been generated, relevant output is
added to the git repo and there is no need to regenerate it.

Install additional data development requirements

	pip install -r requirements-dev.txt

Test data can be generated in data directory as follows

	cd data
	./make_test_data.py

## bash code cell blocks

bash code cell blocks by default return garbled output - see
<https://github.com/takluyver/bash_kernel/issues/115>. The solution is
to set bracketed output off. It suffices to add a code cell at the top
of a notebook that runs the bash kernel:


	```{code-cell} bash
	:"tags": ["remove-output"]
	bind 'set enable-bracketed-paste off'
	```
