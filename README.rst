.. figure:: https://raw.githubusercontent.com/wiki/LIVVkit/LIVVkit/imgs/livvkit.png
   :alt:

LIVVkit Extensions (LEX)
========================

This repository holds a collection of extensions to `LIVVkit <https://livvkit.github.io/Docs/index.html>`_
for validation and custom analyses of ice sheet models and their associated Earth
system models.

LEX was first described in [Evans2018]_

.. [Evans2018] Evans, K.J., J.H. Kennedy, D. Lu, M.M. Forrester, S. Price, J. Fyke,
   A.R. Bennett, M.J. Hoffman, I. Tezaur, C.S. Zender, and M. Vizcaino (In Review).
   LIVVkit 2.1: Automated and extensible ice sheet model validation.
   *Geoscientific Model Development.*

Dependencies
------------
The Python dependencies are listed in `pyproject.toml`, but this suite depends on
a set of reanalysis and observational datasets, which are part of E3SM-Diags, available
on several DOE supported machines, including Perlmutter at NERSC, and Chrysalis at ANL's LCRC.

Currently, LEX is designed to run on Perlmutter, but future work is planned to support other
machines.


Environment setup
-----------------

For setting up an environment to which lex and dependencies will be installed,
conda and Python virtualenv are documented here.

- Python virtualenv
.. code:: bash
    git clone https://github.com/LIVVkit/lex.git
    cd lex
    python -m venv .env
    source .env/bin/activate
    pip install --upgrade pip  # Needed if the system pip version < 21.3
    pip install -e .

This will create a virtual environment at `lex/.env`, and install the LEX package as editable with
all its Python requirements to run.

- Conda environment
.. code:: bash
    git clone https://github.com/LIVVkit/lex.git
    cd lex
    [conda, mamba] create -n lex_env python --file requirements.txt
    [conda, mamba] activate lex_env
    pip install -e .

- Not documented here, but also available for environment management
    - `uv <https://docs.astral.sh/uv/>`_
    - `pixi <https://pixi.sh/latest/>`_


Basic usage
-----------

Within the `lex/config` directory there are templates for ELM r05 and r025 resolutions, as well
as pre-existing configurations for several current runs.

To execute any of these analyses, point ``livv`` (LIVVkit's command line interface)
to any of these extensions via the the ``-V/--validate`` option.

For example, to run the minimal example extension, place the output website in the
`vv_test` directory, and serve the output website you'd run this command:

.. code:: bash

    livv -V config/example/example.yml -o vv_test -s


*Note:* All the extension configurations files assume you are working from the
top level ``lex`` directory. You *can* run any of these extensions from any
directory, but you will need to edit the paths in the YAML configuration files so
that ``livv`` can find the required files.

Likewise, you can also apply these analyses to any new model run [#]_ by adjusting
the paths to point to your model run.

.. [#] This assumes the new data files conform to the format of the included data
   files. That is, an extension that analyses output from the CISM-Albany ice
   sheet model will likely be able to analyze any similar CISM-Albany simulation,
   but likely would *not* be able to analyze output from the PISM ice sheet
   model without "massaging" the PISM files into a CISM-Albany like structure, or
   adjusting the extension. *This is a problem we are working on for future LEX
   releases.*


Running cases on PM-CPU
---------------------------
The `run_livv.sh` script will run all the currently available analyses on pm-cpu for a
particular case, e.g.:
.. code:: bash
    cd lex
    ./run_livv.sh v2.1.r025.IGERA5ELM_MLI-deep_firn_1980_2020

The batch script provided will run all current cases on Perlmutter on a compute node in parallel
.. code:: bash
    cd lex
    sbatch run_lex_pm-cpu.sbatch


Developing a custom extension
-----------------------------

See the `LIVVkit documentation <https://livvkit.github.io/Docs/lex.html>`_ for
details on how to develop an extension. Briefly, a absolute minimum working example
is provided by the ``examples/`` extension, which should be copied to provide the
basis for your new extension. All extensions are required to contain a minimal working
example set of data such that they can be run an executed on any machine.

For extensions that require data for which re-host permission cannot be granted,
they must include documentation on how to acquire and use the data as well as either
a small set of processed data or a set of "fake" example data.


Issues, questions, comments, etc.?
----------------------------------

If you would like to suggest features, request tests, discuss contributions,
report bugs, ask questions, or contact us for any reason, use the
`LIVVkit issue tracker <https://github.com/LIVVkit/LIVVkit/issues>`_.
`LEX issue tracker <https://github.com/LIVVkit/lex/issues>`_.

Want to send us a private message?

**Michael E. Kelleher**
:github: @mkstratos

**Joseph H. Kennedy**
:github: @jhkennedy

**Katherine J. Evans**
:github: @kevans32
