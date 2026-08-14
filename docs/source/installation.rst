Installation
============

REVISE supports Python 3.10 and 3.11. Use a clean environment so the selected
OT solver and scientific stack are unambiguous.

Install a published release
---------------------------

.. code-block:: bash

   python -m pip install revise-svc

Published releases can lag the repository. Check the installed
version before expecting a newly documented source-checkout interface.

Install the current source
--------------------------

.. code-block:: bash

   git clone https://github.com/wuys13/REVISE.git
   cd REVISE
   python -m pip install .

For development and repository tests:

.. code-block:: bash

   python -m pip install -e ".[dev]"

The installed Application command is ``revise-reconstruct``. A source
checkout may run the same request with ``python reconstruct.py``. Sim2Real
Benchmark reproduction is a checkout workflow invoked with
``python reproduce/benchmark_main.py``; it is not an installed console script.

Choose the solver dependency
----------------------------

The base package contains reconstruction, benchmarking, the POT implementation,
clustering, the core scientific stack, and installs POT 0.9.5. The
maintained ``sp-SVC`` and ``sc-SVC-sr`` templates select POT.

Standard sc-SVC default solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default solver is TACCO 0.5.0 for both Global Anchoring
and Local Refinement. Install its extra before running the maintained
templates:

.. code-block:: bash

   python -m pip install "revise-svc[tacco]"

From a source checkout, use:

.. code-block:: bash

   python -m pip install ".[tacco]"

A missing, incompatible, or failed TACCO installation stops the run with an
installation error. If a scientifically different POT run is acceptable,
select ``algorithm.ot_method: pot`` explicitly in the Application YAML. REVISE
never selects POT as an automatic fallback.

Optional capabilities
---------------------

Install only the extra required by the workflow:

.. list-table::
   :header-rows: 1
   :widths: 2 3 4

   * - Capability
     - Source-checkout install
     - Purpose
   * - TACCO OT
     - ``python -m pip install ".[tacco]"``
     - TACCO 0.5.0, required by the default standard sc-SVC route
   * - Pathway analysis
     - ``python -m pip install ".[pathway]"``
     - dependencies used by pathway notebooks
   * - Cell-cell interaction
     - ``python -m pip install ".[cci]"``
     - dependencies used by CCI notebooks
   * - Trajectory analysis
     - ``python -m pip install ".[trajectory]"``
     - dependencies used by trajectory notebooks
   * - SpatialData input
     - ``python -m pip install ".[spatialdata]"``
     - SpatialData/Zarr ST input support

After a matching package version is published, replace ``.`` with
``revise-svc`` in those commands. Installing an extra makes a capability
available; it does not select a reconstruction mode or solver. The CCI extra
does not download a CellPhoneDB database, and no install command downloads
paper data or other research resources.

Verified project environment
----------------------------

The maintained local validation environment is Python 3.10.14 with TACCO
0.5.0, Scanpy 1.11.4, Squidpy 1.6.3, and POT 0.9.5. Project runs use a base
seed of 42 and a single-thread execution environment. These are execution and
validation conditions; the package installer does not configure global thread
environment variables for the user.

The repository also carries release-critical constraints for Python 3.10 and
3.11. They are tested constraint sets, not complete transitive lock files.

Application templates
---------------------

A source checkout exposes the five files under ``configs/application/``. The
installed package carries the same requests under
``revise.application/templates``. Copy a packaged template into your project
before editing it:

.. code-block:: python

   from importlib.resources import as_file, files
   from shutil import copyfile

   resource = files("revise.application").joinpath("templates", "VisiumHD.yaml")
   with as_file(resource) as source:
       copyfile(source, "VisiumHD.yaml")

Then follow :doc:`quickstart`. The package does not distribute the real P1CRC,
P2CRC, or mouse-brain H5AD inputs. Paper data and reproduced results are
available at ``https://zenodo.org/records/17705737``; keep local ``raw_data/``
outside Git.

Build the documentation
-----------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   sphinx-build -W --keep-going -b html docs /tmp/revise-docs-html

The documentation build renders the eleven curated notebook snapshots but never
executes them.
