Installation
============

REVISE supports Python 3.10 and 3.11. Use a clean environment so the selected
OT solver and scientific stack are unambiguous.

Install a published release
---------------------------

.. code-block:: bash

   python -m pip install revise-svc

Published releases can lag the repository. Check the installed version before
expecting a newly documented source-checkout interface.

Install the current source
--------------------------

.. code-block:: bash

   git clone https://github.com/wuys13/REVISE.git
   cd REVISE
   python -m pip install .

For development and repository tests:

.. code-block:: bash

   python -m pip install -e ".[dev]"

The installed Application command is ``revise-reconstruct``. A source checkout
may run the same request with ``python reconstruct.py``. Sim2Real-ST Benchmark
reproduction is a checkout workflow invoked with
``python reproduce/benchmark_main.py``; it is not an installed console script.

Choose the solver dependency
----------------------------

The base package contains reconstruction, benchmarking, POT, clustering, and
the core scientific stack. POT is the default OT implementation. REVISE also
supports TACCO as an alternative OT method; the maintained Xenium cluster-mode
template selects TACCO 0.5.0 explicitly:

.. code-block:: bash

   python -m pip install "revise-svc[tacco]"

From a source checkout, use ``python -m pip install ".[tacco]"``. Set
``algorithm.ot_method`` to ``pot`` or ``tacco`` explicitly when choosing a
different supported solver for an Application request.

Optional capabilities
---------------------

These extras provide additional data reading or downstream bioinformatics
capabilities. Install only the extra required by the workflow:

.. list-table::
   :header-rows: 1
   :widths: 2 3 4

   * - Capability
     - Source-checkout install
     - Purpose
   * - Pathway analysis
     - ``python -m pip install ".[pathway]"``
     - dependencies used by pathway notebooks
   * - Cell-cell interaction
     - ``python -m pip install ".[cci]"``
     - dependencies used by CCI notebooks; databases are separate resources
   * - Trajectory analysis
     - ``python -m pip install ".[trajectory]"``
     - dependencies used by trajectory notebooks
   * - SpatialData input
     - ``python -m pip install ".[spatialdata]"``
     - SpatialData/Zarr ST input support

After a matching package version is published, replace ``.`` with
``revise-svc`` in those commands. Installing an extra makes a capability
available; it does not choose a route or mode, and it does not download paper
data.

.. _application-templates:

Application templates
---------------------

A source checkout exposes exactly three maintained files under
``configs/application/``: ``VisiumHD.yaml``, ``Xenium.yaml``, and
``Visium.yaml``. Copy one to a working directory before editing it. The
installed package carries identical resources under ``revise.application/templates``;
copy a packaged template before editing it:

.. code-block:: python

   from importlib.resources import as_file, files
   from shutil import copyfile

   resource = files("revise.application").joinpath("templates", "Xenium.yaml")
   with as_file(resource) as source:
       copyfile(source, "Xenium.yaml")

Then run ``revise-reconstruct --config <local-copy>.yaml``. A source checkout
may instead use ``python reconstruct.py --config configs/application/<name>.yaml``.
The package does not distribute the real
P1CRC, P2CRC, or mouse-brain H5AD inputs; use the reproduction downloads in
the repository README.

Build the documentation
-----------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   sphinx-build -W --keep-going -b html docs /tmp/revise-docs-html

The documentation build renders the preserved notebook snapshots but never
executes them.
