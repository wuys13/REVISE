Installation
============

Current source contract
-----------------------

The unified CLI and optional groups documented here describe the current
repository. Install that exact code with:

.. code-block:: bash

   git clone https://github.com/wuys13/REVISE.git
   cd REVISE
   python -m pip install .

The published package can be installed with ``python -m pip install
revise-svc``, but releases can lag the repository. Check the installed version
before expecting the current CLI or optional groups.

For development and test tools:

.. code-block:: bash

   python -m pip install -e ".[dev]"

The package contains its engine configuration and Application templates. The
normal installed entry point is ``revise-reconstruct``; Python callers use
``from reconstruct import run_application`` with the same YAML.

Application templates
---------------------

A source checkout exposes ``configs/application/Xenium_T.yaml``,
``configs/application/Xenium_Fib.yaml``, ``configs/application/Xenium_Mono.yaml``,
``configs/application/VisiumHD.yaml``, and ``configs/application/Visium.yaml``.
An installed package carries the same files under
``revise.application/templates``. Copy one into the project with
the standard ``importlib.resources`` interface, then edit the copy:

.. code-block:: python

   from importlib.resources import as_file, files
   from shutil import copyfile

   resource = files("revise.application").joinpath("templates", "Xenium_T.yaml")
   with as_file(resource) as source:
       copyfile(source, "Xenium_T.yaml")

This resource-copy recipe does not add a REVISE CLI or Python API. Run the
copied request with ``revise-reconstruct --config Xenium_T.yaml``. Passing an official bare name or
``configs/application/<name>.yaml`` uses an existing external mirror first and
the package resource only when that file is absent.

Dependency layers
-----------------

The base package contains reconstruction, benchmarking, the POT implementation,
clustering, and the core scientific stack. Additional capabilities are
installed only when needed:

.. list-table::
   :header-rows: 1
   :widths: 2 3 4

   * - Capability
     - Source-checkout install
     - Purpose
   * - Standard sc-SVC default solver
     - ``python -m pip install ".[tacco]"``
     - Installs TACCO 0.5.0, required by the default standard sc-SVC route
   * - Pathway analysis
     - ``python -m pip install ".[pathway]"``
     - Dependencies used by pathway notebooks
   * - Cell-cell interaction analysis
     - ``python -m pip install ".[cci]"``
     - Dependencies used by CCI notebooks
   * - Trajectory analysis
     - ``python -m pip install ".[trajectory]"``
     - Dependencies used by trajectory notebooks
   * - SpatialData input
     - ``python -m pip install ".[spatialdata]"``
     - SpatialData/Zarr input support

After a matching package version is published, replace ``.`` with
``revise-svc`` in those commands. Optional dependency selection and runtime
algorithm selection are separate: installing an extra makes that capability
available but does not activate it. Without the TACCO extra, users who accept a
different reconstruction algorithm must explicitly set
``algorithm.ot_method: pot`` in the application YAML;
REVISE never selects POT as an automatic fallback.

For the maintained Xenium request, install
``python -m pip install "revise-svc[tacco]"``. Missing or incompatible TACCO
fails with directed installation guidance and does not fall back to POT.

The CCI extra does not download a CellPhoneDB database. None of these commands
downloads research data or external analysis resources.

Documentation build
-------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   sphinx-build -W --keep-going -b html docs /tmp/revise-docs-html

Research data
-------------

The paper benchmark/application datasets and reproduced results are available
at ``https://zenodo.org/records/17705737``. The real P1CRC and mouse-brain H5AD
files are not distributed inside the package. Real-data end-to-end testing
remains a separate validation step; software wiring does not establish
biological validity or POT/TACCO numerical equivalence.
