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

The package contains ``revise/revise.yaml``, so installed code can construct
``REVISEPipeline()`` without a checkout-relative configuration path.

Dependency layers
-----------------

The base package contains reconstruction, benchmarking, the default OT
implementation, clustering, and the core scientific stack. Additional
capabilities are installed only when needed:

.. list-table::
   :header-rows: 1
   :widths: 2 3 4

   * - Capability
     - Source-checkout install
     - Purpose
   * - Additional OT implementation
     - ``python -m pip install ".[tacco]"``
     - Adds another selectable OT implementation, such as TACCO
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
available but does not activate it.

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
at ``https://zenodo.org/records/17705737``. Real-data end-to-end testing remains
a separate validation step.
