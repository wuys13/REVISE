Installation
============

Version ``0.1.0rc1`` is a candidate, not a claimed current PyPI release. Its
gate matrix targets Python 3.10 and 3.11. Install either an exact candidate
artifact or the current source checkout.

Candidate Wheel
---------------

The release workflow builds the Wheel once and tests those exact bytes outside
the checkout. Install that file with:

.. code-block:: bash

   python -m pip install /path/to/revise_svc-0.1.0rc1-py3-none-any.whl

The installed command is ``revise-reconstruct``. The source-tree script
``reconstruct.py`` is a compatibility wrapper and is not needed after Wheel
installation.

Install a Wheel extra from those same artifact bytes by appending the extra to
the local Wheel path, for example:

.. code-block:: bash

   python -m pip install "/path/to/revise_svc-0.1.0rc1-py3-none-any.whl[tacco]"

The candidate's PyPI locator and final canonical repository URL are pending
owner confirmation. Do not assume that an unpinned public-package install
produces these candidate bytes.

Source Checkout
---------------

From an existing checkout:

.. code-block:: bash

   python -m pip install .

For development and test tools:

.. code-block:: bash

   python -m pip install -e ".[dev]"

The package contains ``revise/revise.yaml``, so installed code can construct
``REVISEPipeline()`` without a checkout-relative configuration path.

Dependencies and extras
-----------------------

POT, Leiden, and the core scientific stack are base dependencies. From a source
checkout, optional features are separated by domain:

.. code-block:: bash

   python -m pip install ".[tacco]"
   python -m pip install ".[pathway]"
   python -m pip install ".[cci]"
   python -m pip install ".[trajectory]"
   python -m pip install ".[spatialdata]"

The ``tacco`` extra installs the supported TACCO 0.5.0 solver dependency; it
does not select TACCO automatically. Use ``--ot-method tacco`` or the GA/LR
configuration described in :doc:`configuration`. Missing or incompatible
optional dependencies fail when their feature is selected.

The CCI extra does not download the CellPhoneDB database. None of these install
commands downloads research data.

Documentation build
-------------------

Build the user documentation in a separate environment:

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   sphinx-build -W --keep-going -b html docs /tmp/revise-docs-html

Sphinx imports the package from the checkout and mocks heavy scientific
dependencies for API discovery. A successful documentation build proves link,
RST, and import-surface consistency; it is not a scientific run.

Research data
-------------

The official benchmark/application archive identity and redistribution terms
are pending owner confirmation. No numeric archive record is presented as
official for this candidate. Real-data end-to-end testing remains deferred
until the other candidate gates are complete and the owner authorizes the data
download.
