Quick Start
===========

This page uses the canonical installed/source reconstruction interface. Paper
reproduction scripts are listed separately because they require a repository
checkout and are not installed by the Wheel.

Prepare inputs
--------------

For ``--sample-name sample``, ``--st-file st.h5ad``, and
``--sc-ref-file sc_ref.h5ad``, use the flat application layout resolved by the
CLI:

.. code-block:: text

   data/
   |-- sample_st.h5ad
   `-- sc_ref.h5ad

The two resolved paths are ``data/sample_st.h5ad`` and ``data/sc_ref.h5ad``.
Both AnnData inputs require a non-empty ``X``, unique ``obs_names`` and unique
``var_names``. The ST input also requires two spatial coordinate columns in
``obsm["spatial"]``. With the default column arguments, hST requires ``Level1``
in the reference ``obs``; iST and sST require both ``Level1`` and ``Level2``.
The files must share at least one gene.

If the reference has the default ``Patient`` column, its values are compared
to ``--sample-name`` after string normalization and only matching rows are
used. Select another column with ``--patient-key``; if the selected column is
absent, no patient-row filter is applied.

Installed command
-----------------

After installing an exact candidate Wheel, run hST with POT:

.. code-block:: bash

   revise-reconstruct \
     --platform hST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --ot-method pot

Use ``--platform iST`` or ``--platform sST`` for the other two declared
routes. iST also accepts ``--select-ct``, ``--cell-type-col``,
``--sub-cell-type-col``, and ``--ist-mapping mean|random``.

Source compatibility wrapper
----------------------------

From the repository root, the same request can be sent through the source
compatibility wrapper:

.. code-block:: bash

   python reconstruct.py \
     --platform hST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --ot-method pot

``reconstruct.py`` delegates to the same ``revise.cli`` implementation as the
installed command, so route and public output rules are shared.

Preflight
---------

Append ``--dry-run`` to either interface to resolve the route, open H5AD inputs
in backed mode, validate required axes/fields/gene overlap, and check the chosen
solver dependency without running reconstruction:

.. code-block:: bash

   revise-reconstruct \
     --platform hST \
     --sample-name sample \
     --data-root data \
     --st-file st.h5ad \
     --sc-ref-file sc_ref.h5ad \
     --output-root output \
     --ot-method pot \
     --dry-run

Preflight intentionally does not scan every expression value and does not run
preprocessing, GA, LR, finalization, or evaluation. A ready result is therefore
input/dependency evidence, not proof of a completed scientific run.

Select POT or TACCO
-------------------

``--ot-method pot`` sets both GA and LR to POT. After installing the optional
TACCO extra, ``--ot-method tacco`` sets both to TACCO. Omit the option to retain
the merged configuration, including an intentionally mixed GA/LR selection.
TACCO failure is fail-closed and does not fall back to POT.

Output
------

The public result is:

.. code-block:: text

   <output-root>/<sample-name>/<platform>-SVC.h5ad

For the three routes this becomes ``hST-SVC.h5ad``, ``iST-SVC.h5ad``, or
``sST-SVC.h5ad``. The file's ``uns["revise_reconstruction"]`` links to the
canonical run's ``provenance.json``. Inspect that manifest before treating a
run as completed; directory existence alone is not success evidence.

Paper reproduction compatibility
--------------------------------

The following paths exist only in a source checkout:

.. code-block:: bash

   python benchmark_main.py --help
   bash benchmark_main.sh

``benchmark_main.py`` runs one confounding family, which may contain several
leaf runs. The shell launcher delegates multiple families to the bounded,
foreground Python coordinator and returns nonzero when a child case fails.
Existing application reconstruction scripts support historical notebook
layouts, but they are compatibility paths rather than the installed quick
start. Their external data identity and real-data results are pending owner
confirmation.
