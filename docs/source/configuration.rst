Configuration
=============

Application users pass one strict YAML to ``python reconstruct.py --config``
or ``revise-reconstruct --config``. The full top-level schema is
``application``, ``paths``, ``algorithm``, ``inputs``,
``global_anchoring``, ``local_refinement``, ``output``, and ``execution`` plus
``schema_version``. Unknown, missing, misspelled, and route-inapplicable fields
fail before scientific computation.

The package-owned ``revise/revise.yaml`` remains the sole engine authority for
defaults, profiles, routing, expert settings, and locked values. It is not an
application field and is not passed to application ``--config``.

Application and Benchmark selection are deliberately separate. Application
YAML supplies ``application.svc_type``; Benchmark supplies its confounding
factor. Both reach the same engine execution Interface, whose router resolves
profile, task, SVC kind, and strategy. Application runtime/provenance contains
``application_route`` and no Benchmark ``confounding`` field.

Root and path contract
----------------------

``paths.root_dir: .`` means the launch current working directory, not the
application YAML directory. Otherwise it must be an existing absolute
directory. All child input, data-root, and output values must be non-empty
relative paths beneath that root: they cannot be absolute, use ``~``, resolve
to the root itself, or escape it with ``..``.

Input modes are mutually exclusive:

.. list-table::
   :header-rows: 1

   * - Mode
     - Required fields
     - Resolution
   * - ``direct``
     - ``inputs.st.path``, ``inputs.st.format``,
       ``inputs.reference.path``, ``inputs.reference.format``, and
       ``inputs.reference.patient_key``
     - ST ``<root-dir>/<st-path>``; reference
       ``<root-dir>/<reference-path>``
   * - ``legacy_layout``
     - ``inputs.data_root``, ``inputs.st.file``, ``inputs.st.format``,
       ``inputs.reference.file``, ``inputs.reference.format``, and
       ``inputs.reference.patient_key``
     - ST ``<root-dir>/<data-root>/<sample-name>_<st-file>``; reference
       ``<root-dir>/<data-root>/<reference-file>``

Both modes resolve ``output.path`` as ``<root-dir>/<output-path>``. Direct mode
clears the internal data-root/file locators and does not probe
``PM_on_cell.csv``. Legacy layout fixes the optional sc-SVC-sr prior at
``<data-root>/PM_on_cell.csv``; it has no separate application field.

Route-specific fields
---------------------

.. list-table::
   :header-rows: 1

   * - Type
     - Required route fields
     - Maintained template
   * - ``sc-SVC``
     - ``global_anchoring.broad_column``,
       ``local_refinement.subtype_column``, and one concrete
       ``local_refinement.select_cell_type``
     - ``Xenium_T.yaml``, ``Xenium_Fib.yaml``, or ``Xenium_Mono.yaml``; each
       carries a fixed ``select_cell_type`` value
   * - ``sp-SVC``
     - ``global_anchoring.broad_column``; optional non-negative finite
       ``local_refinement.strength``
     - ``VisiumHD.yaml`` with strength ``0.2``
   * - ``sc-SVC-sr``
     - ``global_anchoring.broad_column``; optional non-negative finite
       ``local_refinement.strength``
     - ``Visium.yaml`` with strength ``0.0``

``application.sample_name`` is also the reference-selection value when the
configured ``inputs.reference.patient_key`` column exists. If that column is
absent, the reference is not filtered; this is the maintained Visium-template
case.

GA and LR OT selection
----------------------

REVISE has two OT stages: ``ot.ga.solver`` controls Global Anchoring and
``ot.lr.solver`` controls Local Refinement. ``algorithm`` has only the optional
``ot_method`` field. ``algorithm.ot_method`` controls both GA and LR; omitting
it keeps the selected engine profile authoritative. The value is ``pot`` or
``tacco``; there is no public mixed-stage application override.

The standard sc-SVC profile selects TACCO 0.5.0. Missing, incompatible, or
failed TACCO stops preflight or execution with directed failure. REVISE does
not fall back to POT. Explicitly choosing ``algorithm.ot_method: pot`` is a
solver change, not an automatic fallback.

Assignment and local refinement
-------------------------------

Global Anchoring produces validated posterior ``Q`` and ``argmax(Q)`` labels.
The only public local-refinement option is a non-negative finite strength:

.. code-block:: yaml

   local_refinement:
     strength: 0.2

sp-SVC defaults to ``0.2``. sc-SVC-sr defaults to ``0.0``. Strength zero does
not skip Local Refinement; it only disables assignment-posterior strengthening
of the LR cost. sp-SVC conditions local OT with ``Q``; sc-SVC-sr projects
``Q`` to virtual cells first. Standard sc-SVC uses only ``argmax(Q)`` for
cohort routing and rejects strength. There are no policy, compatibility-mode,
reference-mode, or one-hot-fallback application settings.

The manifest records only ``route``, ``applied``, and ``strength`` under
``local_refinement``; mandatory sc-SVC-sr allocation is separate under
``sr_allocation``. It does not expose solver-event telemetry.

Action and evidence
-------------------

``execution.action`` is ``run`` or ``preflight``. The truth table is:

.. list-table::
   :header-rows: 1

   * - YAML action
     - ``--dry-run``
     - Effective action
   * - ``run``
     - absent
     - ``run``
   * - ``run``
     - present
     - ``preflight``
   * - ``preflight``
     - absent or present
     - ``preflight``

Preflight may write run evidence, including ``preflight.json`` and
``provenance.json``, but does not publish a result H5AD.

Application identity is namespaced under ``application_config`` with
``source_path``, ``source_sha256``, ``declared_root``, ``resolved_root``,
``cwd``, ``resolved_paths``, ``declared_action``, ``effective_action``, and
``dry_run_override``. Top-level ``config_path`` and ``config_hash`` remain the
top-level engine configuration identity and are not overwritten by the
application YAML identity.

For an official bare template name, the CLI first checks
``configs/application/<name>.yaml`` in the launch directory and then reads the
package template. An existing explicit file always wins; arbitrary missing
paths are not guessed.

Output
------

sp-SVC and sc-SVC-sr publish
``<output-root>/<sample-name>/SVC.h5ad``. Standard sc-SVC publishes both
``<output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_spatial.h5ad`` and
``<output-root>/<sample-name>/sc-SVC/<cell-type>/sc_SVC_expr.h5ad``.
