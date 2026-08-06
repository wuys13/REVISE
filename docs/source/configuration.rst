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
directory. All child input and output values must be non-empty
relative paths beneath that root: they cannot be absolute, use ``~``, or escape
it with ``..``.

Application has one input mode: exact direct paths. ST and reference resolve as
``<root-dir>/<st-path>`` and ``<root-dir>/<reference-path>``. An optional
sc-SVC-sr prior is declared as ``inputs.pm_on_cell.path``; omission means no
sidecar search and seeded random quota-slot assignment. Benchmark keeps its
legacy data-root semantics outside this Application schema.

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

Application does not filter the reference by patient. Prepare the reference
before running when it contains multiple cohorts.

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

The request runs by default. ``--dry-run`` (or ``dry_run=True`` in
``run_application``) changes only the effective action to ``preflight``.

Preflight may write run evidence, including ``preflight.json`` and
``provenance.json``, in the invocation's dedicated run directory. A dry run
does not write an H5AD.

Application identity is namespaced under ``application_config`` with
``source_path``, ``source_sha256``, ``declared_root``, ``resolved_root``,
``cwd``, ``resolved_inputs``, ``output_paths``, and ``effective_action``.
Top-level ``config_path`` and ``config_hash`` remain the
top-level engine configuration identity and are not overwritten by the
application YAML identity.

For an official bare template name, the CLI first checks
``configs/application/<name>.yaml`` in the launch directory and then reads the
package template. An existing explicit file always wins; arbitrary missing
paths are not guessed.

Output
------

Every Application invocation allocates exactly one run directory at
``<output.dir>/<output.name>/application__<svc-type>/<timestamp_uuid>/``.
Its run envelope contains ``provenance.json``, ``merged_config.json``, and
``run.log``; completed input validation adds ``preflight.json``. A successful
formal sp-SVC or sc-SVC-sr run writes
``<output.name>.h5ad`` directly into that directory. A formal standard sc-SVC
run writes ``<output.name>_spatial.h5ad`` and
``<output.name>_expr.h5ad`` as sibling files in the same directory. ``--dry-run``
gets an independent run directory with no H5AD, and each repeated formal run
gets a new ``<timestamp_uuid>`` leaf. Application does not create temporary
H5AD files, shared fixed-name copies, or use ``os.replace`` for H5AD output.
