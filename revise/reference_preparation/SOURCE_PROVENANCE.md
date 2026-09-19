# Reference preparation source provenance

This package was integrated from the isolated portable bundle at
`ReviseSTAtlas/preparation/external/revise-reference-preparation` on 2026-09-19.
The Atlas checkout HEAD at integration time was
`fe4ac4a769a0518854e14be6426f1ad0bf327e95`, but the portable bundle was an
untracked directory in that checkout. The per-file hashes below, rather than
that checkout commit, identify the integration snapshot.

The portable bundle adapted Global Anchoring validation and screening from
ReviseSTAtlas while removing its registry, catalog, campaign, persistence and
reconstruction ownership. `scoring.py` preserves the numerical implementation
whose original Atlas source SHA-256 was
`ca295677cb57a49caadcd11709bc9fd56372defc654a845640868afc1c6c6d55`.
`config.py`, `paired.py`, `api.py`, and the reference-only YAML contract were
new in that portable bundle. This REVISE integration adds host-controlled GA
metadata, cross-candidate ST-axis checks, content digests, and verified input
resolution; it does not claim biological reference quality.

## Portable integration snapshot

| Portable path | SHA-256 |
| --- | --- |
| `reference_preparation/__init__.py` | `16ab734159ed46ab725aa4a606b0fdd6b802962aa7d9a075c8a1da1b8deaef5f` |
| `reference_preparation/api.py` | `fa307aa16c27542cbd0c806d87c4478561db584f0c500fb7dadedd469fa2fcc4` |
| `reference_preparation/config.py` | `36be85898a3a007862ac89b58605d293bb71d2ef477d43a0fda558178cfed1b7` |
| `reference_preparation/ga_contract.py` | `32e1b81d042d2c0201e2495c890baf44946417b6a18375b34b07f73ad0514081` |
| `reference_preparation/paired.py` | `630c38c719085fdefa91d0d83c113993f7bc469c99f85c5a616330840ec27498` |
| `reference_preparation/scoring.py` | `ca295677cb57a49caadcd11709bc9fd56372defc654a845640868afc1c6c6d55` |
| `reference_preparation/screening.py` | `4237866d33dbd6dbed4c16f8fadfae1ddd05baa892330d6347747d4546a8b4ab` |
| `tests/test_api.py` | `1be87d0d91f8711741ef4a8f4b8558608ff77ab7f22d550baace0feaa09199e2` |
| `tests/test_config.py` | `5504165d44d9b6fdd88dac2f00f9eb32ca516a2bdff95a4f03bb41ab903e3603` |
| `tests/test_ga_contract.py` | `89748b2d2a1ba76716fcb7d72aa50d8c918b621b937fb70ff76708af87405740` |
| `tests/test_paired.py` | `f4b5646095faff3b08d02ed4b581f78997a1321ebe5f91e954269449079241d7` |
| `tests/test_scoring.py` | `a9a0eee918c39110f5603ac04336f70b0e5129a70818b2ff07f81ad70cdc3461` |
| `tests/test_screening.py` | `65f27fbfada294353a88fece9e9829a33800cead838412f99577982552ca0ef2` |
| `examples/candidates.yaml` | `d19f6addcb7b6b3ce5ffd223deab7d5ef1a3b0b5f7cdab6e15fa83c6632d6ea2` |
| `examples/paired.yaml` | `042b4037db226c3819530f7116b7c70eebbccf272fd40b3bfa6c064ebbd633e6` |
| `examples/reference.yaml` | `7eda8912d10519af6d76cd446b95d153657f22179c052156af8076f99534c0a6` |
| `examples/screen-directory.yaml` | `3dd4fe5f163cdcbc1d69bb59e87fafd5f47e253e073a93d5fd70a65f58576225` |
| `examples/screen-list.yaml` | `4dc0a92da9111d73e41d22a1dc09526cc48ae51843ecde1458c3e01ab0245ef4` |

The six portable test modules are retained in `tests/reference_preparation/`
with only the import namespace changed from `reference_preparation` to
`revise.reference_preparation`. REVISE-specific contracts live in the separate
`test_revise_extensions.py` module.
