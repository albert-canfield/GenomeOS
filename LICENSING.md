# Licensing

GenomeOS is open source and stays open source. Two licences are used, and
which one applies depends on which half of the project a file belongs to.
That split is not bureaucracy: it follows the architecture the project has
had from the start ([docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) §10), where
**BioLang** is an engine meant to spread and **GenomeOS** is the application
built on it.

| Part | Paths | Licence | Why |
|---|---|---|---|
| **BioLang engine** — the language, the intermediate representation, the virtual machine, the standard library and the toolchain | `genomeos/lang/`, `genomeos/ir/`, `genomeos/runtime/`, `genomeos/std/`, `genomeos/coords.py`, `genomeos/bio.py` | **Apache License 2.0** ([LICENSE-APACHE](LICENSE-APACHE)) | a language has to be embeddable. Anyone may put BioLang inside their own tool, open or closed, exactly as they would use Python or Node. Apache 2.0 also grants patent rights explicitly, which MIT does not. |
| **GenomeOS application** — genome decoding, molecules, twins, organisms, cancer, therapeutics, the web interface, the command line, the jobs and storage layers | everything else under `genomeos/`, plus `scripts/` and `tests/` | **GNU Affero General Public License v3.0 or later** ([LICENSE](LICENSE)) | this is the part that took the work. Anyone may use it, study it, change it and publish research with it. What they may not do is take it private: a modified version, **including one offered to others over a network**, must make its source available under the same licence. |
| **Documentation** — `docs/`, `README.md`, and the other Markdown files | | **Creative Commons Attribution 4.0** | quote it, teach from it, translate it; credit the project. |
| **Distilled results** — `data/results/` | | derived from their upstream sources, whose terms apply | each summary names its source; see [ACKNOWLEDGEMENTS.md](ACKNOWLEDGEMENTS.md). |
| **BioLang programs** — `data/demo/`, `data/organisms/` | | Apache License 2.0, as part of the engine's standard material | programs are meant to be copied and adapted. |

`SPDX-License-Identifier` headers on the package entry points state this in
machine-readable form. When BioLang is split into its own package (roadmap
milestone 2.0), the engine paths above move across unchanged and the split
becomes two repositories with these same two licences.

## What this means in practice

**If you are a researcher**, nothing here restricts you. Use GenomeOS, change
it, run it on your own data, publish what you find. Cite it if it helped.

**If you are building a tool on BioLang**, the engine is Apache 2.0. Embed it
and keep your own code under whatever licence you like.

**If you are running a modified GenomeOS as a service**, the AGPL asks one
thing: the people using that service must be able to get your modified
source. This is the "protection": improvements to GenomeOS come back to
GenomeOS rather than disappearing into a private product.

**If the AGPL does not suit your situation**, the copyright holder can grant a
separate commercial licence. GenomeOS is written by one author, so that is a
conversation, not a legal impossibility. Ask.

## The data is a separate question from the code

The licence above governs **this source code**. It does not and cannot grant
rights over the public databases and models GenomeOS talks to. Several of them
are free only for non-commercial or academic use, whatever the code's licence
says:

- **AlphaGenome** — the API and the model weights are non-commercial;
  commercial use goes through Google Cloud self-deployment.
- **NetMHCpan** — academic licence from DTU Health Tech; GenomeOS never ships
  or downloads it.
- **COSMIC identifiers**, and some **cBioPortal** studies — their own terms.

A commercial user must check every source they enable. GenomeOS names the
source and the evidence level of every fact it reports precisely so that check
is possible. Full list with licences:
[ACKNOWLEDGEMENTS.md](ACKNOWLEDGEMENTS.md).

## History

GenomeOS was released under the MIT licence on 2026-09-10 and relicensed on
2026-09-11 as described above. Commits published before the change remain
available under MIT: a licence change is not retroactive, and anyone who
obtained the earlier code keeps the rights it gave them. Everything from the
relicensing commit onward is Apache 2.0 or AGPL-3.0-or-later by the table
above.

## Contributing

Contributions are accepted under the licence of the file being changed
(inbound equals outbound). By contributing you also agree that the copyright
holder may offer your contribution under a separate commercial licence, which
is what keeps dual licensing possible. See
[CONTRIBUTING.md](CONTRIBUTING.md).

GenomeOS produces research hypotheses with their evidence and their gaps. It
is not a medical device and gives no clinical advice. Both licences disclaim
warranty and liability; that disclaimer is meant literally here.
