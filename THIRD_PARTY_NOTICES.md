Third-party notices

This repository contains code derived from or based on MOFChecker.

MOFChecker
----------

- Original project: MOFChecker
- Original author: Kevin Jablonka
- Original license: MIT License
- [Source](https://github.com/lamalab-org/mofchecker)

Copyright (c) 2021 Kevin Jablonka

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

PLATON
------

PLATON is not distributed with this repository.

This repository may provide functionality for calling PLATON
when PLATON has been installed separately by the user.
Users are responsible for obtaining PLATON from the official PLATON
website and for complying with the applicable PLATON license terms.

For PLATON licensing and downloads, please refer to the official PLATON
website:

[https://www.platonsoft.nl/platon/pl030000.html](https://www.platonsoft.nl/platon/pl030000.html)

According to the PLATON distribution terms, PLATON is available free of
charge for academic use. For-profit organisations should contact the
PLATON author/licensor before copying or using the programs, as described
on the official PLATON website.

References: 

- A.L.Spek (2003) J. Appl. Cryst. 36, 7-13 & A.L.Spek (2009) Acta Cryst. D65, 148-155.

BLAKE2 (vendored in `libcif/external/`)
---------------------------------------

`libcif` includes the BLAKE2 reference C implementation and SSE-optimized
variants in `libcif/external/ref/` and `libcif/external/sse/`. The source files
retain their original license headers.

- Original project: BLAKE2 reference source code package
- Original author: Samuel Neves <sneves@dei.uc.pt>
- Original license: tri-licensed under CC0 1.0 Universal, OpenSSL License, or
  Apache License 2.0 (at the user's option)
- [Source](https://github.com/BLAKE2/BLAKE2)
- [Project page](https://blake2.net)

Copyright 2012, Samuel Neves <sneves@dei.uc.pt>.

License terms:

- CC0 1.0 Universal: https://creativecommons.org/publicdomain/zero/1.0
- OpenSSL License: https://www.openssl.org/source/license.html
- Apache License 2.0: https://www.apache.org/licenses/LICENSE-2.0

RDKit (optional libcif dependency)
----------------------------------

`libcif` can be built with optional RDKit support (`LIBCIF_RDK`) for
canonical-order computations and SMILES handling. RDKit is not distributed
with this repository.

- Original project: RDKit
- Original license: BSD-3-Clause
- [Source](https://github.com/rdkit/rdkit)

OpenBabel (optional libcif dependency)
--------------------------------------

`libcif` can be built with optional OpenBabel support (`LIBCIF_OBABEL`) for
xyz-to-SMILES conversion. OpenBabel is not distributed with this repository.

- Original project: Open Babel
- Original license: GNU General Public License v2 (GPL-2.0)
- [Source](https://github.com/openbabel/openbabel)

Spglib (optional libcif dependency)
-----------------------------------

`libcif` can be built with optional Spglib support (`LIBCIF_SPGLIB`) for
symmetry and space-group operations. Spglib is not distributed with this
repository.

- Original project: Spglib
- Original license: BSD-3-Clause
- [Source](https://github.com/spglib/spglib)

Boost (optional libcif dependency)
----------------------------------

When `libcif` is built with RDKit support (`LIBCIF_RDK`), Boost components
(`timer`, `system`, `serialization`, `iostreams`) are also required. Boost is
not distributed with this repository.

- Original project: Boost C++ Libraries
- Original license: Boost Software License 1.0 (BSL-1.0)
- [Source](https://www.boost.org)

QMOF database (test fixtures)
-----------------------------

The CIF files in `test/test_cifs/` and the reference outputs in
`test/fixtures/smoke/` are derived from structures in the QMOF database.

- Original dataset: QMOF — Quantum Metal-Organic Framework Database
- Original authors: Andrew S. Rosen et al.  
- Original license: Creative Commons Attribution 4.0 International (CC-BY-4.0)  
- [Source](https://github.com/Andrew-S-Rosen/QMOF)

References:

- A.S. Rosen, S.M. Iyer, D. Ray, Z. Yao, A. Aspuru-Guzik, L. Gagliardi, J.M. Notestein, R.Q. Snurr. "Machine Learning the Quantum-Chemical Properties of Metal–Organic Frameworks for Accelerated Materials Discovery", Matter, 4, 1578-1597 (2021). DOI: 10.1016/j.matt.2021.02.015.
- A.S. Rosen, V. Fung, P. Huck, C.T. O'Donnell, M.K. Horton, D.G. Truhlar, K.A. Persson, J.M. Notestein, R.Q. Snurr. "High-Throughput Predictions of Metal–Organic Framework Electronic Properties: Theoretical Challenges, Graph Neural Networks, and Data Exploration," npj Comput. Mat., 8, 112 (2022). DOI: 10.1038/s41524-022-00796-6.

The QMOF-derived files are redistributed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) and are used solely as
test fixtures. Any reuse of these files must preserve attribution to the original QMOF authors as required by CC-BY-4.0.
