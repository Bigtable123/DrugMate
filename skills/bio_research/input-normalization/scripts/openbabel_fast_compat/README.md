# Open Babel fast coordinate compatibility driver

This private executable is used only when the installed `obabel --gen3d fast`
fails the adapter's coordinate probe. It links against that same installation;
it does not replace its binaries or shared libraries.

Open Babel 3.1.1's `gen3d` operation builds and minimizes `molCopy`, but its
`fast` return exits before copying those coordinates back to the input molecule.
The compatibility operation preserves the upstream algorithm and adds the
missing `*pmol = molCopy` before the early successful returns. Its force-field
settings, distance-geometry fallback, rotor search during that fallback,
stereochemistry helper and retry count are otherwise unchanged. In particular,
the fallback can still take time and must run under the caller's molecule
timeout.

The normal `fast` branch uses the builder plus 100 conjugate-gradient steps,
with no rotor search. Its existing lack of a final stereochemistry recheck is
also unchanged; this fix does not claim to strengthen that scientific contract.

The operation and helper classes and operation ID have private SBDD names to
avoid colliding with the linked library. `driver.cpp` accepts only the adapter's
argument shape:

```text
-iFORMAT INPUT -osdf -O OUTPUT -h --gen3d fast
```

It uses `OBConversion::Convert` and general options, retaining native operation
ordering, hydrogen addition and multi-record parsing behavior. SMILES requests
normally provide one record per invocation. Whole-file calls for other formats
retain Open Babel's default stop-on-read-error behavior.

## Build

Compile all three source files, using the selected Open Babel installation's
headers and library and an Eigen include directory:

```bash
g++ -std=c++17 -O2 -DHAVE_EIGEN \
  -I"$OB_PREFIX/include/openbabel3" -I"$EIGEN_INCLUDE" \
  driver.cpp gen3d_fast.cpp gen3dstereohelper.cpp \
  -L"$OB_PREFIX/lib" -Wl,-rpath,"$OB_PREFIX/lib" -lopenbabel \
  -o "$OUTPUT"
```

`EIGEN_INCLUDE` contains `Eigen/Core`. `HAVE_EIGEN` is mandatory: omitting it
would silently remove the upstream distance-geometry fallback. Open Babel's
public `distgeom.h` also includes `LBFGS.h`, normally installed in the same
`include/openbabel3` directory. The parent adapter owns dependency discovery,
compilation timeout, runtime environment, binary provenance and the output
directory. No compiled artifact belongs in this source directory.

## Upstream source and license

Vendored from the official `openbabel-3-1-1` tag:

- [gen3d.cpp](https://github.com/openbabel/openbabel/blob/openbabel-3-1-1/src/ops/gen3d.cpp),
  saved as `gen3d_fast.cpp`; original SHA-256
  `acc40af6082f646482718b559f1b9575543290a671ac958155cb899a73e0bca8`.
- [gen3dstereohelper.h](https://github.com/openbabel/openbabel/blob/openbabel-3-1-1/src/stereo/gen3dstereohelper.h);
  original SHA-256
  `116544f9c34be94965b9e2131e17ae2dac67cc177d362472bfc6de4294fa26a4`.
- [gen3dstereohelper.cpp](https://github.com/openbabel/openbabel/blob/openbabel-3-1-1/src/stereo/gen3dstereohelper.cpp);
  original SHA-256
  `af3e3e8f93ee8c0729744950654d2931e75f6535c618fc2e89261786e86571b3`.

Original copyright and license notices remain in each file. The gen3d source
and driver use GPL version 2; the helper permits GPL version 2 or later. The
GPL version 2 text is included in `COPYING`. Local modifications on 2026-09-18
are the coordinate copy fix, private symbol names and the helper include path.
