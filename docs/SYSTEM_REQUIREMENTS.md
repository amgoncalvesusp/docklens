# System requirements

DockLens is a single-process desktop application. It does not use a GPU, does
not spawn worker processes for detection, and holds the whole result set in
memory. The footprint figures below were measured on one host; the minimum
and recommended requirements further down are extrapolations from them and
were not tested on smaller machines.

## Measured footprint

Workload: the 150-pose 2M5D screening set of the BSB 2026 paper (3,616 atoms
per MOL2 file), DS-calibrated profile, Discovery Studio-like view. Host:
16-thread AMD desktop, 64 GB RAM, Windows 11, Python 3.11, DockLens 1.0.0.
Memory is the process peak working set from the Win32 counters
(`K32GetProcessMemoryInfo`), so it includes the interpreter, NumPy/pandas,
Matplotlib and the Qt window.

| Stage | Wall time | Peak working set |
|---|---|---|
| Interpreter + imports | 0.5 s | 62 MB |
| Detection, 150 poses | 15.7 s | 94 MB |
| DS-like view, fingerprints, similarity, 91 pose families | 3.3 s | 94 MB |
| Interface with all four workspaces rendered | ~32 s | 272 MB |

The benchmark script is available from the authors on request.

Detection cost grows with the number of poses and with receptor size; the
fingerprint and similarity stages are quadratic in the number of observations,
which is why analyses above 300 observations switch to a disclosed evenly
spaced training sample (see the README).

## Bare minimum

* 64-bit Windows 10/11 or Linux x86-64.
* 2 GB RAM. The measured peak with the full interface open is 272 MB; 2 GB
  leaves room for the operating system and a browser.
* 200 MB free disk (extrapolated; check the size of the release you install).
* No GPU, no internet connection, no Python installation for the prebuilt
  executables.

## Recommended

* 4 GB RAM or more, so that several projects and a molecular viewer can stay
  open alongside DockLens.
* A display large enough to keep the navigation rail, the analysis field and
  the DockLens Lens inspector visible together. No minimum resolution was
  tested.
* An SSD; project files (`.docklens`) and XLSX exports for large campaigns
  are written atomically and can reach tens of megabytes.

## Running from source

Python 3.9 or newer with the pinned versions in `requirements.txt` (NumPy,
pandas, openpyxl, Matplotlib, PyQt5). No RDKit, OpenBabel or PyMOL is required.
On Linux the Qt platform plugin needs the usual X11/Wayland client libraries
(`libxcb`, `libGL`); for headless use (tests, batch detection) set
`QT_QPA_PLATFORM=offscreen`.

## Containers

No official Docker image is provided. DockLens is a Qt desktop program, and a
container would still require an X or Wayland server on the host, which
reintroduces the installation step the single-file executable removes. The
detection core and the test suite run headless and can be containerized for a
fixed build environment:

```
FROM python:3.11-slim
RUN apt-get update && apt-get install -y libgl1 libxcb1 libxkbcommon0 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt requirements-dev.txt ./
RUN pip install -r requirements.txt -r requirements-dev.txt
COPY . .
ENV QT_QPA_PLATFORM=offscreen
CMD ["pytest", "-q"]
```

This snippet is untested and is offered only as a starting point; it is not
part of the release pipeline.
