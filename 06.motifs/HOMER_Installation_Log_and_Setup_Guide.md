# HOMER Installation Log and Setup Guide

This document records the steps used to install **HOMER** and the **mm10** / **hg38** genome packages in a WSL/Linux environment.

---

## Environment

- OS: WSL/Linux
- Working project directory: `~/task4/task5`
- HOMER installation directory: `~/task4/task5/homer`

---

## 1) Create the HOMER installation directory

```bash
mkdir -p ~/task4/task5/homer
cd ~/task4/task5/homer
```

## 2) Download the HOMER installer script

```bash
curl -O http://homer.ucsd.edu/homer/configureHomer.pl
```

At this point, `configureHomer.pl` is downloaded into the `homer` directory.

## 3) Check whether HOMER is already installed

We checked whether the HOMER executable was already available:

```bash
which findMotifsGenome.pl
```

Initially, this returned nothing, meaning HOMER was not yet installed or not in `PATH`.

## 4) Start HOMER installation

```bash
perl configureHomer.pl -install
```

### Potential issue; if the above command works on first trial then you do not need the following updates on the tools

The environment require system tools such as:

- `gcc`
- `g++`
- `make`
- `zip`
- `unzip`

## Install required system dependencies

Since HOMER requires standard build and archive utilities, we installed them with:

```bash
sudo apt update
sudo apt install -y build-essential zip unzip
```

This provides:

- `gcc`
- `g++`
- `make`
- `zip`
- `unzip`

## 5) Run HOMER installation --- Same as 4

After installing dependencies, we re-ran:

```bash
perl configureHomer.pl -install
```

This successfully downloaded and installed the HOMER software package.
During compilation, several warning messages appeared, but the installer completed successfully.

## 6) Install the mouse genome package (`mm10`)

After HOMER itself was installed, we installed the mouse genome and annotation package:

```bash
perl configureHomer.pl -install mm10
```

This installed:

- Mouse organism package
- `mm10` genome package

The installer reported:

- `Finished Installing mouse`
- `Finished Installing mm10`

## 7) Install the human genome package (`hg38`)

We also need `hg38` for the human-specific vs conserved analyses on native human coordinates:

```bash
perl configureHomer.pl -install hg38
```

This installed:

- Human organism package
- `hg38` genome package

The installer should report:

- `Finished Installing human`
- `Finished Installing hg38`

## 8) Add HOMER to the shell `PATH`

To make HOMER commands available from anywhere, we added the HOMER `bin` directory to `PATH`:

```bash
echo 'export PATH="$PATH: ~/task4/task5/homer/bin"' >> ~/.bashrc # you can replace to your own dir
source ~/.bashrc
```

## 9) Verify the installation

### Check that the executable is now found

```bash
which findMotifsGenome.pl
```

Expected output:

```text
~/findMotifsGenome.pl
```

### Check that `mm10` is installed

```bash
perl ~/task4/task5/homer/configureHomer.pl -list | grep mm10 # you can replace to your own dir
```

Expected output includes something like:

```text
+       mm10    v7.0    mouse genome and annotation for UCSC mm10
```

### Check that `hg38` is installed

```bash
perl ~/task4/task5/homer/configureHomer.pl -list | grep hg38 # you can replace to your own dir
```

Expected output includes something like:

```text
+       hg38    v7.0    human genome and annotation for UCSC hg38
```

---

## Notes

- If installation appears broken, re-run:

```bash
perl configureHomer.pl -install homer
```

- Official references:
  - [HOMER motif documentation](http://homer.ucsd.edu/homer/motif/)
  - [HOMER download page](http://homer.ucsd.edu/homer/download.html)
  - [HOMER2 page](http://homer.ucsd.edu/homer/homer2.html)
