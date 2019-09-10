# ctdna-nextflow

## How to run

#### 1. Requirements

These programs must be installed:

- Docker
- Conda (either miniconda or anaconda)
- Git  

These commands should all produce some output:

```
docker run hello-world
conda -h
git help
```

#### 2. Clone this repository

```
git clone https://github.com/erikwaskiewicz/ctdna-project.git
cd ctdna-project
```

#### 3. Make the nextflow conda environment

```
conda env create -f env/environment.yml
conda activate ctdna-project
```

#### 4. Optional: pull the Docker images 

If you do not do this step, the Docker images will be pulled automatically during the pipeline run

```
docker pull erikwaskiewicz/ctdna-base:latest
etc....
```

#### 5. Test the pipeline

Run the pipeline with the test data

`nextflow run ctdna-snv-caller-comparison.nf -profile test`

Run the unit test script

`pytest --verbose`

