# ROCKE-3D Docker Image

[Docker](https://www.docker.com/)/[Podman](https://podman.io/) image to install and run a containerised ROCKE-3D on Fedora.

## Useful links

* [ROCKE-3D Model Description (doi:10.3847/1538-4365/aa7a06)](https://iopscience.iop.org/article/10.3847/1538-4365/aa7a06/meta#apjsaa7a06s3)
* [ROCKE-3D Webpage](https://simplex.giss.nasa.gov/gcm/ROCKE-3D/)
* [ROCKE-3D Compilers & libraries](https://docs.google.com/document/d/1-I8x1Op215f3m3NTtEo_cP2G-lP329pyEEUAzH6Xhog/view)
* [ROCKE-3D Installation](https://docs.google.com/document/d/1yyI0CDx1wEYbwqRsbvczXpdW2teePZ_NgIePTLFHtNA/edit)
* [ROCKE-3D Tutorial Videos](https://www.youtube.com/playlist?list=PLpMmnV3HS7r3KGXX8hmIBR3grXNu5hfW-)
* [ROCKE-3D Publication Supplements (inc rundecks)](https://portal.nccs.nasa.gov/GISS_modelE/ROCKE-3D/publication-supplements/)
* [Docker build help](https://docs.docker.com/engine/reference/commandline/build/)
* [Docker run help](https://docs.docker.com/engine/reference/commandline/run/)


## Installation & running via published image

* [Install Docker desktop](https://www.docker.com/get-started)
* Ensure Docker desktop is running
* Download published image:

```
docker pull woodwardsh/rocke3d:latest
```

* Run container, noting the mounting of local dir `./ModelE_Support` to container `/home/app/ModelE_Support` for shared storage of model output:

```
docker run -it --rm --volume=${PWD}:/home/app/ModelE_Support woodwardsh/rocke3d:latest

# Options:
# -it       interactive && TTY (starts shell inside container)
# --rm      delete container on exit
# --volume  mount local directory inside container
# -w PATH   sets working directory inside container
```

### Podman

* Replace `docker` with `podman`, and note additional options to fix permissions on mounted volumes (see [podman run](https://docs.podman.io/en/latest/markdown/podman-run.1.html)):

```
podman run -it --rm -v ${PWD}/ModelE_Support:/home/app/ModelE_Support --security-opt label=disable woodwardsh/rocke3d:latest
```


## Installation & running via locally built image

* Clone repo & navigate inside:

```
git clone git@github.com:hannahwoodward/docker-rocke3d.git && cd docker-rocke3d
```

* Build image from Dockerfile (~15 min):

```
docker build -t rocke3d .
```

* Or, if debugging:

```
docker build -t rocke3d . --progress=plain --no-cache
```

* Run locally built container:

```
docker run -it --rm -v ${PWD}/ModelE_Support:/home/app/ModelE_Support rocke3d

# Options:
# -it       interactive && TTY (starts shell inside container)
# --rm      delete container on exit
# -v        mount local directory inside container
# -w PATH   sets working directory inside container
```

### Podman

* Build with similar command, replacing `docker` with `podman`:

```
podman build -t rocke3d .
```

* Run, with additional options to fix permissions on mounted volumes (see [podman run](https://docs.podman.io/en/latest/markdown/podman-run.1.html)):

```
podman run -it --rm -v ${PWD}/ModelE_Support:/home/app/ModelE_Support --security-opt label=disable rocke3d
```

## Testing

* Start container
* Run `sh test-earth.sh` (output written to `ModelE_Support/huge_space/E1oM20_Test`)
* Run `sh test-planet.sh` (uses SOCRATES; output written to `ModelE_Support/huge_space/P1SoM40_Test`)


## Diagnostics & Post-processing

* [ROCKE-3D Diagnostics info](https://simplex.giss.nasa.gov/gcm/doc/UserGuide/diagnostics.html)
* `$MODELDIR/model/mk_diags` directory has been added to PATH, which contains scripts to generate readable netcdf model outputs
  * `scaleacc` for interim/accumulative source files (e.g. `PARTIAL.acc$RUN_ID.nc aij`)
  * `sumfiles` to combine multiple acc files across different time periods
  * Documention can be found in `$MODELDIR/model/mk_diags/conventions.txt`


## Publishing image

```
docker login && docker tag rocke3d woodwardsh/rocke3d && docker push woodwardsh/rocke3d
```


## Model info

* Create a rundeck P2{G,S}{A,x,N}{p,q,o}{F,M}40
  * GISS/SOCRATES radiation
  * Atmosphere of PI Earth
  Ocean {p,q,o}
  * M40 4degx5deg with 40 layers in atmosphere, 13 layers in ocean
* Rundeck start/stop times:
  * YEAR is just an index, MONTH always 1 to 12, HOUR always 0 to 23
    * Calendar system used divides days into 24 model "hours" and years into 12 "months", so with varying orbital/rotation periods will therefore not generally be 3600s or 720hrs, respectively
  * `ISTART=2` way to tell the model that initial conditions will be provided (AIC atmosphere, ground GIC, and if relevant ocean OIC)
  * `IRANDI=X` random number generation seed/adding numerical noise, used for cloud generation
  * `master_yr=1850` tells model to use Earth greenhouse gas concentrations
  * `master_yr=0` runs a transient simulation
* For non-Earth continental configurations:
  * Ensure land mass added at south pole to prevent grid singularity issues
  * Change `OCNDYN.f#L5812`: `J10=4` to `J10=1` (or whatever latitude south pole landmass ends at)
