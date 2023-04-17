# ROCKE-3D Docker Image

[Docker](https://www.docker.com/) image to install and run a containerised ROCKE-3D on Fedora.


## Useful links

* [ROCKE-3D Webpage](https://simplex.giss.nasa.gov/gcm/ROCKE-3D/)
* [ROCKE-3D Compilers & libraries](https://docs.google.com/document/d/1-I8x1Op215f3m3NTtEo_cP2G-lP329pyEEUAzH6Xhog/view)
* [ROCKE-3D Installation](https://docs.google.com/document/d/1yyI0CDx1wEYbwqRsbvczXpdW2teePZ_NgIePTLFHtNA/edit)
* [ROCKE-3D Tutorial Videos](https://www.youtube.com/playlist?list=PLpMmnV3HS7r3KGXX8hmIBR3grXNu5hfW-)
* [Docker build help](https://docs.docker.com/engine/reference/commandline/build/)
* [Docker run help](https://docs.docker.com/engine/reference/commandline/run/)


## Installation & running via published image

* [Install Docker desktop](https://www.docker.com/get-started)
* Ensure Docker desktop is running
* Download published image:

```
docker pull woodwardsh/rocke3d:latest
```

* Run container, noting the mounting of local dir `huge_space` to container `/home/app/ModelE_Support/huge_space` for shared storage of model output:

```
docker run -it --rm --volume=${PWD}/huge_space:/home/app/ModelE_Support/huge_space -w /home/app woodwardsh/rocke3d:latest

# Options:
# -it       interactive && TTY (starts shell inside container)
# --rm      delete container on exit
# --volume  mount local directory inside container
# -w PATH   sets working directory inside container
```


## Installation & running via locally built image

* Clone repo & navigate inside:

```
git clone git@github.com:hannahwoodward/docker-rocke3d.git && cd docker-rocke3d
```

* Build image from Dockerfile (~15 min):

```
docker build -t rocke3d .

# Or, if debugging:

docker build -t rocke3d . --progress=plain --no-cache
```

* Run locally built container

```
docker run -it --rm --volume=${PWD}/huge_space:/home/app/ModelE_Support/huge_space -w /home/app rocke3d

# Options:
# -it       interactive && TTY (starts shell inside container)
# --rm      delete container on exit
# --volume  mount local directory inside container
# -w PATH   sets working directory inside container
```


## Testing

* Start container
* Run `sh test-earth.sh` (output written to `huge_space/E1oM20_Test`)
* Run `sh test-planet.sh` (uses SOCRATES; output written to `huge_space/P1SoM40_test`)


## Publishing image

```
docker login && docker tag rocke3d woodwardsh/rocke3d && docker push woodwardsh/rocke3d
```
