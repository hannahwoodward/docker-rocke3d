# ROCKE-3D Docker Image

[Docker](https://www.docker.com/) image to install and run a containerised ROCKE-3D on Fedora.


## Useful links

* [ROCKE-3D Compilers & libraries](https://docs.google.com/document/d/1-I8x1Op215f3m3NTtEo_cP2G-lP329pyEEUAzH6Xhog/view)
* [ROCKE-3D Installation](https://docs.google.com/document/d/1yyI0CDx1wEYbwqRsbvczXpdW2teePZ_NgIePTLFHtNA/edit)
* [Docker build help](https://docs.docker.com/engine/reference/commandline/build/)
* [Docker run help](https://docs.docker.com/engine/reference/commandline/run/)


## Installation

* [Install Docker desktop](https://www.docker.com/get-started)
* Ensure Docker desktop is running
* Download published image:

```
docker pull woodwardsh/rocke3d:latest
```

## Building image from scratch

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

## Publishing image

```
docker login && docker tag rocke3d woodwardsh/rocke3d && docker push woodwardsh/rocke3d
```

## TODO

[] Copy tests from ROCKE-3D installation doc
[] Add option of mounting external volume for output netcdfs
[] Finish README with running image & use
