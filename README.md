# ROCKE-3D Docker Image

[Docker](https://www.docker.com/) image to install and run a containerised ROCKE-3D on Fedora.


## Useful links

* [ROCKE-3D Compilers & libraries](https://docs.google.com/document/d/1-I8x1Op215f3m3NTtEo_cP2G-lP329pyEEUAzH6Xhog/view)
* [ROCKE-3D Installation](https://docs.google.com/document/d/1yyI0CDx1wEYbwqRsbvczXpdW2teePZ_NgIePTLFHtNA/edit)
* [Docker build help](https://docs.docker.com/engine/reference/commandline/build/)
* [Docker run help](https://docs.docker.com/engine/reference/commandline/run/)


## Building & running image from scratch

* Clone repo & navigate inside:

```
git clone git@github.com:hannahwoodward/docker-rocke3d.git && cd docker-rocke3d
```

* Build image from Dockerfile:

```
docker build -t rocke3d .
```
