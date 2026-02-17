# Contianers {#containers}

MaCh3 containers are stored in https://github.com/mach3-software/MaCh3/pkgs/container/mach3

## Setting up
If you are on cluster then due to security reasons most likely you will have singularity or apptainer, in VERY rare cases docker. If you want to setup on your private machine then use the following.

### Linux Users
Install docker or podman via your distribution package manager, instructions below: YOU ONLY NEED ONE OF THE THREE
docker: [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)

OR docker-desktop: [https://docs.docker.com/desktop/install/linux-install/](https://docs.docker.com/desktop/install/linux-install/)

OR podman: [https://podman.io/docs/installation](https://podman.io/docs/installation)

 docker-desktop might be slightly easier to use for those not as comfortable sysadmining linux.

### Mac Users:
 Install orbstack or docker-desktop. I highly recommend orbstack over docker-desktop:  [https://orbstack.dev/download](https://orbstack.dev/download)

Select the right version for your mac architecture 
### Windows Users:
Install docker-desktop: [https://docs.docker.com/desktop/install/windows-install/](https://docs.docker.com/desktop/install/windows-install/)


# How to use

## Docker
To pull: `docker pull ghcr.io/mach3-software/mach3:alma9latest`<br>
now that you have pulled it use: `docker run -it --name mach3-container ghcr.io/mach3-software/mach3:alma9latest`<br>
Congratulations you just entered MaCh3 container<br>

If encountering any access issues please use:
```
export CR_PAT=<TOKEN>
echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin
```

## Singularity
To pull: `singularity pull docker://ghcr.io/mach3-software/mach3:alma9latest`<br>
now that you pulled it use: `singularity shell mach3_alma9latest.sif`<br>
Congratulations you just entered MaCh3 container<br>
now 
```
cd ${MACH3_INSTALL_DIR}
source bin/setup.MaCh3.sh
```
Your MaCh3 is setup and ready.

If encountering any access issues please use:
```
export SINGULARITY_DOCKER_USERNAME=<USERNAME>
export SINGULARITY_DOCKER_PASSWORD=<TOKEN>
```

Alternatively, you can treat the container as a black box and execute exe inside the container like this.
```
singularity exec <my_container.sif> do_something.exe /path/to/input <arguments>
```


## Words of caution
They are ephemeral: When you exit the main container process, everything is thrown away.
Don't expect changes you make to the container image to persist in different containers from the same image.


## Build experiment specific MaCh3 using Core as base
Docker allow 

#Use MaCh3 container as a base
FROM ghcr.io/mach3-software/mach3:alma9latest AS mach3t2k_build
```
// Clone Your MaCh3
RUN git XXXXXXX.YourMaCh3.git
// cmake
RUN cmake ../
// compile
RUN make && make install
```

For more see https://github.com/KSkwarczynski/OA_Docker/blob/main/MaCh3T2K/Dockerfile

## GPU
TODO!!!

## Acknowledgments
Highly inspired by lessons from Luke Pickering