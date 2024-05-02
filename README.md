# pbta-ancestry

### To reproduce the code in this repository

1. Clone the repository:
```
git clone git@github.com:d3b-center/pbta-ancestry.git
```

2. Pull Docker container:
```
docker pull pgc-images.sbgenomics.com/d3b-bixu/pbta-ancestry:latest
```

3. Start the Docker container; from the `pbta-ancestry` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-ancestry pgc-images.sbgenomics.com/d3b-bixu/pbta-ancestry:latest
```

4. Execute the shell within the docker image; from the `pbta-ancestry` folder, run: 
```
docker exec -ti <CONTAINER_NAME> bash
```