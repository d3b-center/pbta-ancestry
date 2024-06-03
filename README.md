# Genetic ancestry superpopulations show distinct prevalence and outcomes across pediatric central nervous system tumors 

Ryan J. Corbett, Cricket C. Gullickson, Zhuangzhuang Geng, Miguel A. Brown, Bo Zhang, Chuwei Zhong, Nicholas Van Kuren, Antonia Chroni, Christopher Blackden, Ammar S. Naqvi, Alexa Plisiewicz, Sean McHugh, Emmett Drake, Kaitlin Lehmann, Tom B. Davidson, Michael Prados, Phillip B. Storm, Adam C. Resnick, Angela J. Waanders, Sebastian M. Waszak, Sabine Mueller, Jo Lynne Rokita, Cassie Kline

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
docker run --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/pbta-ancestry pgc-images.sbgenomics.com/d3b-bixu/pbta-ancestry:latest
```
NOTE: if running on a MacBook with M1 chip, please incude the argument `--platform linux/amd64`

4. Execute the shell within the docker image; from the `pbta-ancestry` folder, run: 
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Run the download_data.sh shell script to obtain latest data files:

```
bash download_data.sh
```

### Contact

For questions, please submit an issue or send an email to Ryan Corbett ([@rjcorb](https://github.com/rjcorb)): corbettr@chop.edu