This package contains the R code implementing the computational analysis steps, in the form of Rmarkdown documents, that enable readers to reproduce all major figures and results reported in manuscript "Multi-omics reveals clinically relevant proliferative drive associated with mTOR-MYC-OXPHOS activity in chronic lymphocytic leukemia" by Lu and Cannizzaro et al.. In addition, the CLLPDestimate function in the R package can be used to compute CLL-PD from compatible gene expression data provided by users. Instructions can be found in the vignette of the package. 

* [User guide for CLLPDestimate function](https://huber-group-embl.github.io/mofaCLL/CLLPDestimate.html)

* [Walkthrough of the analyses described in our paper](https://huber-group-embl.github.io/mofaCLL/analysisProcedure.html)


## Use docker container for reproducing the analysis results

Excuting the R markdown files in this package requires a large collection of packages, which can be quite cumbersome and may reduce reproducibility when the packages update in the future. Therefore we created a docker container for the system environment that is necessary to run all the analysis steps in the R markdown files. The docker container can be found at Docker Hub: https://hub.docker.com/r/lujunyan1118/mofacll 

### How to use

1. Install docker on your local computer.

2. Run the command below in a terminal.

```
docker pull lujunyan1118/mofacll

docker run -d -p 8787:8787 -e DISABLE_AUTH=true lujunyan1118/mofacll
```

3. Navigate to http://localhost:8787 in a browser and an R Studio Server interface should appear. 

4. In the current directory, you should see a folder called "vignettes_standalone", which contains the R Markdown files for reproducing the analysis results.

**Note:** the last part of the analysis (part6.Rmd) involves a single-cell CyTOF dataset, which requires a relatively large amount of memory. If you experience crashing of the R environment, please consider increasing the memory usage limit of docker.


### More Help

Please see https://hub.docker.com/r/rocker/rstudio and https://github.com/rocker-org/rocker/wiki for the additional documentation of the rocker/rstudio container. 


## Copyright (C) 2021  Junyan Lu

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Email: junyan.lu@embl.de
