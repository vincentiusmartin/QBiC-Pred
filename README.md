# QBiC-Pred (Quantitative predictions of transcription factor binding changes due to sequence variants)

![alt text](http://qbic.gcb.duke.edu/static/images/headerlogo.png)

QBiC-Pred (http://qbic.gcb.duke.edu/) is a web server that uses OLS 6-mer models trained on universal PBM data to predict TF binding changes due to sequence variants, as well as the significance of the changes (which implicitly depends on the quality of the data and models).

## Directory structure
1. **generate_prediction**
This directory contains code to generate 12-mer table used by the website to make prediction
2. **website**
This directory contains all code files used by QBiC-Pred web server. The web server is written using Python Flask for the backend and bootstrap+Ajax for the frontend.
