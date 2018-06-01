Predicting w/ PathoPredictor
==============================

Use docker container [samesense/pathopredictor](https://hub.docker.com/r/samesense/pathopredictor/) to make predictions. The container includes all software, but annotation datasets need to be downloaded. You need to specify a data directory mount when using docker. Your input and output files will live here.

```
docker run -e USER=$USER -e USERID=$UI --user $(id -u) -it -v /local/data/:/opt/pathopredictor/data/ samesense/pathopredictor /bin/bash
```

Put your vcf file(s) in `/local/data/interim/user_preds/`. Usinf the `interim/user_preds/` directory is required, but only use the path to data when envoking docker `-v /local/data/`.

Edit `/opt/pathopredictor/configs/sm_predict_ex_config.json` so that it contains the names (w/ no vcf suffix) of your vcf files.

While in the docker container, navigate to `/opt/pathopredictor/src/rules/`. The prediction pipeline is run with Snakemake:

```
/opt/conda/envs/pathopredictor/bin/snakemake --configfile ../../configs/sm_predict_ex_config.json -s sf_predict.py all_predictions
```

To see the pipeline rules that will be run:

```
/opt/conda/envs/pathopredictor/bin/snakemake --configfile ../../configs/sm_predict_ex_config.json -s sf_predict.py -n all_predictions
```
