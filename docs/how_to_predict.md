Predicting w/ PathoPredictor
==============================

Use docker container [samesense/pathopredictor](https://hub.docker.com/r/samesense/pathopredictor/) to make predictions. The container includes all software, but annotation datasets need to be downloaded.

```
docker run -e USER=$USER -e NB_UID=$UID -e USERID=$UI --user $(id -u) -it -v /mnt/isilon/:/mnt/
```
