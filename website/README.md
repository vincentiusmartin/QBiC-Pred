# QBiC-Pred (Quantitative predictions of transcription factor binding changes due to sequence variants)

## To generate TF mapping:
1. Go to `mapping_generator` folder.
2. Open `mapping_generator.py`. Then edit path in `gene2upbm` and `tfdb`.
3. Run `python3 mapping_generator.py`
4. This will create `mapping_data` folder that should contain:
   * **tflist.txt** contains all transcription factors separated by new line.
   * **pbmlist.txt** contains all pbm filenames used by QBiC.
   * **pbmtohugo.txt** contains mapping from pbm to hgnc name.
       Each line contains mapping in the format of `pbmname:tf1,tf2,...`.
   * **hugotopbm.txt** contains mapping from hgnc name to pbm.
       Each line contains mapping in the format of
       `familyname->tf1:pbm1,pbm2,...;tf2:pbm1,pbm2,...`

## Daemonize web server command

### Apache (used mainly to restart Flask)
1. systemctl reload apache2
2. systemctl restart apache2

### Celery (using systemctl)
Installation is taken from [here](https://thomassileo.name/blog/2012/08/20/how-to-keep-celery-running-with-supervisor/)

- Run supervisord: venv2/bin/supervisord
- Log: supervisorctl tail celeryd
- Restart: supervisorctl restart celeryd

### Daemonize redis
- Since we are using RediSearch, from RediSearch/build: redis-server --loadmodule ./redisearch.so --daemonize yes

## Roadmap (since Feb 2019)
- [x] advanced option
- [ ] accept input 15mer sequences
- [ ] put timestamp on recent jobs and enable deletion on them
