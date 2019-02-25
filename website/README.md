# QBiC-Pred (Quantitative predictions of transcription factor binding changes due to sequence variants)

![alt text](http://qbic.gcb.duke.edu/static/images/headerlogo.png)

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
- Reload: `systemctl reload apache2`
- Restart: `systemctl restart apache2`
- Debugging: `tail /var/log/apache2/error.log`

### Celery (using systemctl)
Installation is taken from [here](https://thomassileo.name/blog/2012/08/20/how-to-keep-celery-running-with-supervisor/)

- Run `supervisord: venv2/bin/supervisord`
- Log: `venv2/bin/supervisorctl tail celeryd`
- Restart: `venv2/bin/supervisorctl restart celeryd`

### Daemonize redis
- Since we are using RediSearch, from RediSearch/build: `redis-server --loadmodule ./redisearch.so --daemonize yes`

## Roadmap (since Feb 2019)
- [x] advanced option for E-score
- [x] accept input 17mer sequences + mutated
- [x] fix bug on the upload
- [x] download, only the filtered
- [ ] loading message/icon when downloaded file is being prepared
- [ ] job name to be the name of the downloaded prediction file
- [ ] acknowledgment and ownership section
- [ ] handle expiry after items are deleted in Redis
- [ ] put timestamp on recent jobs and enable deletion on them
- [ ] email notification
- [ ] remove sharedlist variable from predict if truly unecessary

## Things to discuss
- about the moving image, make something similar to the lab website

## Finished (2018) Roadmap
- [x] loading message
- [x] save user data in redis instead of cookie
- [x] parse chromosome another versions
- [x] delete restbl in /tmp after expired time
- [x] integrate e-score calculation
- [x] search in result DataTable
- [x] input verification
- [x] show user input configuration
- [x] show gene in the output
- [x] family - gene field must match
- [x] things to output -- BEST EFFORT
- [x] match redis session expired date and user recent job list -- FINAGLED W/ CELERY (still need work tho)
