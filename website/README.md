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

## Roadmap
- [x] roadmap1
- [x] roadmap2
- [x] roadmap3
- [ ] roadmap4
