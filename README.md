# TEstrainer
A pipeline to remove multicopy genes from repeat libraries

## Script descriptions
- initial_domain_finder.R - Used to identify common domains found in each family of TE

  - Relies on Repbase and rpstblastn

- strainer.R

  - Removes repeats containing domains yet lacking appropriate domains from repeat library.

    - Plan to add masking step to remove in appropriate portions of chimeric repeats

- TEtrim.py - Used to trim repeats, intended as repeat-specific alternative to trimal

- test.py - script used for development of other python scripts/learning python
