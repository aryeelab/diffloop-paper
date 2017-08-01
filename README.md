<br><br>
# diffloop: a computational framework for identifying and analyzing differential DNA loops from sequencing data

Online resources and supplemental information

### Authors
[Caleb Lareau](mailto:caleblareau@g.harvard.edu) and [Martin Aryee](https://aryee.mgh.harvard.edu/)

### Links
- [Code repository](https://github.com/aryeelab/diffloop-paper)
- [Detailed HTML Vignette](vignette/diffloop_vignette.html), which contains a reproducible walk-through the core functionality of `diffloop`. Cloning and analyzing the data in the repository linked above enables users to similarly follow the analysis structure of the vignette. 
- [Writeup comparing variable](vignette/models_sizeFactors/writeup.pdf) contains a `.pdf` writeup
of assessing differential looping associations under different model specifications. In particular,
we assess specifying negative binomial regression models compared to using limma/voom as well
as examine performing variable differential loop association as a function of binned loop distance. 

### About
This repository contains more detailed use cases and analyses associated with the `diffloop` workflow.
Figure 1 was created using the R/Shiny interface [DNAlandscapeR](https://dnalandscaper.aryeelab.org). 

### Important information
- Install [git large file storage](https://git-lfs.github.com/) before cloning repository to access full data (~ 2 GB total).
- Link to ChIA-PET Pre-processing software: [mango](https://github.com/dphansti/mango)
- A recent version of [diffloop](https://bioconductor.org/packages/release/bioc/html/diffloop.html)
is required to run the code in this repository.
We'd recommend installing from Bioconductor to ensure a recent, stable version--
```
source("https://bioconductor.org/biocLite.R")
biocLite("diffloop")
```

- An older vignette showcasing some additional functionality and analysis not highlighted in this 
repository can be found [here](https://rpubs.com/caleblareau/diffloop_vignette).

<br><br>



