![322084859_1616597628753563_6938419133639811086_n](https://github.com/4n8x/GRN-Vtool/assets/51384420/ca38704f-4465-486f-9cdf-0fcbac36be75)
![GRNVtool_Pipeline](https://drive.google.com/uc?id=1mqSF4mNshepZZEbIg_TMZFZo-ZeifzL0)
## GRN Vtool
https://grnvtool.shinyapps.io/grn-vtool/
### Website Goal

##### GRN-Vtool is a website that allows biologists to run many GRN tools, in our website you we aid for:
- Efficient GRN studies that save biologists from many technical difficulties.
- Presenting the result of different tools in one place.
- Providing biologists with an interactive network where they can click on genes to get extra information will help the biologist to get further insight into the results.
- Making the Search process easier were it will help the biologist to navigate in a large network, while also marking connections for easier analysis.
- Visualize 2 GRN tools results at the same time for comparative analysis.
- Downloading results of the analysis. 


The field of Gene Regulatory Network (GRN) construction plays a pivotal role in unraveling cellular processes and gene pathways. However, conducting GRN analytical studies has posed significant challenges for biologists. Technical expertise is often required for tool installation and configuration, while the reliance on programming command lines hampers flexibility in adjusting tool parameters. Moreover, the absence of a standardized format for visualizing GRNs across different tools complicates comparative analysis. 

GRN VTOOL is a shiny app that offers pre-installed GRN visualization tools, including [DIANE](https://github.com/OceaneCsn/DIANE) and [SeqNet](https://github.com/tgrimes/SeqNet), which enable researchers to overcome technical hurdles. We chose to implement [DIANE](https://github.com/OceaneCsn/DIANE)  and [SeqNet](https://github.com/tgrimes/SeqNet) as they use [GENIE3](https://github.com/aertslab/GENIE3) and [GeneNetWeaver](https://github.com/tschaffter/genenetweaver) and [SynTReN](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-43) for GRN inference. GRN VTOOL presents results from the inference methods in an interactive dashboard with a unified interface.  A notable feature of GRN VTOOL is its capability to facilitate the comparison of 2 GRNs results in one session. These GRNs could result from different inference methods or resulted from using one inference tool on different datasets.

This is helpful for biologists, particularly in patient scenarios, where they can conveniently analyze GRNs before and after treatments or interventions thus eliminating the need for manual comparison. GRN VTOOL empowers biologists to efficiently conduct comparative analyses, gain deeper insights into regulatory dynamics, and accelerate their research, particularly in patient-focused studies. 


![GRNs Comparative Analysis](https://drive.google.com/uc?id=1a9lw88b5vZdin-HovQAl_mjefZOGkf7F)


User can also do preprocessings for the uploaded dataset such as normalization, differential expresion and clustering:

![Normalize](https://drive.google.com/uc?id=1s0ZgzGLLu-7DrPS7hcwdCJIDIQZK-pWn)


![Differential Expression](https://drive.google.com/uc?id=184hfaowqnTBcspPIYIsx_GQNHbYmwzeD)


![Clustering](https://drive.google.com/uc?id=1Nt7I2xqUXbU_25Uptn0pT9cnDxxUdcbr)



User can also view 2 datasets GRN analysis results as follow:

![Normalize2](https://drive.google.com/uc?id=1WEW2zgVwf-zP0fSeZIz5VmUdc1fH76OB)

![Differential Expression2](https://drive.google.com/uc?id=1t3hra2QbPzWjIPuRfli-BoKS0FZKE9kI)

![Clustering2](https://drive.google.com/uc?id=1qVwrNSPagkI9jRnSmpeOWmDZTOLiCULz)


![DIANE2](https://drive.google.com/uc?id=1m29I7ce_2stfQYcDFo6EWeRmjWb_i85n)


![SeqNet2](https://drive.google.com/uc?id=1Mh_2os19ip5SQGGDumVumBpl-j7PLHS4)

User can also downloads the GRN analysis processes he made, for GRN plotting only SeqNet results get downloaded:

![Download](https://drive.google.com/uc?id=1VpJ1MLguWj84H1yN78fyMOfcmJAGtjz8)
#### Datasets usage and preprocessing:
- Dataset used is 	"Molecular plant responses to combined abiotic stresses put a spotlight on unknown and abundant genes", which was done in expression profiling by high throughput sequencing. Our project mainly utilizes [Arabidopsis thaliana organism](https://bio.libretexts.org/Bookshelves/Introductory_and_General_Biology/Biology_(Kimball)/19%3A_The_Diversity_of_Life/19.01%3A_Eukaryotic_Life/19.1.06%3A_Arabidopsis_Thaliana_-_A_Model_Organism), hence this dataset is compatible with our pipeline.
  
- You can download the datasets that are valid for usage in our website from the [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146206).

The dataset needs to be preprocessed by removing all columns except those under raw counts from kallisto, you can use a preprocessed dataset from Datasets folder.

#### Launching instructions
You can use GRN-Vtool by entering the website link and registering an account in our wesbite to login, this account can be used to login to our website later.

