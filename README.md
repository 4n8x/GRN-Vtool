![322084859_1616597628753563_6938419133639811086_n](https://github.com/4n8x/GRN-Vtool/assets/51384420/ca38704f-4465-486f-9cdf-0fcbac36be75)
## GRN Vtool
https://grnvtool.shinyapps.io/grn-vtool/
### Website Goal

##### GRN-Vtool is a website that allows biologists to run many GRN tools, in our website you we aid for:
- Efficient GRN studies that save biologists from many technical difficulties.
- Presenting the result of different tools in one place where the similarities and differences are automatically highlighted helps biologists to conduct efficient comparative analysis.
- Providing biologists with an interactive network where they can click on genes to get extra information will help the biologist to get further insight into the results.
- Making the Search process easier were it will help the biologist to navigate in a large network.
- Partial network extraction will aid them in productive analysis. 


The field of Gene Regulatory Network (GRN) construction plays a pivotal role in unraveling cellular processes and gene pathways. However, conducting GRN analytical studies has posed significant challenges for biologists. Technical expertise is often required for tool installation and configuration, while the reliance on programming command lines hampers flexibility in adjusting tool parameters. Moreover, the absence of a standardized format for visualizing GRNs across different tools complicates comparative analysis. 

GRN VTOOLS is a web-based platform that offers pre-installed GRN visualization tools, including DIANE and SeqNet, which enable researchers to overcome technical hurdles. We chose to implement DIANE  and SeqNet as they use GENIE3 and GeneNetWeaver, SynTReN, and Rogers for GRN inferenc. GENIE3 is a random forest ensemble machine learning procedure that was among the best performers of the DREAM challenges. GRN VTOOLS presents results from the inference methods in an interactive dashboard with a unified interface.  A notable feature of GRN VTOOL is its capability to facilitate the comparison of two GRNs in one session. These GRNs could result from different inference methods or resulted from using one inference tool on different datasets.

This is helpful for biologists, particularly in patient scenarios, where they can conveniently analyze GRNs before and after treatments or interventions thus eliminating the need for manual comparison. GRN VTOOL empowers biologists to efficiently conduct comparative analyses, gain deeper insights into regulatory dynamics, and accelerate their research, particularly in patient-focused studies. 

 

#### SeqNet is the second method we chose for GRN inference, it uses Twelve co-expression methods 

<img width="468" alt="SeqNetTable" src="https://user-images.githubusercontent.com/51384420/236009398-70b4a36b-1d10-4825-b889-6b83b9c5c57d.png">


#### You can use Two tools on the Website to generate your network [SeqNet](https://github.com/tgrimes/SeqNet),[DIANE](https://github.com/OceaneCsn/DIANE).

#### Datasets usage and preprocessing:
- Dataset used is 	"Molecular plant responses to combined abiotic stresses put a spotlight on unknown and abundant genes", which was done in expression profiling by high throughput sequencing. Our project mainly utilizes Arabidopsis thaliana organism, hence this dataset is compatible with our pipeline.[2]
  
- You can download the datasets that are valid for usage in our website from the [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146206).

The dataset needs to be preprocessed by removing all columns except those under raw counts from kallisto, including "target_id" column.

#### Launching instructions
You can use GRN-Vtool by entering the website and following the steps:
- Register an account in our wesbite to login, this account can be used to login to our website later.
- Upload your Dataset from Data Tab > Upload Data, you can upload up to 2 datasets.
- Perform Normalization, Differential Expression, Clustering if needed, the 2nd file uploaded will be assumed to have those already done to the dataset.
- Generate the network from Generate Network Tab after choosing a tool (SeqNet or DIANE, or both).
- You can download the results to your local machine as csv format, you can only download the results after performing the entire pipeline for processing.


 #### References
 [1]
 [2] https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146206
