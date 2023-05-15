## GRN Vtool
https://grn-vtool.shinyapps.io/grn-vtool/
### Website Goal

##### GRN-Vtool is a website that allows biologists to run many GRN tools, in our website you we aid for:
- Efficient GRN studies that save biologists from many technical difficulties.
- Presenting the result of different tools in one-place where the similarities and differences are automatically highlighted helps biologists to conduct efficient comparative analysis.
- Providing the biologists with an interactive network where they can click on genes to get extra information will help the biologist to get further insight on the results.
- Making the Search process easier were it will help the biologist to navigate in large network.
- Partial network extraction will aid them for productive analysis. 

### Technology used
#### In GRN-Vtool this we chose to use DIANE as the first tool to infer GRN from gene expression data. 
DIANE is an R package, is a simple exploratory visualization, that allows the user to observe the normalized expression levels of several genes of interest. DIANE enables gene expression profiles clustering using the statistical framework for inferring mixture models through an Expectation-Maximization (EM) algorithm[21],[22].In DIANE, the package chosen for GRN reconstruction is GENIE3, a machine learning procedure that was among the best performers of the DREAM challenges. GENIE3 uses Random Forests[23], which is a machine learning method based on the inference of a collection of regression trees. It has the advantage of being a non-parametric procedure, requiring very few modelling or biological priors, while being able to capture interactions and high order combinatorics between regulators. 

 

#### SeqNet is the second method we chose for GRN inference, it uses Twelve co-expression methods 
are considered: ARACNE [24], BC3NET[25],C3NET [26], CLR [27],DWLasso[28],GENIE3 [29],Glasso[30],MR-NET [31],Pearson correlation, Shrinkage[32],Silencer[33],and Spearman correlation. ARACNE, CLR, and MRNET were run using the default settings provided by the minet R package [34], with "spearman" being used as the entropy estimator. In Figure you can see a summary of the main SeqNet functions.
<img width="468" alt="SeqNetTable" src="https://user-images.githubusercontent.com/51384420/236009398-70b4a36b-1d10-4825-b889-6b83b9c5c57d.png">


#### You can use Two tools in the Website to generate your network [SeqNet](https://github.com/tgrimes/SeqNet),[DIANE](https://github.com/OceaneCsn/DIANE).


###Lunching instructions
You can use GRN-Vtool by entring the website and folow the steps:
- Upload your normlized dataset from Data Tab>Upload Data
- Choose your appropriate tool to use for generating the network from Tools Tab
- Generate the network from Generate Network Tab
- You can download the output in your local machine
- You can reguster to save you genareted networks for further usage 
 
