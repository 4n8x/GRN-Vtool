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

### Technology used
#### In GRN-Vtool we chose to use DIANE as the first tool to infer GRN from gene expression data. 
DIANE is an R package, is a simple exploratory visualization, that allows the user to observe the normalized expression levels of several genes of interest. DIANE enables gene expression profile clustering using the statistical framework for inferring mixture models through an Expectation-Maximization (EM) algorithm. In DIANE, the package chosen for GRN reconstruction is GENIE3, a machine-learning procedure among the best DREAM challenges performers. GENIE3 uses Random Forests[23], which is a machine-learning method based on the inference of a collection of regression trees. It has the advantage of being a non-parametric procedure, requiring very few modeling or biological priors while being able to capture interactions and high-order combinatorics between regulators. 

 

#### SeqNet is the second method we chose for GRN inference, it uses Twelve co-expression methods 

<img width="468" alt="SeqNetTable" src="https://user-images.githubusercontent.com/51384420/236009398-70b4a36b-1d10-4825-b889-6b83b9c5c57d.png">


#### You can use Two tools on the Website to generate your network [SeqNet](https://github.com/tgrimes/SeqNet),[DIANE](https://github.com/OceaneCsn/DIANE).


### Lunching instructions
You can use GRN-Vtool by entering the website and following the steps:
- Upload your normalized dataset from Data Tab>Upload Data
- Choose your appropriate tool to use for generating the network from the Tools Tab
- Generate the network from Generate Network Tab
- You can download the output from your local machine
- You can register to save your generated networks for further usage 
 
