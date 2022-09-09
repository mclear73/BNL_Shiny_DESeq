## RNA-UI   


### Description 
This Shiny application allows for the analysis of Raw count data, differential expression analysis, and some post-DEG analyses.

The original UI design and framework was based off of the similar shiny application, [GENAVi](https://github.com/alpreyes/GENAVi). Here is a link to the [manuscript](https://link.springer.com/article/10.1186/s12864-019-6073-7). 

GENAVi was limited to Human or Mouse data only. Thus, while I appreciated parts of their analysis framework, it was not applicable to the broad range of data sets that can be encountered outside of human genomics.

**Thus far, RNA-UI has been adapted to integrate the following:**

1) Alternative DE analysis through edgeR  
2) Added a sample selection option so that a subset of data can be visualized  
3) Addition of 95% confidence interval clustering on PCA  
4) Substantially modified enrichment analysis for custom datasets  
5) KEGG Mapping  
6) Keyword and gene set extractions from annotated data for in-depth analysis of DEG lists and subsequent generation of figures  

**This code is still under development and should only be used as a reference. **

### To Deploy
- Clone the repository into a local directory where you plan on running your analysis
- Open your favorite R IDE such as RStudio
- Open the server.R file
- Press run app
- Follow the instructions in your R IDE for downloading the required packages
- Depending on your IDE, one of the following will happen:
	- RStudio will open a UI window that you can interact with
	- R will open a webpage that is hosted locally. You can interact with the app here