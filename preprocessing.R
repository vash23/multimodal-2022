library(Seurat)
library(SeuratDisk)
library(Matrix)
library(rhdf5)

####     convert h5 to seurat format
#axis0 - features names
#axis1 - cell ID
#block0_values - count

#read RNA values
rna<-h5read("train_cite_inputs.h5","train_cite_inputs/axis0")
cell.ID<-h5read("train_cite_inputs.h5","train_cite_inputs/axis1")
cite.rna<-h5read("train_cite_inputs.h5","train_cite_inputs/block0_values")

#create sparse matrix
sparse.cite.rna<-Matrix(cite.rna, sparse = TRUE)

# remove dense matrix to save on memory
rm(cite.rna)

#add cell ID and RNA 
rownames(sparse.cite.rna)<-rna
colnames(sparse.cite.rna)<-cell.ID

#read protein values
protein<-h5read("train_cite_targets.h5","train_cite_targets/axis0")
cell.ID2<-h5read("train_cite_targets.h5","train_cite_targets/axis1")
cite.protein<-h5read("train_cite_targets.h5","train_cite_targets/block0_values")

#create sparse matrix
sparse.cite.protein<-Matrix(cite.protein, sparse = TRUE)

# remove dense matrix to save on memory
rm(cite.protein)

#check if cell IDs are equal
all.equal(cell.ID,cell.ID2)
rm(cell.ID2)

#add cell ID and protein
rownames(sparse.cite.protein)<-protein
colnames(sparse.cite.protein)<-cell.ID

#read test RNA values
test.rna<-h5read("test_cite_inputs.h5","test_cite_inputs/axis0")
cell.ID<-h5read("test_cite_inputs.h5","test_cite_inputs/axis1")
cite.test.rna<-h5read("test_cite_inputs.h5","test_cite_inputs/block0_values")

#create sparse matrix
sparse.cite.test.rna<-Matrix(cite.test.rna, sparse = TRUE)

#check if RNAs on test and train data sets are the same
all.equal(rna,test.rna)
rm(test.rna)

# remove dense matrix to save on memory
rm(cite.test.rna)

#add cell ID and test.rna
rownames(sparse.cite.test.rna)<-rna
colnames(sparse.cite.test.rna)<-cell.ID

#create Seurat object
cite <- CreateSeuratObject(counts = sparse.cite.rna)

#create new assays to store protein
protein.assay <- CreateAssayObject(counts = sparse.cite.protein)

#add the assay to Seurat object
cite[["protein"]] <- protein.assay

#verify assays
Assays(cite)

#remove sparse data to save on memory
rm(sparse.cite.rna)
rm(sparse.cite.protein)

rownames(cite[["protein"]])

#EDA
# default assay is RNA
DefaultAssay(cite) <- "RNA"
DefaultAssay(cite)

#cite <- NormalizeData(cite)
cite <- FindVariableFeatures(cite)
cite <- ScaleData(cite)
cite <- RunPCA(cite, verbose = FALSE)
cite <- FindNeighbors(cite, dims = 1:30)
cite <- FindClusters(cite, resolution = 0.8, verbose = FALSE)
cite <- RunUMAP(cite, dims = 1:30)
DimPlot(cite, label = TRUE)

#save to sseurat h5
SaveH5Seurat(cite, overwrite = TRUE)
