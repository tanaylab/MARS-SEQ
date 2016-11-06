# scrdb design

## general outlines (2/11/16)

 - The purpose of this project is to create an infrastructure for holding all of the single cell RNA data in one place, 
  and provide an interface for loading selected parts of the data to R, process them, and save for later use. (Similar to what is done with Misha)
 - This DB will include both data created in this lab, as well as data imported from other labs, including scRNA data produced using protocols other than Mars-SEQ 
 which of course have different biases and other issues).
 
### Use cases 
1. **Upload:** User has acquired new data, she needs to run a script that will upload the raw data to the file system in the right format.
The script will upload the data to a specific folder, and also add a record in the root directory of the DB, adding it to the list of datasets.
2. User wants to **extract** data from the DB. 
	- She initializes an *scrdb metadata* object by supplying the home directory. The object reads the *experiments* file from the root and returns the list of available datasets.
	Alternatively, the user supplies a subfolder. In this case, the object goes over all sub-directories and searches for metadata files.
	- The user can now review all the cells/batches under this folder and ask for a subset (by experiment ID, batch , minimal UMI count, subject ID etc.)
	- The subset is loaded to memory and a new *data* object is created, alog with a *metadata* fields describing how the data was acquired
	 (which directories, what filter was used, date of creation, user...)
3. **Data normalization** 
	- The *data* class should  have a normalization method that would allow to chose from one of the built-in normalizations (down-sample, divide by colSum...) 
	or apply a user defined normalization. The method of normalization should be recorded in the metadata of this object. 
4. **Basic analysis pipeline:** The user can automatically create basic statistics and plots: Gene correlations, cell correlations, 
5. **Saving results**
	- After filtering and normalizing data, the user should use a *save to file* method that will store the filtered and normalized objects, along with the metadata, 
	including documentation on what filtering and normalizations were used. 

### File system and data representation 
- **Home directory:** /net/mraid14/export/tgdata/db/tgdb/scrdb. 
- This folder will contain a sub directory per species/build (mm9, hg19...)
- Each experiment (i.e. data from a specific paper, a single sequencing run etc.) will be located in a separate subdirectory. Related experiments should be grouped under one parent folder. For example:  
	- /net/mraid14/export/tgdata/db/tgdb/mm9/scrdb/embryo/e8.5_20160101 
	- /net/mraid14/export/tgdata/db/tgdb/mm9/scrdb/embryo/20160101
- Each subdirectory will contain: **scdb.meta** file, describing the experiment; **cells.txt** file, describing each cell; **data.Rda** file, containing the gene expression counts of the cells.
1. **scdb_meta.r**: A text file. contains the following fields: 
	
Name | Type | Mandatory? | Example | Value 
--- | :---: | --- | --- |---
desc | c | yes | "Mouse embryonic cells day 8.5..." |  Short description of the data.
lab | c | yes | "Tanay" | Lab name from which the data was taken (usually surname of the PI) 
date_created | c | yes | "01012016" | The date the data was created. For internal data - usually date of sequencing. For external data - date of paper publishing/GEO record creation  
date_added | c | yes | "01012016" | The date of adding data to the DB. Format DDMMYYYY
paper | c | no | "Paul et al. (2015) Cell" | Paper name 
geo | list | no | list=("GSE72857", "GSE72858") | GEO accession numbers. 
url | c | no | "" |  URL from which the dataset was downloaded. 

2. **cells.txt**: asd
3. **data.Rda**:
	
	
	
### Other notes

---

## Tasks: 
1. Create a list of all available UMI tables and their metadata, including downloaded from the web.
2. Formally describe the format of the cells data frame, what are the obligatory fields, what are variable?
3. Write a script for the upload use case.

  
---

## Design

- save a unique name and a short display name for the experiment. The unique name can be extracted from the directory name. 
