# CCSynergy
---

## Authors:
Sayed-Rzgar Hosseini and Xiaobo Zhou

---


## Abstract
Combination therapy is a promising strategy for confronting the complexity of cancer. However, experimental exploration of the vast space of potential drug combinations is costly and unfeasible. Therefore, computational methods for predicting drug synergy are much-needed for narrowing down this space, especially when examining new cellular contexts. Here, we thus introduce CCSynergy, a flexible, context-aware and integrative deep learning framework that we have established to unleash the potential of the Chemical Checker extended drug bioactivity profiles for the purpose of drug synergy prediction. We have shown that CCSynergy enables predictions of superior accuracy, remarkable robustness and improved context-generalizability as compared to the state-of-the-art methods in the field. Having established the potential of CCSynergy for generating experimentally validated predictions, we next exhaustively explored the untested drug combination space. This resulted in a compendium of potentially synergistic drug combinations on hundreds of cancer cell lines, which can guide future experimental screens.

---
### Requirements:
- Python 3.7.0
- Tensorflow 2.1.0
- Tensorboard 2.1.0
- Keras 2.3.1
- Numpy 1.21.6 
- Scipy 1.4.1
- Pandas 1.0.1
- Scikit-Learn 0.22.1

---

### Input files:
CCSynergy uses the following files as inputs, which can be found in the /Data folder in this repository.


#### 1. Drug synergy data
We have included both Merck.csv and Sanger.csv drug synergy dataset, which were used in this study (folder: /Data/Synergy_Data).
The rows in these files correspond to a unique (drug pairs + Cell line) triplet, and the columns represent drug and cell line names, their indices, and the corresponding synergy scores.
Index1 and Index2 indicate the row numbers corresponding to the drug1 and drug2 in the given Drug_Single_Info file, which is located in the /Data/Drug_Representation/ folder. Furthermore, Cindex and Tindex indicate the cell line and tissue indices respectively.


#### 2. Drug representation data
CCSynergy uses Chemical Checker drug signatures of type 2 to represent a given drug. We have downloaded the 25 CC spaces for each drug in the Merck and Sanger dataset, which are represented as 25 .csv files and are located in their corresponding subfolders in the /Data/Drug_Representation/ folder. The rows in these 25 .csv files correspond to a given drug, which is specified by the same row number in the corresonding Drug_Single_Info file. The 128 columns represent the signature vector of length 128 for each drug. 


#### 3. Cell line representation data
We have introduced 5 variants of CCSynergy each using one of the five different cell line representation methods, which we have stored as five .csv files in /Data/Cell_Representation folder. Each column in the .csv files correspond to a distinct cell line, which is represented as a vector of length 100.   


#### 4. Cross validation indexing files
CCSynergy requires three single-column input files specifying the testing, training and validation indices. 
We have provided an example in the /Data/Cross_Validation folder, which is generated within our CV1 scheme for the Merck dataset.
In general, these files can be generated using the CV1_splitter or CV2_splitter functions in the Cross_Validation.R program in this repository.
These functions take DIR (full path to the working directory), Dataset name ("Sanger" or "Merck"),seed value (for reproducibility), fold number and tissue index (only in CV2 scheme) as inputs and  generate the three .csv index files in the Cross_Validation folder.

---

### Parameters:
CCSynergy program requires the following parameters:
1. Variables (defined by the user)
- DIR (String): This specifies the full path to the working directory (Note: make sure that the /Data folder (with its four sub-folders) exists in the working directory)
- DATAset (String): This specifies the synergy dataset to be used, which is either "Merck" or "Sanger"
- CELLindex (Integer): The index of cell-line representation method that is an integer ranging between 1 and 5.
- DRUGindex (Integer): The index of Chemical Checker drug signature that is an integer ranging between 0 and 24. Each index corresponds to one of the 25 CC spaces: A1 (0) ... A5 (4) ... B5 (9) ... E5 (24).
- TASK (Integer): An integer specifying the nature of the learning task: regression (1) or classification (2). 
- FOLDnumber (Integer): The testing fold number, which is an integer between 1 and 5
- CVtype (Integer): CV1 (1) or CV2 (2)
  
2. Static parameters (already determined based on hyper-parameter optimization)
- n1 (Integer): The number of neurons in the first hidden layer (Default=2000)
- n2 (Integer): The number of neurons in the second hidden layer (Default=1000)
- n3 (Integer): The number of neurons in the third hidden layer (Default=500)
- batch (Integer): Batch size (Default=128)
- lr (Floating number): Learning rate (Default=0.0001)
- num (Integer): The size of input vector (Default=356) 

---

### Output files:
The main output file is simply a .csv file with two columns: 1) Real values and 2)Predicted values. This file can then be used in downstream analyses to calculate various metrics. Note that the file name is specified in line 98 of the CCSynergy code.
Furthermore, the model is also saved as a folder in the working directory, and its name is specified in line 90 of the CCSynergy code. 

---

### Examples:

#### Example 1:
In this example, we use the Merck drug synergy dataset for a regression task within the CV1 scheme to get predictions for the testing fold (i.e., FOLD_TEST=1) based on CCSynergy III (i.e., CELLindex=3) and by using C4 CC signatures as drug features (i.e., DRUGindex=13). 

Step1:
We need to generate the three index files in the Cross_Validation sub-folder by running the following R code.

We assume that the working directory is "/home/User/CCSynergy" and we set the seed as 28.

```Rscript
> source("/home/User/CCSynergy/Cross_Validation.R")
> CV1_Splitter("/home/User/CCSynergy","Merck",28,1)
```

Step2:
We then need to run the CCSynergy program as follows:
```shell
$ python3 /home/User/CCSynergy/CCSynergy.py "/home/User/CCSynergy" "Merck" 3 13 1 1 1
```

The output will be stored as a .csv file named as Results_Merck_CCSynergy3_C4_Fold1_CV1.csv, which will be located in the working directory.


#### Example 2:
In this example, we use the Sanger drug synergy dataset for a classification task within the CV2 scheme to get predictions for the testing fold (i.e., FOLD_TEST=4) and breast tissue (i.e., TINDEX=2) based on CCSynergy IV (i.e., CELLindex=5) and by using E3 CC signatures as drug features (i.e., DRUGindex=22). 

Step1:
We need to generate the three index files in the Cross_Validation sub-folder by running the following R code.

We assume that the working directory is "/home/User/CCSynergy" and we set the seed as 66.

```Rscript
> source("/home/User/CCSynergy/Cross_Validation.R")
> CV2_Splitter("/home/User/CCSynergy","Sanger",66,4,2)
```

Step2:
We then need to run the CCSynergy program as follows:
```shell
$ python3 /home/User/CCSynergy/CCSynergy.py "/home/User/CCSynergy" "Sanger" 5 22 2 4 2

```

The output will be stored as a .csv file named as Results_Sanger_CCSynergy5_E3_Fold4_CV2.csv, which will be located in the working directory.


---

## Citation

---

## License

- **[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)**
