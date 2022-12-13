import pandas as pd
import numpy as np
import os
import random
import sys
import math
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.layers import Dense, Dropout, Input, Activation,Flatten
from sklearn.metrics import mean_squared_error
import scipy.stats as stats
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from tensorflow.keras import models, layers
from sklearn.metrics import r2_score,mean_absolute_error
import tensorflow as tf
from tensorflow import keras
from tensorflow.python.keras import backend as K
from tensorflow.python.ops import array_ops
from tensorflow.python.ops import math_ops
from sklearn.metrics import mean_squared_error


def Prepare_Features(DIR,DATAset,DRUGindex,num):
    FeatureName=['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5','D1','D2','D3','D4','D5','E1','E2','E3','E4','E5']
    data=pd.read_csv(str(DIR)+"/Data/Synergy_Data/"+str(DATAset)+".csv")
    cell_features=pd.read_csv(str(DIR)+"/Data/Cell_Representation/Cell"+str(num)+".csv")
    drug_features=np.array(pd.read_csv(str(DIR)+"/Data/Drug_Representation/"+str(DATAset)+"/"+FeatureName[DRUGindex]+".csv"))
    DD=np.array(data).shape[0]
    features=np.zeros(((2*DD),356))
    data.Cell_Line=data.Cell_Line.astype(str)
    DF=drug_features[:,1:]
    for i in range(DD):
        d1=data.Index1[i]
        d2=data.Index2[i]
        cell=str(data.Cell_Line[i]).strip() 
        features[i,:]=np.concatenate((DF[d1-1,:],DF[d2-1,:],np.array(cell_features[cell])),axis=0)#A-B orientation
        features[(i+DD),:]=np.concatenate((DF[d2-1,:],DF[d1-1,:],np.array(cell_features[cell])),axis=0)#B-A orientation
        print(i)
    return features


def Synergy_Scores(DIR,DATAset):
    data=pd.read_csv(str(DIR)+"/Data/Synergy_Data/"+str(DATAset)+".csv")
    DD=np.array(data).shape[0]
    responses=np.zeros(((2*DD)))
    Synergy=data['Synergy']
    for i in range(DD):
        responses[i]=Synergy[i]#A-B orientation
        responses[(i+DD)]=Synergy[i]#B-A orientation      
    return responses


def DNN(inputLength,n1=2000,n2=1000,n3=500,lr=0.0001,TASK=1):
    model = models.Sequential()
    model.add(layers.Dense(n1,kernel_initializer="he_normal", input_shape=[inputLength]))
    model.add(layers.Dropout(0.5))
    model.add(layers.Dense(n2, activation='relu',kernel_initializer="he_normal"))
    model.add(layers.Dropout(0.3))
    model.add(layers.Dense(n3, activation='tanh',kernel_initializer="he_normal"))
    if TASK==1: #Regression
       model.add(layers.Dense(1))
       model.compile( optimizer=keras.optimizers.Adam(learning_rate=float(lr),beta_1=0.9, beta_2=0.999, amsgrad=False), loss='mean_squared_error',metrics=['mse', 'mae'])
    else: #Classification
       model.add(layers.Dense(1,activation = "sigmoid")) 
       model.compile( optimizer=keras.optimizers.Adam(learning_rate=float(lr),beta_1=0.9, beta_2=0.999, amsgrad=False), loss = "binary_crossentropy", metrics=["accuracy"])
    return model


def CCSynergy(DIR,DATAset,CELLindex,DRUGindex,n1,n2,n3,batch,lr,num,TASK,FOLDnumber,CVtype):
    #Reading synergy data and indices
    data=pd.read_csv(str(DIR)+"/Data/Synergy_Data/"+str(DATAset)+".csv")    
    Training=pd.read_csv(str(DIR)+"/Data/Cross_Validation/Training_index.csv")## Needs to be modified depending on Cross validation scheme, synergy dataset and the testing fold
    Testing=pd.read_csv(str(DIR)+"/Data/Cross_Validation/Testing_index.csv")## Needs to be modified depending on Cross validation scheme, synergy dataset and the testing fold
    Validation=pd.read_csv(str(DIR)+"/Data/Cross_Validation/Validation_index.csv")## Needs to be modified depending on Cross validation scheme, synergy dataset and the testing fold
    FeatureName=['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','C1','C2','C3','C4','C5','D1','D2','D3','D4','D5','E1','E2','E3','E4','E5']
    #Features and synergy scores  
    X=Prepare_Features(DIR,DATAset,DRUGindex,CELLindex)
    Y=Synergy_Scores(DIR,DATAset)
    #Train, Test and Validation data
    X_train=X[Training.Index]
    Y_train=Y[Training.Index] 
    Y_train=np.array(Y_train)        
    X_test=X[Testing.Index]
    Y_test=Y[Testing.Index]    
    X_val=X[Validation.Index]
    Y_val=Y[Validation.Index]
    #Model training and prediction
    DNN_model=DNN(num,n1,n2,n3,lr,TASK)
    cb_check = ModelCheckpoint((str(DIR)+'/Results_'+str(DATAset)+'_CCSynergy'+str(CELLindex)+'_'+str(FeatureName[DRUGindex])+'_Fold'+str(FOLDnumber)+'_CV'+str(CVtype)), verbose=1, monitor='val_loss',save_best_only=True, mode='auto')
    DNN_model.fit(x=X_train,y=Y_train,batch_size=batch,epochs = 1000,shuffle=True,validation_data = (X_val,Y_val),callbacks=[EarlyStopping(monitor='val_loss', mode='auto', patience = 10),cb_check] )
    DNN_model = tf.keras.models.load_model((str(DIR)+'/Results_'+str(DATAset)+'_CCSynergy'+str(CELLindex)+'_'+str(FeatureName[DRUGindex])+'_Fold'+str(FOLDnumber)+'_CV'+str(CVtype)))
    predict_test=DNN_model.predict(X_test)
    #Final output
    Final=pd.DataFrame()
    Final['Real']=Y_test     
    Final['Prediction']=predict_test
    Final.to_csv((str(DIR)+'/Results_'+str(DATAset)+'_CCSynergy'+str(CELLindex)+'_'+str(FeatureName[DRUGindex])+'_Fold'+str(FOLDnumber)+'_CV'+str(CVtype)+'.csv'))## The final output file 

   
def main():
    ## Parameters and their default values
    DIR="" #Working directory (Full Path)
    DATAset="Merck" #Dataset to be used (e.g. "Merck" or "Sanger")
    CELLindex=0    #Index for cell line representation method (1-5)
    DRUGindex=0    #Index for drug representation method (0-24)
    TASK=1  #Regression (1) or Classification (2) task
    FOLDnumber=0 #an integer ranginging between 1 and 5
    CVtype=0 #CV1 (1) or CV2 (2)
    n1=2000 #number of neurons in the first layer
    n2=1000 #number of neurons in the second layer
    n3=500  #number of neurons in the third layer
    lr=0.0001 #learning rate
    batch=128 #batch size
    num=356   #size of the input vector
    ## Parameters to be specified by the user
    DIR=sys.argv[1]
    DATAset=sys.argv[2]
    CELLindex=sys.argv[3]
    DRUGindex=sys.argv[4]
    TASK=sys.argv[5]
    FOLDnumber=sys.argv[6]
    CVtype=sys.argv[7]
    ## Running the software program
    CCSynergy(str(DIR),str(DATAset),int(CELLindex),int(DRUGindex),int(n1),int(n2),int(n3),int(batch),float(lr),int(num),int(TASK),int(FOLDnumber),int(CVtype))


main()


