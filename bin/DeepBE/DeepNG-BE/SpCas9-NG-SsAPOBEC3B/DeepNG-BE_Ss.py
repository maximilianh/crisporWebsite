import tensorflow as tf
import numpy as np
import pandas as pd
import itertools

header=['target + PAM','feature']

dataset_ = pd.read_csv('DeepNG-BE_input_example.csv',header=None,names=header)
final_model =  tf.keras.models.load_model('DeepNG-BE_Ss_model',compile=False)

start=6
length=7

def preprocess_seq(data):
    print("Start preprocessing the sequence done 2d")
    global length
    global start

    DATA_X = np.zeros((len(data),length,4), dtype=int)
    print(np.shape(data), len(data), length)
    for l in range(len(data)):
        for i in range(length):

            try: data[l][i+start]
            except: print(data[l], i+start, length, len(data))

            if data[l][i+start]in "Aa":    DATA_X[l, i, 0] = 1
            elif data[l][i+start] in "Cc": DATA_X[l, i, 1] = 1
            elif data[l][i+start] in "Gg": DATA_X[l, i, 2] = 1
            elif data[l][i+start] in "Tt": DATA_X[l, i, 3] = 1
            elif data[l][i+start] in "Xx": pass
            else:
                print("Non-ATGC character " + data[l])
                print(i)
                print(data[l][i+start])
                sys.exit()
        #loop end: i
    #loop end: l
    print("Preprocessing the sequence done")
    return DATA_X

seq_length = 30

def preprocess_seq_pam(data):
    print("Start preprocessing the sequence done 2d")
    global seq_length

    DATA_X = np.zeros((len(data),seq_length,4), dtype=int)
    print(np.shape(data), len(data), seq_length)
    for l in range(len(data)):
        for i in range(seq_length):

            try: data[l][i]
            except: print(data[l], i, seq_length, len(data))

            if data[l][i]in "Aa":    DATA_X[l, i, 0] = 1
            elif data[l][i] in "Cc": DATA_X[l, i, 1] = 1
            elif data[l][i] in "Gg": DATA_X[l, i, 2] = 1
            elif data[l][i] in "Tt": DATA_X[l, i, 3] = 1
            elif data[l][i] in "Xx": pass
            else:
                print("Non-ATGC character " + data[l])
                print(i)
                print(data[l][i])
                sys.exit()
        #loop end: i
    #loop end: l
    print("Preprocessing the sequence done")
    return DATA_X


dataset_seq_masked = preprocess_seq(dataset_['target + PAM'])
dataset_seq_masked = pd.Series(list(dataset_seq_masked),name='seq')

dataset_seq_masked_pam = preprocess_seq_pam(dataset_['target + PAM'])
dataset_seq_masked_pam = pd.Series(list(dataset_seq_masked_pam ),name='seq_pam')

dataset_all = pd.concat([dataset_,dataset_seq_masked,dataset_seq_masked_pam],axis=1)

X_test_seq = np.stack(dataset_all['seq'])
X_test_seq_pam = np.stack(dataset_all['seq_pam'])
hyperparameter_prediction = final_model.predict([X_test_seq_pam,X_test_seq_pam], batch_size=128)
hyperparameter_prediction = pd.DataFrame(hyperparameter_prediction)

alist = []
for x in range(length):
    b = []
    for i in range(1,2**x+1):
        for k in range(1,2**(6-x)+1):
            b.append(2**(7-x)*i-k-1)
    alist.append(b)

for i in range(len(hyperparameter_prediction)):
    eff = hyperparameter_prediction.iloc[i].sum()
    motif = dataset_all.iloc[0,i][start:start+length]
    non = []
    for k in range(len(motif)):
        if motif[k] in 'Cc':
            non.append(2**(6-k))
    on = []
    for j in range(len(non)):
        iter = list(itertools.combinations(non,j+1))
        for k in iter:
            on.append(sum(k))
    on = [x - 1 for x in on]
    off = list(set(list(range(127)))-set(on))
    for m in off:
        hyperparameter_prediction[hyperparameter_prediction==hyperparameter_prediction.iloc[i,m]] = 0
    patt = hyperparameter_prediction.iloc[i].sum()

    hyperparameter_prediction.iloc[i] = (hyperparameter_prediction.iloc[i]/patt)*eff

prediction = hyperparameter_prediction.stack().reset_index().iloc[:,2]
seq = pd.DataFrame(np.repeat(dataset_all['target + PAM'].reset_index(drop=True).values,127,axis=0))
pam = pd.DataFrame(np.repeat(dataset_all['feature'].reset_index(drop=True).values,127,axis=0))
hyperparameter_prediction=pd.concat([seq,seq,pam,prediction], axis=1)
hyperparameter_prediction.columns=['target + PAM','edited output','PAM_variant','Prediction score']

for i in range(int(len(hyperparameter_prediction)/127)):
    for k in range(len(alist)):
        for j in alist[k]:
            hyperparameter_prediction.iloc[i+j,1] = hyperparameter_prediction.iloc[i+j,1][:k+6] + "t" + hyperparameter_prediction.iloc[i+j,1][k+1+6:]
hyperparameter_prediction = hyperparameter_prediction[hyperparameter_prediction['Prediction score'] != 0]
hyperparameter_prediction.to_excel('prediction_result.xlsx' ,engine='openpyxl')
