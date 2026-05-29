import tensorflow as tf
import numpy as np
import pandas as pd
import itertools

header=['target + PAM','feature']

dataset_ = pd.read_csv('DeepBE_input_example.csv',header=None,names=header)
final_model =  tf.keras.models.load_model('DeepBE_YE1_model',compile=False)

model0 = tf.keras.models.load_model('PAM_variant_SpCas9_model.h5',compile=False)
model1 = tf.keras.models.load_model('PAM_variant_VRQR_model.h5',compile=False)
model2 = tf.keras.models.load_model('PAM_variant_NG_model.h5',compile=False)
model3 = tf.keras.models.load_model('PAM_variant_NRRH_model.h5',compile=False)
model4 = tf.keras.models.load_model('PAM_variant_NRTH_model.h5',compile=False)
model5 = tf.keras.models.load_model('PAM_variant_NRCH_model.h5',compile=False)
model6 = tf.keras.models.load_model('PAM_variant_SpG_model.h5',compile=False)
model7 = tf.keras.models.load_model('PAM_variant_SpRY_model.h5',compile=False)
model8 = tf.keras.models.load_model('PAM_variant_Sc++_model.h5',compile=False)

model_list = [model0, model1, model2, model3, model4, model5, model6, model7, model8]

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

start_proportion=4
length_proportion=20

def preprocess_seq_proportion(data):
    print("Start preprocessing the sequence done 2d")
    global length_proportion
    global start_proportion

    DATA_X = np.zeros((len(data),length_proportion,4), dtype=int)
    print(np.shape(data), len(data), length_proportion)
    for l in range(len(data)):
        for i in range(length_proportion):

            try: data[l][i+start_proportion]
            except: print(data[l], i+start_proportion, length_proportion, len(data))

            if data[l][i+start_proportion]in "Aa":    DATA_X[l, i, 0] = 1
            elif data[l][i+start_proportion] in "Cc": DATA_X[l, i, 1] = 1
            elif data[l][i+start_proportion] in "Gg": DATA_X[l, i, 2] = 1
            elif data[l][i+start_proportion] in "Tt": DATA_X[l, i, 3] = 1
            elif data[l][i+start_proportion] in "Xx": pass
            else:
                print("Non-ATGC character " + data[l])
                print(i)
                print(data[l][i+start_proportion])
                sys.exit()
        #loop end: i
    #loop end: l
    print("Preprocessing the sequence done")
    return DATA_X


dataset_seq_masked = preprocess_seq(dataset_['target + PAM'])
dataset_seq_masked = pd.Series(list(dataset_seq_masked),name='seq')

dataset_seq_masked_pam = preprocess_seq_pam(dataset_['target + PAM'])
dataset_seq_masked_pam = pd.Series(list(dataset_seq_masked_pam ),name='seq_pam')

dataset_seq_masked_proportion = preprocess_seq_proportion(dataset_['target + PAM'])
dataset_seq_masked_proportion = pd.Series(list(dataset_seq_masked_proportion ),name='seq_proportion')

dataset_all = pd.concat([dataset_,dataset_seq_masked,dataset_seq_masked_pam,dataset_seq_masked_proportion],axis=1)

X_test_seq_predict = np.stack(dataset_all['seq'])
X_test_seq_pam_predict = np.stack(dataset_all['seq_pam'])
X_test_seq_proportion_predict = np.stack(dataset_all['seq_proportion'])

for i in range(len(model_list)):
    hyperparameter_prediction = model_list[0].predict(X_test_seq_pam_predict, batch_size=128)
    hyperparameter_prediction = pd.DataFrame(hyperparameter_prediction)
    hyperparameter_prediction.columns = ['Prediction'+str(i)]
    dataset_all = pd.concat([dataset_all.reset_index(drop=True),hyperparameter_prediction.reset_index(drop=True)],axis=1)

for i in range(len(dataset_all)):
    dataset_all.iloc[i,np.r_[5:dataset_all.iloc[i,1]+5, dataset_all.iloc[i,1]+6:14]] = 0

dataset_all = pd.concat([dataset_all,dataset_all.iloc[:,5:14].sum(axis=1)], axis=1)
dataset_all = dataset_all.rename(columns={0: 'Sum'})

X_test_seq = np.stack(dataset_all['seq'])
X_test_seq_pam = np.stack(dataset_all['seq_pam'])
X_test_seq_proportion = np.stack(dataset_all['seq_proportion'])
X_test_PAM_predict = np.stack(dataset_all['Sum'])

hyperparameter_prediction = final_model.predict([[X_test_seq,X_test_PAM_predict],X_test_seq_proportion], batch_size=128)
hyperparameter_prediction = pd.DataFrame(hyperparameter_prediction)

alist = []
short_length = 5
short_start = 7
for x in range(short_length):
    b = []
    for i in range(1,2**x+1):
        for k in range(1,2**(4-x)+1):
            b.append(2**(5-x)*i-k-1)
    alist.append(b)

for i in range(len(hyperparameter_prediction)):
    eff = hyperparameter_prediction.iloc[i].sum()
    motif = dataset_all.iloc[0,i][short_start:short_start+short_length]
    non = []
    for k in range(len(motif)):
        if motif[k] in 'Cc':
            non.append(2**(4-k))
    on = []
    for j in range(len(non)):
        iter = list(itertools.combinations(non,j+1))
        for k in iter:
            on.append(sum(k))
    on = [x - 1 for x in on]
    off = list(set(list(range(31)))-set(on))
    for m in off:
        hyperparameter_prediction[hyperparameter_prediction==hyperparameter_prediction.iloc[i,m]] = 0
    patt = hyperparameter_prediction.iloc[i].sum()

    hyperparameter_prediction.iloc[i] = (hyperparameter_prediction.iloc[i]/patt)*eff

prediction = hyperparameter_prediction.stack().reset_index().iloc[:,2]
seq = pd.DataFrame(np.repeat(dataset_all['target + PAM'].reset_index(drop=True).values,31,axis=0))
pam = pd.DataFrame(np.repeat(dataset_all['feature'].reset_index(drop=True).values,31,axis=0))
hyperparameter_prediction=pd.concat([seq,seq,pam,prediction], axis=1)
hyperparameter_prediction.columns=['target + PAM','edited output','PAM_variant','Prediction score']

for i in range(int(len(hyperparameter_prediction)/31)):
    for k in range(len(alist)):
        for j in alist[k]:
            hyperparameter_prediction.iloc[i+j,1] = hyperparameter_prediction.iloc[i+j,1][:k+7] + "t" + hyperparameter_prediction.iloc[i+j,1][k+1+7:]
hyperparameter_prediction = hyperparameter_prediction[hyperparameter_prediction['Prediction score'] != 0]
hyperparameter_prediction.to_excel('prediction_result.xlsx' ,engine='openpyxl')
