import tensorflow as tf
import numpy as np
import pandas as pd

header=['target + PAM','feature']

dataset_ = pd.read_csv('PAM_variant_input_example.csv',header=None,names=header)
final_model =  tf.keras.models.load_model('PAM_variant_NRCH_model.h5',compile=False)

def preprocess_seq(data,length):
    print("Start preprocessing the sequence done 2d")

    DATA_X = np.zeros((len(data),length,4), dtype=int)
    print(np.shape(data), len(data), length)
    for l in range(len(data)):
        for i in range(length):

            try: data.iloc[l][i]
            except: print(data.iloc[l],l, i, length, len(data))

            if data.iloc[l][i]in "Aa":    DATA_X[l, i, 0] = 1
            elif data.iloc[l][i] in "Cc": DATA_X[l, i, 1] = 1
            elif data.iloc[l][i] in "Gg": DATA_X[l, i, 2] = 1
            elif data.iloc[l][i] in "Tt": DATA_X[l, i, 3] = 1
            elif data.iloc[l][i] in "Xx": pass
            else:
                print("Non-ATGC character " + data[l])
                print(i)
                print(data.iloc[l][i])
                sys.exit()
        #loop end: i
    #loop end: l
    print("Preprocessing the sequence done")
    return DATA_X

dataset_seq_masked = preprocess_seq(dataset_['target + PAM'],30)

dataset_seq_masked = pd.Series(list(dataset_seq_masked),name='seq')

dataset_all = pd.concat([dataset_,dataset_seq_masked],axis=1)

X_test_seq = np.stack(dataset_all['seq'])
hyperparameter_prediction = final_model.predict(X_test_seq, batch_size=128)
hyperparameter_prediction = pd.DataFrame(hyperparameter_prediction)

hyperparameter_prediction=pd.concat([dataset_all['target + PAM'].reset_index(drop=True),dataset_all['feature'].reset_index(drop=True),hyperparameter_prediction.reset_index(drop=True)],axis=1,ignore_index=True)
hyperparameter_prediction.columns=['target + PAM','PAM_variant','Prediction score']
hyperparameter_prediction.to_excel('prediction_result.xlsx' ,engine='openpyxl')
