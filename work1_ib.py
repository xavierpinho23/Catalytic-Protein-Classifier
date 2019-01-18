import pandas as pd
import os
import numpy as np
import urllib.request as ur
from sklearn import cross_validation
from sklearn import svm
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc,roc_auc_score
from sklearn.multiclass import OneVsRestClassifier
from sklearn import metrics

### 1 Work - Introduction to Bioinformatics    
### Xavier Pinho & Jorge Melo, University of Coimbra, 2018/2019

home = os.getcwd()

link = "http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640%20go:3824&fil=&force=no&format=fasta"
f = ur.urlopen(link)
reader = f.read().decode("utf-8")

names = []
seq = []
seq1 = []
dicts =[]
count = 0
list1 = "https://github.com/xavierpinho23/Music_DataBase".split(">")

for i in list1:
	if count == 0:
		pass
	else:
		list2 = i.split("|")
		if len(list2) < 2:
			pass
		else:
			names.append(list2[1])
			lista3 = i.split("\n")
			seq.append(lista3[1:])
	count = count + 1

for k in seq:
    string = ""
    for l in k:
        string = string + l
    seq1.append(string)

#aminoacids     
letters = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","Y","Z","X"]

frequency_list = []

for j in range(len(names)):
    aa = seq1[j]
    dic = {}
    for i in letras:
        if i not in dic:
            dic[i] = 0
    for k in aa:
        if k not in dic:
            dic[k] = 1
        else:
            dic[k] = dic[k] + 1
    frequency_list.append(dic)
    
#create dataframe with proteins with catalytic activity
df = pd.DataFrame(frequency_list)
df.insert(0, 'names', names)
#create target column with the value 1 (1 = catalytic) 
df['target'] = 1

#export dataframe to csv file
#filename = "path/table1.csv"
#pd.DataFrame.to_csv(df,  filename,header = True,)

#############################################################################
new_link = "http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=proteome:UP000005640&fil=&force=no&format=fasta"
f_novo = ur.urlopen(new_link)
new_reader = f_novo.read().decode("utf-8")

new_names = []
new_seq = []
new_seq1 = []
new_dicts =[]
new_count = 0
new_list1 = new_reader.split(">")

for i in new_list1:
	if new_count == 0:
		pass
	else:
		new_list2 = i.split("|")
		if len(new_list2) < 2:
			pass
		else:
			new_names.append(new_list2[1])
			new_list3 = i.split("\n")
			new_seq.append(new_list3[1:])
	new_count = new_count + 1

for k in new_seq:
    string_novo = ""
    for l in k:
        string_novo = string_novo + l
    new_seq1.append(string_novo)
    
frequency_list_new = []

for j in range(len(new_names)):
    aa = new_seq1[j]
    dic_new = {}
    for i in letras:
        if i not in dic_new:
            dic_new[i] = 0
    for k in aa:
        if k not in dic_new:
            dic_new[k] = 1
        else:
            dic_new[k] = dic_new[k] + 1
    frequency_list_new.append(dic_new)

#create dataframe with all the proteins
df_novo = pd.DataFrame(frequency_list_new)
df_novo.insert(0, 'names', new_names)
#create target column with value 0 for non-catalytic proteins 
df_novo['target'] = 0

for i in range(len(new_names)):
    if new_names[i] in names:
        #if the protein it's on the 1st dataframe then it has catalytic activity
        df_novo.set_value(col = 'target', index = i, value = 1)

#export dataframe to a csv file
#filename_new = "path/tabel_new.csv"
#pd.DataFrame.to_csv(df_novo, filename_new, header = True)

#df_novo['target'].value_counts()
#0    58486
#1    14613
#data non-balanced

#split dataset into train and test
X_train, X_test, y_train, y_test = cross_validation.train_test_split(df_novo.values[:,1:-1], df_novo.target, test_size=0.4, random_state=0)
print("X_train size:" + str(X_train.shape))
print("X_test  size:" + str(X_test.shape))
print("y_train size:" + str(y_train.shape))
print("y_test  size:" + str(y_test.shape))

#ignore the 1st column and convert object to float
X_train = X_train[:,:].astype(float)
y_train = y_train[:]
X_test  = X_test[:,:].astype(float)
y_test  = y_test[:]

#SVM - Support Vector Machines
clf = svm.SVC(kernel='linear', C=1, probability = True).fit(X_train, y_train)
acc = clf.score(X_test, y_test)
print("Accuracy: " + str(acc))

#ROC e AUC
probability = clf.predict_proba(X_test)
pds = probability[:,1]
fpr,tpr, threshold = metrics.roc_curve(y_test, pds)
roc_auc = metrics.auc(fpr,tpr)

#ROC plot
plt.figure()
plt.plot(fpr, tpr,'r', label = 'AUC - %0.2f' % roc_auc)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.show()

#Data balancing
#SVM - Support Vector Machines
clf1 = svm.SVC(kernel = 'linear',C=1, probability = True, class_weight = 'balanced').fit(X_train,y_train)
acc1 = clf1.score(X_test,y_test)
print("Accuracy: " + str(acc1))

#ROC e AUC
probability1 = clf1.predict_proba(X_test)
pds1 = probability1[:,1]
fpr1,tpr1, threshold1 = metrics.roc_curve(y_test, pds1)
roc_auc1 = metrics.auc(fpr1,tpr1)

#ROC plot
plt.figure()
plt.plot(fpr1, tpr1,'r', label = 'AUC - %0.2f' % roc_auc1)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.show()

########################################################################
#python exe1.py > out.txt
output_txt =  str(home) + "/output.txt"
tabel = df_novo.head()
f = open(output_txt, 'w')
f.write("Tabela Final:" + "\n")
f.write(str(tabel))
f.write("\n" + "Accuracy: " + str(acc))
f.write("\n" + "AUC: " + str(roc_auc))
f.write("\n" + "Accuracy with balanced data: " + str(acc1))
f.write("\n" + "AUC: " + str(roc_auc1))
f.close()
