#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:19:20 2023

@author: mingyi
"""
import os
from math import floor
from math import ceil
import random
import fnmatch
import numpy as np
import subprocess
from math import log
from sklearn.ensemble import AdaBoostClassifier
from joblib import dump,load

my_env = os.environ.copy()
my_env['PATH']=os.getenv('PATH')
my_env['DATAPATH']=os.getenv('DATAPATH')


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def Seq_folder(path):
    """ find all the sequences file in one folder and creat a dictionary,
    which keys are the short name of sequences
    values are the path of seuqences files
    ex. ("/home/user/Desktop/")
    """
    Seq_folder = dict()
    Seq_list=find("*.seq",path)
    for i in Seq_list:
        Seq_name='.'.join(i.strip().split('.')[:-1])
        Seq_name=os.path.basename(Seq_name)
        Seq_folder[Seq_name]=i
    return Seq_folder

def Seq_family_list(path):
    """ find all folder's name that includes Sequence file and creat a list,
    which will be used in gernating the configure file for Decoy Finder,
    values are the path of seuqences files
    ex. ("/home/user/Desktop/")
    """
    Seq_F = list()
    Seq_list=find("*.seq",path)
    for i in Seq_list:
        Seq_name='.'.join(i.strip().split('.')[:-1])
        Seq_name=os.path.dirname(Seq_name).split('/')[-1]
        if Seq_name not in Seq_F:
            Seq_F.append(Seq_name)
    return Seq_F

class DecoyFinder_training:
    '''this class will generate all the required file for training, which both include decoy from other family and shuffle'''
    def __init__(self,file_path,output_path):
        family_name_list = Seq_family_list(file_path)

        self.svm_data_set(family_name_list,file_path,output_path,300,1)
        self.svm_data_set(family_name_list,file_path,output_path,300,2)
        self.svm_data_set(family_name_list,file_path,output_path,300,3)
        self.svm_shuffle_id_set(family_name_list,file_path,output_path,300)
        self.svm_shuffle_id_set(family_name_list,file_path,output_path,300,shuffle = 1)
        self.svm_shuffle_id_set(family_name_list,file_path,output_path,300,shuffle = 2)
        self.svm_shuffle_id_set(family_name_list,file_path,output_path,300,shuffle = 3)
    def shuffle_sequence(self,input_sequence = None, output_dir= None,shuffle_precentage=None):
        
        Sequence_path ='.'.join(input_sequence.strip().split('.')[:-1])
        Seq_name=os.path.basename(Sequence_path)
        A=open(input_sequence,'r').readlines()
        Sequence_list = []
        fire_secure = 0
        output_folder = 0
        for line in A:
            if ';' in line:
                continue
            elif Seq_name not in line:
                
                Oneline=line.split()
                Oneline=[i for i in Oneline[0]]
                if '1' in Oneline:
                    Oneline.remove('1')
                Sequence_list.extend(Oneline)
                
                #print(Oneline[0])
            
        if shuffle_precentage:
            #
            random_list=[]
            position_list = []
            shuffle_number= floor(shuffle_precentage*len(Sequence_list))
            for i in range(shuffle_number):
                random_position=random.randint(1,len(Sequence_list))-1
                while random_position in position_list:
                    random_position=random.randint(1,len(Sequence_list))-1
                random_list.append(Sequence_list[random_position])
                position_list.append(random_position)
                Sequence_list[random_position]=0
            random.shuffle(position_list)
            for i in range(len(random_list)):
                Sequence_list[position_list[i]]=random_list[i]
        else:
            random.shuffle(Sequence_list)
        try:
            os.makedirs(output_dir+Seq_name+"_random_0/")
            output_folder = output_dir+Seq_name+"_random_0/"
            sfile = open(output_folder+Seq_name+"_random_0"+".seq",'w')
            sfile.write(';\n')
        # The first non-comment line contains the sequence name
            sfile.write(Seq_name+"_random_0"+'\n')
            # Break up the sequence into 20 nucleotide chunks and write
            sequence_len = 0
            for i in range(len(Sequence_list)):
                sfile.write(Sequence_list[i])
                sequence_len +=1
                if sequence_len == 20:
                    sfile.write('\n')
                    sequence_len = 0
                if i == len(Sequence_list)-1:
                        sfile.write('\n')
            # Seq files end with a 1 character
            sfile.write('1\n')
            sfile.close()
            fire_secure = 1
        except:
            number = 1
        if fire_secure == 0:
            for i in range(0,100):
                if os.path.exists(output_dir+Seq_name+"_random_"+str(i)+"/"):
                    number = i+1
            count = 0
            while number != 0 and not os.path.exists(output_dir+Seq_name+"_random_"+str(number)+"/") and count == 0:
                try:
                    code = str(number)
                    os.makedirs(output_dir+Seq_name+"_random_"+code+"/")
                    output_folder = output_dir+Seq_name+"_random_"+code+"/"
                    sfile = open(output_folder+Seq_name+"_random_"+str(number)+".seq",'w')
                    sfile.write(';\n')
            
                    # The first non-comment line contains the sequence name
                    sfile.write(Seq_name+"_random_"+str(number)+'\n')
                    
                    # Eliminate gaps in the aligned sequence
                    
                    # Break up the sequence into 20 nucleotide chunks and write
                    sequence_len = 0
                    for i in range(len(Sequence_list)):
                        sfile.write(Sequence_list[i])
                        sequence_len +=1
                        if sequence_len == 20:
                            sfile.write('\n')
                            sequence_len = 0
                        if i == len(Sequence_list)-1:
                            sfile.write('\n')
                    # Seq files end with a 1 character
                    sfile.write('1\n')
                    sfile.close()
                    count =+1
                except:
                    number +=1
        return output_folder
    def write_file(self,Seq_folder,output_path,non_folder = None,family_name =None,processer = None ):
        """
        write_file:
                  This function can search for all the sequence file in the folder
                  and then create a config file for TurboFold
       Argument:
       input_path:the path to sequence file folder
       
       output_path:the path to output file folder
       
       non_h: a list of nonhomologouse sequence
       ex. (A,"/home/user/Desktop/,B")
       """
     #put all the sequences path in a list and nonhomolougus will have two lists
     # one for name, NH, which will be used later to write# $NH line
     #one for path, which will be used to write Seq = path...
        Seq_list =list()
        NH=list()
        NH_S=list()
        NH_K=list()
        for i in Seq_folder.values():
            Seq_list.append(i)
        
        if non_folder:
           for i in non_folder.values():
                NH_S.append(i)
           for i in non_folder.keys():
                NH.append(i)
           file_path = open(output_path+"conf.txt","w")
           HM_list=list()
           NH_list=list()
           I=0
           O=0
           #this loop will run the maximum time when number of sequences are
           #smaller than 8, or run 8 times when the number of sequences are bigger than 8
           #it will randomly pick 8 homologous and 8 nonhomologous sequences
           for i in range (ceil(len(Seq_list))):
               if I <25:
                   S = Seq_list[random.randint(0,len(Seq_list)-1)]
                   while S in HM_list:
                       S = Seq_list[random.randint(0,len(Seq_list)-1)]
               else:
                   continue
               HM_list.append(S)
               I +=1
           for i in range(ceil(len(NH_S))):
               if O <25:
                   
                   C = random.randint(0,len(NH_S)-1)
                   N = NH_S[C]
                   N2 = NH[C]
                   while N in NH_list:
                       C = random.randint(0,len(NH_S)-1)
                       N = NH_S[C]
                       N2 = NH[C]
               else:
                   continue
               O +=1
               NH_list.append(N)
               NH_K.append(N2)
               
           
           
           count=0
           Numb_seq = len(HM_list)+len(NH_list)
           
           for i in range (ceil(Numb_seq)):
               number=i+1
               if i <= (len(HM_list)-1):
                   Seq = "Seq" + str(number)
                   file_path.write(Seq + " = " + HM_list[i]+'\n')
                   Seq_name = os.path.basename(HM_list[i].rstrip("\n").split("/")[-1])[:-4]
                   Ct  ="Ct"+ str(number)
                   file_path.write(Ct +" = "+ output_path+Seq_name+".ct"+'\n')
                   Save = "Save"+str(number)
                   file_path.write(Save + " = "+output_path+Seq_name+".pfs"+'\n')
               else:
                   while count <= (len(NH_list)-1):
                       Seq = "Seq" + str(number)
                       file_path.write(Seq + " = " + NH_list[count]+'\n')
                       Seq_name = os.path.basename(NH_list[count].rstrip("\n").split("/")[-1])[:-4]
                       Ct  ="Ct"+ str(number)
                       file_path.write(Ct +" = "+ output_path+Seq_name+".cst"+'\n')
                       Save = "Save"+str(number)
                       file_path.write(Save + " = "+output_path+Seq_name+".pfs"+'\n')
                       number += 1
                       count +=1
        #this part is for only homologus sequences, and it will write out all the
           #sequences in the folder.
        else:
           Numb_seq = len(Seq_list)
        #create an empty config file
           file_path = open(output_path+"conf.txt","w")
        #write out the config file accounting to the number of sequece file found
           for i in range (Numb_seq):
                number=str(i+1)
                Seq = "Seq" + number
                file_path.write(Seq + " = " + Seq_list[i]+'\n')
            
           for i in range (Numb_seq):
                number=str(i+1)
                Seq_name = os.path.basename(Seq_list[i].rstrip("\n").split("/")[-1])[:-4]
                Ct  ="Ct"+ number
                file_path.write(Ct +" = "+ output_path+Seq_name+".ct"+'\n')
                
           for i in range (Numb_seq):
                Seq_name = os.path.basename(Seq_list[i].rstrip("\n").split("/")[-1])[:-4]
                number=str(i+1)
                Save = "Save"+number
                file_path.write(Save + " = "+output_path+Seq_name+".pfs"+'\n')
           file_path.write('StartingSaveFiles = {')
           for i in range (Numb_seq):
               Seq_name = os.path.basename(Seq_list[i].rstrip("\n").split("/")[-1])[:-4]
               file_path.write(output_path+Seq_name+"_s.pfs")
           file_path.write('}\n')
        Numb = str(Numb_seq)
        if processer:
            file_path.write('Processors ='+str(processer)+'\n')
        file_path.write("SequenceNumber ="+Numb+'\n')
        file_path.write("OutAln = " + output_path +"output.aln"+'\n')
        file_path.write("Mode = MEA"+'\n')
       
        file_path.close()
        
    def svm_shuffle_id_set(self,shuffle_name_list=None,file_path=None,out_path=None,pair_number=None,shuffle = None,loop_number= None):
        '''generate shuffle sequence based on pairwise identity'''
        seq_path=dict()
        precentage_for_sequence = dict()
        for i in range(len(shuffle_name_list)):
           seq_name = shuffle_name_list[i]
           seq_path[seq_name] = Seq_folder(file_path+seq_name+"/")
           list_seq = [a for a in seq_path[seq_name]]
           list_value = [a for a in seq_path[seq_name].values()]
           seq_pick_dict = dict()
           for m in range(30):
               seq_pick_dict[list_seq[m]]=list_value[m]
           if not os.path.exists(file_path+seq_name+"_30_homologus/"):
               os.makedirs(file_path+seq_name+"_30_homologus/")
           if not os.path.exists(file_path+seq_name+"_30_homologus/conf.txt"):
               self.write_file(seq_pick_dict,file_path+seq_name+"_30_homologus/",family_name = seq_name+"_30")
           if not os.path.exists(file_path+seq_name+"_30_homologus/total_id.txt"):
               AB = Sequence_alignment(file_path+seq_name+"_30_homologus/conf.txt")
           id_open=open(file_path+seq_name+"_30_homologus/total_id.txt","r").readlines()
           id_precentage = 1-(float(id_open[0].split("\t")[1])/100)
           precentage_for_sequence[seq_name]=id_precentage
        self.svm_shuffle_data_set(shuffle_name_list,file_path,out_path,pair_number,precentage_for_sequence,shuffle)

    def svm_shuffle_data_set(self,family_name_list=None,file_path=None,out_path = None,pick_number=None,shuffle_precent = None,shuffle_pick = None,set_pick_homo=None,Seq_taken_before = None):
        """shuffle_precent = 0 ~ 1 will be a dictionary for different RNA family
        this function will write all the configure file with decoy from shuffle"""    
        seq_path=dict()
        Family_taken_all = []
        Shuffle_dict = dict()
        
        if shuffle_precent:
            shuffle_precentage = shuffle_precent
        else:
            shuffle_precentage = random.random()
        if shuffle_pick:
            shuffle_pick_number = shuffle_pick
        else:
            shuffle_pick_number = 0
        for i in range(len(family_name_list)):
           seq_name = family_name_list[i]
           seq_path[seq_name] = Seq_folder(file_path+seq_name+"/")
           list_seq = [a for a in seq_path[seq_name]]
           for m in range(len(list_seq)):
               Shuffle_dict[list_seq[m]]=0
        
        if pick_number:
            random_number = pick_number
        else:
            random_number = 100
        if Seq_taken_before:
            Seq_been_taken=Seq_taken_before
        else:
            Seq_been_taken = dict()
            #list match the sequence used
            for i in range(len(family_name_list)):
                homo = family_name_list[i]
                Seq_been_taken[homo]=[]
        #select R-family
        for i in range(random_number):
        # need at least one family to pair
            if len(Family_taken_all) == len(family_name_list):
                print("sequence number is not enough for pairing")
                raise
            Seq_for_aligment = dict()
            fire_secure = 0
            
            random_pick_family = random.randint(0,len(family_name_list)-1)
            family_homo = family_name_list[random_pick_family]
            while family_homo in Family_taken_all:
                random_pick_family = random.randint(0,len(family_name_list)-1)
                family_homo = family_name_list[random_pick_family]
            
            homo_seq_name = [i for i in seq_path[family_homo]]
            
            
            homo_seq_list = [i for i in seq_path[family_homo].values()]
            
            homo_seq_taken = dict()
            #pick one for shuffle
            path_for_shuffle = out_path+"shuffle_sequence/"+family_homo+"/"
            shuffle_folder = []
            for i in range(shuffle_pick_number):
                random_pick_shuffle = random.randint(0,len(homo_seq_name)-1)
                shuffle_seq_picked_name = homo_seq_name[random_pick_shuffle]
                shuffle_seq_picked_path = homo_seq_list[random_pick_shuffle]
                if type(shuffle_precentage) == dict:
                    precentage = shuffle_precentage[family_homo]
                    
                    self.shuffle_sequence(shuffle_seq_picked_path,path_for_shuffle,precentage)
                else:
                    shuffle_precentage = random.random()
                    
                    while shuffle_precentage < 0.3 or shuffle_precentage > 0.7:
                        shuffle_precentage = random.random()
                        if 0.7>shuffle_precentage>0.3:
                            break
                        
                    self.shuffle_sequence(shuffle_seq_picked_path,path_for_shuffle,shuffle_precentage)
                    
                    
                Folder_random = Seq_folder(path_for_shuffle+shuffle_seq_picked_name+"_random_"+str(Shuffle_dict[shuffle_seq_picked_name])+"/")
                shuffle_folder.append(Folder_random)            
                Folder_key = [a for a in Folder_random.keys()]
                Folder_value = [b for b in Folder_random.values()]
                Seq_for_aligment[Folder_key[0]] = Folder_value[0]            
                Shuffle_dict[shuffle_seq_picked_name] += 1
            if set_pick_homo:
                set_pick_number=set_pick_homo
            else:
                homo_number = random.randint(5,20) 
                set_pick_number= homo_number
            for m in range(set_pick_number):
                seq_tried = []
                seq_picking_number = random.randint(0,len(homo_seq_name)-1)
                seq_picked = homo_seq_name[seq_picking_number]
                while seq_picked in Seq_been_taken[family_homo]:
                    seq_picking_number = random.randint(0,len(homo_seq_name)-1)
                    seq_picked = homo_seq_name[seq_picking_number]
                    seq_tried.append(seq_picked)
                    if len(seq_tried) == len(homo_seq_name):
                        if not family_homo in Family_taken_all:
                            Family_taken_all.append(family_homo)
                        break
                homo_seq_taken[homo_seq_name[seq_picking_number]]=homo_seq_list[seq_picking_number]
                Seq_for_aligment[homo_seq_name[seq_picking_number]]=homo_seq_list[seq_picking_number]
                Seq_been_taken[family_homo].append(homo_seq_name[seq_picking_number])
            if len(homo_seq_taken) == 0:
                pass
            
            path_for_aligment = out_path+"data_pairing_shuffle/"+family_homo
            number=0
            try:
                os.makedirs(path_for_aligment+"_0/")
                fire_secure = 1
                self.write_file(Seq_for_aligment,path_for_aligment+"_0/",family_name = "homo_"+family_homo+"_1_shuffled")
                Sfile = open(path_for_aligment+"_0/contamin_seq.txt","w")
                contamin_for_write = []
                for p in range(len(shuffle_folder)):
                    A = [i for i in shuffle_folder[p].values()]
                    contamin_for_write.append(A)
                for l in range(len(contamin_for_write)):
                    contamin_seq_true_name_file = open(contamin_for_write[l][0],"r").readlines()
                    Sfile.write(contamin_seq_true_name_file[1]+"\n")
                Sfile.close()
                
            except:
                number = 1
            if fire_secure == 0:
                count = 0
                for i in range(0,500):
                    if os.path.exists(path_for_aligment+"_"+str(i)+"/"):
                        number = i+1
                       
                    while number != 0 and not os.path.exists(path_for_aligment+"_"+str(number)+"/") and count == 0:
                        try:
                            code = str(number)
                            os.makedirs(path_for_aligment+"_"+code+"/")
                            output_folder = path_for_aligment+"_"+code+"/"
                            count =+1
                            self.write_file(Seq_for_aligment,output_folder,family_name ="homo_"+family_homo+"_1_shuffled")
                            
                            Sfile = open(output_folder+"contamin_seq.txt","w")
                            contamin_for_write = []
                            for p in len(shuffle_folder):
                                A = [i for i in shuffle_folder[p].keys()]
                                contamin_for_write.append(A)
                            for l in range(len(contamin_for_write)):
                                contamin_seq_true_name_file = open(contamin_for_write[l],"r").readlines()
                            Sfile.close()
                            
                        except:
                            number +=1
        return Seq_been_taken

    

    def svm_data_set(self,family_name_list=None,file_path=None,output_path = None,pick_number=None,set_pick_cont=None,set_pick_homo=None,Seq_taken_before = None):
        '''this function will write all the configure file with decoy from other families'''    
        seq_path=dict()
        Family_taken_all = []
        for i in range(len(family_name_list)):
           seq_name = family_name_list[i]
           seq_path[seq_name] = Seq_folder(file_path+seq_name+"/")
        if pick_number:
            random_number = pick_number
        else:
            random_number = 100
        if Seq_taken_before:
            Seq_been_taken=Seq_taken_before
        else:
            Seq_been_taken = dict()
            #list match the sequence used
            for i in range(len(family_name_list)):
                homo = family_name_list[i]
                for k in range(len(family_name_list)):
                    contam = family_name_list[k]
                    if contam == homo:
                        continue
                    else:
                        Seq_been_taken[homo,contam]=[]
                        Seq_been_taken[contam,homo]=[]
            #select R-family
            for i in range(random_number):
            # need at least two family to pair
                if len(Family_taken_all) == len(family_name_list)-1:
                    print("sequence number is not enough for pairing")
                    raise
                Seq_for_aligment = dict()
                fire_secure = 0
                
                random_pick_family = random.randint(0,len(family_name_list)-1)
                family_homo = family_name_list[random_pick_family]
                while family_homo in Family_taken_all:
                    random_pick_family = random.randint(0,len(family_name_list)-1)
                    family_homo = family_name_list[random_pick_family]
                random_pick_family_2 = random.randint(0,len(family_name_list)-1)
                while random_pick_family_2 == random_pick_family:
                    random_pick_family_2 = random.randint(0,len(family_name_list)-1)
                    family_cont = family_name_list[random_pick_family_2]
                    if family_cont != family_homo:
                        break
                family_cont = family_name_list[random_pick_family_2] 
                while family_cont in Family_taken_all:
                    random_pick_family_2 = random.randint(0,len(family_name_list)-1)
                    try_family = []
                    while random_pick_family_2 == random_pick_family:
                        random_pick_family_2 = random.randint(0,len(family_name_list)-1)
                        if family_name_list[random_pick_family_2] not in try_family:
                            try_family = family_name_list[random_pick_family_2]
                        if len(try_family) == len(family_name_list):
                            pass
                    
                    family_cont = family_name_list[random_pick_family_2]
                
                homo_seq_name = [i for i in seq_path[family_homo]]
                contamin_seq_name = [i for i in seq_path[family_cont]]
                homo_seq_list = [i for i in seq_path[family_homo].values()]
                contamin_seq_list = [i for i in seq_path[family_cont].values()]
                homo_seq_taken = dict()
                contamin_seq_taken = dict()
                if set_pick_homo:
                    set_pick_number=set_pick_homo
                else:
                    homo_number = random.randint(5,20) 
                    set_pick_number= homo_number
                for m in range(set_pick_number):
                    seq_tried = []
                    seq_picking_number = random.randint(0,len(homo_seq_name)-1)
                    seq_picked = homo_seq_name[seq_picking_number]
                    while seq_picked in Seq_been_taken[family_homo,family_cont]:
                        seq_picking_number = random.randint(0,len(homo_seq_name)-1)
                        seq_picked = homo_seq_name[seq_picking_number]
                        seq_tried.append(seq_picked)
                        if len(seq_tried) == len(homo_seq_name):
                            if not family_homo in Family_taken_all:
                                Family_taken_all.append(family_homo)
                            break
                    homo_seq_taken[homo_seq_name[seq_picking_number]]=homo_seq_list[seq_picking_number]
                    Seq_for_aligment[homo_seq_name[seq_picking_number]]=homo_seq_list[seq_picking_number]
                    Seq_been_taken[family_homo,family_cont].append(homo_seq_name[seq_picking_number])
                    Seq_been_taken[family_cont,family_homo].append(homo_seq_name[seq_picking_number])
                if len(homo_seq_taken) == 0:
                    pass
                if set_pick_cont:
                    set_pick_cont=set_pick_cont
                else:
                    set_pick_cont=1
                for n in range(set_pick_cont):
                    seq_tried_cont = []
                    seq_picking_number = random.randint(0,len(contamin_seq_name)-1)
                    seq_picked = contamin_seq_name[seq_picking_number]
                    while seq_picked in Seq_been_taken[family_homo,family_cont]:
                        seq_picking_number = random.randint(0,len(contamin_seq_name)-1)
                        seq_picked = contamin_seq_name[seq_picking_number]
                        seq_tried_cont.append(seq_picked)
                        if len(seq_tried_cont) == len(homo_seq_name):
                            if not family_cont in Family_taken_all:
                                Family_taken_all.append(family_cont)
                            break
                
                    contamin_seq_taken[contamin_seq_name[seq_picking_number]]=contamin_seq_list[seq_picking_number]
                    Seq_for_aligment[contamin_seq_name[seq_picking_number]]=contamin_seq_list[seq_picking_number]
                    Seq_been_taken[family_homo,family_cont].append(contamin_seq_name[seq_picking_number])
                    Seq_been_taken[family_cont,family_homo].append(contamin_seq_name[seq_picking_number])
                if len(contamin_seq_taken) == 0:
                    pass
                path_for_aligment = output_path+"data_pairing_cont/"+"homo_"+family_homo+"_cont_"+family_cont
                number=0
                try:
                    os.makedirs(path_for_aligment+"_0/")
                    fire_secure = 1
                    self.write_file(Seq_for_aligment,path_for_aligment+"_0/",family_name = "homo_"+family_homo+"_cont_"+family_cont)
                    Sfile = open(path_for_aligment+"_0/contamin_seq.txt","w")
                    contamin_for_write = [i for i in contamin_seq_taken.keys()]
                    for l in range(len(contamin_for_write)):
                        contamin_seq_true_name_file = open(contamin_seq_taken[contamin_for_write[l]],"r").readlines()
                        Sfile.write(contamin_seq_true_name_file[1]+"\n")
                    Sfile.close()
                    
                except:
                    number = 1
                if fire_secure == 0:
                    count = 0
                    for i in range(0,500):
                        if os.path.exists(path_for_aligment+"_"+str(i)+"/"):
                            number = i+1
                           
                        while number != 0 and not os.path.exists(path_for_aligment+"_"+str(number)+"/") and count == 0:
                            try:
                                code = str(number)
                                os.makedirs(path_for_aligment+"_"+code+"/")
                                output_folder = path_for_aligment+"_"+code+"/"
                                count =+1
                                self.write_file(Seq_for_aligment,output_folder,family_name ="homo_"+family_homo+"_cont_"+family_cont)
                                
                                Sfile = open(output_folder+"contamin_seq.txt","w")
                                contamin_for_write = [i for i in contamin_seq_taken.keys()]
                                for l in range(len(contamin_for_write)):
                                    contamin_seq_true_name_file = open(contamin_seq_taken[contamin_for_write[l]],"r").readlines()
                                    Sfile.write(contamin_seq_true_name_file[1]+"\n")
                                Sfile.close()
                                
                            except:
                                number +=1
            return Seq_been_taken
                                

class Sequence_alignment:
    """This class can use config file to create the dictionary for Sequence 
    name and match_score for Sequence
    (only part for pairwise identity was used, rest are not used)
    Argument:
             conf_file: the path for the config file
    Returns: Sequence name dictionary
             base pair probabilites
             Sequence aligment
             match_score for two sequence
             match_score for randomly picked sequence"""
    
    def __init__(self, conf_file=None, rand_number = None,smp = None):
        self.Seqname_dict=None#name dictionary
        self.probabilities=None#base pair probailites
        self.NH = None#list of nonhomologous sequences
        self.algn=None#Sequence aligment
        self.match_score=None#match score for two sequence
        self.rand_s=None#match_score for randomly picked sequence
        self.EX = None#export data path
        self.HS =None#export data path of histogram
        self.RS = None
        self.FS = None
#        self.features= None
        self.family_name = None
        self.FG = None
        self.SP = None
       # self.FH = None
        self.match_score_bigger_than_one = None
        self.alpha_1 = None
        self.alpha_2 = None
        
        self.seq_id = None
        if conf_file:
            self.seq_id = dict()
            if smp:
                self.ParseTFconfig(conf_file,smp)
            else:
                self.ParseTFconfig(conf_file)
            if self.NH:
                if rand_number:
                    self.rand_match(rand_number,self.NH)
                else:
                    self.rand_match(2000,self.NH)
            else:
                if rand_number:
                    self.rand_match(rand_number,self.NH)
                else:
                    self.rand_match(2000,self.NH)
            
            self.export_data()
    def ParseTFconfig(self, conf_file,smp=None):
        TFconf=open(conf_file,"r")
        self.Seqname_dict=dict()
        self.probabilities=dict()
        for line in TFconf:
            words=line.split('=')
            file_type=words[0].strip()
            file_path=words[1].strip()
            Seq_name='.'.join(line.strip().split('.')[:-1])
            Seq_name=os.path.basename(Seq_name)

            
         
            if file_type[:3] == 'Seq':
                if os.path.exists(file_path):
                    self.Seqname(line)
                    #update name dictionarysum(npoj.pp_matrix[:,base1])
                  
                else:
                    if file_type == 'SequenceNumber':
                        continue
                    else:
                    #error
                        print ("no seqence file")
            elif file_type == 'OutAln':
                # import alignment
                self.read_algn(file_path)
            elif file_type[:4] == 'Save':
                #import pair probabilities
                pp_file = (".").join(file_path.split(".")[:-1])+".pp"
                if not os.path.exists(file_path):
                    if not os.path.exists(pp_file):
                        if smp:
                            subprocess.call(["TurboFold-smp",conf_file],env=my_env)
                        else:   
                            subprocess.call(["TurboFold",conf_file],env=my_env)
                    else:
                        self.probabilities[self.Seqname_dict[Seq_name]]=PP_Obj(pp_file)
                else:
                    self.probabilities[self.Seqname_dict[Seq_name]]=PP_Obj(file_path)
#                elif os.path.exists(pp_file):
#                    self.probabilities[self.Seqname_dict[Seq_name]]=PP_Obj(pp_file)
            # # $NH in config file can determine which sequences are nonhomologus
            # and then put all the nonhomologus sequences into a list with the actual name
            elif file_type == '# $NH':
                NHL = file_path.split()
                self.NH = list()
                for i in NHL:
                    self.NH.append(self.Seqname_dict[i])

            elif file_type == '# $EX':
                # # $EX = path in config file can determine export data's saving path
                # or it will save to Desktop
                self.EX = file_path
                  #the savingfile should have the sequence name the same as the 
                #real sequence name in config file, Or this will not work
                #due to the key errors of dicitonary 
            elif file_type == '# $HS':
                 
                 self.HS = file_path
            elif file_type == '# $RS':
                 
                 self.RS = file_path
            elif file_type == '# $FS':
                 
                 self.FS = file_path
            elif file_type == '# $SP':
         
                 self.SP = file_path
            elif file_type == '# $FG':
                
                 self.FG = file_path
                 self.family_name = Seq_name
            #elif file_type == '# $FH':
                
                # self.FH = file_path
                 
        
    def read_algn(self,filename,name_dict=None):
            """This function is used to read the output sequence alignment, which is
        in ClustalW format, then create a dictionary to record name of HM sequences 
        and the sequenceself.EX
            Argument:
                 conf_file: the path for the config file
            Returns:
                 algn: a dictionary of homlogous sequences
        
        """
        
            A=open(filename,'r').readlines()
        
            self.algn=dict()
        
            for line in A:
                Oneline=line.split()
        
            #print(Oneline[0])
                if len(Oneline)!=2:
                    continue
            
                if name_dict:
                    name = name_dict(Oneline[0])
                    
                else:
                    name= Oneline[0]
                    
                if name in self.algn.keys(): 
                    self.algn[name] = self.algn[name] + Oneline[1]
                else:
                    self.algn[name] = Oneline[1]
    
    #        self.name_dict_r = dict()
    #        for key, value in self.name_dict:
    #            self.name_dict_r[value]=key
                
            return self.algn

        
   
            
    def Seqname(self,line):
        """ Using config file's infromation to update the dictionary 
         of sequence name and short name for sequence file
         Argument:
             conf_file: the path for the config file
         Returns:
             Seqname dictionary: keys[short name, such as 5s2_1] 
                                 values[sequence name, such as _____ARTHROB._GLOBIFORMIS_ATCC8010_(GRAM___ACTINOBACTER]
        """
        
        
#        MST=open(conf_file,'r').readlines()
        import re
    

#        Seq_name=os.path.basename(line.rstrip("\n").split(".")[0])
        Seq_name = '.'.join(line.strip().split('.')[:-1])
        Seq_name=os.path.basename(Seq_name)
        #get the short name for sequence, which help to identify the sequence file
        Seq_path= line.rstrip("\n").split("=")[1].strip()
        #get the sequence file path from the config file and open the sequence file
        #from sequence file, we can get the real name of the sequence
        
        Seq_file=open(Seq_path,'r')
        try:
            line = Seq_file.readline()    
        except UnicodeDecodeError:
            print(line)
            
        while line:
        #identify the first line of the sequence file, which is not what we want
            if line[0]!=';':
                break
            else:
                line = Seq_file.readline()    
        
                
        Seq_name2=line
        Seq_name2 = Seq_name2.replace('+','_')
        Seq_name2 = Seq_name2.replace(',','_')
        Seq_name2 = re.sub('\s','_',Seq_name2)
        Seq_name2 = Seq_name2[:-1]
                    
            
        self.Seqname_dict[Seq_name]=Seq_name2

    

    def random_PairSeq(self,NH=[]):
        """ This function can randomly pick two sequences for random match score
            calculation from the aligment dictionary
            Argument: NH_seq is an optional input, which is an list of nonhomologous
                      sequence
            Returns: two random sequences
        """
        #put all the sequence name from alignment dicitonary into a list, seq_list
        seq_list = [i for i in self.algn.keys()]
        #
        if len(NH)>0:
            seq1 = seq_list[random.randint(0,len(seq_list)-1)]
            
            while seq1 in NH:
                seq1 = seq_list[random.randint(0,len(seq_list)-1)]
                
            seq2 = NH[random.randint(0,len(NH)-1)]
            while seq1 == seq2:
                seq1 = random.randint(0,len(seq_list)-1)
                seq1 = seq_list[seq1]
            
        else:
            for i in range (2):
                # range can determine how many sequence to pick
                # for match socre calculation, we only need two
                # so seq1we use range (2) here
                seq1 = seq_list[random.randint(0,len(seq_list)-1)]#pick sequence randomly, len(seq_list) need to -1, if use 0 to len(seq_list). 
                seq2 = seq_list[random.randint(0,len(seq_list)-1)]#if use 1 to len(seq_list), len(seq_list) do not need to -1, or it will never count the last sequence into the calculation
                while seq1==seq2:
                    # if seq1 is same as seq2, then this will randomly pick another seq2
                    seq2 = seq_list[random.randint(0,len(seq_list)-1)]
        
        if seq1 < seq2:
            #this part is make sure that the same pair is not picked twice
            #the seq1 and seq2 are now only numbers, such as seq1 =1, seq2 =2
            #so we can use if loop to make sure there will not a situation like
            # 1,2 and 2,1 pair appeared, since they will only be 1,2 in set or list.
            return seq1, seq2
        else:
            return seq2, seq1
                
    def match_soc(self,seq1,seq2):
        """match_soc
        this function will use two sequences to calculate the match socres, by using
        the equation: i = a nucleotide at position i, k the same 
                    P(i,k)= a1 * square root of (i pair up * k pair up + i pair down * k pair down)+ a2* square root(i unpair *k unpair)+a3
        α1= 1.0, α2 = 0.8, α3 = 0.5 were used as defaults for TurboFold II, which
        are the optimal parameters.
        
        Arguments:
            seq1 {str}:
                the sequence randomly picked by the random pick function, or
                manually given
                
            seq2 {str}:
                 the sequnece randomly picked by the random pick function, or
                manually given, make sure not the same as the sequence 1
            
        Returns: rand_s[keys: name of two sequences]=(value: match score of positions)
            match_score {dictioanry [float]}: the match socre will return as a
            dictionary with the two sequences as keys and match socre of each 
            position as the values.
        
                
        """
        # Initialize match score vector with zeros
        match_score= [0 for i in range(len(self.algn[seq1]))]
        alpha_1= [0 for i in range(len(self.algn[seq1]))]
        alpha_2= [0 for i in range(len(self.algn[seq1]))]
        alpha_1_seq1_up= [0 for i in range(len(self.algn[seq1]))]
        alpha_1_seq2_up= [0 for i in range(len(self.algn[seq1]))]
        alpha_1_seq1_down= [0 for i in range(len(self.algn[seq1]))]
        alpha_1_seq2_down= [0 for i in range(len(self.algn[seq1]))]
        alpha_2_seq1_unpair= [0 for i in range(len(self.algn[seq1]))]
        alpha_2_seq2_unpair= [0 for i in range(len(self.algn[seq1]))]
        alpha_1_pair=[0 for i in range(len(self.algn[seq1]))]
        alpha_2_unpair=[0 for i in range(len(self.algn[seq1]))]
        #n_position = []
        # get the pairwise by function get_pairwis
        Pairwise_algn, Pair_i = self.get_pairwise(seq1,seq2)
        Position_seq_1 = [ k for k in self.algn[seq1]]
        Position_seq_2 = [ k for k in self.algn[seq2]]
        total_id = 0
        match_id = 0
        for i in range(len(Position_seq_1)):
            P1 = Position_seq_1[i]
            P2 = Position_seq_2[i]
            if P1 != "-":
                if P2 !="-":
                    if P1 == P2:
                        total_id+=1
                        match_id+=1                        
                    elif P1 != P2:
                        total_id+=1
        id_precentage = match_id/total_id*100
        pair_name = str(seq1)+"_"+str(seq2)            
        self.seq_id[pair_name]=id_precentage
        # use i to get all the position of sequence 
        for i in range(len(self.algn[seq1])):
            if Pairwise_algn[seq1][i] != 0:
                a = Pairwise_algn[seq1][i]-1
               #n_position.append(Pairwise_algn[seq1][i])
        #for #i in range(len(n_position)):
                #i = n_position[i]
                try:
                        #make sure the we are puting the aligned base into calculation 
                        #match_score[Pair_i[seq1][i]]= 1.0*sqrt(
                    try:
                        self.probabilities[seq1].prob_up(Pairwise_algn[seq2][a])
                    except:
                        print(i)
                        print("seq1")
                    try:
                        self.probabilities[seq2].prob_up(a)
                    except:
                        print("seq2")
                        print(i)
                    
                    M1 = sqrt(self.probabilities[seq1].prob_up(Pairwise_algn[seq2][a]-1)*self.probabilities[seq2].prob_up(a))
                    M1 += sqrt(self.probabilities[seq1].prob_down(Pairwise_algn[seq2][a]-1)*self.probabilities[seq2].prob_down(a))
                    M1 = 1.0*M1
                    alpha_1_seq1_up[i]=self.probabilities[seq1].prob_up(Pairwise_algn[seq2][a]-1)
                    alpha_1_seq2_up[i]=self.probabilities[seq2].prob_up(a)
                    alpha_1_seq1_down[i]=self.probabilities[seq1].prob_down(Pairwise_algn[seq2][a]-1)
                    alpha_1_seq2_down[i]=self.probabilities[seq2].prob_down(a)
                    alpha_2_seq1_unpair[i]=self.probabilities[seq1].unpair(Pairwise_algn[seq2][a]-1)
                    alpha_2_seq2_unpair[i]=self.probabilities[seq2].unpair(a)
                    alpha_1_pair[i]=M1
                    A2 = sqrt(self.probabilities[seq1].unpair(Pairwise_algn[seq2][a]-1)*self.probabilities[seq2].unpair(a))
                    alpha_2_unpair[i]=A2                    
                    if M1 != 0 or A2 != 0:
                        alpha_1[i] = M1
                        alpha_2 [i] = A2
                    M1 += 0.8*sqrt(self.probabilities[seq1].unpair(Pairwise_algn[seq2][a]-1)*self.probabilities[seq2].unpair(a))
                    
                    if M1 != 0:
                        match_score[i] = M1
                        match_score[i]+=0.5
                        
                except KeyError:
                    print('missing seq')

       
        for i in [_ for _ in range(len(alpha_1))][::-1]:
          if Pairwise_algn[seq1][i] == 0:
              if Pairwise_algn[seq2][i] == 0:
                  del alpha_1[i] 
                  del alpha_2[i]
        

#            
        return match_score, alpha_1, alpha_2

    def rand_match(self, count,NH=None):
        """
        """
        self.rand_s = dict()
        self.features = dict()
        self.alpha_1 = dict()
        self.alpha_2 = dict()
        seq_list = [i for i in self.algn.keys()]
        if NH:
            rand_m = set()
            #initialize a dictionary for match score
           
            
            Max_count= len(NH)*(len(seq_list)-len(NH))
            if count > Max_count:
                count = Max_count
            while len(rand_m) < count:
                rand_m.update([self.random_PairSeq(NH)])
            for seq in rand_m:
                match_score, alpha_1, alpha_2 = self.match_soc(seq[0],seq[1])
                self.rand_s[seq]= match_score
                self.alpha_1[seq] = alpha_1
                self.alpha_2[seq] = alpha_2
#                self.features[seq]=self.get_features(seq[0],seq[1])
        else:
            N= len(seq_list)
            Max_count = int(N*(N-1)/2)
            if count >Max_count:
                count = Max_count
            #initialize a set for random picked sequences
            rand_m = set()
            #initialize a dictionary for match score
            
            #the while loop will pick 100 randomed sequence pairs
            #and put the sequence into the set
            while len(rand_m) < count:
                rand_m.update([self.random_PairSeq()])
            # for loop will calculate all the sequence pairs and put them into dicitonary
            for seq in rand_m:
                match_score, alpha_1, alpha_2 = self.match_soc(seq[0],seq[1])
                self.rand_s[seq]= match_score
                self.alpha_1[seq] = alpha_1
                self.alpha_2[seq] = alpha_2
#                self.features[seq]=self.get_features(seq[0],seq[1])
        

# Step 4: Write a function that will return any pairwise sequence alignment from
#         the multiple sequence alignment.  The function should be of the form 
#         paiSeq1_c.append(j)rwise = get_pairwise(MSA, seq1, seq2).  The output pairwise should 
#         be a dictionary such that pairwise[seq1][position1] returns the 
#         nucleotide position in sequence 2 that aligns with the nucleotide at
#         position1 in sequence 1.

    def get_pairwise(self, seq1, seq2):
        """ seq1 is the name of frist sequence,seq 2 is the name of second sequence. 
        The function randomly pick two sequence from
        the aliment file and compare the alignment.
        Arguments:
            seq1 {str}:
                the sequence randomly picked by the random pick function, or
                manually given
                
            seq2 {str}:
                 the sequnece randomly picked by the random pick function, or
                manually given, make sure not the same as the sequence 1
            
        Returns: Pairwise: a dicitonary have the sequence 1 's aligned 
                           base pair's position in seuqence 2. so as the sequence 2
                 Pair_i:   a dicitonary have the map of sequence 1 and seuqence 2.
                           No - in the sqeuence.
        """
        #initialize for Pairwise and Pair_i
        Pairwise=dict()
        Pair_i=dict()
        #Error_Check, will report if the seq1 and seq2 are exist or not
        try:
            Seq1_a=self.algn[seq1]
        except:
            print("no seq1")
            return
        try:
            Seq2_a=self.algn[seq2]
        except:
            print("no seq2")
            return
        #get rid of the space at the end so that we can have the true sequence
        Seq1_b= (Seq1_a.strip())
        Seq2_b=(Seq2_a.strip())
        #initialize a list for poisitions of basepair in sequence
        j=0
        k=0
        Seq1_c=list()
        Seq2_c=list()
        Seq1_m=list()
        Seq2_m=list()
        #the for loop will go thourgh all the position of the sequence
        #and put the aligned base pair into the dictionary
        #so as the poisiton information for each sequence to Pair_i
        for i in range(len(Seq1_b)):
            if Seq1_b[i]!='-':
                j+=1
               
            if Seq2_b[i]!='-':
                k+=1
                
            if Seq1_b[i] =='-':
                   Seq1_c.append(0)
            else:
                if Seq2_b[i]!='-':
                    Seq1_c.append(k)
                    
                else:
                    Seq1_c.append(0)
                    
    
                Seq1_m.append(i)
                
            if Seq2_b[i] =='-':
                    Seq2_c.append(0)
                    
            else:
                if Seq1_b[i]!='-':
                    Seq2_c.append(j)
                   
                else:
                    Seq2_c.append(0)
                   
                Seq2_m.append(i)
    #        else:
    #            Seq1_c.append(0)
    #        if Seq1_b[i] == Seq2_b[i] !='-':
    #            Seq2_c.append(j)
    #        else:
    #            Seq2_c.append(0)
        #print (Seq1_cPair_)
        #print (Seq2_c)
    #    Pairwise[seq1,seq2] = Seq2_c
        Pairwise[seq2] = Seq2_c
        Pairwise[seq1] = Seq1_c
        Pair_i[seq2] = Seq2_m
        Pair_i[seq1] = Seq1_m
        
        return Pairwise, Pair_i   
    def export_data(self,file_path=None):
        X = [i for i in self.rand_s.values()]
        #Y = [i for i in self.alpha_1.values()]
        #Z = [i for i in self.alpha_2.values()]
        #Y = Y.extend
        #Y = np.array(Y,dtype =float)
        #Z = np.array(Z,dtype =float)
        #Z_value, bins_H, bins_V = heatmap_count(self.alpha_1.values,self.alpha_2.values,0,1.0,20)
        #bins_H = X , bins_V =Y
        id_value = [i for i in self.seq_id.values()]
        Total_average_id = 0
        for i in range(len(id_value)):
            Total_average_id += id_value[i]
        Total_average_id = Total_average_id/len(id_value)
        K = self.export_his(X)
        if file_path:
#            safe_save(file_path+"Match_data.csv",X,delimiter=",")
            np.savetxt(file_path+"his.csv",K,delimiter=",")
            #safe_save(file_path+"alpha_1.csv",Y,delimiter=",")
            #safe_save(file_path+"alpha_2.csv",Z,delimiter=",")
            family_name = "Test"
            plot_match_score_histogram(K[0], K[1], family_name,file_path+family_name+".png")
            #heatmap_plot(bins_H,bins_V,Z_values, family_name, file_path+family_name+"_Heatmap.png")
            self.export_rands(file_path+"rand.txt",self.rand_s)
#            self.export_features(self.features,file_path+"features.txt")
        else:
#            if self.EX:
#              safe_save(self.EX,X,delimiter=",")
                #A_path = '.'.join((self.EX).strip().split('.')[:-1])
                #safe_save(A_path+"_alpha_1.csv",Y,delimiter=",")
                #safe_save(A_path+"_alpha_2.csv",Z,delimiter=",")
#            else:
               # safe_save("/home/mzhu/Desktop/Match_data_alpha_1.csv",Y,delimiter=",")
                #safe_save("/home/mzhu/Desktop/Match_data_alpha_2.csv",Z,delimiter=",")
#               safe_save("/home/mzhu/Desktop/Match_data.csv",X,delimiter=",")
            if self.HS:
                np.savetxt(self.HS,K,delimiter=",")
            else:
                np.savetxt("/home/mzhu/Desktop/his.csv",K,delimiter=",")
#            if self.RS:
#                self.export_rands(self.RS,self.rand_s)
#            else:
#                self.export_rands("/home/mzhu/Desktop/rand.txt",self.rand_s)
            if self.FS:
#               
                ID_path = os.path.dirname(self.FS)
                with open(ID_path+"/total_id.txt",'w') as f:
                    f.write("Total_average_pair_wise_identity"+"\t"+str(Total_average_id))
                f.close()
            else:
#                self.export_features(self.features,"/home/mzhu/Desktop/features.txt")
                with open("/home/mzhu/Desktop/total_id.txt",'w') as f:
                    f.write("Total_average_pair_wise_identity"+"\t"+str(Total_average_id))
                f.close()

    def export_his(self,data=None,lower = None, upper=None,bins=None):
        if lower:
            lower = lower
        else:
            lower=0
        if upper:
            upper = upper
        else:
            upper = 1.6
        if bins:
            bins= bins
        else:
            bins = 32

        bin_width = (upper - lower)/bins
    
        bin_points = []
        counts = []
        precent = []
        span = upper - lower
        AlMatch = 0
        # Find the boundary for each bin and initialize the counts list
        for i in range(bins):
            bin_points.append(lower+ i*bin_width)
            counts.append(0)
        # Add the last bin boundary
        bin_points.append(lower+ bins*bin_width)
        
        # Distribute the counts to the bins
        for item in data:
            # Calculate the bin number
            for i in item:
                index = floor((i-lower)/span*bins)
            
            # Add a count to the bin
                if i == upper:
                    counts[index-1] += 1
                else:
                    counts[index] += 1
        
        # Find the center of each bin
        bins = [(bin_points[i]+bin_points[i+1])/2 for i in range(len(bin_points)-1)]    
        
        for i in range(len(counts)):
            if i != 0:
                AlMatch += counts[i]
        for i in range(len(bins)):
            if i == 0:
                precent.append(0)
            else:
                precent.append(counts[i]/AlMatch)
                
        return bins,precent,counts
        
        
       
        
    def export_rands(self,path = None,rand_s = None):
        if not os.path.exists(path):
            file = open(path,"w")
        else:
            file = open("/home/mzhu/Desktop/randmatchsave.txt","w")
        if rand_s:
            N = [i for i in rand_s.keys()]
            V = [i for i in rand_s.values()]
            for i in range(len(N)):
                name = str(N[i])
                value = str(V[i])
                file.write(name +" = "+value+"\n")
            file.close()
        else:
            print ("no data")
            
    def match_his(self,match_soc = None, lower= None, upper = None, bins = None):
        if lower:
            lower = lower
        else:
            lower=0
        if upper:
            upper = upper
        else:
            upper = 1.6
        if bins:
            bins= bins
        else:
            bins = 32

        
        bin_width = (upper - lower)/bins
    
        bin_points = []
        counts = []
        precent = []
        span = upper - lower
        
        # Find the boundary for each bin and initialize the counts list
        for i in range(bins):
            bin_points.append(lower+ i*bin_width)
            counts.append(0)
        # Add the last bin boundary
        bin_points.append(lower+ bins*bin_width)
        
        # Distribute the counts to the bins
        for item in match_soc:
            index = floor((item-lower)/span*bins)
        
        # Add a count to the bin
            if item == upper:
                counts[index-1] += 1
            else:
                counts[index] += 1
        
        # Find the center of each bin
        bins = [(bin_points[i]+bin_points[i+1])/2 for i in range(len(bin_points)-1)]    
            
        
        for i in range(len(counts)):
            if i == 0:
                AlMatch = 0
            else:
                AlMatch += counts[i]
        for i in range(len(bins)):
            if bins[i] == bins[0]:
                precent.append(0)
            else:
                if AlMatch == 0:
                    print(list(self.rand_s.keys())[list(self.rand_s.values()).index(match_soc)])
                else:
                    precent.append(counts[i]/AlMatch)
                
        return bins,precent,counts
        
            
    def export_test_his (self,file_path=None,lower=None,upper=None,bins=None):
        X = [i for i in self.rand_s.values()]
        X = np.array(X)
        if lower:
            lower = lower
        else:
            lower=0
        if upper:
            upper = upper
        else:
            upper = 1.6
        if bins:
            bins= bins
        else:
            bins = 32
        
        K = self.export_his(X,lower,upper,bins)
        
        num = str(bins)
        path = file_path + "test_his"+ num +".cvs"
        np.savetxt(path,K,delimiter=",")
        

    
def safe_save(path, array, **kwargs):
    if not os.path.exists(path):
        np.savetxt(path, array, **kwargs)
    else:
        print('File already exists!')
        raise


    

def count_bins(data, lower, upper, bins):
    """count_bins
    Generate counts for histogram output
    
    Arguments:
        data {iterable numeric type}:
            A iterable object that produces numeric types.  Can be an array.
        
        lower {numeric}:
            A numeric type that specifies the lower bounds for the bins
            
        upper {numeric}:
            A numeric type that specifies the upper bounds for the bins
            
        bins {int}:
            A integer specifying the number of bins to put the counts in
            
    Returns:
        counts {list{int}}:
            A list of counts for the bins
            
        bins {list{float}}:
            A list of the center point of each bin            
    """
    bin_width = (upper - lower)/bins

    bin_points = []
    counts = []
    span = upper - lower
    
    # Find the boundary for each bin and initialize the counts list
    for i in range(bins):
        bin_points.append(lower+ i*bin_width)
        counts.append(0)
    # Add the last bin boundary
    bin_points.append(lower+ bins*bin_width)
    
    # Distribute the counts to the bins
    for i in data:
        # Calculate the bin number
            index = floor((i-lower)/span*bins)
        
        # Add a count to the bin
            if i == upper:
                counts[index-1] += 1
            else:
                counts[index] += 1
    
    # Find the center of each bin
    bins = [(bin_points[i]+bin_points[i+1])/2 for i in range(len(bin_points)-1)]    
    
    return counts, bins


class DecoyFinder_feature_model_trainning:
    """This class can use config file to create the dictionary for Sequence in
    configure file of TurboFold
    all feature will be stored in self.hyper
    Argument:
             conf_file: the path for the config file, and create a Sequence name dictionary with features, includes KLscore, Z_score,
    aligment shann, individual shann,default shann,difference shann to identify the decoy sequence
    Returns:  decoy sequence """
    def __init__(self, conf_file=None,smp = None):
        self.Seqname_dict=dict()#name dictionary
        self.probabilities=dict()#base pair probailites
        self.algn=dict()#Sequence aligment
        self.algn_to_Seqname_file = dict()# a dictionary of number of file to name of seqeunces to get alignment information
        self.alpha_1 = dict()#PS Score
        self.alpha_2 = dict()#US Score
        self.Z_value = dict()#Z_value dictionary will be used to calculate KLscore
        self.heatmap = dict()#heatmap data just call the sequence name, it will show heatmap data for each sequence
        self.hyper = []#hyper feature return
        
        
        #load conf file
        if conf_file:
            if smp:
                self.ParseTFconfig(conf_file,smp)
            else:
                self.ParseTFconfig(conf_file)
            
            self.all_match()#get all the PS and US Score and calculate Z value
            self.decoy_identify()#add all the feature
    
    def ParseTFconfig(self, conf_file,smp=None):
        TFconf=open(conf_file,"r")
        
        for line in TFconf:
            words=line.split('=')
            if len(words) < 2:
                continue
            else:
                file_type=words[0].strip()
                file_path=words[1].strip()


#            if file_type[:3] == 'Seq':
#                if os.path.exists(file_path):
#                    if file_type[3] in self.Seqname_dict.keys():
#                        self.Seqname_dict[file_type[3]].update({'Sequence_file':file_path})
#                    else:
#                        self.Seqname_dict[file_type[3]] = {'Sequence_file':file_path}
#                else:
#                    if file_type == 'SequenceNumber':
#                        continue
#                    else:
#                    #error
#                        print ("no seqence file")
                if file_type[:2] == 'Ct':
                    Ct_number = file_type.split('=')[0]
                    Ct_number = Ct_number.split('Ct')[1]
                    if os.path.exists(file_path):
                        if Ct_number in self.Seqname_dict.keys():
                            self.Seqname_dict[Ct_number].update({'CT_file':file_path})
                        else:
                            self.Seqname_dict[Ct_number] = {'CT_file':file_path}
                    else:
                        if smp:
                                subprocess.call(["TurboFold-smp",conf_file],env=my_env)
                                
                                if Ct_number in self.Seqname_dict.keys():
                                    self.Seqname_dict[Ct_number].update({'CT_file':file_path})
                                else:
                                    self.Seqname_dict[Ct_number] = {'CT_file':file_path}
                        else:   
                                subprocess.call(["TurboFold",conf_file],env=my_env)
                                
                                if Ct_number in self.Seqname_dict.keys():
                                    self.Seqname_dict[Ct_number].update({'CT_file':file_path})
                                else:
                                    self.Seqname_dict[Ct_number] = {'CT_file':file_path}
                elif file_type == 'OutAln':
                    # import alignment
                    if os.path.exists(file_path):
                        self.read_algn(file_path)
                        self.algn_path = file_path
                    else:
                        if smp:
                                subprocess.call(["TurboFold-smp",conf_file],env=my_env)
                                if os.path.exists(file_path):
                                    self.read_algn(file_path)
                                    self.algn_path = file_path
                        else:   
                                subprocess.call(["TurboFold",conf_file],env=my_env)
                                if os.path.exists(file_path):
                                    self.read_algn(file_path)
                                    self.algn_path = file_path
                elif file_type[:4] == 'Save' and file_type != 'SaveFiles' :
                    #import pair probabilities
                    Save_number = file_type.split('=')[0]
                    Save_number = Save_number.split('Save')[1]
                    pp_file = (".").join(file_path.split(".")[:-1])+".pp"
                    if os.path.exists(file_path):
                        
                        if Save_number in self.Seqname_dict.keys():
                            self.Seqname_dict[Save_number].update({'Save_file':file_path})
                        else:
                            self.Seqname_dict[Save_number] = {'Save_file':file_path}
                    if not os.path.exists(file_path):
                        if not os.path.exists(pp_file):
                            if smp:
                                subprocess.call(["TurboFold-smp",conf_file],env=my_env)
                                
                                if Save_number in self.Seqname_dict.keys():
                                    self.Seqname_dict[Save_number].update({'Save_file':file_path})
                                else:
                                    self.Seqname_dict[Save_number] = {'Save_file':file_path}
                            else:   
                                subprocess.call(["TurboFold",conf_file],env=my_env)
                                
                                if Save_number in self.Seqname_dict.keys():
                                    self.Seqname_dict[Save_number].update({'Save_file':file_path})
                                else:
                                    self.Seqname_dict[Save_number] = {'Save_file':file_path}
                        else:
                            self.probabilities[Save_number]=PP_Obj(pp_file)
                            
                            if Save_number in self.Seqname_dict.keys():
                                self.Seqname_dict[Save_number].update({'Save_file':file_path})
                            else:
                                self.Seqname_dict[Save_number] = {'Save_file':file_path}
                    else:
                        self.probabilities[Save_number]=PP_Obj(file_path)
                        
                        if Save_number in self.Seqname_dict.keys():
                            self.Seqname_dict[Save_number].update({'Save_file':file_path})
                        else:
                            self.Seqname_dict[Save_number] = {'Save_file':file_path}
        
        
        #            elif file_type == 'InSeq':
        #                file_path = file_path.replace('{','')
        #                file_path = file_path.replace('}','')
        #                file_path_list = file_path.split(';')
        #                for seqs_number in range(len(file_path_list)):
        #                    if not str(seqs_number+1) in self.Seqname_dict.keys():
        #                        self.Seqname_dict[str(seqs_number+1)] = {'Sequence_file':file_path_list[seqs_number].rstrip()}
        #                    else:
        #                        self.Seqname_dict[str(seqs_number+1)].update({'Sequence_file':file_path_list[seqs_number].rstrip()})
                elif file_type == 'OutCT':
                    file_path = file_path.replace('{','')
                    file_path = file_path.replace('}','')
                    file_path_list = file_path.split(';')
                    try:
                        file_path_list.remove('')
                    except:
                        pass
                    for seqs_number in range(len(file_path_list)):
                        if os.path.exists(file_path_list[seqs_number]):
                            if not str(seqs_number+1) in self.Seqname_dict.keys():
                                self.Seqname_dict[str(seqs_number+1)] = {'CT_file':file_path_list[seqs_number].rstrip()}
                            else:
                                self.Seqname_dict[str(seqs_number+1)].update({'CT_file':file_path_list[seqs_number].rstrip()})
                        else:
                            if smp:
                                subprocess.call(["TurboFold-smp",conf_file],env=my_env)
                                if not str(seqs_number+1) in self.Seqname_dict.keys():
                                    self.Seqname_dict[str(seqs_number+1)] = {'CT_file':file_path_list[seqs_number].rstrip()}
                                else:
                                    self.Seqname_dict[str(seqs_number+1)].update({'CT_file':file_path_list[seqs_number].rstrip()})
                            else:   
                                subprocess.call(["TurboFold",conf_file],env=my_env)
                                if not str(seqs_number+1) in self.Seqname_dict.keys():
                                    self.Seqname_dict[str(seqs_number+1)] = {'CT_file':file_path_list[seqs_number].rstrip()}
                                else:
                                    self.Seqname_dict[str(seqs_number+1)].update({'CT_file':file_path_list[seqs_number].rstrip()})
                
                elif file_type == 'SaveFiles':
                    file_path = file_path.replace('{','')
                    file_path = file_path.replace('}','')
                    file_path_list = file_path.split(';')
                    try:
                        file_path_list.remove('')
                    except:
                        pass
                    for seqs_number in range(len(file_path_list)):
                        if not str(seqs_number+1) in self.Seqname_dict.keys():
                            self.Seqname_dict[str(seqs_number+1)] = {'Save_file':file_path_list[seqs_number].rstrip()}
                        else:
                            self.Seqname_dict[str(seqs_number+1)].update({'Save_file':file_path_list[seqs_number].rstrip()})
                        path_file = file_path_list[seqs_number]
                        
                        if not os.path.exists(path_file):
                            pp_file = (".").join(path_file.split(".")[:-1])+".pp"
                            if not os.path.exists(pp_file):
                                if smp:
                                    subprocess.call(["TurboFold-smp",conf_file],env=my_env)
                                else:   
                                    subprocess.call(["TurboFold",conf_file],env=my_env)
                            else:
                                self.probabilities[str(seqs_number+1)]=PP_Obj(pp_file)
                        else:
                            self.probabilities[str(seqs_number+1)]=PP_Obj(path_file)    
                
                elif file_type == 'StartingSaveFiles':
                    file_path = file_path.replace('{','')
                    file_path = file_path.replace('}','')
                    file_path_list = file_path.split(';')
                    try:
                        file_path_list.remove('')
                    except:
                        pass
                    for seqs_number in range(len(file_path_list)):
                        if os.path.exists(file_path_list[seqs_number]):
                            if not str(seqs_number+1) in self.Seqname_dict.keys():
                                self.Seqname_dict[str(seqs_number+1)] = {'StartingSaveFiles':file_path_list[seqs_number].rstrip()}
                            else:
                                self.Seqname_dict[str(seqs_number+1)].update({'StartingSaveFiles':file_path_list[seqs_number].rstrip()})
                        else:
                            if smp:
                                subprocess.call(["TurboFold-smp",conf_file],env=my_env)
                                if not str(seqs_number+1) in self.Seqname_dict.keys():
                                    self.Seqname_dict[str(seqs_number+1)] = {'StartingSaveFiles':file_path_list[seqs_number].rstrip()}
                                else:
                                    self.Seqname_dict[str(seqs_number+1)].update({'StartingSaveFiles':file_path_list[seqs_number].rstrip()})
                            else:   
                                subprocess.call(["TurboFold",conf_file],env=my_env)
                                if not str(seqs_number+1) in self.Seqname_dict.keys():
                                    self.Seqname_dict[str(seqs_number+1)] = {'StartingSaveFiles':file_path_list[seqs_number].rstrip()}
                                else:
                                    self.Seqname_dict[str(seqs_number+1)].update({'StartingSaveFiles':file_path_list[seqs_number].rstrip()})
        for num in range(len([ i for i in self.algn.keys()])):
            name = [ i for i in self.algn.keys()][num]
            num_dic = str(num+1)
            self.algn_to_Seqname_file[name] = num_dic
    
    
    def read_algn(self,filename,name_dict=None):
            """This function is used to read the output sequence alignment, which is
        in ClustalW format, then create a dictionary to record name of HM sequences 
        and the sequenceself.EX
            Argument:
                 conf_file: the path for the config file
            Returns:
                 algn: a dictionary of homlogous sequences
        
        """
            
            A=open(filename,'r').readlines()
            
            
            for line in A:
                Oneline=line.split()
                
                
                if len(Oneline)!=2:
                    continue
                
                if name_dict:
                    name = name_dict(Oneline[0])
                
                
                else:
                    name= Oneline[0]
                
                if name in self.algn.keys(): 
                    self.algn[name] = self.algn[name] + Oneline[1]
                else:
                    self.algn[name] = Oneline[1]
            
            
            
            return self.algn
    
    
    
    
    def match_soc(self,seq1,seq2):
        """match_soc 
        this function will use two sequences to calculate the pairing socres and unpairing scores, by using
        the equation: i = a nucleotide at position i, k the same 
                    PS(i,k) = square root of (i pair up * k pair up + i pair down * k pair down)
                    US(i,k) = square root(i unpair *k unpair)
        
        Arguments:
            seq1 {str}:
                the sequence is manually given
            
            seq2 {str}:
                 the sequnece is manually given, make sure not the same as the sequence 1
        
        Returns: paring socre[keys: name of sequences]=(value: paring socre of positions)
            unparing socre[keys: name of sequences]=(value: unparing socre of positions) 
        
        
        """
        # Initialize match score vector with zeros. self.algn is dicitonary of alignment of sequences
        alpha_1= [0 for i in range(len(self.algn[seq1]))]
        alpha_2= [0 for i in range(len(self.algn[seq1]))]
        
        
        Pairwise_algn, Pair_i = self.get_pairwise(seq1,seq2)
        Position_seq_1 = [ k for k in self.algn[seq1]]
        Position_seq_2 = [ k for k in self.algn[seq2]]
        total_id = 0
        match_id = 0
        for i in range(len(Position_seq_1)):
            P1 = Position_seq_1[i]
            P2 = Position_seq_2[i]
            if P1 != "-":
                if P2 !="-":
                    if P1 == P2:
                        total_id+=1
                        match_id+=1                        
                    elif P1 != P2:
                        total_id+=1
        
        
        # use i to get all the position of sequence 
        for i in range(len(self.algn[seq1])):
            if Pairwise_algn[seq1][i] != 0:
                a = Pairwise_algn[seq1][i]-1
                
                try:
                    
                    try:
                        self.probabilities[self.algn_to_Seqname_file[seq1]].prob_up(Pairwise_algn[seq2][a])
                    except:
                        print(i)
                        print(seq1)
                    try:
                        self.probabilities[self.algn_to_Seqname_file[seq2]].prob_up(a)
                    except:
                        print(seq2)
                        
                        print(i)
                    #M1 is pairing score
                    M1 = sqrt(self.probabilities[self.algn_to_Seqname_file[seq1]].prob_up(Pairwise_algn[seq2][a]-1)*self.probabilities[self.algn_to_Seqname_file[seq2]].prob_up(a))
                    M1 += sqrt(self.probabilities[self.algn_to_Seqname_file[seq1]].prob_down(Pairwise_algn[seq2][a]-1)*self.probabilities[self.algn_to_Seqname_file[seq2]].prob_down(a))
                    
                    #A2 is unpairing score
                    A2 = sqrt(self.probabilities[self.algn_to_Seqname_file[seq1]].unpair(Pairwise_algn[seq2][a]-1)*self.probabilities[self.algn_to_Seqname_file[seq2]].unpair(a))

#                    if M1 != 0 or A2 != 0:
                    alpha_1[i] = M1
                    alpha_2 [i] = A2
                
                
                
                except KeyError:
                    print('missing seq')

#            
        return alpha_1, alpha_2
    
    def all_match(self):
        """
        Get all the pairing and unpairing score for each sequence
        """
        
        #get sequence name in alignment
        seq_list = [i for i in self.algn.keys()]
        N= len(seq_list)
        #run through all the sequences
        for i in range(N):
            sequence_1 = seq_list[i]
            #run through all the sequences other than i
            for m in range(N):
                if m != i:
                    alpha_1, alpha_2 = self.match_soc(sequence_1,seq_list[m])
                    if sequence_1 in self.alpha_1.keys():
                        #add all the pairing and unpairing score for aligned position
                        for p in range(len(self.alpha_1[sequence_1])):
                            self.alpha_1[sequence_1][p] = self.alpha_1[sequence_1][p] + alpha_1[p]
                            self.alpha_2[sequence_1][p] = self.alpha_2[sequence_1][p] + alpha_2[p]
                    else:
                        self.alpha_1[sequence_1] = alpha_1
                        self.alpha_2[sequence_1] = alpha_2
            # normalized all the score with number of sequences in calculation -1  because sequence i is not calculated score with itself
            for p in range(len(self.alpha_1[sequence_1])):
                self.alpha_1[sequence_1][p] = self.alpha_1[sequence_1][p]/(N-1)
                self.alpha_2[sequence_1][p] = self.alpha_2[sequence_1][p]/(N-1)
            #bins the score into 30 bins
            span = 1
            bin_width = 1/30
            bins = 30
            bin_points=[]
            upper = 1
            lower = 0
            Total_count = 0
            for a in range(bins+1):
                bin_points.append(lower+ a*bin_width)
            Z_value = [[0 for n in range(bins+1)] for j in range(bins+1)]
            
            #created a heatmap for PS and US
            X = [i for i in self.alpha_1[sequence_1]]
            Y = [i for i in self.alpha_2[sequence_1]]
            for k, m  in zip(X,Y):    
                index_X = floor((k-lower)/span*bins)
                index_Y = floor((m-lower)/span*bins)
                if k == upper:
                    if m ==upper:
                        Z_value[index_X-1][index_Y-1] += 1
                        Total_count += 1
                    else:
                        Z_value[index_X-1][index_Y] += 1
                        Total_count += 1
                else:
                    if m ==upper:
                        Z_value[index_X][index_Y-1] += 1
                        Total_count += 1
                    else:
                        Z_value[index_X][index_Y] += 1
                        Total_count += 1
            Z_value_final = []
            Z_value_heatmap = [[0 for n in range(bins+1)] for j in range(bins+1)]
            for i in range(len(Z_value)):
                for j in range(len(Z_value)):
                    Z_value_heatmap[i][j] = Z_value[i][j]/Total_count*100
            self.heatmap[sequence_1] = Z_value_heatmap
            #finished the heatmap
            for o in range(len(Z_value)-1):
                for p in range(len(Z_value)-1-o):
                    Z_value_final.append(Z_value[o][p])
            
            self.Z_value[sequence_1] = Z_value_final
    
    
    def get_pairwise(self, seq1, seq2):
        """ seq1 is the name of frist sequence,seq 2 is the name of second sequence. 
        The function randomly pick two sequence from
        the aliment file and compare the alignment.
        Arguments:
            seq1 {str}:
                the sequence randomly picked by the random pick function, or
                manually given
            
            seq2 {str}:
                 the sequnece randomly picked by the random pick function, or
                manually given, make sure not the same as the sequence 1
        
        Returns: Pairwise: a dicitonary have the sequence 1 's aligned 
                           base pair's position in seuqence 2. so as the sequence 2
                 Pair_i:   a dicitonary have the map of sequence 1 and seuqence 2.
                           No - in the sqeuence.
        """
        #initialize for Pairwise and Pair_i
        Pairwise=dict()
        Pair_i=dict()
        #Error_Check, will report if the seq1 and seq2 are exist or not
        try:
            Seq1_a=self.algn[seq1]
        except:
            print("no seq1")
            return
        try:
            Seq2_a=self.algn[seq2]
        except:
            print("no seq2")
            return
        #get rid of the space at the end so that we can have the true sequence
        Seq1_b= (Seq1_a.strip())
        Seq2_b=(Seq2_a.strip())
        #initialize a list for poisitions of basepair in sequence
        j=0
        k=0
        Seq1_c=list()
        Seq2_c=list()
        Seq1_m=list()
        Seq2_m=list()
        #the for loop will go thourgh all the position of the sequence
        #and put the aligned base pair into the dictionary
        #so as the poisiton information for each sequence to Pair_i
        for i in range(len(Seq1_b)):
            if Seq1_b[i]!='-':
                j+=1
            
            if Seq2_b[i]!='-':
                k+=1
            
            if Seq1_b[i] =='-':
                   Seq1_c.append(0)
            else:
                if Seq2_b[i]!='-':
                    Seq1_c.append(k)
                
                else:
                    Seq1_c.append(0)
                
                
                Seq1_m.append(i)
            
            if Seq2_b[i] =='-':
                    Seq2_c.append(0)
            
            else:
                if Seq1_b[i]!='-':
                    Seq2_c.append(j)                   
                else:
                    Seq2_c.append(0)
                
                Seq2_m.append(i)
        
        Pairwise[seq2] = Seq2_c
        Pairwise[seq1] = Seq1_c
        Pair_i[seq2] = Seq2_m
        Pair_i[seq1] = Seq1_m
        
        return Pairwise, Pair_i 
    def enf2_feature(self):
            self.efn2 =dict()
            Seq_name = [i for i in self.algn.keys()]
            for i in range(len(self.algn.keys())):
                ct = self.Seqname_dict[str(i+1)]['CT_file']
                subprocess.call(["efn2",ct,'foo.efn2'],env=my_env)
                e_file = open ('foo.efn2','r')
                os.remove('foo.efn2')
                for line in e_file:
                    energy = line.split('=')[-1]
                    energy = energy.split('±')[0]
                    self.efn2[Seq_name[i]] = float(energy)
    def decoy_identify(self):
            '''add all the feature together,and identify the decoy sequences'''
            Seq_name = [i for i in self.algn.keys()]
            self.enf2_feature()
            
            hyper_feature=[]
            for i in range(len(Seq_name)):
                feature_list = []
                target_count = self.Z_value[Seq_name[i]]
                sequence_length = sum(target_count)
                pseduco_m=[0 for k in target_count]
                counts_rest = [0 for k in target_count]
                pseduco_rest = [0 for k in target_count]
                for m in range(len(Seq_name)):
                    if m !=i:
                        for counts in range(len(target_count)):
                            counts_rest[counts] += self.Z_value[Seq_name[m]][counts]/(len(self.algn.keys())-1)
                sequence_length_rest = sum(counts_rest)
                for counts in range(len(target_count)):
                    pseduco_m[counts] += (target_count[counts]+1/sequence_length)/sequence_length
                    pseduco_rest[counts] += (counts_rest[counts]+1/sequence_length_rest)/sequence_length_rest
                KL_Score = 0
                for k in range(len(pseduco_m)):
                    KL_Score += pseduco_m[k]*np.log(pseduco_m[k]/pseduco_rest[k])
                
                feature_list.append(KL_Score)
                
                energy_list = []
                energy_total = 0
                for m in range(len(Seq_name)):
                    if i != m:
                        energy_total += float(self.efn2[Seq_name[m]])
                        energy_list.append(self.efn2[Seq_name[m]])
                efn2_mean = energy_total/(len(self.algn.keys())-1)
                efn2_standard_distrubution = np.std(energy_list)
                z_score = (self.efn2[Seq_name[i]]-efn2_mean)/efn2_standard_distrubution
                feature_list.append(z_score)
                
                pfs_file = self.Seqname_dict[str(i+1)]['Save_file']
                subprocess.call(["ProbabilityPlot",pfs_file,'foo.pp','--text'],env=my_env)
                
                pp_data = open('foo.pp',"r").readlines()
                shannon_entropy = 0
                for line in pp_data:
                    na = line.split("\t")
                    if len(na) == 1:
                        sequence_length = float(na[0].rstrip())
                    if len(na) == 3:
                        if na[0] != "i":
                            probability_log = -float(na[2].rstrip())
                            probability = 10**(probability_log)
                            shannon_entropy += -(probability*probability_log)
                shannon_entropy = shannon_entropy/sequence_length
                os.remove('foo.pp')
                feature_list.append(shannon_entropy)
                
                aligment_entropy = aligment_shannon_entropy(self.algn_path)
                feature_list.append(aligment_entropy[Seq_name[i]])
                pfs_file_s = self.Seqname_dict[str(i+1)]['StartingSaveFiles']
                subprocess.call(["ProbabilityPlot",pfs_file_s,'foo.pp','--text'],env=my_env)
                
                pp_data = open('foo.pp',"r").readlines()
                shannon_entropy_s = 0
                for line in pp_data:
                    na = line.split("\t")
                    if len(na) == 1:
                        sequence_length_s = float(na[0].rstrip())
                    if len(na) == 3:
                        if na[0] != "i":
                            probability_log = -float(na[2].rstrip())
                            probability = 10**(probability_log)
                            shannon_entropy_s += -(probability*probability_log)
                shannon_entropy_s = shannon_entropy_s/sequence_length_s
                os.remove('foo.pp')
                feature_list.append(shannon_entropy_s-shannon_entropy)
                feature_list.append(shannon_entropy_s)
                hyper_feature.append(feature_list)
            
            self.hyper = hyper_feature

def aligment_shannon_entropy(algn_file = None):
    '''input: aligment file
       calculate the aligment shannon entropy
    '''
    
    out_aln_dict = dict()
    out_entropy = dict()
    
    
#    Position_probability=list()
    out_aln_open = open(algn_file,"r").readlines()
    for line in out_aln_open:
                Oneline=line.split()
        
            #print(Oneline[0])
                if len(Oneline)!=2:
                    continue
            
                name= Oneline[0]
                
                nucleotide = [l for l in Oneline[1].split()[0]]
                if name in out_aln_dict.keys():
                    out_aln_dict[name] += nucleotide
                else:
                    out_aln_dict[name] = nucleotide
                    
    nucleotide_list = [i for i in out_aln_dict.values()]
    nucleotide_name = [i for i in out_aln_dict.keys()]
    Position_probability = [ 0 for i in nucleotide_list[0]]
    
    for i in range(len(nucleotide_list[0])):
        A = 0
        U = 0
        C = 0
        G = 0 
        gap = 0
        for b in range(len(nucleotide_list)):
            if nucleotide_list[b][i] == 'A':
                A +=1
            elif nucleotide_list[b][i] == 'U':
                U +=1
            elif nucleotide_list[b][i] == 'C':
                C +=1
            elif nucleotide_list[b][i] == 'G':
                G +=1
            elif nucleotide_list[b][i] == '-':
                gap +=1
                
        Position_probability[i] = ({'A':A/len(nucleotide_list),'U':U/len(nucleotide_list),'C':C/len(nucleotide_list),'G':G/len(nucleotide_list),'-':gap/len(nucleotide_list)})
            
        
    for i in nucleotide_name:     
        list_of_sequence = out_aln_dict[i]
        out_entropy_map = [ 0 for i in nucleotide_list[0]]
        for b in range(len(list_of_sequence)):
            try:
                out_entropy_map[b] += -(Position_probability[b][list_of_sequence[b]]*(log(Position_probability[b][list_of_sequence[b]])))
            except:
                out_entropy_map[b] += 0
        out_entro = sum(out_entropy_map)/len(out_entropy_map)
        out_entropy[i]=out_entro 

        
    return out_entropy


def generate_features(conf_folder = None):
    '''Input: directory of all configure file 
        this function will generate all the feature for training '''
    conf_file = find('conf.txt',conf_folder)
    '''list with all feature information'''
    hyper_feature_list = []
    '''list with homologus and decoy information'''
    identify = []
#    undone_list=[]
    '''get thro all configure files and add the feature to a list'''
    for i in conf_file:
        conf_file_dir = os.path.dirname(i)
        contamin_file = find('contamin_seq.txt',conf_file_dir)[0]
        contamin_Seq_name_file = open(contamin_file,'r').readlines()
        contamin_Seq_name = []
        if len(contamin_Seq_name_file) != 0:
            for con_seq in contamin_Seq_name_file:
                if con_seq != '\n':
                    contamin_Seq_name.append(con_seq.rstrip())
#        try:            
        AB = DecoyFinder_feature_model_trainning(i)
        seq_name_list = [ i for i in AB.algn.keys()]
        for seq_name in seq_name_list:
            if not seq_name in contamin_Seq_name:
                identify.append(0)
            else:
                identify.append(1)
        hyper_feature_list.extend(AB.hyper)
#        except:
#            undone_list.append(i)
    return hyper_feature_list,identify

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='DecoyFinder train model',description='This script will generate configure file for DecoyFinder to train and output the trained AdaBoost model')
    parser.add_argument('sequence_file', help='sequence file directory')
    parser.add_argument('output', help='configure file output directory')
    parser.add_argument('Model_output', help='AdaBoost trained file output directory')
    args = parser.parse_args()
    '''generate all the required file'''
    DecoyFinder_training(args.sequence_file,args.output_path)
    '''train AdaBoost model'''
    Trian_feature_list,Trian_identity = generate_features(args.output_path)
    Adaclassifer_model = AdaBoostClassifier(n_estimators=500,algorithm='SAMME')
    Adaclassifer_model.fit(Trian_feature_list,Trian_identity)
    dump(AdaBoostClassifier,args.Model_output)