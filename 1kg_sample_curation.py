#!/usr/local/bin/python3
# this file is to revise HG00100 and HG00101 sample mixed in database
from pymongo import MongoClient
import json
import pandas as pd
import numpy as np
import math

def get_var(sample):
    content=[]
    with open('./imports/' +sample + '.txt','r') as f:
        lines = f.readlines()
        for line in lines:
            #line=line.replace(',','\t')

            if sample in line:
                line = line.replace('\n', '')
                column = line.replace('#','').split('\t')
                #print(column)
            if '#' not in line:
                line = line.replace('\n', '')
                #print(line.split('\t'))
                content.append(line.split('\t'))
    df=pd.DataFrame(content, columns=column)
    return df



def bio_cur(old,new):

    i = 0
    while i < new_var.shape[0]:
        new_var.at[i,'start'] = math.trunc(int(new_var.at[i,'POS'])-1)
        new_var.at[i, 'end'] = new_var.at[i, 'ID'].split(':')[3].split('-')[1]
        new_var.at[i, 'chr'] = new_var.at[i, 'ID'].split(':')[2].split('r')[1]
        if 'LOSS' in str(new_var.at[i, 'ID'].split(':')[1]):
            new_var.at[i, 'var'] = 'DEL'
        if 'GAIN' in new_var.at[i, 'ID'].split(':')[1]:
            new_var.at[i, 'var'] = 'DUP'
        i = i + 1
    return new_var

old='HG00101'
new='HG00100'
new_var=get_var(new)


client = MongoClient()
db = client['progenetix']
bios = db['biosamples'].find({'id':'onekgbs-'+old})
for bio_document in bios:
    bio_copy = bio_document
    bio_copy['id']='onekgbs-'+new
    bio_copy['individual_id']='onekgind-'+new
    bio_copy['info']['callset_ids'] = ['onekgcs-' + new]
    del bio_copy['_id']
#db['biosamples'].insert_one(bio_copy)

cs = db['callsets'].find({'id':'onekgcs-'+old})
for cs_document in cs:
    cs_copy = cs_document
    cs_copy['id']='onekgcs-'+new
    cs_copy['individual_id']='onekgind-'+new
    cs_copy['biosample_id'] = 'onekgbs-' + new

    del cs_copy['_id']
#db['callsets'].insert_one(cs_copy)

ind = db['individuals'].find({'id':'onekgind-'+old})
for ind_document in ind:
    ind_copy = ind_document
    ind_copy['id']='onekgind-'+new
    ind_copy['sex']={
		"id" : "PATO:0020002",
		"label" : "female genotypic sex"
	}


    del ind_copy['_id']
#db['individuals'].insert_one(ind_copy)
new_var=get_var(new)
new_var1=bio_cur(old,new)

i=0
while i<new_var1.shape[0]:
    variant_internal_id=str(new_var1.at[i,'chr'])+':'+str(math.trunc(new_var1.at[i,'start']))+'-'+str(new_var1.at[i,'end'])+':'+str(new_var1.at[i,'var'])
    vars = db['variants'].find_one({"biosample_id":'onekgbs-'+old})
    condition={'$and':[{"variant_internal_id":variant_internal_id},{"biosample_id":'onekgbs-'+old}]}

    var_copy = vars
    var_copy['biosample_id'] = 'onekgbs-' + new
    var_copy['individual_id'] = 'onekgind-' + new
    var_copy['callset_id'] = 'onekgcs-' + new
    var_copy["variant_internal_id"] = variant_internal_id
    var_copy['info'] = {
        "var_length": int(int(new_var1.at[i, 'end']) - int(new_var1.at[i, 'start'])),
        "cn_count": int(new_var1.at[i, new].split(':')[2])
        #"cnv_value": float(new_var1.at[i, new].split(':')[1])
    }
    var_copy['location']['interval'] = {
        "start": {
            "type": "Number",
            "value": int(new_var1.at[i, 'start'])
        },
        "end": {
            "type": "Number",
            "value": int(new_var1.at[i, 'end'])
        }
    }
    if 'DEL' in str(new_var1.at[i, 'var']):
        var_copy['relative_copy_class'] = "partial loss"
        var_copy['variant_state'] = {"id": "EFO:0030067", "label": "copy number loss"}
    else:
        var_copy['relative_copy_class'] = "low-level gain"
        var_copy['variant_state'] = {
            "id": "EFO:0030070",
            "label": "copy number gain"
        }
    del var_copy['_id']
    var={}
    var['start'] = int(new_var1.at[i, 'start'])
    var['end'] = int(new_var1.at[i, 'end'])
    var['variant_type']=str(new_var1.at[i, 'var'])
    var['reference_name']=str(new_var1.at[i,'chr'])
    var['info']={
        "var_length": int(int(new_var1.at[i, 'end']) - int(new_var1.at[i, 'start'])),
        "cnv_value": int(new_var1.at[i, new].split(':')[2])
    }
    var['biosample_id'] = 'onekgbs-' + new
    var['individual_id'] = 'onekgind-' + new
    var['callset_id'] = 'onekgcs-' + new
    var["digset"] = variant_internal_id

    #id=db['variants'].insert_one(var_copy).inserted_id
    id = db['variants'].insert_one(var).inserted_id
    db['variants'].update_one({'_id':id},{'$set':{'id':'onekgvar-'+str(id)}})


    
    i=i+1
