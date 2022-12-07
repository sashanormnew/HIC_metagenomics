#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys 
import os



def with_assembly_info(plasflow,PlDB_all,VirV,MobR,assembly,Name,MobT):
    MobT = MobT.loc[MobT["rep_type(s)"]!="-"]
    MobT["sample_id"] = [i.split()[0] for i in MobT["sample_id"]]
    MobT = MobT [["sample_id","rep_type(s)","predicted_mobility"]]
    MobT = MobT.rename(columns = {"sample_id":"contig_name"})

    #print(MobT.head())
    assembly = assembly.rename(columns = {0:"assembly_name",1:"contig_name",2:"assembly_id",3:"assembly_tax"})
    assembly["assembly_tax"] = ["_".join(i.split()) for i in assembly["assembly_tax"]]
    
    
    MobR["NODE"] = MobR["contig_id"]
    #MobR["contig_id"] = [i.split()[0] for i in MobR["contig_id"]]
    #print(MobR)
    
    #PlDB_all = PlDB_all.loc[(PlDB_all[5]>=80) & (PlDB_all[6]>=60)]
    PlDB_all=PlDB_all.drop_duplicates(subset=[0],inplace=False)
    PlDB_all = PlDB_all.loc[(PlDB_all[5]>=90) & (PlDB_all[6]>=80)] 
    
    plasflow_pl = plasflow.loc[plasflow.label.str.contains("plasmid")][["contig_name","contig_length","label"]]
    
    MobR["molecule_type"] = ["plasmid" if i!="-" else "-" for i in MobR["mash_neighbor_identification"]]
    
    plasflow_pl = plasflow_pl[["contig_name","label"]]
    plasflow_pl = plasflow_pl.rename(columns = {'label':'Plasflow'})
    
    PlDB_all = PlDB_all.rename(columns = {0:'contig_name'})
    PlDB_all["NCBI_id0.8cov0.6"]=PlDB_all[8]
    PlDB_all = PlDB_all[["contig_name","NCBI_id0.8cov0.6"]]
    
    Vir=VirV.rename(columns = {'Contig name':'contig_name'})
    MobR = MobR.rename(columns = {'contig_id':'contig_name'})
    
    Vir = Vir[["contig_name","Prediction","Length","Circular"]]
    Vir = Vir.rename(columns = {'Prediction':'ViralVerify'})
    
    MobR = MobR[["contig_name","molecule_type","NODE","mash_neighbor_identification"]]
    MobR = MobR.rename(columns = {'molecule_type':'Mob_recon'})
    a= pd.merge(pd.merge(pd.merge(PlDB_all,Vir,how='outer')
                         ,MobR,how='outer'),plasflow_pl,how='outer')
    #a1 = pd.merge(a,assembly,how='outer')  
    a1 = pd.merge(pd.merge(a,assembly,how='outer'),MobT,how='outer')
    
    a1 = a1[["contig_name","NODE","NCBI_id0.8cov0.6","Plasflow","Mob_recon","ViralVerify","Circular","Length","mash_neighbor_identification","assembly_name","assembly_tax","rep_type(s)","predicted_mobility"]]
   
    nn = []

    a1["rep_type(s)"] = a1["rep_type(s)"].fillna("Nope")
    for n in range(a1.shape[0]):
        l = list(a1.iloc[n,[1,2,3,4]])
        pl = " ".join([str(k) for k in l]).count("plasmid")
        Pl = " ".join([str(k) for k in l]).count("Plasmid")
        pl_gen = Pl+pl
        if Pl>0 or pl_gen>2 or list(a1.iloc[n,[-2]])[0]!="Nope":
            nn.append(n)
    
    a11 = a1.iloc[nn]
    
    a11 = a11.fillna("Nope")
    
    a11["assembly_Length"]= [int(i.split("_")[3]) if i!='Nope' else 0 for i in a11["assembly_name"]]
    a11["assembly_name"] =["_".join(i.split("_")[0:2]) if i!='Nope' else i for i in a11["assembly_name"]]
    
    a11["assembly_contig_Length"] = np.where(a11["assembly_name"]!="Nope", a11["assembly_Length"], a11["Length"])
    a11["contig_ass_name"] = np.where(a11["assembly_name"]!="Nope", a11["assembly_name"], a11["contig_name"])
    
    
    a11["Plasmid_annotation"] = ["_".join(k.split()[1:3]) for k in a11["NCBI_id0.8cov0.6"]]
    a11["ViralVerify"] = ["_".join(k.split()) for k in a11["ViralVerify"]] 
    pl= []
    for i in range(a11.shape[0]):
        n1,n2,n3 = a11.iloc[i,[0,6,9]]
        if n3.startswith("RN") or "+" in n2:
            pl.append("Plasmid_full_assembly")
        else:
            pl.append("Plamid_contig") 
    a11["Plasmid_status"] = pl
    
    a11 = a11.drop(columns = ["Plasflow","Circular","NCBI_id0.8cov0.6","assembly_Length","mash_neighbor_identification"])
    return a11
    
    


# In[120]:


def without_assembly_info(plasflow,PlDB_all,VirV,MobR,Name):
    MobR["NODE"] = ["".join(i.split()[1].split(":")[1:]) for i in MobR["contig_id"]]
    MobR["contig_id"] = [i.split()[0] for i in MobR["contig_id"]]
    
    PlDB_all = PlDB_all.loc[(PlDB_all[5]>=80) & (PlDB_all[6]>=60)]
    PlDB_all=PlDB_all.drop_duplicates(subset=[0],inplace=False)
    PlDB_all
    
    plasflow_pl = plasflow.loc[plasflow.label.str.contains("plasmid")][["contig_name","contig_length","label"]]
    MobR["molecule_type"] = ["plasmid" if i!="-" else "-" for i in MobR["mash_neighbor_identification"]]
    plasflow_pl = plasflow_pl[["contig_name","label"]]
    plasflow_pl = plasflow_pl.rename(columns = {'label':'Plasflow'})
    PlDB_all = PlDB_all.rename(columns = {0:'contig_name'})
    PlDB_all["NCBI_id0.8cov0.6"]=PlDB_all[8]
    PlDB_all = PlDB_all[["contig_name","NCBI_id0.8cov0.6"]]
    Vir=VirV.rename(columns = {'Contig name':'contig_name'})
    MobR = MobR.rename(columns = {'contig_id':'contig_name'})
    Vir = Vir[["contig_name","Prediction","Length","Circular"]]
    Vir = Vir.rename(columns = {'Prediction':'ViralVerify'})
    MobR = MobR[["contig_name","molecule_type","NODE","mash_neighbor_identification"]]
    MobR = MobR.rename(columns = {'molecule_type':'Mob_recon'})
    
    
    a= pd.merge(pd.merge(pd.merge(PlDB_all,Vir,how='outer')
                         ,MobR,how='outer'),plasflow_pl,how='outer')
    
    a1 = a[["contig_name","NODE","NCBI_id0.8cov0.6","Plasflow","Mob_recon","ViralVerify","Circular","Length"]]
    
    nn = []
    for n in range(a1.shape[0]):
        l = list(a1.iloc[n,[1,2,3,4]])
        pl = " ".join([str(k) for k in l]).count("plasmid")
        Pl = " ".join([str(k) for k in l]).count("Plasmid")
        pl_gen = Pl+pl
        if pl_gen>2:
            nn.append(n)
    
    a11 = a1.iloc[nn]
    a11 = a11.fillna("Nope")
    
    a11 = a11.loc[a11["NCBI_id0.8cov0.6"].str.contains("[^Uncultured]")]
    #a11["Name"] = ["_".join(k.split()[1:3]) for k in a11["NCBI_id0.8cov0.6"]]
    
    pl= []
    for i in range(a11.shape[0]):
        n1,n2 = a11.iloc[i,[0,6]]
        if "+" in n2:
            #pl.append("_".join(n3.split()))
            pl.append("Plasmid_full_assembly")
        else:
            #print(n1,n2)
            pl.append("Plamid_contig")
    a11["Plasmid_status"] = pl
    a11 = a11[["contig_name","NODE","Plasflow","Mob_recon","ViralVerify","Length","Plasmid_status","Name"]]
    a11 = a11.drop(columns = ["ViralVerify"])
    
    return a11
    


# In[ ]:


if os.stat(sys.argv[5]).st_size != 0:
    plasflow = pd.read_table(sys.argv[1])
    PlDB_all = pd.read_table(sys.argv[2],header=None)
    VirV= pd.read_csv(sys.argv[3],sep =",")
    MobR = pd.read_csv(sys.argv[4],sep="\t")
    assembly = pd.read_table(sys.argv[5],header=None)
    Name = sys.argv[6]
    MobT = pd.read_table(sys.argv[7]) 
    res = with_assembly_info(plasflow,PlDB_all,VirV,MobR,assembly,Name,MobT)
    name = f'{Name}_PLAS_OUT_mobr_VrV.txt'
    print(res.head())
    res.to_csv(name, header=None, index=None, sep='\t', mode='w')
    #print(res)

else:
    plasflow = pd.read_table(sys.argv[1])
    PlDB_all = pd.read_table(sys.argv[2],header=None)
    VirV= pd.read_csv(sys.argv[3],sep =",")
    MobR = pd.read_csv(sys.argv[4],sep="\t")
    Name = sys.argv[6]
    res = without_assembly_info(plasflow,PlDB_all,VirV,MobR,Name)
    name = f'{Name}_PLAS_OUT_mobr_VrV.txt'
    res.to_csv(name, header=None, index=None, sep='\t', mode='w')
    

