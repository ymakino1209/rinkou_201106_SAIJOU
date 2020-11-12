#!/usr/bin/env python
# coding: utf-8

# # 発展課題（Python使用）
# - 各遺伝子について、exon領域を “start-end,start-end,…”と出力する処理を実装。NIPBLで例示。
# - 遺伝子のタイプごとにexon数をカウントし、分布をヒストグラムで表示する。
# 
# - (とりあえずここまで)
# 
# - それぞれの遺伝子ファイルについて、「何らかの方法で」遺伝子の重複を除いた修正ファイルを作成する。
# - ２つの遺伝子ファイル間で共通する遺伝子、共通でない遺伝子がいくつ存在するかカウントする。
# - ２つの遺伝子ファイル間で共通する遺伝子について、 refFlat.hg38.txtにgene typeの列を追加したファイルを作成する。

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# ### 各遺伝子について、exon領域を “start-end,start-end,…”と出力する処理を実装。NIPBLで例示。

# In[2]:


UCSC = pd.read_table("refFlat.hg38.txt", skiprows=1, names=('geneName', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds'))
UCSC
# そのまま読み込むと、はじめのカラムが#geneNameという名前になってしまう。
# そこで、header部分は読み込まず、namesを使って自分で付けることにした。


# In[3]:


# "geneName" 列の値が NIPBLと一致する行をdfとして取得。
# .copy()にしないと、attribute errorとなる。
NIPBL_UCSC_df = UCSC.query("geneName == 'NIPBL'").copy()
NIPBL_UCSC_df

# cdsEndの異なる二種類が表示される


# In[4]:


# exonStartsごとに、","で分割してみる。仮の"SplitStars"に値を入力。
NIPBL_UCSC_df["SplitStarts"] = NIPBL_UCSC_df["exonStarts"].str.split(",")
NIPBL_UCSC_df["SplitStarts"] 


# In[5]:


# exonEndsごとに、","で分割してみる。仮の"SplitEnds"に値を入力。
NIPBL_UCSC_df["SplitEnds"] = NIPBL_UCSC_df["exonEnds"].str.split(",")
NIPBL_UCSC_df["SplitEnds"] 


# In[6]:


# exon1の開始点は.str.get(0)で取得できる
NIPBL_UCSC_df["exon1"] = NIPBL_UCSC_df["SplitStarts"].str.get(0)
NIPBL_UCSC_df["exon1"]


# In[7]:


# exon2の開始点は.str.get(1)で取得できる
NIPBL_UCSC_df["exon2"] = NIPBL_UCSC_df["SplitStarts"].str.get(1)
NIPBL_UCSC_df["exon2"]


# In[8]:


# exonCountsを取り出す方法
print(NIPBL_UCSC_df.at[19085, 'exonCount'])


# In[9]:


# for文で、exonCountsの数だけstart-endを表示する

for i in range((NIPBL_UCSC_df.at[19086, 'exonCount'])): #19086行のほうが多いのでこちらを使用
    print ("exon"+ str (i+1))
    print(NIPBL_UCSC_df["SplitStarts"].str.get(i) + "-" + NIPBL_UCSC_df["SplitEnds"].str.get(i))


# - 最低限こなしたが
# - 19085行と19086行に分けたほうがいいかも知れない
# - 出力を","でつないで、一行で出すべきだったかも

# ### 遺伝子のタイプごとにexon数をカウントし、分布をヒストグラムで表示する。

# In[10]:


# typeのあるデータEnsemblを読み込む
Ensembl = pd.read_table("refFlat.GRCh38.gene.name.txt", names=('geneName', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'geneType', 'type'))
# Ensemblにはheaderがないので、names = でheaderを指定
Ensembl


# In[11]:


Ensembl_s = Ensembl.sort_values(['type']) # typeごとにソートする
df = Ensembl_s[["type", "exonCount"]] # 必要な部分のみ取り出す
df


# In[12]:


df.hist(by="type", figsize=(60, 40))


# - 最低限こなしたが
# - 見づらい
# - 42種類多すぎ

# - lncRNA
# - miRNA
# - protein_coding
# - の３つくらいでいいかも

# ### それぞれの遺伝子ファイルについて、「何らかの方法で」遺伝子の重複を除いた修正ファイルを作成する。

# - 重複の原因はisoformなど。
# - 最も良いのは、代表的な転写産物を残すこと（最も広く発現している、最も発現量が多い）だが、それは知る由もない。
# - 次善の策として、最も長い転写産物が機能を持っていそう、、、ということで、「長いものを残す」。

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[13]:


# "geneName" 列の値が NIPBLと一致する行を表示
NIPBL_Ensembl_df = Ensembl.query('geneName == "NIPBL"')
NIPBL_Ensembl_df
# Ensemblは代表的なIsoformを抜き出しているので、１行のみ

