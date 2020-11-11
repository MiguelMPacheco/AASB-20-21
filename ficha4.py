# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 18:29:17 2020

@author: Dias
"""


def valida(seq):
    #seq = input("Insira uma sequencia: ").upper()
    dna = "ACGT"
    res = True
    pos = 0
    while pos < len(seq) and res:
        if seq.upper()[pos] not in dna:
            res = False
        else:
            pos += 1
    return res

def ficheirotxt(FileHandle):
    
    file = open(FileHandle, "r") 
    lines = file.readlines() 
    seq = ""
    for line in lines:
        seq += line.replace("\n","") 
    file.close() 
    return seq

def ficheirofasta(FileHandle):
    
    file = open(FileHandle, "r") 
    fasta_file = file.read() 
    fasta_list = fasta_file.split("\n") 
    seq_list = fasta_list[1:] 
    seq = "".join(seq_list) 
    return seq

#%%



def complemento_inverso(seq):
    if valida(seq):
        return (seq[::-1].upper().replace("A", "t").replace("T","a").replace("C","g").replace("G","c").upper())
    else:
        raise TypeError("Sequencia invalida")

def transcricao(seq):
    if valida(seq):
        return seq.upper().replace("T", "U")
    else:
        raise TypeError("Sequencia invalida")


def codoes(seq):
    
    if valida(seq):
        resultado = ""
        seq2 = list(seq.upper())
        
        while len(seq2) % 3 != 0:
            seq2.pop(-1)
        
        seq_nova = "".join(seq2)
        
     
        
        for letter in seq_nova:
            if (len(resultado) + 1) % 4 == 0:
                resultado += " " + letter
            else:
                resultado += letter
        
            
        return resultado
    
    else:
        raise TypeError("Sequencia invalida")

def traducao(seq):
    
    if valida(seq):
        
        gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
        
        
        seq_correct = codoes(seq).split()
        resultado = ""
        for codao in seq_correct:
            if codao in gencode.keys():
                resultado += gencode[codao]
        
        return resultado
    
    else:
        raise TypeError("Sequencia invalida")





#%%



def contar_bases(seq):
    
    dna = "ATCG"
    dic = {}
    if valida(seq):
        for letter in seq.upper():
            if letter not in dic:
                dic[letter] = 1
            else:    
                dic[letter] +=1
        
        for letter2 in dna:
            if letter2 not in dic:
                dic[letter2] = 0
                
        return dic
    
        
    else:
        raise TypeError("Sequencia invalida")
            



#%%
def reading_frames(seq):
    lista = []
    for i in range(3):
        lista.append(''.join(traducao(seq[i:])))
    for i in range (3):
        lista.append(traducao(complemento_inverso(seq)[i:]))
    return lista



#%%
def get_proteinas():
    prota = input().upper()
    print(*sorted(sorted(prota.split()), key = len, reverse = True))








